/*
Copyright (C) 2013-2014 Regents of the University of Minnesota.
This file is part of InMAP.

InMAP is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

InMAP is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with InMAP.  If not, see <http://www.gnu.org/licenses/>.
*/

package main

import (
	"bufio"
	"encoding/gob"
	"encoding/json"
	"flag"
	"fmt"
	"io/ioutil"
	//	"log"
	"math"
	"os"
	"path/filepath"
	"runtime"
	//	"sort"
	"strconv"
	"strings"
	"sync"
	"time"

	//	"text/tabwriter"

	//"github.com/ctessum/geomconv"
	//"github.com/ctessum/geomop"
	//"github.com/ctessum/inmap/lib.inmap"
	//"github.com/ctessum/shapefile"
	//"github.com/jonas-p/go-shp"
	//"github.com/patrick-higgins/rtreego"
	//"github.com/twpayne/gogeom/geom"

	"code.google.com/p/lvd.go/cdf"
	//	"devM/science"
	"devM/sparse"
)

var configFile = flag.String("config", "none", "Path to configuration file")

const version = "BasicBeta1ConcOnly"

type dataHolder struct {
	dims        []string
	Description string
	Units       string
	data        *sparse.DenseArray
}

// Chemical mass conversions
const (
	// grams per mole
	mwNOx = 46.0055
	mwN   = 14.0067
	mwNO3 = 62.00501
	mwNH3 = 17.03056
	mwNH4 = 18.03851
	mwS   = 32.0655
	mwSO2 = 64.0644
	mwSO4 = 96.0632
	// ratios
	NOxToN = mwN / mwNOx
	NtoNO3 = mwNO3 / mwN
	SOxToS = mwSO2 / mwS
	StoSO4 = mwS / mwSO4
	NH3ToN = mwN / mwNH3
	NtoNH4 = mwNH4 / mwN
)

//const tolerance = 0.0003 // tolerance for convergence
const tolerance = 0.91    // tolerance for convergence
const checkPeriod = 3600. // seconds, how often to check for convergence
const daysPerSecond = 1. / 3600. / 24.
const topLayerToCalc = 47 // The top layer to do calculations for

// These are the names of pollutants accepted as emissions (μg/s)
var EmisNames = []string{"VOC", "NOx", "NH3", "SOx", "PM2_5"}

// These are the names of pollutants within the model
var polNames = []string{"gOrg", "pOrg", // gaseous and particulate organic matter
	"PM2_5",      // PM2.5
	"gNH", "pNH", // gaseous and particulate N in ammonia
	"gS", "pS", // gaseous and particulate S in sulfur
	"gNO", "pNO"} // gaseous and particulate N in nitrate

var OutputVariables = []string{"VOC", "SOA", "PrimaryPM2_5", "NH3", "pNH4",
	"SOx", "pSO4", "NOx", "pNO3", "TotalPM2_5"}

//"Total deaths", "White deaths", "Non-white deaths",
//"High income deaths", "Low income deaths",
//"High income white deaths", "Low income non-white deaths"}

// Indicies of individual pollutants in arrays.
const (
	igOrg, ipOrg, iPM2_5, igNH, ipNH, igS, ipS, igNO, ipNO = 0, 1, 2, 3, 4, 5, 6, 7, 8
)

type InMAPdata struct {
	Data          []*Cell // One data holder for each grid cell
	Dt            float64 // seconds
	Nlayers       int     // number of model layers
	LayerStart    []int   // start index of each layer (inclusive)
	LayerEnd      []int   // end index of each layer (exclusive)
	NumIterations int     // number of itereations to limit run to
	westBoundary  []*Cell // boundary cells
	eastBoundary  []*Cell // boundary cells
	northBoundary []*Cell // boundary cells
	southBoundary []*Cell // boundary cells
	topBoundary   []*Cell // boundary cells; assume bottom boundary is the same as lowest layer
}

//  Set the time step using the Courant–Friedrichs–Lewy (CFL) condition.
func (d *InMAPdata) setTstepCFL() {
	const Cmax = 1.
	val := 0.
	for _, c := range d.Data {
		thisval := max(c.UPlusSpeed/c.Dx, c.UMinusSpeed/c.Dx,
			c.VPlusSpeed/c.Dy, c.VMinusSpeed/c.Dy,
			c.WPlusSpeed/c.Dz, c.WMinusSpeed/c.Dz)
		if thisval > val {
			val = thisval
		}
	}
	d.Dt = Cmax / math.Pow(3., 0.5) / val // seconds
}

// Run air quality model. Emissions are assumed to be in units
// of μg/s, and must only include the pollutants listed in "EmisNames".
// Output is in the form of map[pollutant][layer][row]concentration,
// in units of μg/m3.
func (d *InMAPdata) Run(emissions map[string][]float64) (
	outputConc map[string][][]float64) {

	startTime := time.Now()
	timeStepTime := time.Now()

	// Emissions: all except PM2.5 go to gas phase
	for pol, arr := range emissions {
		switch pol {
		case "VOC":
			d.addEmisFlux(arr, 1., igOrg)
		case "NOx":
			d.addEmisFlux(arr, NOxToN, igNO)
		case "NH3":
			d.addEmisFlux(arr, NH3ToN, igNH)
		case "SOx":
			d.addEmisFlux(arr, SOxToS, igS)
		case "PM2_5":
			d.addEmisFlux(arr, 1., iPM2_5)
		default:
			panic(fmt.Sprintf("Unknown emissions pollutant %v.", pol))
		}
	}

	oldSum := make([]float64, len(polNames))
	iteration := 0
	nDaysRun := 0.
	timeSinceLastCheck := 0.
	nprocs := runtime.GOMAXPROCS(0) // number of processors
	funcChan := make([]chan func(*Cell, *InMAPdata), nprocs)
	var wg sync.WaitGroup

	for procNum := 0; procNum < nprocs; procNum++ {
		funcChan[procNum] = make(chan func(*Cell, *InMAPdata), 1)
		// Start thread for concurrent computations
		go d.doScience(nprocs, procNum, funcChan[procNum], &wg)
	}

	// make list of science functions to run at each timestep
	scienceFuncs := []func(c *Cell, d *InMAPdata){
		func(c *Cell, d *InMAPdata) { c.addEmissionsFlux(d) },
		func(c *Cell, d *InMAPdata) {
			c.UpwindAdvection(d.Dt)
			c.Mixing(d.Dt)
			c.Chemistry(d)
			c.DryDeposition(d)
			c.WetDeposition(d.Dt)
		}}

	for { // Run main calculation loop until pollutant concentrations stabilize

		// Send all of the science functions to the concurrent
		// processors for calculating
		wg.Add(len(scienceFuncs) * nprocs)
		for _, function := range scienceFuncs {
			for pp := 0; pp < nprocs; pp++ {
				funcChan[pp] <- function
			}
		}

		// do some things while waiting for the science to finish
		iteration++
		nDaysRun += d.Dt * daysPerSecond
		fmt.Printf("马上。。。Iteration %-4d  walltime=%6.3gh  Δwalltime=%4.2gs  "+
			"timestep=%2.0fs  day=%.3g\n",
			iteration, time.Since(startTime).Hours(),
			time.Since(timeStepTime).Seconds(), d.Dt, nDaysRun)
		timeStepTime = time.Now()
		timeSinceLastCheck += d.Dt

		// bug check:
		bugRow := 8610
		fmt.Println("SO2 g", d.Data[bugRow].Cf[5], "SO2 p", d.Data[bugRow].Cf[6])
		fmt.Println("PM2_5", d.Data[bugRow].Cf[2])

		// Occasionally, check to see if the pollutant concentrations have converged
		if timeSinceLastCheck >= checkPeriod {
			wg.Wait() // Wait for the science to finish, only when we need to check
			// for convergence.
			timeToQuit := true
			timeSinceLastCheck = 0.
			for ii, pol := range polNames {
				var sum float64
				for _, c := range d.Data {
					sum += c.Cf[ii]
				}
				if !checkConvergence(sum, oldSum[ii], pol) {
					timeToQuit = false
				}
				if d.NumIterations > 0 && iteration > d.NumIterations {
					timeToQuit = true
				}
				oldSum[ii] = sum
			}
			if timeToQuit {
				break // leave calculation loop because we're finished
			}
		}
	}
	// Prepare output data
	outputConc = make(map[string][][]float64)
	for _, name := range OutputVariables {
		outputConc[name] = make([][]float64, d.Nlayers)
		for k := 0; k < d.Nlayers; k++ {
			outputConc[name][k] = d.toArrayBasic(name, k)
		}
	}
	return
}

// Calculate emissions flux given emissions array in units of μg/s
// and a scale for molecular mass conversion.
func (d *InMAPdata) addEmisFlux(arr []float64, scale float64, iPol int) {
	for row, val := range arr {
		fluxScale := 1. / d.Data[row].Dx / d.Data[row].Dy /
			d.Data[row].Dz // μg/s /m/m/m = μg/m3/s
		d.Data[row].emisFlux[iPol] = val * scale * fluxScale
	}
	return
}

// Carry out the atmospheric chemistry and physics calculations
func (d *InMAPdata) doScience(nprocs, procNum int,
	funcChan chan func(*Cell, *InMAPdata), wg *sync.WaitGroup) {
	var c *Cell
	for f := range funcChan {
		for ii := procNum; ii < len(d.Data); ii += nprocs {
			c = d.Data[ii]
			c.lock.Lock() // Lock the cell to avoid race conditions
			if c.Layer <= topLayerToCalc {
				f(c, d) // run function
			}
			c.lock.Unlock() // Unlock the cell: we're done editing it
		}
		wg.Done()
	}
}

// Convert the concentration data into a regular array
func (d *InMAPdata) toArrayBasic(pol string, layer int) []float64 {
	o := make([]float64, d.LayerEnd[layer]-d.LayerStart[layer])
	for i, c := range d.Data[d.LayerStart[layer]:d.LayerEnd[layer]] {
		c.lock.RLock()
		switch pol {
		case "VOC":
			o[i] = c.Cf[igOrg]
		case "SOA":
			o[i] = c.Cf[ipOrg]
		case "PrimaryPM2_5":
			o[i] = c.Cf[iPM2_5]
		case "TotalPM2_5":
			o[i] = c.Cf[iPM2_5] + c.Cf[ipOrg] + c.Cf[ipNH] + c.Cf[ipS] + c.Cf[ipNO]
		case "NH3":
			o[i] = c.Cf[igNH] / NH3ToN
		case "pNH4":
			o[i] = c.Cf[ipNH] * NtoNH4
		case "SOx":
			o[i] = c.Cf[igS] / SOxToS
		case "pSO4":
			o[i] = c.Cf[ipS] * StoSO4
		case "NOx":
			o[i] = c.Cf[igNO] / NOxToN
		case "pNO3":
			o[i] = c.Cf[ipNO] * NtoNO3
		case "VOCemissions":
			o[i] = c.emisFlux[igOrg]
		case "NOxemissions":
			o[i] = c.emisFlux[igNO]
		case "NH3emissions":
			o[i] = c.emisFlux[igNH]
		case "SOxemissions":
			o[i] = c.emisFlux[igS]
		case "PM2_5emissions":
			o[i] = c.emisFlux[iPM2_5]
		case "UPlusSpeed":
			o[i] = c.UPlusSpeed
		case "UMinusSpeed":
			o[i] = c.UMinusSpeed
		case "VPlusSpeed":
			o[i] = c.VPlusSpeed
		case "VMinusSpeed":
			o[i] = c.VMinusSpeed
		case "WPlusSpeed":
			o[i] = c.WPlusSpeed
		case "WMinusSpeed":
			o[i] = c.WMinusSpeed
		case "Organicpartitioning":
			o[i] = c.OrgPartitioning
		case "Sulfurpartitioning":
			o[i] = c.SPartitioning
		case "Nitratepartitioning":
			o[i] = c.NOPartitioning
		case "Ammoniapartitioning":
			o[i] = c.NHPartitioning
		case "Particlewetdeposition":
			o[i] = c.ParticleWetDep
		case "SO2wetdeposition":
			o[i] = c.SO2WetDep
		case "Non-SO2gaswetdeposition":
			o[i] = c.OtherGasWetDep
		case "Kxxyy":
			o[i] = c.Kxxyy
		case "Kzz":
			o[i] = c.Kzz
		case "M2u":
			o[i] = c.M2u
		case "M2d":
			o[i] = c.M2d
		case "PblTopLayer":
			o[i] = c.PblTopLayer
		default:
			panic(fmt.Sprintf("Unknown variable %v.", pol))
		}
		c.lock.RUnlock()
	}
	return o
}

// Data for a single grid cell
type Cell struct {
	//	Geom                           geom.T       // Cell geometry
	//	WebMapGeom                     geom.T       // Cell geometry in web map (mercator) coordinate system
	UPlusSpeed, UMinusSpeed        float64      // [m/s]
	VPlusSpeed, VMinusSpeed        float64      // [m/s]
	WPlusSpeed, WMinusSpeed        float64      // [m/s]
	OrgPartitioning, SPartitioning float64      // gaseous fraction
	NOPartitioning, NHPartitioning float64      // gaseous fraction
	ParticleWetDep, SO2WetDep      float64      // wet deposition rate [1/s]
	OtherGasWetDep                 float64      // wet deposition rate [1/s]
	ParticleDryDep, NH3DryDep      float64      // Dry deposition velocities [m/s]
	SO2DryDep, VOCDryDep           float64      // Dry deposition velocities [m/s]
	NOxDryDep                      float64      // Dry deposition velocities [m/s]
	SO2oxidation                   float64      // SO2 oxidation to SO4 by HO and H2O2 [1/s]
	Kzz                            float64      // Grid center vertical diffusivity after applying convective fraction [m2/s]
	KzzAbove, KzzBelow             []float64    // horizontal diffusivity [m2/s] (staggered grid)
	Kxxyy                          float64      // Grid center horizontal diffusivity [m2/s]
	KyySouth, KyyNorth             []float64    // horizontal diffusivity [m2/s] (staggered grid)
	KxxWest, KxxEast               []float64    // horizontal diffusivity at [m2/s] (staggered grid)
	M2u                            float64      // ACM2 upward mixing (Pleim 2007) [1/s]
	M2d                            float64      // ACM2 downward mixing (Pleim 2007) [1/s]
	TotalPop, WhitePop             float64      // Population [people/grid cell]
	TotalPoor, WhitePoor           float64      // Poor population [people/grid cell]
	AllCauseMortality              float64      // Mortalities per 100,000 people per year
	PblTopLayer                    float64      // k index of boundary layer top
	Dx, Dy, Dz                     float64      // grid size [meters]
	Volume                         float64      // [cubic meters]
	Row                            int          // master cell index
	Ci                             []float64    // concentrations at beginning of time step [μg/m3]
	Cf                             []float64    // concentrations at end of time step [μg/m3]
	emisFlux                       []float64    //  emissions [μg/m3/s]
	West                           []*Cell      // Neighbors to the East
	East                           []*Cell      // Neighbors to the West
	South                          []*Cell      // Neighbors to the South
	North                          []*Cell      // Neighbors to the North
	Below                          []*Cell      // Neighbors below
	Above                          []*Cell      // Neighbors above
	GroundLevel                    []*Cell      // Neighbors at ground level
	WestFrac, EastFrac             []float64    // Fraction of cell covered by each neighbor (adds up to 1).
	NorthFrac, SouthFrac           []float64    // Fraction of cell covered by each neighbor (adds up to 1).
	AboveFrac, BelowFrac           []float64    // Fraction of cell covered by each neighbor (adds up to 1).
	GroundLevelFrac                []float64    // Fraction of cell above to each ground level cell (adds up to 1).
	IWest                          []int        // Row indexes of neighbors to the East
	IEast                          []int        // Row indexes of neighbors to the West
	ISouth                         []int        // Row indexes of neighbors to the South
	INorth                         []int        // Row indexes of neighbors to the north
	IBelow                         []int        // Row indexes of neighbors below
	IAbove                         []int        // Row indexes of neighbors above
	IGroundLevel                   []int        // Row indexes of neighbors at ground level
	DxPlusHalf                     []float64    // Distance between centers of cell and East [m]
	DxMinusHalf                    []float64    // Distance between centers of cell and West [m]
	DyPlusHalf                     []float64    // Distance between centers of cell and North [m]
	DyMinusHalf                    []float64    // Distance between centers of cell and South [m]
	DzPlusHalf                     []float64    // Distance between centers of cell and Above [m]
	DzMinusHalf                    []float64    // Distance between centers of cell and Below [m]
	Layer                          int          // layer index of grid cell
	Temperature                    float64      // Average temperature, K
	WindSpeed                      float64      // RMS wind speed, [m/s]
	S1                             float64      // stability parameter [?]
	SClass                         float64      // stability class: "0=Unstable; 1=Stable
	lock                           sync.RWMutex // Avoid cell being written by one subroutine and read by another at the same time.
}

func (c *Cell) prepare() {
	c.Volume = c.Dx * c.Dy * c.Dz
	c.Ci = make([]float64, len(polNames))
	c.Cf = make([]float64, len(polNames))
	c.emisFlux = make([]float64, len(polNames))
}

func (c *Cell) makecopy() *Cell {
	c2 := new(Cell)
	c2.Dx, c2.Dy, c2.Dz = c.Dx, c.Dy, c.Dz
	c2.Kxxyy = c.Kxxyy
	c2.prepare()
	return c2
}

// Calculate center-to-center cell distance,
// fractions of grid cell covered by each neighbor
// and harmonic mean staggered-grid diffusivities.
func (cell *Cell) neighborInfo() {
	cell.DxPlusHalf = make([]float64, len(cell.East))
	cell.EastFrac = make([]float64, len(cell.East))
	cell.KxxEast = make([]float64, len(cell.East))
	for i, c := range cell.East {
		cell.DxPlusHalf[i] = (cell.Dx + c.Dx) / 2.
		cell.EastFrac[i] = min(c.Dy/cell.Dy, 1.)
		cell.KxxEast[i] = harmonicMean(cell.Kxxyy, c.Kxxyy)
	}
	cell.DxMinusHalf = make([]float64, len(cell.West))
	cell.WestFrac = make([]float64, len(cell.West))
	cell.KxxWest = make([]float64, len(cell.West))
	for i, c := range cell.West {
		cell.DxMinusHalf[i] = (cell.Dx + c.Dx) / 2.
		cell.WestFrac[i] = min(c.Dy/cell.Dy, 1.)
		cell.KxxWest[i] = harmonicMean(cell.Kxxyy, c.Kxxyy)
	}
	cell.DyPlusHalf = make([]float64, len(cell.North))
	cell.NorthFrac = make([]float64, len(cell.North))
	cell.KyyNorth = make([]float64, len(cell.North))
	for i, c := range cell.North {
		cell.DyPlusHalf[i] = (cell.Dy + c.Dy) / 2.
		cell.NorthFrac[i] = min(c.Dx/cell.Dx, 1.)
		cell.KyyNorth[i] = harmonicMean(cell.Kxxyy, c.Kxxyy)
	}
	cell.DyMinusHalf = make([]float64, len(cell.South))
	cell.SouthFrac = make([]float64, len(cell.South))
	cell.KyySouth = make([]float64, len(cell.South))
	for i, c := range cell.South {
		cell.DyMinusHalf[i] = (cell.Dy + c.Dy) / 2.
		cell.SouthFrac[i] = min(c.Dx/cell.Dx, 1.)
		cell.KyySouth[i] = harmonicMean(cell.Kxxyy, c.Kxxyy)
	}
	cell.DzPlusHalf = make([]float64, len(cell.Above))
	cell.AboveFrac = make([]float64, len(cell.Above))
	cell.KzzAbove = make([]float64, len(cell.Above))
	for i, c := range cell.Above {
		cell.DzPlusHalf[i] = (cell.Dz + c.Dz) / 2.
		cell.AboveFrac[i] = min((c.Dx*c.Dy)/(cell.Dx*cell.Dy), 1.)
		cell.KzzAbove[i] = harmonicMean(cell.Kzz, c.Kzz)
	}
	cell.DzMinusHalf = make([]float64, len(cell.Below))
	cell.BelowFrac = make([]float64, len(cell.Below))
	cell.KzzBelow = make([]float64, len(cell.Below))
	for i, c := range cell.Below {
		cell.DzMinusHalf[i] = (cell.Dz + c.Dz) / 2.
		cell.BelowFrac[i] = min((c.Dx*c.Dy)/(cell.Dx*cell.Dy), 1.)
		cell.KzzBelow[i] = harmonicMean(cell.Kzz, c.Kzz)
	}
	cell.GroundLevelFrac = make([]float64, len(cell.GroundLevel))
	for i, c := range cell.GroundLevel {
		cell.GroundLevelFrac[i] = min((c.Dx*c.Dy)/(cell.Dx*cell.Dy), 1.)
	}
}

// Add in emissions flux to each cell at every time step, also
// set initial concentrations to final concentrations from previous
// time step, and set old velocities to velocities from previous time
// step.
func (c *Cell) addEmissionsFlux(d *InMAPdata) {
	for i, _ := range polNames {
		c.Cf[i] += c.emisFlux[i] * d.Dt
		c.Ci[i] = c.Cf[i]
	}
}

type configData struct {
	// Path to location of baseline meteorology and pollutant data,
	// where [layer] is a stand-in for the model layer number. The files
	// should be in Gob format (http://golang.org/pkg/encoding/gob/).
	// Can include environment variables.
	InMAPdataTemplate string

	NumLayers     int // Number of vertical layers to use in the model
	NumProcessors int // Number of processors to use for calculations

	// Paths to emissions shapefiles.
	// Can be elevated or ground level; elevated files need to have columns
	// labeled "height", "diam", "temp", and "velocity" containing stack
	// information in units of m, m, K, and m/s, respectively.
	// Emissions will be allocated from the geometries in the shape file
	// to the InMAP computational grid, but the mapping projection of the
	// shapefile must be the same as the projection InMAP uses.
	// Can include environment variables.
	EmissionsShapefiles []string

	// Path to desired output file location, where [layer] is a stand-in
	// for the model layer number. Can include environment variables.
	OutputTemplate string

	// If true, output data for all model layers. If false, only output
	// the lowest layer.
	OutputAllLayers bool

	// Number of iterations to calculate. If < 1, convergence
	// is automatically calculated.
	NumIterations int

	// Port for hosting web page. If HTTPport is `8080`, then the GUI
	// would be viewed by visiting `localhost:8080` in a web browser.
	// If HTTPport is "", then the web server doesn't run.
	HTTPport string
}

func main() {
	flag.Parse()
	if *configFile == "" {
		fmt.Println("Need to specify configuration file as in " +
			"`inmap -config=configFile.json`")
		os.Exit(1)
	}
	config := readConfigFile(*configFile)

	fmt.Println("\n",
		"-------------------------------------------------------\n",
		"                       Welcome!\n",
		"  Global (In)tervention (M)odel for (A)ir (P)ollution  \n",
		"               Version "+version+"                \n",
		"                   Copyright 2013-2014                 \n",
		"         Regents of the University of Minnesota        \n",
		"-------------------------------------------------------\n")

	runtime.GOMAXPROCS(config.NumProcessors)

	fmt.Println("Reading input data...")
	// InMAPdataTemplate is path to netcdf output file
	d := initInMAPdataGlobalBasic(config.InMAPdataTemplate,
		config.NumLayers, config.NumIterations, config.HTTPport)

	// Input units = tons/year; output units = μg/s
	const massConv = 907184740000.       // μg per short ton
	const timeConv = 3600. * 8760.       // seconds per year
	const emisConv = massConv / timeConv // convert tons/year to μg/s

	// need to fill:
	// emissions[inmap.EmisNames[j]][row] += val * weightFactor
	// where value is in ug/s and weightFactor is 1 (for now)

	emissions := make(map[string][]float64)
	for _, pol := range EmisNames {
		if _, ok := emissions[pol]; !ok {
			emissions[pol] = make([]float64, len(d.Data))
		}
	}

	// Next, add emissions, how about to a group of cells in China

	emissions["SOx"][8610] += 50000.000 * emisConv
	emissions["SOx"][8754] += 50000.000 * emisConv
	emissions["SOx"][8611] += 50000.000 * emisConv
	emissions["SOx"][8755] += 50000.000 * emisConv
	emissions["SOx"][8901] += 50000.000 * emisConv
	emissions["SOx"][9045] += 50000.000 * emisConv
	emissions["SOx"][9046] += 50000.000 * emisConv
	emissions["SOx"][9190] += 50000.000 * emisConv
	emissions["SOx"][9047] += 50000.000 * emisConv
	emissions["SOx"][9191] += 50000.000 * emisConv
	emissions["SOx"][47922] += 50000.000 * emisConv
	emissions["SOx"][48066] += 50000.000 * emisConv
	emissions["SOx"][47923] += 50000.000 * emisConv
	emissions["SOx"][48067] += 50000.000 * emisConv
	emissions["SOx"][48213] += 50000.000 * emisConv
	emissions["SOx"][48357] += 50000.000 * emisConv
	emissions["SOx"][48358] += 50000.000 * emisConv
	emissions["SOx"][48502] += 50000.000 * emisConv
	emissions["SOx"][48359] += 50000.000 * emisConv
	emissions["SOx"][48503] += 50000.000 * emisConv
	emissions["SOx"][100338] += 50000.000 * emisConv
	emissions["SOx"][100482] += 50000.000 * emisConv
	emissions["SOx"][100339] += 50000.000 * emisConv
	emissions["SOx"][100483] += 50000.000 * emisConv
	emissions["SOx"][100629] += 50000.000 * emisConv
	emissions["SOx"][100773] += 50000.000 * emisConv
	emissions["SOx"][100774] += 50000.000 * emisConv
	emissions["SOx"][100918] += 50000.000 * emisConv
	emissions["SOx"][100775] += 50000.000 * emisConv
	emissions["SOx"][100919] += 50000.000 * emisConv

	emissions["PM2_5"][8610] += 50000.000 * emisConv
	emissions["PM2_5"][8611] += 50000.000 * emisConv

	/*
		// Add in emissions shapefiles
		// Load emissions into rtree for fast searching
		emisTree := rtreego.NewTree(25, 50)
		for _, fname := range config.EmissionsShapefiles {
			fmt.Println("Loading emissions shapefile:\n", fname)
			fname = strings.Replace(fname, ".shp", "", -1)
			f1, err := os.Open(fname + ".shp")
			if err != nil {
				//panic(err)
				continue
			}
			shp, err := shapefile.OpenShapefile(f1)
			if err != nil {
				panic(err)
			}
			f2, err := os.Open(fname + ".dbf")
			if err != nil {
				panic(err)
			}
			dbf, err := shapefile.OpenDBFFile(f2)
			if err != nil {
				panic(err)
			}
			for i := 0; i < int(dbf.DBFFileHeader.NumRecords); i++ {
				sRec, err := shp.NextRecord()
				if err != nil {
					panic(err)
				}
				fields, err := dbf.NextRecord()
				if err != nil {
					panic(err)
				}
				e := new(emisRecord)
				e.emis = make([]float64, len(inmap.EmisNames))
				e.g = sRec.Geometry
				for ii, pol := range inmap.EmisNames {
					// Input units = tons/year; output units = μg/s
					const massConv = 907184740000.       // μg per short ton
					const timeConv = 3600. * 8760.       // seconds per year
					const emisConv = massConv / timeConv // convert tons/year to μg/s
					if iii, ok := dbf.FieldIndicies[pol]; ok {
						switch fields[iii].(type) {
						case float64:
							e.emis[ii] += fields[iii].(float64) * emisConv
							if math.IsNaN(e.emis[ii]) {
								e.emis[ii] = 0.
							}
						case int:
							e.emis[ii] += float64(fields[iii].(int)) * emisConv
						}
					}
				}
				if iii, ok := dbf.FieldIndicies["height"]; ok {
					e.height = fields[iii].(float64) // stack height [m]
					if math.IsNaN(e.height) {
						e.height = 0.
					}
				}
				if iii, ok := dbf.FieldIndicies["diam"]; ok {
					e.diam = fields[iii].(float64) // stack diameter [m]
					if math.IsNaN(e.diam) {
						e.diam = 0.
					}
				}
				if iii, ok := dbf.FieldIndicies["temp"]; ok {
					e.temp = fields[iii].(float64) // stack temperature [K]
					if math.IsNaN(e.temp) {
						e.temp = 0.
					}
				}
				if iii, ok := dbf.FieldIndicies["velocity"]; ok {
					e.velocity = fields[iii].(float64) // stack velocity [m/s]
					if math.IsNaN(e.velocity) {
						e.velocity = 0.
					}
				}
				e.bounds, err = geomconv.GeomToRect(e.g)
				if err != nil {
					panic(err)
				}
				emisTree.Insert(e)
			}
			f1.Close()
			f2.Close()
		}

		fmt.Println("Allocating emissions to grid cells...")
		// allocate emissions to appropriate grid cells
		for i := d.LayerStart[0]; i < d.LayerEnd[0]; i++ {
			cell := d.Data[i]
			bounds, err := geomconv.GeomToRect(cell.Geom)
			if err != nil {
				panic(err)
			}
			for _, eTemp := range emisTree.SearchIntersect(bounds) {
				e := eTemp.(*emisRecord)
				var intersection geom.T
				switch e.g.(type) {
				case geom.Point:
					in, err := geomop.Within(e.g, cell.Geom)
					if err != nil {
						panic(err)
					}
					if in {
						intersection = e.g
					} else {
						continue
					}
				default:
					intersection, err = geomop.Construct(e.g, cell.Geom,
						geomop.INTERSECTION)
					if err != nil {
						panic(err)
					}
				}
				if intersection == nil {
					continue
				}
				var weightFactor float64 // fraction of geometry in grid cell
				switch e.g.(type) {
				case geom.Polygon, geom.MultiPolygon:
					weightFactor = geomop.Area(intersection) / geomop.Area(e.g)
				case geom.LineString, geom.MultiLineString:
					weightFactor = geomop.Length(intersection) / geomop.Length(e.g)
				case geom.Point:
					weightFactor = 1.
				default:
					panic(geomop.UnsupportedGeometryError{intersection})
				}
				var plumeRow int
				if e.height > 0. { // calculate plume rise
					plumeRow, err = d.CalcPlumeRise(
						e.height, e.diam, e.temp, e.velocity, i)
					if err != nil {
						panic(err)
					}
				} else {
					plumeRow = i
				}
				for j, val := range e.emis {
					emissions[inmap.EmisNames[j]][plumeRow] += val * weightFactor
				}
			}
		}
	*/

	for pol, arr := range emissions {
		sum := 0.
		for _, val := range arr {
			sum += val
		}
		fmt.Printf("%v, %g ug/s\n", pol, sum)
	}

	// Run model
	finalConc := d.Run(emissions) //, config.OutputAllLayers)

	shapeNC := []int{15, 91, 144}
	OutNCA := make([]*sparse.DenseArray, 4)

	ncNames := []string{"SOx", "pSO4", "TotalPM2_5", "PrimaryPM2_5"}
	for i, name := range ncNames {
		OutNCA[i] = sparse.ZerosDense(shapeNC...)
		for z := 0; z < shapeNC[0]; z++ {
			indexNC := 0
			for y := 0; y < shapeNC[1]; y++ {
				for x := 0; x < shapeNC[2]; x++ {
					OutNCA[i].Set(finalConc[name][z][indexNC], z, y, x)
					indexNC++
				}
			}
		}
	}

	// write out data to file
	outputFile := "basicOutputV1.nc"
	fmt.Printf("Writing out data to %v...\n", outputFile)
	h := cdf.NewHeader(
		[]string{"x", "y", "z"},
		[]int{shapeNC[2], shapeNC[1], shapeNC[0]})
	h.AddAttribute("", "comment", "First Output from globalInmapBasic")

	dataOut := map[string]dataHolder{
		"SOx": dataHolder{[]string{"z", "y", "x"},
			"Average ∆concentration of SO2",
			"ug m-3", OutNCA[0]},
		"pSO4": dataHolder{[]string{"z", "y", "x"},
			"Average ∆concentration of sulfate",
			"ug m-3", OutNCA[1]},
		"TotalPM2_5": dataHolder{[]string{"z", "y", "x"},
			"Average ∆concentration of total PM2.5",
			"ug m-3", OutNCA[2]},
		"PrimaryPM2_5": dataHolder{[]string{"z", "y", "x"},
			"Average ∆concentration of primary PM2.5",
			"ug m-3", OutNCA[3]}}

	for name, d := range dataOut {
		h.AddVariable(name, d.dims, []float32{0})
		h.AddAttribute(name, "description", d.Description)
		h.AddAttribute(name, "units", d.Units)
	}
	h.Define()
	ff, err := os.Create(outputFile)
	if err != nil {
		panic(err)
	}
	f, err := cdf.Create(ff, h) // writes the header to ff
	if err != nil {
		panic(err)
	}
	for name, d := range dataOut {
		writeNCF(f, name, d.data)
	}
	err = cdf.UpdateNumRecs(ff)
	if err != nil {
		panic(err)
	}
	ff.Close()

	/*
		writeOutput(finalConc, d, config.OutputTemplate, config.OutputAllLayers)

		fmt.Println("\nIntake fraction results:")
		breathingRate := 15. // [m³/day]
		iF := d.IntakeFraction(breathingRate)
		// Write iF to stdout
		w := tabwriter.NewWriter(os.Stdout, 0, 8, 1, '\t', 0)
		var popList []string
		for _, m := range iF {
			for p := range m {
				popList = append(popList, p)
			}
			break
		}
		sort.Strings(popList)
		fmt.Fprintln(w, strings.Join(append([]string{"pol"}, popList...), "\t"))
		for pol, m := range iF {
			temp := make([]string, len(popList))
			for i, pop := range popList {
				temp[i] = fmt.Sprintf("%.3g", m[pop])
			}
			fmt.Fprintln(w, strings.Join(append([]string{pol}, temp...), "\t"))
		}
		w.Flush()
	*/

	fmt.Println("\n",
		"------------------------------------\n",
		"      Global InMAP Completed!\n",
		"------------------------------------\n")
} // end of main

// Dev's functions/types for basic InMAP run.
func initInMAPdataGlobalBasic(filetemplate string, nLayers int, numIterations int, httpPort string) *InMAPdata {
	inputData := make([][]*Cell, nLayers)
	d := new(InMAPdata)
	d.NumIterations = numIterations
	d.Nlayers = nLayers
	d.LayerStart = make([]int, nLayers)
	d.LayerEnd = make([]int, nLayers)
	// actually, going to make basic variable grids first.
	var wg sync.WaitGroup
	wg.Add(nLayers)
	for k := 0; k < nLayers; k++ {
		go func(k int) {
			filename := strings.Replace(filetemplate, "[layer]",
				fmt.Sprintf("%v", k), -1)
			f, err := os.Open(filename)
			if err != nil {
				fmt.Println(err.Error())
				os.Exit(1)
			}
			g := gob.NewDecoder(f)
			g.Decode(&inputData[k])
			d.LayerStart[k] = 0
			d.LayerEnd[k] = len(inputData[k])
			f.Close()
			wg.Done()
		}(k)
	}
	wg.Wait()
	ncells := 0
	// Adjust so beginning of layer is at end of previous layer
	for k := 0; k < nLayers; k++ {
		d.LayerStart[k] += ncells
		d.LayerEnd[k] += ncells
		ncells += len(inputData[k])
	}
	// set up data holders
	d.Data = make([]*Cell, ncells)
	for _, indata := range inputData {
		for _, c := range indata {
			c.prepare()
			d.Data[c.Row] = c
		}
	}
	d.westBoundary = make([]*Cell, 0, 200)
	d.eastBoundary = make([]*Cell, 0, 200)
	d.southBoundary = make([]*Cell, 0, 200)
	d.northBoundary = make([]*Cell, 0, 200)
	d.topBoundary = make([]*Cell, 0, 200)
	nprocs := runtime.GOMAXPROCS(0)
	wg.Add(nprocs)
	for procNum := 0; procNum < nprocs; procNum++ {
		go func(procNum int) {
			for ii := procNum; ii < len(d.Data); ii += nprocs {
				cell := d.Data[ii]
				// Link cells to neighbors and/or boundaries.
				if len(cell.IWest) == 0 {
					c := cell.makecopy()
					cell.West = []*Cell{c}
					d.westBoundary = append(d.westBoundary, c)
				} else {
					cell.West = make([]*Cell, len(cell.IWest))
					for i, row := range cell.IWest {
						cell.West[i] = d.Data[row]
					}
					cell.IWest = nil
				}
				if len(cell.IEast) == 0 {
					c := cell.makecopy()
					cell.East = []*Cell{c}
					d.eastBoundary = append(d.eastBoundary, c)
				} else {
					cell.East = make([]*Cell, len(cell.IEast))
					for i, row := range cell.IEast {
						cell.East[i] = d.Data[row]
					}
					cell.IEast = nil
				}
				if len(cell.ISouth) == 0 {
					c := cell.makecopy()
					cell.South = []*Cell{c}
					d.southBoundary = append(d.southBoundary, c)
				} else {
					cell.South = make([]*Cell, len(cell.ISouth))
					for i, row := range cell.ISouth {
						cell.South[i] = d.Data[row]
					}
					cell.ISouth = nil
				}
				if len(cell.INorth) == 0 {
					c := cell.makecopy()
					cell.North = []*Cell{c}
					d.northBoundary = append(d.northBoundary, c)
				} else {
					cell.North = make([]*Cell, len(cell.INorth))
					for i, row := range cell.INorth {
						cell.North[i] = d.Data[row]
					}
					cell.INorth = nil
				}
				if len(cell.IAbove) == 0 {
					c := cell.makecopy()
					cell.Above = []*Cell{c}
					d.topBoundary = append(d.topBoundary, c)
				} else {
					cell.Above = make([]*Cell, len(cell.IAbove))
					for i, row := range cell.IAbove {
						cell.Above[i] = d.Data[row]
					}
					cell.IAbove = nil
				}
				if cell.Layer != 0 {
					cell.Below = make([]*Cell, len(cell.IBelow))
					cell.GroundLevel = make([]*Cell, len(cell.IGroundLevel))
					for i, row := range cell.IBelow {
						cell.Below[i] = d.Data[row]
					}
					for i, row := range cell.IGroundLevel {
						cell.GroundLevel[i] = d.Data[row]
					}
					cell.IBelow = nil
					cell.IGroundLevel = nil
				} else { // assume bottom boundary is the same as lowest layer.
					cell.Below = []*Cell{d.Data[cell.Row]}
					cell.GroundLevel = []*Cell{d.Data[cell.Row]}
				}
				cell.neighborInfo()
			}
			wg.Done()
		}(procNum)
	}
	wg.Wait()

	d.setTstepCFL() // Set time step
	//d.setTstepRuleOfThumb() // Set time step
	// go d.WebServer(httpPort)

	return d
}

// other needed functions:
func checkConvergence(newSum, oldSum float64, Var string) bool {
	bias := (newSum - oldSum) / oldSum
	fmt.Printf("%v: total mass difference = %3.2g%% from last check.\n",
		Var, bias*100)
	if math.Abs(bias) > tolerance || math.IsInf(bias, 0) {
		return false
	} else {
		return true
	}
}

func harmonicMean(a, b float64) float64 {
	return 2. * a * b / (a + b)
}

func writeNCF(f *cdf.File, Var string, data *sparse.DenseArray) {
	data32 := make([]float32, len(data.Elements))
	for i, e := range data.Elements {
		data32[i] = float32(e)
	}
	end := f.Header.Lengths(Var)
	start := make([]int, len(end))
	w := f.Writer(Var, start, end)
	_, err := w.Write(data32)
	if err != nil {
		panic(err)
	}
}

// science.go functions pasted here:

// boundary layer and based on Wilson (2004) for above the boundary layer.
// Also calculate horizontal mixing.
func (c *Cell) Mixing(Δt float64) {
	for ii, _ := range c.Cf {
		// Pleim (2007) Equation 10.
		if c.Layer < f2i(c.PblTopLayer) { // Convective mixing
			for i, g := range c.GroundLevel { // Upward convection
				c.Cf[ii] += g.M2u * g.Ci[ii] * Δt * c.GroundLevelFrac[i]
			}
			for i, a := range c.Above { // Balancing downward mixing
				c.Cf[ii] += (a.M2d*a.Ci[ii]*a.Dz/c.Dz - c.M2d*c.Ci[ii]) *
					Δt * c.AboveFrac[i]
			}
		}
		for i, a := range c.Above { // Mixing with above
			c.Cf[ii] += 1. / c.Dz * (c.KzzAbove[i] * (a.Ci[ii] - c.Ci[ii]) /
				c.DzPlusHalf[i]) * Δt * c.AboveFrac[i]
		}
		for i, b := range c.Below { // Mixing with below
			c.Cf[ii] += 1. / c.Dz * (c.KzzBelow[i] * (b.Ci[ii] - c.Ci[ii]) /
				c.DzMinusHalf[i]) * Δt * c.BelowFrac[i]
		}
		// Horizontal mixing
		for i, w := range c.West { // Mixing with West
			c.Cf[ii] += 1. / c.Dx * (c.KxxWest[i] *
				(w.Ci[ii] - c.Ci[ii]) / c.DxMinusHalf[i]) * Δt * c.WestFrac[i]
		}
		for i, e := range c.East { // Mixing with East
			c.Cf[ii] += 1. / c.Dx * (c.KxxEast[i] *
				(e.Ci[ii] - c.Ci[ii]) / c.DxPlusHalf[i]) * Δt * c.EastFrac[i]
		}
		for i, s := range c.South { // Mixing with South
			c.Cf[ii] += 1. / c.Dy * (c.KyySouth[i] *
				(s.Ci[ii] - c.Ci[ii]) / c.DyMinusHalf[i]) * Δt * c.SouthFrac[i]
		}
		for i, n := range c.North { // Mixing with North
			c.Cf[ii] += 1. / c.Dy * (c.KyyNorth[i] *
				(n.Ci[ii] - c.Ci[ii]) / c.DyPlusHalf[i]) * Δt * c.NorthFrac[i]
		}
	}
}

// Calculates advective flux in West and East directions
// using upwind flux-form spatial approximation for δ(uq)/δx.
// Returns mass flux per unit area per unit time.
func (c *Cell) westEastFlux(ii int) float64 {
	var flux float64
	for i, w := range c.West {
		flux += (w.UPlusSpeed*w.Ci[ii] -
			c.Ci[ii]*c.UMinusSpeed) * c.WestFrac[i]
	}
	for i, e := range c.East {
		flux += (e.UMinusSpeed*e.Ci[ii] -
			c.Ci[ii]*c.UPlusSpeed) * c.EastFrac[i]
	}
	return flux
}

// Calculates advective flux in South and North directions
// using upwind flux-form spatial approximation for δ(uq)/δx.
// Returns mass flux per unit area per unit time.
func (c *Cell) southNorthFlux(ii int) float64 {
	var flux float64
	for i, s := range c.South {
		flux += (s.VPlusSpeed*s.Ci[ii] -
			c.Ci[ii]*c.VMinusSpeed) * c.SouthFrac[i]
	}
	for i, n := range c.North {
		flux += (n.VMinusSpeed*n.Ci[ii] -
			c.Ci[ii]*c.VPlusSpeed) * c.NorthFrac[i]
	}
	return flux
}

// Calculates advective flux in Below and Above directions
// using upwind flux-form spatial approximation for δ(uq)/δx.
// Returns mass flux per unit area per unit time.
func (c *Cell) belowAboveFlux(ii int) float64 {
	var flux float64
	if c.Layer != 0 { // Can't advect downwards from bottom cell
		for i, b := range c.Below {
			flux += (b.WPlusSpeed*b.Ci[ii] -
				c.Ci[ii]*c.WMinusSpeed) * c.BelowFrac[i]
		}
	}
	for i, a := range c.Above {
		flux += (a.WMinusSpeed*a.Ci[ii] -
			c.Ci[ii]*c.WPlusSpeed) * c.AboveFrac[i]
	}
	return flux
}

// Calculates advection in the cell based
// on the upwind differences scheme.
func (c *Cell) UpwindAdvection(Δt float64) {
	for ii, _ := range c.Cf {
		c.Cf[ii] += c.westEastFlux(ii) / c.Dx * Δt
		c.Cf[ii] += c.southNorthFlux(ii) / c.Dy * Δt
		c.Cf[ii] += c.belowAboveFlux(ii) / c.Dz * Δt
	}
}

// Calculates the secondary formation of PM2.5.
// Explicitely calculates formation of particulate sulfate
// from gaseous and aqueous SO2.
// Partitions organic matter ("gOrg" and "pOrg"), the
// nitrogen in nitrate ("gNO and pNO"), and the nitrogen in ammonia ("gNH" and
// "pNH) between gaseous and particulate phase
// based on the spatially explicit partioning present in the baseline data.
func (c *Cell) Chemistry(d *InMAPdata) {

	// All SO4 forms particles, so sulfur particle formation is limited by the
	// SO2 -> SO4 reaction.
	ΔS := c.SO2oxidation * c.Cf[igS] * d.Dt
	c.Cf[igS] -= ΔS
	c.Cf[ipS] += ΔS

	// VOC/SOA partitioning
	totalOrg := c.Cf[igOrg] + c.Cf[ipOrg]
	c.Cf[igOrg] = totalOrg * c.OrgPartitioning
	c.Cf[ipOrg] = totalOrg * (1 - c.OrgPartitioning)

	// NH3 / NH4 partitioning
	totalNH := c.Cf[igNH] + c.Cf[ipNH]
	c.Cf[igNH] = totalNH * c.NHPartitioning
	c.Cf[ipNH] = totalNH * (1 - c.NHPartitioning)

	// NOx / pN0 partitioning
	totalNO := c.Cf[igNO] + c.Cf[ipNO]
	c.Cf[igNO] = totalNO * c.NOPartitioning
	c.Cf[ipNO] = totalNO * (1 - c.NOPartitioning)

}

// Calculates particle removal by dry deposition
func (c *Cell) DryDeposition(d *InMAPdata) {
	if c.Layer == 0 {
		fac := 1. / c.Dz * d.Dt
		noxfac := math.Max(1-c.NOxDryDep*fac, 0.0)
		so2fac := math.Max(1-c.SO2DryDep*fac, 0.0)
		vocfac := math.Max(1-c.VOCDryDep*fac, 0.0)
		nh3fac := math.Max(1-c.NH3DryDep*fac, 0.0)
		pm25fac := math.Max(1-c.ParticleDryDep*fac, 0.0)
		c.Cf[igOrg] *= vocfac
		c.Cf[ipOrg] *= pm25fac
		c.Cf[iPM2_5] *= pm25fac
		c.Cf[igNH] *= nh3fac
		c.Cf[ipNH] *= pm25fac
		c.Cf[igS] *= so2fac
		c.Cf[ipS] *= pm25fac
		c.Cf[igNO] *= noxfac
		c.Cf[ipNO] *= pm25fac
	}
}

// Calculates particle removal by wet deposition
func (c *Cell) WetDeposition(Δt float64) {
	particleFrac := math.Max(1.-c.ParticleWetDep*Δt, 0.0)
	SO2Frac := math.Max(1.-c.SO2WetDep*Δt, 0.0)
	otherGasFrac := math.Max(1-c.OtherGasWetDep*Δt, 0.0)
	c.Cf[igOrg] *= otherGasFrac
	c.Cf[ipOrg] *= particleFrac
	c.Cf[iPM2_5] *= particleFrac
	c.Cf[igNH] *= otherGasFrac
	c.Cf[ipNH] *= particleFrac
	c.Cf[igS] *= SO2Frac
	c.Cf[ipS] *= particleFrac
	c.Cf[igNO] *= otherGasFrac
	c.Cf[ipNO] *= particleFrac
}

// convert float to int (rounding)
func f2i(f float64) int {
	return int(f + 0.5)
}

func max(vals ...float64) float64 {
	m := 0.
	for _, v := range vals {
		if v > m {
			m = v
		}
	}
	return m
}
func min(v1, v2 float64) float64 {
	if v1 < v2 {
		return v1
	} else {
		return v2
	}
}

// ----------------------------------------------------------------------------
// original inmap functions/types below here

/*
type emisRecord struct {
	g        geom.T
	emis     []float64
	bounds   *rtreego.Rect
	height   float64 // stack height [m]
	diam     float64 // stack diameter [m]
	temp     float64 // stack temperature [K]
	velocity float64 // stack velocity [m/s]
}

func (e emisRecord) Bounds() *rtreego.Rect {
	return e.bounds
}

// write data out to shapefile
func writeOutput(results map[string][][]float64, d *inmap.InMAPdata,
	outFileTemplate string, writeAllLayers bool) {

	// Projection definition. This may need to be changed for a different
	// spatial domain.
	const proj4 = `PROJCS["Lambert_Conformal_Conic",GEOGCS["GCS_unnamed ellipse",DATUM["D_unknown",SPHEROID["Unknown",6370997,0]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Conformal_Conic"],PARAMETER["standard_parallel_1",33],PARAMETER["standard_parallel_2",45],PARAMETER["latitude_of_origin",40],PARAMETER["central_meridian",-97],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]`

	vars := make([]string, 0, len(results))
	for v := range results {
		vars = append(vars, v)
	}
	sort.Strings(vars)
	fields := make([]shp.Field, len(vars))
	for i, v := range vars {
		fields[i] = shp.FloatField(v, 14, 8)
	}

	var nlayers int
	if writeAllLayers {
		nlayers = d.Nlayers
	} else {
		nlayers = 1
	}
	row := 0
	for k := 0; k < nlayers; k++ {

		filename := strings.Replace(outFileTemplate, "[layer]",
			fmt.Sprintf("%v", k), -1)
		// remove extension and replace it with .shp
		extIndex := strings.LastIndex(filename, ".")
		if extIndex == -1 {
			extIndex = len(filename)
		}
		filename = filename[0:extIndex] + ".shp"
		shape, err := shp.Create(filename, shp.POLYGON)
		if err != nil {
			log.Fatal(err)
		}
		shape.SetFields(fields)

		numRowsInLayer := len(results[vars[0]][k])
		for i := 0; i < numRowsInLayer; i++ {
			s, err := geomconv.Geom2Shp(d.Data[row].Geom)
			if err != nil {
				panic(err)
			}
			shape.Write(s)
			for j, v := range vars {
				shape.WriteAttribute(i, j, results[v][k][i])
			}
			row++
		}
		shape.Close()

		// Create .prj file
		f, err := os.Create(filename[0:extIndex] + ".prj")
		if err != nil {
			log.Fatal(err)
		}
		fmt.Fprint(f, proj4)
		f.Close()
	}
}
*/

func s2i(s string) int {
	i, err := strconv.ParseInt(s, 0, 64)
	if err != nil {
		panic(err)
	}
	return int(i)
}
func s2f(s string) float64 {
	f, err := strconv.ParseFloat(s, 64)
	if err != nil {
		panic(err)
	}
	return f
}

// readConfigFile reads and parses a json configuration file.
// See below for the required variables.
func readConfigFile(filename string) (config *configData) {
	// Open the configuration file
	var (
		file  *os.File
		bytes []byte
		err   error
	)
	file, err = os.Open(filename)
	if err != nil {
		fmt.Printf("The configuration file you have specified, %v, does not "+
			"appear to exist. Please check the file name and location and "+
			"try again.\n", filename)
		os.Exit(1)
	}
	reader := bufio.NewReader(file)
	bytes, err = ioutil.ReadAll(reader)
	if err != nil {
		panic(err)
	}

	config = new(configData)
	err = json.Unmarshal(bytes, config)
	if err != nil {
		fmt.Printf(
			"There has been an error parsing the configuration file.\n"+
				"Please ensure that the file is in valid JSON format\n"+
				"(you can check for errors at http://jsonlint.com/)\n"+
				"and try again!\n\n%v\n\n", err.Error())
		os.Exit(1)
	}

	config.InMAPdataTemplate = os.ExpandEnv(config.InMAPdataTemplate)
	config.OutputTemplate = os.ExpandEnv(config.OutputTemplate)

	for i := 0; i < len(config.EmissionsShapefiles); i++ {
		config.EmissionsShapefiles[i] =
			os.ExpandEnv(config.EmissionsShapefiles[i])
	}

	if config.OutputTemplate == "" {
		fmt.Println("You need to specify an output template in the " +
			"configuration file(for example: " +
			"\"OutputTemplate\":\"output_[layer].geojson\"")
		os.Exit(1)
	}

	outdir := filepath.Dir(config.OutputTemplate)
	err = os.MkdirAll(outdir, os.ModePerm)
	if err != nil {
		fmt.Println(err.Error())
		os.Exit(1)
	}
	return
}
