package main

import (
	"bufio"
	"encoding/json"
	"flag"
	"fmt"
	"io/ioutil"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"runtime/pprof"
	"strings"
	"sync"
	"time"

	"bitbucket.org/ctessum/atmos/acm2"
	"bitbucket.org/ctessum/atmos/emep"
	"bitbucket.org/ctessum/atmos/gocart"
	"bitbucket.org/ctessum/atmos/seinfeld"
	"bitbucket.org/ctessum/atmos/wesely1989"
	"code.google.com/p/lvd.go/cdf"
	//	"bitbucket.org/ctessum/sparse"
	"devM/sparse"
)

type ConfigInfo struct {
	//	Wrfout              string  // Location of WRF output files. [DATE] is a wild card for the simulation date.
	GeosMet             string  // Location of GEOS met files. (GEOSFP.[DATE].A3dyn.2x25.nc)
	GeosMetI            string  // Location of GEOS met files (GEOSFP.[DATE].I3.2x25.nc).
	GeosMetA            string  // Location of GEOS met files (GEOSFP.[DATE].A1.2x25.nc).
	GeosMetAcld         string  // Location of GEOS met files (GEOSFP.[DATE].A3cld.2x25.nc).
	GeosMetAmstE        string  // Location of GEOS met files (GEOSFP.[DATE].A3mstE.2x25.nc).
	GeosChem            string  // Location of GEOS chem files. (gc_output.[DATE].[HOUR]0000.nc).
	GeosApBp            string  // Location of GEOS file containing Ap and Bp for layer height calculations
	GeosFixed           string  // Location of GEOS file containing fixed information (land vs. water)
	MetTimeStep         int     // Timestep for met data, for GEOS it is 3-hrs, WRF is 1-hr
	OutputDir           string  // Directory to put the output files in
	OutputFilePrefix    string  // name for output files
	StartDate           string  // Format = "YYYYMMDD"
	EndDate             string  // Format = "YYYYMMDD"
	Nprocs              int     // number of processors to use
	VariableGrid_x_o    float64 // lower left of output grid, x
	VariableGrid_y_o    float64 // lower left of output grid, y
	VariableGrid_dx     float64 // m
	VariableGrid_dy     float64 // m
	Xnests              []int   // Nesting multiples in the X direction
	Ynests              []int   // Nesting multiples in the Y direction
	HiResLayers         int     // number of layers to do in high resolution (layers above this will be lowest resolution.
	CtmGrid_x_o         float64 // lower left of Chemical Transport Model (CTM) grid, x
	CtmGrid_y_o         float64 // lower left of grid, y
	CtmGrid_dx          float64 // m
	CtmGrid_dy          float64 // m
	CtmGrid_nx          int
	CtmGrid_ny          int
	GridProj            string   // projection info for CTM grid; Proj4 format
	PopDensityCutoff    float64  // limit for people per unit area in the grid cell
	PopCutoff           float64  // limit for total number of people in the grid cell
	BboxOffset          float64  // A number significantly less than the smallest grid size but not small enough to be confused with zero.
	CensusFile          string   // Path to census shapefile
	CensusPopColumns    []string // Shapefile fields containing populations for multiple demographics
	PopGridColumn       string   // Name of field in shapefile to be used for determining variable grid resolution
	MortalityRateFile   string   // Path to the mortality rate shapefile
	MortalityRateColumn string   // Name of field in mortality rate shapefiel containing the mortality rate.
}

const (
	//	wrfFormat    = "2006-01-02_15_04_05"
	geosFormat   = "20060102"
	inDateFormat = "20060102" // reference date for go time.Parse
	tolerance    = 1.e-10     // tolerance for comparing floats

	// physical constants
	MWa      = 28.97   // g/mol, molar mass of air
	mwN      = 46.0055 // g/mol, molar mass of nitrogen
	mwS      = 32.0655 // g/mol, molar mass of sulfur
	mwNH4    = 18.03851
	mwSO4    = 96.0632
	mwNO3    = 62.00501
	g        = 9.80665 // m/s2
	rr       = 287.058 // (J /kg K), specific gas constant for dry air
	rC       = 8.31446 // (m3 Pa / k mol)
	κ        = 0.41    // Von Kármán constant
	atmPerPa = 9.86923267e-6
	avNum    = 6.02214e23
)

var (
	start       time.Time
	end         time.Time
	current     time.Time
	numTsteps   float64 // this is total hours / met time step (8 time steps per day for GEOS)
	numTsteps1h float64 // this is total hours
	configFile  *string = flag.String("config", "none", "Path to configuration file")
	config              = new(ConfigInfo)
)

var (
	/*
		// RACM VOC species and molecular weights (g/mol);
		// Only includes anthropogenic precursors to SOA from
		// anthropogenic (aSOA) and biogenic (bSOA) sources as
		// an Ahmadov et al. (2012)
		// assume condensable vapor from SOA has molar mass of 70
		aVOC = map[string]float64{"hc5": 72, "hc8": 114,
					"olt": 42, "oli": 68, "tol": 92, "xyl": 106, "csl": 108,
					"cvasoa1": 70, "cvasoa2": 70, "cvasoa3": 70, "cvasoa4": 70}
				bVOC = map[string]float64{"iso": 68, "api": 136, "sesq": 84.2,
					"lim": 136, "cvbsoa1": 70, "cvbsoa2": 70,
					"cvbsoa3": 70, "cvbsoa4": 70}
			// VBS SOA species (anthropogenic only)
			aSOA = map[string]float64{"asoa1i": 1, "asoa1j": 1, "asoa2i": 1,
				"asoa2j": 1, "asoa3i": 1, "asoa3j": 1, "asoa4i": 1, "asoa4j": 1}
			// VBS SOA species (biogenic only)
			bSOA = map[string]float64{"bsoa1i": 1, "bsoa1j": 1, "bsoa2i": 1,
				"bsoa2j": 1, "bsoa3i": 1, "bsoa3j": 1, "bsoa4i": 1, "bsoa4j": 1}
	*/

	// GEOS aVOC and bVOC -- currently incomplete list simply for coding development purposes
	aVOC = map[string]float64{"IJ-AVG-S__BENZ": 12, "IJ-AVG-S__TOLU": 12, "IJ-AVG-S__XYLE": 12}
	bVOC = map[string]float64{"IJ-AVG-S__ISOP": 12, "IJ-AVG-S__LIMO": 136.23}

	aSOA = map[string]float64{"IJ-AVG-S__ASOA1": 150, "IJ-AVG-S__ASOA2": 150, "IJ-AVG-S__ASOA3": 150}
	bSOA = map[string]float64{"IJ-AVG-S__ISOA1": 150, "IJ-AVG-S__ISOA2": 150, "IJ-AVG-S__ISOA3": 150}

	// NO and NO2 species and molecular weights, multiplied by their
	// nitrogen fractions
	// NOx = map[string]float64{"IJ-AVG-S__NO": 30 * 30 / mwN, "IJ-AVG-S__NO2": 46 * 46 / mwN}
	NOx = map[string]float64{"IJ-AVG-S__NO": 46, "IJ-AVG-S__NO2": 46}
	// mass of N in NO
	//NO = map[string]float64{"no": 1.}
	NO = map[string]float64{"IJ-AVG-S__NO": 46.}
	// mass of N in  NO2
	NO2 = map[string]float64{"IJ-AVG-S__NO2": 46.}
	//NO2 = map[string]float64{"no2": 1.}
	// MADE particulate NO species, nitrogen fraction
	// pNO = map[string]float64{"IJ-AVG-S__NIT": mwNO3 / mwN, "no3aj": mwNO3 / mwN}
	pNO = map[string]float64{"IJ-AVG-S__NIT": 62, "IJ-AVG-S__NITs": 62}
	// RACM SOx species and molecular weights
	//SOx = map[string]float64{"so2": 64 * 64 / mwS, "sulf": 98 * 98 / mwS}
	SOx = map[string]float64{"IJ-AVG-S__SO2": 64, "IJ-AVG-S__DMS": 62}
	// MADE particulate Sulfur species; sulfur fraction
	//pS  = map[string]float64{"so4ai": mwSO4 / mwS, "so4aj": mwSO4 / mwS}
	pS = map[string]float64{"IJ-AVG-S__SO4": 96, "IJ-AVG-S__SO4s": 96}
	//NH3 = map[string]float64{"nh3": 17.03056 * 17.03056 / mwN}
	NH3 = map[string]float64{"IJ-AVG-S__NH3": 17.0}
	// MADE particulate ammonia species, nitrogen fraction
	//pNH = map[string]float64{"nh4ai": mwNH4 / mwN, "nh4aj": mwNH4 / mwN}
	pNH = map[string]float64{"IJ-AVG-S__NH4": 18}

	// totalPM25 = map[string]float64{"PM2_5_DRY": 1.}
)

var (
	// GEOS ap and bp values for calculating layer heights
	// Must update manually if changing to a new model with different
	// vertical grid setup
	// switching everything to chem levels!
	ap     = sparse.ZerosDense(48)
	bp     = sparse.ZerosDense(48)
	FRLAND *sparse.DenseArray
	cLev   = map[int]int{0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11, 12: 12, 13: 13, 14: 14, 15: 15,
		16: 16, 17: 17, 18: 18, 19: 19, 20: 20, 21: 21, 22: 22, 23: 23, 24: 24, 25: 25, 26: 26, 27: 27, 28: 28, 29: 29, 30: 30,
		31: 31, 32: 32, 33: 33, 34: 34, 35: 35,
		36: 36, 37: 36,
		38: 37, 39: 37,
		40: 38, 41: 38,
		42: 39, 43: 39,
		44: 40, 45: 40, 46: 40, 47: 40,
		48: 41, 49: 41, 50: 41, 51: 41,
		52: 42, 53: 42, 54: 42, 55: 42,
		56: 43, 57: 43, 58: 43, 59: 43,
		60: 44, 61: 44, 62: 44, 63: 44,
		64: 45, 65: 45, 66: 45, 67: 45,
		68: 46, 69: 46, 70: 46, 71: 46}
	chemVert = 47 // GEOS-Chem has 47 vertical layers
	metVert  = 72 // GEOS-Met variables are on 72 vertical layers
)

func init() {
	var err error

	flag.Parse()
	if *configFile == "" {
		fmt.Println("Need to specify configuration file as in " +
			"`aim -config=configFile.json`")
		os.Exit(1)
		// simple run: "go run geos2inmapBeta.1.go -config=configGEOSbeta.1.json"
	}
	ReadConfigFile(*configFile)

	start, err = time.Parse(inDateFormat, config.StartDate)
	if err != nil {
		panic(err)
	}
	end, err = time.Parse(inDateFormat, config.EndDate)
	if err != nil {
		panic(err)
	}
	end = end.AddDate(0, 0, 1) // add 1 day to the end
	numTsteps = end.Sub(start).Hours() / float64(config.MetTimeStep)
	numTsteps1h = end.Sub(start).Hours()

	runtime.GOMAXPROCS(config.Nprocs)

	getApBp()
	getFixed()
}

var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to file")

func main() {

	flag.Parse()
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			panic(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	fmt.Println("numTsteps as int:", int(numTsteps))
	fmt.Println("numTsteps as float64:", numTsteps)

	fileNames := []string{config.GeosMet, config.GeosMetI, config.GeosMetA, config.GeosMetAcld, config.GeosMetAmstE, config.GeosChem}

	// make a list of all needed variables by file type
	varListMet1 := []string{}
	varListMet2 := []string{}
	varListMet3 := []string{}
	varListMet4 := []string{}
	varListMet5 := []string{}
	varListChem := []string{}

	varListMet1 = append(varListMet1, "U")
	varListMet1 = append(varListMet1, "V")
	varListMet1 = append(varListMet1, "OMEGA")
	varListMet2 = append(varListMet2, "T")
	varListMet2 = append(varListMet2, "PS")
	varListMet3 = append(varListMet3, "USTAR")
	varListMet3 = append(varListMet3, "HFLUX")
	varListMet3 = append(varListMet3, "Z0M")
	varListMet3 = append(varListMet3, "PBLH")
	varListMet3 = append(varListMet3, "PRECTOT")
	varListMet3 = append(varListMet3, "PARDF")
	varListMet3 = append(varListMet3, "PARDR")
	varListMet4 = append(varListMet4, "QL")
	varListMet4 = append(varListMet4, "CLOUD")
	varListMet5 = append(varListMet5, "PFLCU")
	varListMet5 = append(varListMet5, "PFLLSAN")

	// get all the chemistry variables from their MAPS
	varListChem = append(varListChem, "TIME-SER__AIRDEN") // air density (molec / cm3)
	varListChem = append(varListChem, "TIME-SER__OH")     // OH density (molec / cm3)
	varListChem = append(varListChem, "IJ-AVG-S__H2O2")
	for k := range aVOC {
		varListChem = append(varListChem, k)
	}
	for k := range aSOA {
		varListChem = append(varListChem, k)
	}
	for k := range bVOC {
		varListChem = append(varListChem, k)
	}
	for k := range bSOA {
		varListChem = append(varListChem, k)
	}
	for k := range NOx {
		varListChem = append(varListChem, k)
	}
	for k := range pNO {
		varListChem = append(varListChem, k)
	}
	for k := range SOx {
		varListChem = append(varListChem, k)
	}
	for k := range pS {
		varListChem = append(varListChem, k)
	}
	for k := range NH3 {
		varListChem = append(varListChem, k)
	}
	for k := range pNH {
		varListChem = append(varListChem, k)
	}

	varAllList := []string{}
	varAllList = append(varAllList, varListMet1...)
	varAllList = append(varAllList, varListMet2...)
	varAllList = append(varAllList, varListMet3...)
	varAllList = append(varAllList, varListMet4...)
	varAllList = append(varAllList, varListMet5...)
	varAllList = append(varAllList, varListChem...)

	// calculate wind speed and direction
	windDirectionChan := make(chan []*sparse.DenseArray)
	windDirectionChanOut := make(chan []*sparse.DenseArray)
	go calcWindDirection(windDirectionChan, windDirectionChanOut, varAllList)

	uvwAvgChan := make(chan []*sparse.DenseArray)
	uvwAvgChanOut := make(chan *sparse.DenseArray)
	go calcWindSpeed(uvwAvgChan, uvwAvgChanOut, varAllList)

	// get average pblh height
	pblhChan := make(chan []*sparse.DenseArray)
	avgPblhOut := make(chan *sparse.DenseArray)
	go average(pblhChan, avgPblhOut, "PBLH", varAllList)

	// get average temps
	tempChan := make(chan []*sparse.DenseArray)
	avgTempOut := make(chan *sparse.DenseArray)
	go average(tempChan, avgTempOut, "T", varAllList)

	// get average Air Density
	adChan := make(chan []*sparse.DenseArray)
	adOut := make(chan *sparse.DenseArray)
	go average(adChan, adOut, "TIME-SER__AIRDEN", varAllList)

	// calculate partitioning
	aVOCaSOAchan := make(chan []*sparse.DenseArray)
	aSOAchanOut := make(chan *sparse.DenseArray)
	go calcPartitioning(aVOCaSOAchan, aSOAchanOut, aVOC, aSOA, varAllList)

	bVOCbSOAchan := make(chan []*sparse.DenseArray)
	bSOAchanOut := make(chan *sparse.DenseArray)
	go calcPartitioning(bVOCbSOAchan, bSOAchanOut, bVOC, bSOA, varAllList)

	// NEED to double check all MW assumptions
	// NEED bVOC, bSOA

	NOxPNOxchan := make(chan []*sparse.DenseArray)
	NOxPNOxchanOut := make(chan *sparse.DenseArray)
	go calcPartitioning(NOxPNOxchan, NOxPNOxchanOut, NOx, pNO, varAllList)

	SOxPSOxchan := make(chan []*sparse.DenseArray)
	SOxPSOxchanOut := make(chan *sparse.DenseArray)
	go calcPartitioning(SOxPSOxchan, SOxPSOxchanOut, SOx, pS, varAllList)

	NH3NH4chan := make(chan []*sparse.DenseArray)
	NH3NH4chanOut := make(chan *sparse.DenseArray)
	go calcPartitioning(NH3NH4chan, NH3NH4chanOut, NH3, pNH, varAllList)

	NONO2chan := make(chan []*sparse.DenseArray)
	NONO2chanOut := make(chan *sparse.DenseArray)
	go calcPartitioning(NONO2chan, NONO2chanOut, NO, NO2, varAllList)

	avgPMchan := make(chan []*sparse.DenseArray)
	avgPMchanOut := make(chan *sparse.DenseArray)
	go averagePMSum(avgPMchan, avgPMchanOut, varAllList, pS, pNH, pNO, aSOA, bSOA)

	wetDepchan := make(chan []*sparse.DenseArray)
	wetDepchanOut := make(chan []*sparse.DenseArray)
	go calcWetDeposition(wetDepchan, wetDepchanOut, varAllList)

	// Calculate stability for plume rise, vertical mixing,
	// and chemical reaction rates
	SMCchan := make(chan []*sparse.DenseArray)
	SMCchanOut := make(chan []*sparse.DenseArray)
	go StabilityMixingChemistry(SMCchan, SMCchanOut, varAllList)

	// READ CHEM AND MET FILES
	iterateGeosTimeSteps("GEOS, met and chem iteration: ",
		fileNames, varListMet1, varListMet2, varListMet3, varListMet4, varListMet5, varListChem,
		windDirectionChan, pblhChan, tempChan, uvwAvgChan, adChan,
		aVOCaSOAchan, bVOCbSOAchan, NOxPNOxchan, SOxPSOxchan, NH3NH4chan, NONO2chan, avgPMchan, SMCchan, wetDepchan)

	// met results
	windDirectionChan <- nil
	windDirData := <-windDirectionChanOut
	uPlusSpeed := windDirData[0]
	uMinusSpeed := windDirData[1]
	vPlusSpeed := windDirData[2]
	vMinusSpeed := windDirData[3]
	wPlusSpeed := windDirData[4]
	wMinusSpeed := windDirData[5]

	uvwAvgChan <- nil
	windSpeed := <-uvwAvgChanOut

	pblhChan <- nil
	pblh := <-avgPblhOut

	tempChan <- nil
	temperature := <-avgTempOut

	adChan <- nil
	alt := <-adOut
	alt = densityKGM3(alt)
	for i, _ := range alt.Elements {
		alt.Elements[i] = 1 / alt.Elements[i]
	}

	// partitioning results
	aVOCaSOAchan <- nil
	soaPart := <-aSOAchanOut
	sogAVG := <-aSOAchanOut
	soaAVG := <-aSOAchanOut

	bVOCbSOAchan <- nil
	bsoaPart := <-bSOAchanOut
	bsogAVG := <-bSOAchanOut
	bsoaAVG := <-bSOAchanOut

	NOxPNOxchan <- nil
	NOPartitioning := <-NOxPNOxchanOut
	gNO := <-NOxPNOxchanOut
	pNO := <-NOxPNOxchanOut

	SOxPSOxchan <- nil
	SPartitioning := <-SOxPSOxchanOut
	gS := <-SOxPSOxchanOut
	pS := <-SOxPSOxchanOut

	NH3NH4chan <- nil
	NHPartitioning := <-NH3NH4chanOut
	gNH := <-NH3NH4chanOut
	pNH := <-NH3NH4chanOut

	NONO2chan <- nil
	NO_NO2partitioning := <-NONO2chanOut
	<-NONO2chanOut
	<-NONO2chanOut

	avgPMchan <- nil
	totalpm25 := <-avgPMchanOut

	// SMC and layer height results
	SMCchan <- nil
	SMCoutArray := <-SMCchanOut
	particleDryDep := SMCoutArray[0]
	layerHeights := SMCoutArray[1]
	Dz := SMCoutArray[2]
	SO2oxidation := SMCoutArray[3]
	S1 := SMCoutArray[4]
	Sclass := SMCoutArray[5]
	Kzz := SMCoutArray[6]
	M2u := SMCoutArray[7]
	M2d := SMCoutArray[8]
	Kxxyy := SMCoutArray[9]
	SO2DryDep := SMCoutArray[10]
	NOxDryDep := SMCoutArray[11]
	NH3DryDep := SMCoutArray[12]
	VOCDryDep := SMCoutArray[13]

	// wet deposition results
	wetDepchan <- nil
	wetDepOutArray := <-wetDepchanOut
	particleWetDep := wetDepOutArray[0]
	SO2WetDep := wetDepOutArray[1]
	otherGasWetDep := wetDepOutArray[2]

	// write out data to file
	outputFile := filepath.Join(config.OutputDir, config.OutputFilePrefix+".ncf")
	fmt.Printf("Writing out data to %v...\n", outputFile)
	h := cdf.NewHeader(
		[]string{"x", "y", "z", "zStagger", "zChem"},
		[]int{uPlusSpeed.Shape[2], uPlusSpeed.Shape[1], uPlusSpeed.Shape[0],
			uPlusSpeed.Shape[0] + 1, soaAVG.Shape[0]})
	h.AddAttribute("", "comment", "Meteorology and baseline chemistry data file")

	//	uPlusSpeed = avgToChemLevel(uPlusSpeed)

	data := map[string]dataHolder{
		"LayerHeights": dataHolder{[]string{"zStagger", "y", "x"},
			"Height at edge of layer", "m", layerHeights},
		"Dz": dataHolder{[]string{"z", "y", "x"},
			"Vertical grid size", "m", Dz},
		"UPlusSpeed": dataHolder{[]string{"zChem", "y", "x"},
			"Average speed of wind going in +U direction", "m/s", uPlusSpeed},
		"UMinusSpeed": dataHolder{[]string{"z", "y", "x"},
			"Average speed of wind going in -U direction", "m/s", uMinusSpeed},
		"VPlusSpeed": dataHolder{[]string{"z", "y", "x"},
			"Average speed of wind going in +V direction", "m/s", vPlusSpeed},
		"VMinusSpeed": dataHolder{[]string{"z", "y", "x"},
			"Average speed of wind going in -V direction", "m/s", vMinusSpeed},
		"WPlusSpeed": dataHolder{[]string{"z", "y", "x"},
			"Average speed of wind going in +W direction", "m/s", wPlusSpeed},
		"WMinusSpeed": dataHolder{[]string{"z", "y", "x"},
			"Average speed of wind going in -W direction", "m/s", wMinusSpeed},
		"WindSpeed": dataHolder{[]string{"z", "y", "x"},
			"RMS wind speed", "m s-1", windSpeed},
		"Pblh": dataHolder{[]string{"y", "x"},
			"Planetary boundary layer height", "m", pblh},
		"ParticleDryDep": dataHolder{[]string{"z", "y", "x"},
			"Dry deposition velocity for particles", "m s-1", particleDryDep},
		"SO2DryDep": dataHolder{[]string{"z", "y", "x"},
			"Dry deposition velocity for SO2", "m s-1", SO2DryDep},
		"NOxDryDep": dataHolder{[]string{"z", "y", "x"},
			"Dry deposition velocity for NOx", "m s-1", NOxDryDep},
		"NH3DryDep": dataHolder{[]string{"z", "y", "x"},
			"Dry deposition velocity for NH3", "m s-1", NH3DryDep},
		"VOCDryDep": dataHolder{[]string{"z", "y", "x"},
			"Dry deposition velocity for VOCs", "m s-1", VOCDryDep},
		"ParticleWetDep": dataHolder{[]string{"z", "y", "x"},
			"Wet deposition rate constant for fine particles",
			"s-1", particleWetDep},
		"SO2WetDep": dataHolder{[]string{"z", "y", "x"},
			"Wet deposition rate constant for SO2 gas", "s-1", SO2WetDep},
		"OtherGasWetDep": dataHolder{[]string{"z", "y", "x"},
			"Wet deposition rate constant for other gases", "s-1", otherGasWetDep},
		"NOPartitioning": dataHolder{[]string{"zChem", "y", "x"},
			"Mass fraction of N from NOx in particle {vs. gas} phase", "fraction",
			NOPartitioning},
		"gNO": dataHolder{[]string{"zChem", "y", "x"},
			"Average concentration of nitrogen fraction of gaseous NOx", "ug m-3",
			gNO},
		"pNO": dataHolder{[]string{"zChem", "y", "x"},
			"Average concentration of nitrogen fraction of particulate NO3",
			"ug m-3", pNO},
		"NHPartitioning": dataHolder{[]string{"z", "y", "x"},
			"Mass fraction of N from NH3 in particle {vs. gas} phase", "fraction",
			NHPartitioning},
		"gNH": dataHolder{[]string{"z", "y", "x"},
			"Average concentration of nitrogen fraction of gaseous ammonia",
			"ug m-3", gNH},
		"pNH": dataHolder{[]string{"z", "y", "x"},
			"Average concentration of nitrogen fraction of particulate ammonium",
			"ug m-3", pNH},
		"NO_NO2partitioning": dataHolder{[]string{"z", "y", "x"},
			"Mass fraction of N in NOx that exists as NO.", "fraction",
			NO_NO2partitioning},
		"SO2oxidation": dataHolder{[]string{"z", "y", "x"},
			"Rate of SO2 oxidation to SO4 by hydroxyl radical and H2O2",
			"s-1", SO2oxidation},
		"SPartitioning": dataHolder{[]string{"z", "y", "x"},
			"Mass fraction of S from SOx in particle {vs. gas} phase", "fraction",
			SPartitioning},
		"gS": dataHolder{[]string{"z", "y", "x"},
			"Average concentration of sulfur fraction of gaseous SOx", "ug m-3",
			gS},
		"pS": dataHolder{[]string{"z", "y", "x"},
			"Average concentration of sulfur fraction of particulate sulfate",
			"ug m-3", pS},
		"aOrgPartitioning": dataHolder{[]string{"zChem", "y", "x"},
			"Mass fraction of anthropogenic organic matter in particle {vs. gas} phase",
			"fraction", soaPart},
		"aVOC": dataHolder{[]string{"zChem", "y", "x"},
			"Average anthropogenic VOC concentration", "ug m-3", sogAVG},
		"aSOA": dataHolder{[]string{"zChem", "y", "x"},
			"Average anthropogenic secondary organic aerosol concentration", "ug m-3", soaAVG},
		"bOrgPartitioning": dataHolder{[]string{"z", "y", "x"},
			"Mass fraction of biogenic organic matter in particle {vs. gas} phase",
			"fraction", bsoaPart},
		"bVOC": dataHolder{[]string{"z", "y", "x"},
			"Average biogenic VOC concentration", "ug m-3", bsogAVG},
		"bSOA": dataHolder{[]string{"z", "y", "x"},
			"Average biogenic secondary organic aerosol concentration", "ug m-3", bsoaAVG},
		"Kzz": dataHolder{[]string{"zStagger", "y", "x"},
			"Vertical turbulent diffusivity", "m2 s-1", Kzz},
		"M2u": dataHolder{[]string{"z", "y", "x"},
			"ACM2 nonlocal upward mixing {Pleim 2007}", "s-1", M2u},
		"M2d": dataHolder{[]string{"z", "y", "x"},
			"ACM2 nonlocal downward mixing {Pleim 2007}", "s-1", M2d},
		"Kxxyy": dataHolder{[]string{"z", "y", "x"},
			"Horizontal eddy diffusion coefficient", "m2 s-1", Kxxyy},
		"Temperature": dataHolder{[]string{"z", "y", "x"},
			"Average Temperature", "K", temperature},
		"S1": dataHolder{[]string{"z", "y", "x"},
			"Stability parameter", "?", S1},
		"Sclass": dataHolder{[]string{"z", "y", "x"},
			"Stability parameter", "0=Unstable; 1=Stable", Sclass},
		"alt": dataHolder{[]string{"z", "y", "x"},
			"Inverse density", "m3 kg-1", alt},
		"TotalPM25": dataHolder{[]string{"z", "y", "x"},
			"Total PM2.5 concentration", "ug m-3", totalpm25}}

	testB := data["UPlusSpeed"].dims
	fmt.Println("data[UPlusSpeed].dims = ", testB)

	for name, d := range data {
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
	for name, d := range data {
		writeNCF(f, name, d.data)
	}
	err = cdf.UpdateNumRecs(ff)
	if err != nil {
		panic(err)
	}
	ff.Close()
	variableGridBasic(data)

} // End of main

func getOrderMap(iMap map[string]float64, varAllList []string) map[string]int {
	// Note, when itereating over maps, you do not get the same order
	// each time. So this function returns a second map with the order of the string in varAllList.
	order := make(map[string]int)
	for pol := range iMap {
		for i, polC := range varAllList {
			if polC == pol {
				order[pol] = i
			}
		}
	}
	return order
}

func getOrder(pol string, varAllList []string) int {
	for i, polC := range varAllList {
		if polC == pol {
			return i
		}
	}
	msg := "getOrder did not find " + pol
	panic(msg)
	//return order
}

func layerPresAvg(ps float64, k int) float64 {
	// units of hPa
	lp := ps*bp.Elements[k] + ap.Elements[k]
	lp += ps*bp.Elements[k+1] + ap.Elements[k+1]
	lp /= 2
	return lp
}

func unStag3d(in *sparse.DenseArray) *sparse.DenseArray {
	out := sparse.ZerosDense(in.Shape[0]-1, in.Shape[1], in.Shape[2])
	for i := 0; i < in.Shape[2]; i++ {
		for j := 0; j < in.Shape[1]; j++ {
			for k := 0; k < in.Shape[0]-1; k++ {
				out.Set((in.Get(k, j, i)+in.Get(k+1, j, i))/2., k, j, i)
			}
		}
	}
	return out
}

func densityKGM3(ad *sparse.DenseArray) *sparse.DenseArray {
	out := sparse.ZerosDense(ad.Shape...)
	for i, val := range ad.Elements {
		out.Elements[i] += val * (MWa / avNum) * 1000.
	}
	return out
}
func findQrain(pflcuA, pfllsanA, cf, Δz, airDen *sparse.DenseArray) *sparse.DenseArray {
	qrain := sparse.ZerosDense(airDen.Shape...)
	for i := 0; i < len(qrain.Elements); i++ {
		if cf.Elements[i] > 0.001 {
			qrain.Elements[i] += ((pflcuA.Elements[i] + pfllsanA.Elements[i]) * 60 * 60 * float64(config.MetTimeStep)) /
				(Δz.Elements[i] * airDen.Elements[i])
			//(cf.Elements[i] * Δz.Elements[i] * airDen.Elements[i])
		}
	}
	return qrain
}

func calcWetDeposition(wetDepchan, wetDepchanOut chan []*sparse.DenseArray, varAllList []string) {
	var wdParticle, wdSO2, wdOtherGas *sparse.DenseArray
	var airDen, pflcuA, pfllsanA, Δz, qrain *sparse.DenseArray

	airDenI := getOrder("TIME-SER__AIRDEN", varAllList) // Chem
	tI := getOrder("T", varAllList)
	psI := getOrder("PS", varAllList)
	cfI := getOrder("CLOUD", varAllList)
	pflcuI := getOrder("PFLCU", varAllList)
	pfllsanI := getOrder("PFLLSAN", varAllList)

	firstData := true
	for {
		data := <-wetDepchan
		if data == nil {
			outArray := make([]*sparse.DenseArray, 3)
			outArray[0] = arrayAverage(wdParticle)
			outArray[1] = arrayAverage(wdSO2)
			outArray[2] = arrayAverage(wdOtherGas)
			wetDepchanOut <- outArray
			return
		}
		// cloudFrac := <-cloudFracChan // frac
		// alt := <-altChan             // m3/kg
		if firstData {
			wdParticle = sparse.ZerosDense(data[tI].Shape...) // units = 1/s
			wdSO2 = sparse.ZerosDense(data[tI].Shape...)      // units = 1/s
			wdOtherGas = sparse.ZerosDense(data[tI].Shape...) // units = 1/s
			airDen = sparse.ZerosDense(data[tI].Shape...)
			pflcuA = sparse.ZerosDense(data[tI].Shape...)
			pfllsanA = sparse.ZerosDense(data[tI].Shape...)
			Δz = sparse.ZerosDense(data[tI].Shape...)
			qrain = sparse.ZerosDense(data[tI].Shape...)
			firstData = false
		}

		_, Δz = calcLayerHeights(data[tI], data[psI])

		pflcuA = unStag3d(data[pflcuI])
		pfllsanA = unStag3d(data[pfllsanI])

		airDen = densityKGM3(data[airDenI]) // convert airDen to kg/m3

		qrain = findQrain(pflcuA, pfllsanA, data[cfI], Δz, airDen)

		for i := 0; i < len(qrain.Elements); i++ {
			wdp, wds, wdo := emep.WetDeposition(data[cfI].Elements[i],
				qrain.Elements[i], airDen.Elements[i], Δz.Elements[i])
			wdParticle.Elements[i] += wdp
			wdSO2.Elements[i] += wds
			wdOtherGas.Elements[i] += wdo
		}
	}

}

func StabilityMixingChemistry(inChan, outChan chan []*sparse.DenseArray, varAllList []string) {
	var u, h, hflux, ρ, gocartObk, zo, p float64
	var kPblTop int
	var particleDryDep *sparse.DenseArray
	var SO2DryDep *sparse.DenseArray
	var NOxDryDep *sparse.DenseArray
	var NH3DryDep *sparse.DenseArray
	var VOCDryDep *sparse.DenseArray
	var layerHeightsAvg, DzAvg *sparse.DenseArray
	var SO2oxidation *sparse.DenseArray
	var S1 *sparse.DenseArray
	var Sclass *sparse.DenseArray
	var Kzz *sparse.DenseArray
	var M2d *sparse.DenseArray
	var M2u *sparse.DenseArray
	var Kyy *sparse.DenseArray

	airDenI := getOrder("TIME-SER__AIRDEN", varAllList) // Chem
	ohI := getOrder("TIME-SER__OH", varAllList)         // Chem
	h2o2I := getOrder("IJ-AVG-S__H2O2", varAllList)     // Chem
	tI := getOrder("T", varAllList)                     // I3
	ustarI := getOrder("USTAR", varAllList)             // A1
	hfluxI := getOrder("HFLUX", varAllList)             // A1
	z0I := getOrder("Z0M", varAllList)                  // A1
	psI := getOrder("PS", varAllList)                   // I3
	pblhI := getOrder("PBLH", varAllList)               // A1
	qlI := getOrder("QL", varAllList)                   // A3cld
	precipFluxI := getOrder("PRECTOT", varAllList)      // A1
	pardfI := getOrder("PARDF", varAllList)             // A1
	pardrI := getOrder("PARDR", varAllList)             // A1

	const (
		po    = 101300. // Pa, reference pressure
		kappa = 0.2854  // related to von karman's constant
		Cp    = 1006.   // m2/s2-K; specific heat of air
	)

	firstData := true
	for {
		data := <-inChan
		if data == nil {
			layerHeightsAvg = arrayAverage(layerHeightsAvg)
			// Check for mass balance in convection coefficients
			fmt.Println("M2u and M2d are not passing the mass balance test!")
			/*
				for k := 0; k < M2u.Shape[0]-2; k++ {
					for j := 0; j < M2u.Shape[1]; j++ {
						for i := 0; i < M2u.Shape[2]; i++ {
							z := layerHeightsAvg.Get(k, j, i)
							zabove := layerHeightsAvg.Get(k+1, j, i)
							z2above := layerHeightsAvg.Get(k+2, j, i)
							Δzratio := (z2above - zabove) / (zabove - z)
							m2u := M2u.Get(k, j, i)
							val := m2u - M2d.Get(k, j, i) +
								M2d.Get(k+1, j, i)*Δzratio
							if math.Abs(val/m2u) > 1.e-8 {
								panic(fmt.Errorf("M2u and M2d don't match: "+
									"(k,j,i)=(%v,%v,%v); val=%v; m2u=%v; "+
									"m2d=%v, m2dAbove=%v",
									k, j, i, val, m2u, M2d.Get(k, j, i),
									M2d.Get(k+1, j, i)))
							}
						}
					}
				} */
			// convert Kzz to unstaggered grid
			KzzUnstaggered := sparse.ZerosDense(SO2oxidation.Shape...)
			for j := 0; j < KzzUnstaggered.Shape[1]; j++ {
				for i := 0; i < KzzUnstaggered.Shape[2]; i++ {
					for k := 0; k < KzzUnstaggered.Shape[0]; k++ {
						KzzUnstaggered.Set(
							(Kzz.Get(k, j, i)+Kzz.Get(k+1, j, i))/2.,
							k, j, i)
					}
				}
			}

			outArrays := make([]*sparse.DenseArray, 14)
			outArrays[0] = arrayAverage(particleDryDep)
			outArrays[1] = layerHeightsAvg
			outArrays[2] = arrayAverage(DzAvg)
			outArrays[3] = arrayAverage(SO2oxidation)
			outArrays[4] = arrayAverage(S1)
			outArrays[5] = arrayAverage(Sclass)
			outArrays[6] = arrayAverage(KzzUnstaggered)
			outArrays[7] = arrayAverage(M2u)
			outArrays[8] = arrayAverage(M2d)
			outArrays[9] = arrayAverage(Kyy)
			outArrays[10] = arrayAverage(SO2DryDep)
			outArrays[11] = arrayAverage(NOxDryDep)
			outArrays[12] = arrayAverage(NH3DryDep)
			outArrays[13] = arrayAverage(VOCDryDep)
			outChan <- outArrays
			return
		}
		if firstData {
			layerHeightsAvg = sparse.ZerosDense(data[tI].Shape[0]+1, data[tI].Shape[1], data[tI].Shape[2])
			DzAvg = sparse.ZerosDense(data[tI].Shape...)

			S1 = sparse.ZerosDense(data[tI].Shape...)
			Sclass = sparse.ZerosDense(data[tI].Shape...)
			particleDryDep = sparse.ZerosDense(data[tI].Shape...) // units = m/s
			SO2DryDep = sparse.ZerosDense(data[tI].Shape...)      // units = m/s
			NOxDryDep = sparse.ZerosDense(data[tI].Shape...)      // units = m/s
			NH3DryDep = sparse.ZerosDense(data[tI].Shape...)      // units = m/s
			VOCDryDep = sparse.ZerosDense(data[tI].Shape...)      // units = m/s
			Kzz = sparse.ZerosDense(layerHeightsAvg.Shape...)     // units = m2/s
			M2u = sparse.ZerosDense(data[tI].Shape...)            // units = 1/s
			M2d = sparse.ZerosDense(data[tI].Shape...)            // units = 1/s
			SO2oxidation = sparse.ZerosDense(data[tI].Shape...)   // units = 1/s
			Kyy = sparse.ZerosDense(data[tI].Shape...)            // units = m2/s

			firstData = false
		}

		layerHeights, Dz := calcLayerHeights(data[tI], data[psI])
		layerHeightsAvg.AddDense(layerHeights)
		DzAvg.AddDense(Dz)

		for i := 0; i < layerHeights.Shape[2]; i++ {
			for j := 0; j < layerHeights.Shape[1]; j++ {

				// 2-D:

				for k := 0; k < layerHeights.Shape[0]; k++ {
					if layerHeights.Get(k, j, i) >= data[pblhI].Get(j, i) {
						kPblTop = k
						break
					}
				}
				// Calculate boundary layer average temperature (K)
				To := 0.
				for k := 0; k < layerHeights.Shape[0]; k++ {
					if k == kPblTop {
						To /= float64(k)
						break
					}
					To += data[tI].Get(k, j, i)
				}

				// Calculate convective mixing rate
				u = data[ustarI].Get(j, i) // friction velocity
				h = layerHeights.Get(kPblTop, j, i)
				hflux = data[hfluxI].Get(j, i)                         // heat flux [W m-2]
				ρ = data[airDenI].Get(0, j, i) * (MWa / avNum) * 1000. // density [kg/m3]
				L := acm2.ObukhovLen(hflux, ρ, To, u)                  // Monin-Obukhov length [m]
				fconv := acm2.ConvectiveFraction(L, h)
				m2u := acm2.M2u(layerHeights.Get(1, j, i),
					layerHeights.Get(2, j, i), h, L, u, fconv)

				// Calculate dry deposition
				p = data[psI].Get(j, i)
				p = layerPresAvg(p, 0) * 100 // convert hPa to Pa
				// z: [m] surface layer; assumed to be 10% of boundary layer.
				// with a max of 100 meters (simplification here...)
				z := h / 10.
				if z > 100. {
					z = 100
				}
				gocartObk = gocart.ObhukovLen(hflux, ρ, To, u)

				zo = data[z0I].Get(j, i) // NOTE: GEOS roughness larger than WRF -- larger grid size? Incorporates topography?
				// Why is WRF look up roughness so small?!? all land types < meter? Could that be right?
				//USGSz0[lu] // roughness length [m]

				weselyLU := wesely1989.Range // middle value (range land) as GEOS does not have USGS lu categories
				if FRLAND.Get(j, i) < 0.51 {
					weselyLU = wesely1989.Water
				}

				const dParticle = 0.3e-6 // [m], Seinfeld & Pandis fig 8.11
				const ρparticle = 1830.  // [kg/m3] Jacobson (2005) Ex. 13.5
				const Θsurface = 0.      // surface slope [rad]; Assume surface is flat.

				// Dry deposition for Particles:
				particleDryDep.AddVal(
					gocart.ParticleDryDep(gocartObk, u, To, h, zo, dParticle/2., ρparticle, p), 0, j, i)

				// Dry deposition for gases

				// This is not the best way to tell what season it is.
				//var iSeasonP seinfeld.SeasonalCategory // for particles
				var iSeasonG wesely1989.SeasonCategory // for gases
				switch {
				case To > 273.+20.:
					//iSeasonP = seinfeld.Midsummer
					iSeasonG = wesely1989.Midsummer
				case To <= 273.+20 && To > 273.+10.:
					//iSeasonP = seinfeld.Autumn
					iSeasonG = wesely1989.Autumn
				case To <= 273.+10 && To > 273.+0.:
					//iSeasonP = seinfeld.LateAutumn
					iSeasonG = wesely1989.LateAutumn
				default:
					//iSeasonP = seinfeld.Winter
					iSeasonG = wesely1989.Winter
				}

				const dew = false // don't know if there's dew.
				rain := data[precipFluxI].Get(j, i) > 1.e-9

				G := data[pardfI].Get(j, i) + data[pardrI].Get(j, i) // + glw.Get(j, i) // irradiation [W/m2]

				// Dry deposition for gases:

				SO2DryDep.AddVal(
					seinfeld.DryDepGas(z, zo, u, L, To, ρ,
						G, Θsurface,
						wesely1989.So2Data, iSeasonG,
						weselyLU, rain, dew, true, false), 0, j, i)
				NOxDryDep.AddVal(
					seinfeld.DryDepGas(z, zo, u, L, To, ρ,
						G, Θsurface,
						wesely1989.No2Data, iSeasonG,
						weselyLU, rain, dew, false, false), 0, j, i)
				NH3DryDep.AddVal(
					seinfeld.DryDepGas(z, zo, u, L, To, ρ,
						G, Θsurface,
						wesely1989.Nh3Data, iSeasonG,
						weselyLU, rain, dew, false, false), 0, j, i)
				VOCDryDep.AddVal(
					seinfeld.DryDepGas(z, zo, u, L, To, ρ,
						G, Θsurface,
						wesely1989.OraData, iSeasonG,
						weselyLU, rain, dew, false, false), 0, j, i)

				// 3-D:
				for k := 0; k < data[airDenI].Shape[0]; k++ {
					pL := layerPresAvg(data[psI].Get(j, i), k) * 100

					// Stability
					var dtheta_dz = 0. // potential temperature gradient
					Tval := data[tI].Get(k, j, i) * math.Pow(po/pL, kappa)
					if k > 0 {
						pM1 := layerPresAvg(data[psI].Get(j, i), k-1) * 100
						TvalM1 := data[tI].Get(k-1, j, i) * math.Pow(po/pM1, kappa)
						dtheta_dz = (Tval - TvalM1) /
							(layerHeights.Get(k, j, i) -
								layerHeights.Get(k-1, j, i)) // K/m
					}
					// Stability parameter
					pressureCorrection := math.Pow(pL/po, kappa)
					s1 := dtheta_dz / data[tI].Get(k, j, i) * pressureCorrection // not sure why
					// there is a pressureCorrection here (following wrf2inmap.go, which has 2 pressureCorrections).
					S1.AddVal(s1, k, j, i)
					// Stability class
					if dtheta_dz < 0.005 {
						Sclass.AddVal(0., k, j, i)
					} else {
						Sclass.AddVal(1., k, j, i)
					}

					// Mixing
					z := layerHeights.Get(k, j, i)
					zabove := layerHeights.Get(k+1, j, i)
					zcenter := (layerHeights.Get(k, j, i) +
						layerHeights.Get(k+1, j, i)) / 2
					Δz := zabove - z

					const freeAtmKzz = 3. // [m2 s-1]
					if k >= kPblTop {     // free atmosphere (unstaggered grid)
						Kzz.AddVal(freeAtmKzz, k, j, i)
						Kyy.AddVal(freeAtmKzz, k, j, i)
						if k == data[tI].Shape[0]-1 { // Top Layer
							Kzz.AddVal(freeAtmKzz, k+1, j, i)
						}
					} else { // Boundary layer (unstaggered grid)
						Kzz.AddVal(acm2.Kzz(z, h, L, u, fconv), k, j, i)
						M2d.AddVal(acm2.M2d(m2u, z, Δz, h), k, j, i)
						M2u.AddVal(m2u, k, j, i)
						kmyy := acm2.CalculateKm(zcenter, h, L, u)
						Kyy.AddVal(kmyy, k, j, i)

						//////////////////////////
						//                                              m2d := acm2.M2d(m2u, z, Δz, h)
						//                                      z2 := LayerHeights.Get(k+1, j, i)
						//                                      Δz2 := LayerHeights.Get(k+1, j, i) - z2
						//                                              m2d2 := acm2.M2d(m2u, z2, Δz2, h)

						/////////////////////////
					}

					// Gas phase sulfur chemistry
					// const cm3perm3 = 100. * 100. * 100.
					// const molarMassAir = 28.97 / 1000.             // kg/mol
					// const airFactor = molarMassAir / Na * cm3perm3 // kg/molec.* cm3/m3

					M := data[airDenI].Get(k, j, i) // molec. air / cm3

					hoConc := data[ohI].Get(k, j, i) // molec. HO / cm3

					// SO2 oxidation rate (Stockwell 1997, Table 2d)
					const kinf = 1.5e-12
					ko := 3.e-31 * math.Pow(data[tI].Get(k, j, i)/300., -3.3)
					SO2rate := (ko * M / (1 + ko*M/kinf)) * math.Pow(0.6,
						1./(1+math.Pow(math.Log10(ko*M/kinf), 2.))) // cm3/molec/s
					kSO2 := SO2rate * hoConc

					// Aqueous phase sulfur chemistry
					qCloudVal := data[qlI].Get(k, j, i)
					if qCloudVal > 0. {
						const pH = 3.5 // doesn't really matter for SO2
						qCloudVal *=
							data[airDenI].Get(k, j, i) * (MWa / avNum) // convert to volume frac. (double checked units here)
						kSO2 += seinfeld.SulfurH2O2aqueousOxidationRate(
							data[h2o2I].Get(k, j, i), pH, data[tI].Get(k, j, i), pL*atmPerPa,
							qCloudVal)
					}
					SO2oxidation.AddVal(kSO2, k, j, i) // 1/s
				}

			}
		}
	}
}

// marginal partitioning
// calcPartitioning(aVOCaSOAchan, aSOAchanOut, aVOC, aSOA, varAllList)
func calcPartitioning(polchan chan []*sparse.DenseArray, outchan chan *sparse.DenseArray, gMap, pMap map[string]float64, varAllList []string) {
	var gas, particle, gasdata, particledata, oldgas, oldparticle, partitioning *sparse.DenseArray
	var gOrder, pOrder map[string]int
	gOrder = getOrderMap(gMap, varAllList)
	pOrder = getOrderMap(pMap, varAllList)
	airDenName := "TIME-SER__AIRDEN"
	airDenI := getOrder(airDenName, varAllList)

	firstData := true
	for {
		data := <-polchan
		if data == nil {
			outchan <- arrayAverage(partitioning)
			outchan <- arrayAverage(gas)
			outchan <- arrayAverage(particle)
			return
		}

		if firstData {
			fmt.Println(gOrder)
			fmt.Println(pOrder)

			// In the first time step, just copy the arrays to the
			// old arrays; don't do any calculations.
			partitioning = sparse.ZerosDense(data[airDenI].Shape...)
			gas = sparse.ZerosDense(data[airDenI].Shape...)
			particle = sparse.ZerosDense(data[airDenI].Shape...)
			gasdata = sparse.ZerosDense(data[airDenI].Shape...)
			particledata = sparse.ZerosDense(data[airDenI].Shape...)

			// now sum tot-gas mass and tot-particle mass
			oldgas = sparse.ZerosDense(data[airDenI].Shape...)
			oldparticle = sparse.ZerosDense(data[airDenI].Shape...)
			for pol, factor := range gMap {
				for i, val := range data[gOrder[pol]].Elements {
					// convert ppb to μg/m3 -- using air density from GEOS chem
					oldgas.Elements[i] += val * factor * (1. / 1000) * (data[airDenI].Elements[i] * 1.0e6 / avNum)
				}
			}
			for pol, factor := range pMap {
				for i, val := range data[pOrder[pol]].Elements {
					// still to do:
					// convert ????? to ug/m3, listed as ppb, what does that mean for particles?
					// maybe the same thing as gas, they just give a different MW that converts it.
					// oldparticle.Elements[i] += val * factor / alt.Elements[i]
					oldparticle.Elements[i] += val * factor * (1. / 1000) * (data[airDenI].Elements[i] * 1.0e6 / avNum)
				}
			}
			firstData = false
			continue
		}

		gasdata.Scale(0)
		particledata.Scale(0)
		for pol, factor := range gMap { // add new gas mass
			for i, val := range data[gOrder[pol]].Elements {
				// convert ppm to μg/m3 -- no alt in GEOS chem, will need to update
				// gasdata.Elements[i] += val * factor / MWa * 1000. / alt.Elements[i]
				gasdata.Elements[i] += val * factor * (1. / 1000) * (data[airDenI].Elements[i] * 1.0e6 / avNum)
				gas.Elements[i] += val * factor * (1. / 1000) * (data[airDenI].Elements[i] * 1.0e6 / avNum)
			}
		}
		for pol, factor := range pMap { // add new particles
			for i, val := range data[pOrder[pol]].Elements {
				// convert ppm to μg/m3 -- no alt in GEOS chem, will need to update
				// particle.Elements[i] += val * factor * (1. / 1000/) * molesPm3a[i]
				// particledata.Elements[i] += val * factor / alt.Elements[i]
				particledata.Elements[i] += val * factor * (1. / 1000) * (data[airDenI].Elements[i] * 1.0e6 / avNum)
				particle.Elements[i] += val * factor * (1. / 1000) * (data[airDenI].Elements[i] * 1.0e6 / avNum)
			}
		}
		for i, particleval := range particledata.Elements {
			particlechange := particleval - oldparticle.Elements[i]
			totalchange := particlechange + (gasdata.Elements[i] -
				oldgas.Elements[i])
			// Calculate the marginal partitioning coefficient, which is the
			// change in particle concentration divided by the change in overall
			// concentration. Force the coefficient to be between zero and
			// one.
			partitioning.Elements[i] +=
				math.Min(math.Max(particlechange/totalchange, 0), 1)
		}
		oldgas = gasdata.Copy()
		oldparticle = particledata.Copy()
	}
}

func average(datachan chan []*sparse.DenseArray, outChan chan *sparse.DenseArray, name string, varAllList []string) {
	var avgdata *sparse.DenseArray

	nameI := getOrder(name, varAllList)

	firstData := true
	for {
		data := <-datachan
		if data == nil {
			for i, val := range avgdata.Elements {
				avgdata.Elements[i] = val / numTsteps
			}
			outChan <- avgdata
			return
		}
		if firstData {
			avgdata = sparse.ZerosDense(data[nameI].Shape...)
			firstData = false
		}
		avgdata.AddDense(data[nameI])
	}
}

func averagePMSum(datachan chan []*sparse.DenseArray, outChan chan *sparse.DenseArray, varAllList []string, maps ...map[string]float64) {
	var avgdata *sparse.DenseArray
	var order []int
	var mw []float64

	airDenI := getOrder("TIME-SER__AIRDEN", varAllList)

	for _, m := range maps {
		for n, mwn := range m {
			order = append(order, getOrder(n, varAllList))
			mw = append(mw, mwn)
		}
	}

	firstData := true
	for {
		data := <-datachan
		if data == nil {
			for i, val := range avgdata.Elements {
				avgdata.Elements[i] = val / numTsteps
			}
			outChan <- avgdata
			return
		}
		if firstData {
			avgdata = sparse.ZerosDense(data[order[0]].Shape...)
			firstData = false
		}
		for i, I := range order {
			for ii, val := range data[I].Elements {
				val *= (mw[i] * (1. / 1000) * (data[airDenI].Elements[ii] * 1.0e6 / avNum))
				avgdata.Elements[ii] += val
			}
		}
	}

}

func getFixed() {
	fixedFile := config.GeosFixed
	f := new(cdfFile)
	var err error
	f.ff, err = os.Open(fixedFile)
	if err != nil {
		panic(err)
	}
	f.f, err = cdf.Open(f.ff)
	if err != nil {
		panic(err)
	}
	// --------- based loosely on readNCF:
	pol := "FRLAND"
	dims := f.f.Header.Lengths(pol)
	if len(dims) == 0 {
		panic(fmt.Sprintf("Variable %v not in fixed file.", pol))
	}
	dims = dims[1:]
	FRLAND = sparse.ZerosDense(dims[0], dims[1])
	nread := 1
	for _, dim := range dims {
		nread *= dim
	}
	r := f.f.Reader(pol, nil, nil)
	buf := r.Zero(nread)
	//      f.ncfLock.Lock()
	n, err := r.Read(buf)
	//      f.ncfLock.Unlock()
	if err != nil {
		panic(err)
	}
	fmt.Println(n)
	fmt.Println("in getFixed, dims =", dims)
	for i, val := range buf.([]float32) {
		FRLAND.Elements[i] = float64(val)
	}

	//pol = "Bp"
	//rb := f.f.Reader(pol, nil, nil)
	//bufb := rb.Zero(nread)
	//nb, err := rb.Read(bufb)
	//if err != nil {
	//	panic(err)
	//}
	//fmt.Println(nb)
	//for i, val := range bufb.([]float32) {
	//	bp.Elements[i] = float64(val)
	//}
	f.ff.Close()

}

func getApBp() {
	d := config.StartDate
	h := "03"
	cfilename := strings.Replace(config.GeosChem, "[DATE]", d, -1)
	cfilename = strings.Replace(cfilename, "[HOUR]", h, -1)
	layerFile := cfilename // config.GeosApBp

	f := new(cdfFile)
	var err error
	f.ff, err = os.Open(layerFile)
	if err != nil {
		panic(err)
	}
	f.f, err = cdf.Open(f.ff)
	if err != nil {
		panic(err)
	}
	// --------- based loosely on readNCF:
	pol := "Ap"
	dims := f.f.Header.Lengths(pol)
	if len(dims) == 0 {
		panic(fmt.Sprintf("Variable %v not in ApBp file.", pol))
	}
	// fmt.Println("ap dims =", dims)
	nread := dims[0]
	r := f.f.Reader(pol, nil, nil)
	buf := r.Zero(nread)
	//      f.ncfLock.Lock()
	n, err := r.Read(buf)
	//      f.ncfLock.Unlock()
	if err != nil {
		panic(err)
	}
	fmt.Println(n)
	fmt.Println("in getApBp, dims =", dims)
	for i, val := range buf.([]float32) {
		ap.Elements[i] = float64(val)
	}

	pol = "Bp"
	rb := f.f.Reader(pol, nil, nil)
	bufb := rb.Zero(nread)
	nb, err := rb.Read(bufb)
	if err != nil {
		panic(err)
	}
	fmt.Println(nb)
	for i, val := range bufb.([]float32) {
		bp.Elements[i] = float64(val)
	}
	f.ff.Close()
}

// calcLayerHeights calculates the heights above the ground
// of the layers (in meters). Using the hypsometric equation.
//
func calcLayerHeights(tL, spL *sparse.DenseArray) (*sparse.DenseArray, *sparse.DenseArray) {
	var layerHeights *sparse.DenseArray
	// var layerHeightsI *sparse.DenseArray
	var Dz *sparse.DenseArray
	var p *sparse.DenseArray

	layerHeights = sparse.ZerosDense(tL.Shape[0]+1, tL.Shape[1], tL.Shape[2])
	// layerHeightsI = sparse.ZerosDense(tL.Shape[0]+1, tL.Shape[1], tL.Shape[2])
	p = sparse.ZerosDense(tL.Shape[0]+1, tL.Shape[1], tL.Shape[2])
	Dz = sparse.ZerosDense(tL.Shape...)

	for k := 0; k < tL.Shape[0]+1; k++ {
		for j := 0; j < tL.Shape[1]; j++ {
			for i := 0; i < tL.Shape[2]; i++ {
				p.Set((spL.Get(j, i)*bp.Elements[k])+ap.Elements[k], k, j, i) // spL in units of hPa
				if k > 0 {
					t := tL.Get(k-1, j, i)                                                    // tL in units on K
					h := -1 * math.Log(float64(p.Get(k, j, i)/p.Get(k-1, j, i))) * rr * t / g // in meters
					// layerHeightsI.Set(h+layerHeightsI.Get(k-1, j, i), k, j, i)
					layerHeights.Set(h+layerHeights.Get(k-1, j, i), k, j, i)
					Dz.Set(h, k-1, j, i)
				}
				if k < 1 {
					h := float64(0)
					layerHeights.Set(h, k, j, i)
				}
			}
		}
	}
	return layerHeights, Dz
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

type cdfReaderFunc func([]*cdfFile, int, *sync.WaitGroup)

type cdfFile struct {
	f       *cdf.File
	ff      *os.File
	ncfLock sync.Mutex
}

func readChemNCF(pol string, f *cdfFile) (data *sparse.DenseArray) {
	dims := f.f.Header.Lengths(pol)
	if len(dims) == 0 {
		panic(fmt.Sprintf("Variable %v not in file.", pol))
	}
	nread := 1
	for _, dim := range dims {
		nread *= dim
	}
	//	start, end := make([]int, len(dims)+1), make([]int, len(dims)+1)
	//	start[0], end[0] = hour, hour+1
	r := f.f.Reader(pol, nil, nil)
	buf := r.Zero(nread)
	f.ncfLock.Lock()
	_, err := r.Read(buf)
	f.ncfLock.Unlock()
	if err != nil {
		panic(err)
	}
	fmt.Println("in readChemNCF, dims =", dims)
	data = sparse.ZerosDense(dims...)
	for i, val := range buf.([]float32) {
		data.Elements[i] = float64(val)
	}
	return

}

func readNCF(pol string, f *cdfFile, hour int) (data *sparse.DenseArray) {
	dims := f.f.Header.Lengths(pol)
	if len(dims) == 0 {
		panic(fmt.Sprintf("Variable %v not in file.", pol))
	}
	dims = dims[1:]
	nread := 1
	for _, dim := range dims {
		nread *= dim
	}
	start, end := make([]int, len(dims)+1), make([]int, len(dims)+1)
	start[0], end[0] = hour, hour+1
	r := f.f.Reader(pol, start, end)
	buf := r.Zero(nread)
	f.ncfLock.Lock()
	_, err := r.Read(buf)
	f.ncfLock.Unlock()
	if err != nil {
		panic(err)
	}
	fmt.Println("in readNCF, dims =", dims)
	data = sparse.ZerosDense(dims...)
	for i, val := range buf.([]float32) {
		data.Elements[i] = float64(val)
	}
	return
}

type dataHolder struct {
	dims        []string
	Description string
	Units       string
	data        *sparse.DenseArray
}

func processFile(datechan chan string, fileNames []string, finishchan chan int,
	funcs ...cdfReaderFunc) {

	totFiles := len(fileNames)
	f := make([]*cdfFile, totFiles)

	var err error
	for d := range datechan {

		for i := 0; i < totFiles; i++ {
			f[i] = new(cdfFile)
			filename := strings.Replace(fileNames[i], "[DATE]", d, -1)

			f[i].ff, err = os.Open(filename)
			if err != nil {
				panic(err)
			}
			f[i].f, err = cdf.Open(f[i].ff)
			if err != nil {
				panic(err)
			}
		}
		for hour := 0; hour < (24 / config.MetTimeStep); hour++ {
			fmt.Println("hour =", hour)
			var wg sync.WaitGroup
			wg.Add(len(funcs))
			for _, fn := range funcs {
				go fn(f, hour, &wg)
			}
			wg.Wait()
		}
		for i := 0; i < totFiles; i++ {
			f[i].ff.Close()
		}
	}
	finishchan <- 0
}

func readArrayVar(Var []string, fNum []int, tStep []int, datachans ...chan []*sparse.DenseArray) cdfReaderFunc {
	return func(f []*cdfFile, hour int, wg *sync.WaitGroup) {
		defer wg.Done()
		fmt.Println("got to readArrayVar function, hour =", hour)

		vLen := len(Var)
		data := make([]*sparse.DenseArray, vLen)
		if (2 * len(Var)) != (len(fNum) + len(tStep)) {
			panic("readArrayVar has unequal length inputs")
		}

		for i, VarI := range Var {
			if tStep[i] != config.MetTimeStep {
				// Here is where I add 1 hour or 3 hour time step issues.
				if tStep[i] != 1 {
					panic("time step other than 3 or 1 not implemented")
				}
				data1hr := make([]*sparse.DenseArray, 3)
				hour1hr := hour * 3
				for j := 0; j < 3; j++ {
					data1hr[j] = readNCF(VarI, f[fNum[i]], hour1hr+j)
				}
				data1hr[0].AddDense(data1hr[1])
				data1hr[0].AddDense(data1hr[2])
				data1hr[0].Scale(1.0 / 3.0)
				data[i] = data1hr[0]
			} else {
				data[i] = readNCF(VarI, f[fNum[i]], hour)
			}
		}
		for _, datachanI := range datachans {
			datachanI <- data
		}
	}
}

func StagToStagChemLevel(in3d *sparse.DenseArray) *sparse.DenseArray {
	var out *sparse.DenseArray
	if in3d.Shape[0] != metVert+1 {
		panic("not correct vertical layers on met data incoming for stag to stag mapping")
	}

	out = sparse.ZerosDense(chemVert+1, in3d.Shape[1], in3d.Shape[2])
	clev := 0

	for k := 0; k < in3d.Shape[0]; k++ {
		switch {
		case k == in3d.Shape[0]-1:
			clev = chemVert
		case cLev[k] == k:
			clev = k
		case cLev[k] != k && cLev[k] != cLev[k-1]:
			clev = cLev[k]
		default:
			continue
		}
		for j := 0; j < in3d.Shape[1]; j++ {
			for i := 0; i < in3d.Shape[2]; i++ {
				out.AddVal(in3d.Get(k, j, i), clev, j, i)
			}
		}
	}
	return out
}

func avgToChemLevel(in3d *sparse.DenseArray) *sparse.DenseArray {
	var out *sparse.DenseArray
	if in3d.Shape[0] != metVert {
		panic("not correct vertical layers on met data incoming")
	}

	out = sparse.ZerosDense(chemVert, in3d.Shape[1], in3d.Shape[2])
	clev := 0
	levCount := make([]float64, chemVert)
	for k := 0; k < metVert; k++ {
		clev = cLev[k]
		levCount[clev]++
	}
	for k := 0; k < metVert; k++ {
		clev = cLev[k]
		for j := 0; j < in3d.Shape[1]; j++ {
			for i := 0; i < in3d.Shape[2]; i++ {
				out.AddVal(in3d.Get(k, j, i)/levCount[clev], clev, j, i)
			}
		}
	}
	return out
}

func iterateGeosTimeSteps(msg string, fileNames []string,
	metVarNames1, metVarNames2, metVarNames3, metVarNames4, metVarNames5, chemVarNames []string,
	datachans ...chan []*sparse.DenseArray) {

	varLens := []int{len(metVarNames1), len(metVarNames2), len(metVarNames3), len(metVarNames4), len(metVarNames5), len(chemVarNames)}
	totFiles := len(varLens)
	totVars := 0
	for _, i := range varLens {
		totVars += i
	}
	data := make([]*sparse.DenseArray, totVars)
	dataM1 := make([]*sparse.DenseArray, totVars)

	f := make([]*cdfFile, totFiles-1)
	c := new(cdfFile)
	var err error

	hours := []string{"00", "03", "06", "09", "12", "15", "18", "21"}

	firstHour := true
	delta, _ := time.ParseDuration("24h")
	for now := start; now.Before(end); now = now.Add(delta) {
		d := now.Format(geosFormat)
		fmt.Println(msg + d + "...")
		for i := 0; i < totFiles-1; i++ {
			if varLens[i] > 0 {
				filename := strings.Replace(fileNames[i], "[DATE]", d, -1)
				f[i] = new(cdfFile)
				f[i].ff, err = os.Open(filename)
				if err != nil {
					panic(err)
				}
				f[i].f, err = cdf.Open(f[i].ff)
				if err != nil {
					panic(err)
				}
			}
		}

		// loop time, and keep opening and closing chem files
		for hourI := 0; hourI < 8; hourI++ { // base 3-hr time step, adjust for chem files (open/close)
			// and met file with hourly outputs (average to 3 hr)
			polNum := 0
			if varLens[0] > 0 {
				for _, polName := range metVarNames1 { // A3 file has times-avg values at 90, 270, 450 ...
					// all files set timing of A3 here
					data[polNum] = avgToChemLevel(readNCF(polName, f[0], hourI))
					polNum++
				}
			}
			if varLens[1] > 0 { // I3 file has time-inst values at 0, 180, 360 ...
				if firstHour {
					for _, polName := range metVarNames2 {
						if polName == "PS" {
							data[polNum] = readNCF(polName, f[1], hourI)
						} else {
							data[polNum] = avgToChemLevel(readNCF(polName, f[1], hourI))
						}
						dataM1[polNum] = data[polNum].Copy()
						polNum++
					}
				} else {
					for _, polName := range metVarNames2 {
						var dataHold *sparse.DenseArray
						if polName == "PS" {
							data[polNum] = readNCF(polName, f[1], hourI)
						} else {
							data[polNum] = avgToChemLevel(readNCF(polName, f[1], hourI))
						}
						dataHold = data[polNum].Copy()
						data[polNum].AddDense(dataM1[polNum])
						data[polNum].Scale(1.0 / 2.0)
						dataM1[polNum] = dataHold.Copy()
						polNum++
					}
				}
			}
			if varLens[2] > 0 {
				for _, polName := range metVarNames3 { // A1 file here has time-avg values at 30, 90, 150 ...
					data1hr := make([]*sparse.DenseArray, 3)
					hour1hr := hourI * 3
					for j := 0; j < 3; j++ {
						data1hr[j] = readNCF(polName, f[2], hour1hr+j)
					}
					data1hr[0].AddDense(data1hr[1])
					data1hr[0].AddDense(data1hr[2])
					data1hr[0].Scale(1.0 / 3.0)
					data[polNum] = data1hr[0]
					polNum++
				}
			}
			if varLens[3] > 0 {
				for _, polName := range metVarNames4 { // A3 file has times-avg values at 90, 270, 450 ...
					data[polNum] = avgToChemLevel(readNCF(polName, f[3], hourI))
					polNum++
				}
			}
			if varLens[4] > 0 {
				for _, polName := range metVarNames5 { // A3 file has times-avg values at 90, 270, 450 ...
					// this file has flux values at stag vertical levels
					data[polNum] = StagToStagChemLevel(readNCF(polName, f[4], hourI))
					polNum++
				}
			}
			if varLens[totFiles-1] > 0 {
				// chem file has: "instantaneous species concentrations in mixing ratio (and aerosols in mole ratio)
				// at 0, 3, 6, 9, 12, 15, 18, and 21 UTC for each day."
				// gc_output.[DATE].[HOUR]0000.nc
				cfilename := strings.Replace(fileNames[totFiles-1], "[DATE]", d, -1)
				cfilename = strings.Replace(cfilename, "[HOUR]", hours[hourI], -1)

				c.ff, err = os.Open(cfilename)
				if err != nil {
					panic(err)
				}
				c.f, err = cdf.Open(c.ff)
				if err != nil {
					panic(err)
				}

				if firstHour {
					for _, polName := range chemVarNames {
						data[polNum] = readChemNCF(polName, c)
						dataM1[polNum] = data[polNum].Copy()
						polNum++
					}
				} else {
					for _, polName := range chemVarNames {
						var dataHold *sparse.DenseArray
						data[polNum] = readChemNCF(polName, c)
						dataHold = data[polNum].Copy()
						data[polNum].AddDense(dataM1[polNum])
						data[polNum].Scale(1.0 / 2.0)
						dataM1[polNum] = dataHold.Copy()
						polNum++
					}
				}
				c.ff.Close()
			}
			for _, datachanI := range datachans {
				datachanI <- data
			}
			firstHour = false
		} // end of hour loop

		for i := 0; i < totFiles-1; i++ {
			if varLens[i] > 0 {
				f[i].ff.Close()
			}
		}
	}
}

func iterateTimeSteps(msg string, fileNames []string, funcs ...cdfReaderFunc) {

	datechan := make(chan string)
	finishchan := make(chan int)
	fmt.Println("runtime.GOMAXPROCS(0)", runtime.GOMAXPROCS(0))
	for i := 0; i < runtime.GOMAXPROCS(0); i++ {
		go processFile(datechan, fileNames, finishchan, funcs...)
	}
	delta, _ := time.ParseDuration("24h")
	for now := start; now.Before(end); now = now.Add(delta) {
		d := now.Format(geosFormat)
		fmt.Println(msg + d + "...")
		//		file := strings.Replace(fileName, "[DATE]", d, -1)
		datechan <- d
	}
	close(datechan)
	for i := 0; i < runtime.GOMAXPROCS(0); i++ {
		<-finishchan
	}
}

// Calculate RMS wind speed
func calcWindSpeed(uvwIn chan []*sparse.DenseArray, outChan chan *sparse.DenseArray, varAllList []string) {
	var speed *sparse.DenseArray

	uI := getOrder("U", varAllList)
	vI := getOrder("V", varAllList)
	wwI := getOrder("OMEGA", varAllList)
	tI := getOrder("T", varAllList)
	spI := getOrder("PS", varAllList)

	firstData := true
	var dims []int
	for {
		data := <-uvwIn
		if data == nil {
			//			fmt.Println("in windSpeed, u == nil")
			for i, val := range speed.Elements {
				speed.Elements[i] = val / numTsteps
			}
			outChan <- speed
			return
		}
		u := data[uI]
		v := data[vI]
		ww := data[wwI]
		temp := data[tI]
		sp := data[spI]
		if firstData {
			// get unstaggered grid sizes
			dims = make([]int, len(u.Shape))
			for i, ulen := range u.Shape {
				vlen := v.Shape[i]
				wlen := ww.Shape[i]
				dims[i] = minInt(ulen, vlen, wlen)
				// dims[i] = minInt(ulen, vlen) // , wlen)
			}
			speed = sparse.ZerosDense(dims...)
			firstData = false
		}
		//		fmt.Println("in windSpeed, dims = ", dims)
		//		fmt.Println("u,v at 222:", u.Get(2, 2, 2), v.Get(2, 2, 2))
		w := omega2w(ww, sp, temp)
		for k := 0; k < dims[0]; k++ {
			for j := 0; j < dims[1]; j++ {
				for i := 0; i < dims[2]; i++ {
					//					fmt.Println("i,i+1 for u.Get i =", i)
					ucenter := math.Abs(u.Get(k, j, i)) // GEOS output not stag so don't need to average
					vcenter := math.Abs(v.Get(k, j, i))
					wcenter := math.Abs(w.Get(k, j, i)) // +
					//						math.Abs(w.Get(k+1, j, i))) / 2.
					s := math.Pow(math.Pow(ucenter, 2.)+
						math.Pow(vcenter, 2.)+math.Pow(wcenter, 2.), 0.5)
					// s := math.Pow(math.Pow(ucenter, 2.)+math.Pow(vcenter, 2.), 0.5)
					speed.AddVal(s, k, j, i)
				}
			}
		}
	}
	return
}

func omega2w(ww, sp, temp *sparse.DenseArray) (w *sparse.DenseArray) {
	p := sparse.ZerosDense(ww.Shape[0]+1, ww.Shape[1], ww.Shape[2])
	Dz := sparse.ZerosDense(ww.Shape...)
	w = sparse.ZerosDense(ww.Shape...)

	for k := 0; k < ww.Shape[0]; k++ {
		for j := 0; j < ww.Shape[1]; j++ {
			for i := 0; i < ww.Shape[2]; i++ {
				p.Set((sp.Get(j, i)*bp.Elements[k])+ap.Elements[k], k, j, i) // spL in units of hPa
				p.Set((sp.Get(j, i)*bp.Elements[k+1])+ap.Elements[k+1], k+1, j, i)
				t := temp.Get(k, j, i)                                                    // tL in units on K
				h := -1 * math.Log(float64(p.Get(k+1, j, i)/p.Get(k, j, i))) * rr * t / g // in meters
				Dz.Set(h, k, j, i)
				//if k == 2 && j == 2 && i == 2 {
				//	fmt.Println("p,p(k+1),Dz,t,g=", p.Get(2, 2, 2), p.Get(3, 2, 2), Dz.Get(2, 2, 2), temp.Get(2, 2, 2), g)
				//	fmt.Println("h,dz 2,2,2=", h, Dz.Get(k, j, i))
				//}
				www := ww.Get(k, j, i) /
					((p.Get(k+1, j, i) - p.Get(k, j, i)) * 100 / Dz.Get(k, j, i)) // * 100 (hPa -> Pa)
				w.Set(www, k, j, i)
			}
		}
	}
	//	fmt.Println("omega2w 2,2,2=", w.Get(2,2,2))
	//	fmt.Println("p,p(k+1),Dz,t,g=", p.Get(2,2,2), p.Get(3,2,2),Dz.Get(2,2,2),temp.Get(2,2,2),g)
	return
}

// Calculate average wind directions and speeds
// func calcWindDirection(uChan, vChan, wChan chan *sparse.DenseArray) {
// func calcWindDirection(uChan, vChan, omegaChan, spWchan, tWchan chan *sparse.DenseArray) {
func calcWindDirection(allChan chan []*sparse.DenseArray, outChan chan []*sparse.DenseArray, varAllList []string) {
	var uPlusSpeed *sparse.DenseArray
	var uMinusSpeed *sparse.DenseArray
	var vPlusSpeed *sparse.DenseArray
	var vMinusSpeed *sparse.DenseArray
	var wPlusSpeed *sparse.DenseArray
	var wMinusSpeed *sparse.DenseArray

	uI := getOrder("U", varAllList)
	vI := getOrder("V", varAllList)
	wwI := getOrder("OMEGA", varAllList)
	tI := getOrder("T", varAllList)
	spI := getOrder("PS", varAllList)

	firstData := true
	var dims []int
	for {
		//	fmt.Println("start of for loop in calcWindDir")
		data := <-allChan
		if data == nil {
			fmt.Println("in data == nil")
			fdata := make([]*sparse.DenseArray, 6)
			fdata[0] = arrayAverage(uPlusSpeed)
			fdata[1] = arrayAverage(uMinusSpeed)
			fdata[2] = arrayAverage(vPlusSpeed)
			fdata[3] = arrayAverage(vMinusSpeed)
			fdata[4] = arrayAverage(wPlusSpeed)
			fdata[5] = arrayAverage(wMinusSpeed)
			outChan <- fdata
			return
		}
		u := data[uI]
		v := data[vI]
		ww := data[wwI]
		temp := data[tI]
		sp := data[spI]
		if firstData {
			// get unstaggered grid sizes
			dims = make([]int, len(u.Shape))
			for i, ulen := range u.Shape {
				vlen := v.Shape[i]
				dims[i] = minInt(ulen, vlen)
			}
			uPlusSpeed = sparse.ZerosDense(dims...)
			uMinusSpeed = sparse.ZerosDense(dims...)
			vPlusSpeed = sparse.ZerosDense(dims...)
			vMinusSpeed = sparse.ZerosDense(dims...)
			wPlusSpeed = sparse.ZerosDense(dims...)
			wMinusSpeed = sparse.ZerosDense(dims...)
			firstData = false
		}
		//	fmt.Println("in calcWindDirection, dims =", dims)
		//	fmt.Println("u, w at 2,2,2 =", u.Get(2,2,2), ww.Get(2,2,2))
		w := omega2w(ww, sp, temp)
		//	fmt.Println("in calcWindDirection, finished omega2w")
		for k := 0; k < dims[0]; k++ {
			for j := 0; j < dims[1]; j++ {
				for i := 0; i < dims[2]; i++ {
					ucenter := u.Get(k, j, i) // + u.Get(k, j, i+1)) / 2.
					vcenter := v.Get(k, j, i) // + v.Get(k, j+1, i)) / 2.
					wcenter := w.Get(k, j, i) // + w.Get(k+1, j, i)) / 2.
					if ucenter > 0 {
						uPlusSpeed.AddVal(ucenter, k, j, i)
					} else {
						uMinusSpeed.AddVal(-ucenter, k, j, i)
					}
					if vcenter > 0 {
						vPlusSpeed.AddVal(vcenter, k, j, i)
					} else {
						vMinusSpeed.AddVal(-vcenter, k, j, i)
					}
					if wcenter > 0 {
						wPlusSpeed.AddVal(wcenter, k, j, i)
					} else {
						wMinusSpeed.AddVal(-wcenter, k, j, i)
					}
				}
			}
		}
	}
	return
}

func arrayAverage(s *sparse.DenseArray) *sparse.DenseArray {
	for i, val := range s.Elements {
		s.Elements[i] = val / numTsteps
	}
	return s
}

// Reads and parse a json configuration file.
func ReadConfigFile(filename string) {
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

	err = json.Unmarshal(bytes, &config) // this command stores config file data into the config struct created above
	if err != nil {
		fmt.Printf(
			"There has been an error parsing the configuration file.\n"+
				"Please ensure that the file is in valid JSON format\n"+
				"(you can check for errors at http://jsonlint.com/)\n"+
				"and try again!\n\n%v\n\n", err.Error())
		os.Exit(1)
	}

	err = os.MkdirAll(config.OutputDir, os.ModePerm)
	if err != nil {
		fmt.Println(err.Error())
		os.Exit(1)
	}
	return
}

func minInt(vals ...int) int {
	minval := vals[0]
	for _, val := range vals {
		if val < minval {
			minval = val
		}
	}
	return minval
}
