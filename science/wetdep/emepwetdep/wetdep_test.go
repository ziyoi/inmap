/*
Copyright © 2017 the InMAP authors.
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

package emepwetdep

import (
	"fmt"
	"testing"

	"github.com/spatialmodel/inmap"
)

func indices() (SO2, OtherGas, PM25) {
	return SO2{5}, OtherGas{3, 7, 0}, PM25{1, 2, 4, 6, 8}
}

// physical constants
const (
	// Molar masses [grams per mole]
	mwNOx = 46.0055
	mwN   = 14.0067 // g/mol, molar mass of nitrogen
	mwNO3 = 62.00501
	mwNH3 = 17.03056
	mwNH4 = 18.03851
	mwS   = 32.0655 // g/mol, molar mass of sulfur
	mwSO2 = 64.0644
	mwSO4 = 96.0632

	// Chemical mass conversions [ratios]
	NOxToN = mwN / mwNOx
	NtoNO3 = mwNO3 / mwN
	SOxToS = mwSO2 / mwS
	StoSO4 = mwS / mwSO4
	NH3ToN = mwN / mwNH3
	NtoNH4 = mwNH4 / mwN
)

// Indicies of individual pollutants in arrays.
const (
	igOrg int = iota
	ipOrg
	iPM2_5
	igNH
	ipNH
	igS
	ipS
	igNO
	ipNO
)

// emisConv lists the accepted names for emissions species, the array
// indices they correspond to, and the
// factors needed to convert [μg/s] of emitted species to [μg/s] of
// model species.
var emisConv = map[string]struct {
	i    int
	conv float64
}{
	"VOC":   {i: igOrg, conv: 1},
	"NOx":   {i: igNO, conv: NOxToN},
	"NH3":   {i: igNH, conv: NH3ToN},
	"SOx":   {i: igS, conv: SOxToS},
	"PM2_5": {i: iPM2_5, conv: 1},
}

// AddEmisFlux adds emissions flux to Cell c based on the given
// pollutant name and amount in units of μg/s. The units of
// the resulting flux are μg/m3/s.
func AddEmisFlux(c *inmap.Cell, name string, val float64) error {
	fluxScale := 1. / c.Dx / c.Dy / c.Dz // μg/s /m/m/m = μg/m3/s
	conv, ok := emisConv[name]
	if !ok {
		return fmt.Errorf("simplechem: '%s' is not a valid emissions species; valid options are VOC, NOx, NH3, SOx, and PM2_5", name)
	}
	c.EmisFlux[conv.i] += val * conv.conv * fluxScale
	return nil
}

func TestWetDeposition(t *testing.T) {
	cfg, ctmdata, pop, popIndices, mr := inmap.VarGridTestData()
	emis := inmap.NewEmissions()

	mutator, err := inmap.PopulationMutator(cfg, popIndices)
	if err != nil {
		t.Error(err)
	}
	d := &inmap.InMAP{
		InitFuncs: []inmap.DomainManipulator{
			cfg.RegularGrid(ctmdata, pop, popIndices, mr, emis, AddEmisFlux),
			cfg.MutateGrid(mutator, ctmdata, pop, mr, emis, AddEmisFlux, nil),
			inmap.SetTimestepCFL(),
		},
		RunFuncs: []inmap.DomainManipulator{
			inmap.Calculations(inmap.AddEmissionsFlux()),
			inmap.Calculations(WetDeposition(indices)),
			inmap.SteadyStateConvergenceCheck(1, cfg.PopGridColumn, nil),
		},
	}
	if err := d.Init(); err != nil {
		t.Error(err)
	}
	for _, c := range d.Cells() {
		for i := range c.Ci {
			c.Cf[i] = 1 // set concentrations to 1
		}
	}
	if err := d.Run(); err != nil {
		t.Error(err)
	}

	for _, c := range d.Cells() {
		for ii, cc := range c.Cf {
			if cc > 1 || cc <= 0.99 {
				t.Errorf("ground-level cell %v pollutant %d should equal be between 0.99 and 1 but is %g", c, ii, cc)
			}
		}
	}
}
