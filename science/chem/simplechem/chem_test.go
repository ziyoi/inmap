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

package simplechem

import (
	"math"
	"testing"

	"github.com/ctessum/geom"
	"github.com/spatialmodel/inmap"
)

const E = 1000000. // emissions

// Test whether mass is conserved during chemical reactions.
func TestChemistry(t *testing.T) {
	const (
		testTolerance = 1.e-8
	)
	cfg, ctmdata, pop, popIndices, mr := inmap.VarGridTestData()
	emis := inmap.NewEmissions()
	emis.Add(&inmap.EmisRecord{
		SOx:  E,
		NOx:  E,
		PM25: E,
		VOC:  E,
		NH3:  E,
		Geom: geom.Point{X: -3999, Y: -3999.},
	}) // ground level emissions

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
			inmap.Calculations(Chemistry()),
			inmap.SteadyStateConvergenceCheck(1, cfg.PopGridColumn, nil),
		},
	}
	if err := d.Init(); err != nil {
		t.Error(err)
	}
	if err := d.Run(); err != nil {
		t.Error(err)
	}

	c := d.Cells()[0]
	sum := 0.
	sum += c.Cf[igOrg] + c.Cf[ipOrg]
	sum += (c.Cf[igNO] + c.Cf[ipNO]) / inmap.NOxToN
	sum += (c.Cf[igNH] + c.Cf[ipNH]) / inmap.NH3ToN
	sum += (c.Cf[igS] + c.Cf[ipS]) / inmap.SOxToS
	sum += c.Cf[iPM2_5]
	sum *= c.Volume

	if c.Cf[ipOrg] == 0 || c.Cf[ipS] == 0 || c.Cf[ipNH] == 0 || c.Cf[ipNO] == 0 {
		t.Error("chemistry appears not to have occured")
	}
	if different(sum, 5*E*d.Dt, testTolerance) {
		t.Error("different")
	}
}

func different(a, b, tolerance float64) bool {
	if 2*math.Abs(a-b)/math.Abs(a+b) > tolerance || math.IsNaN(a) || math.IsNaN(b) {
		return true
	}
	return false
}
