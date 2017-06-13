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
	"testing"

	"github.com/spatialmodel/inmap"
)

func indices() (SO2, OtherGas, PM25) {
	return SO2{5}, OtherGas{3, 7, 0}, PM25{1, 2, 4, 6, 8}
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
			cfg.RegularGrid(ctmdata, pop, popIndices, mr, emis),
			cfg.MutateGrid(mutator, ctmdata, pop, mr, emis, nil),
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
