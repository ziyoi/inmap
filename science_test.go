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

package inmap

import "testing"

// Test whether convective mixing coefficients are balanced in
// a way that conserves mass
func TestConvectiveMixing(t *testing.T) {
	const testTolerance = 1.e-8

	cfg, ctmdata, pop, popIndices, mr := VarGridTestData()
	emis := NewEmissions()

	mutator, err := PopulationMutator(cfg, popIndices)
	if err != nil {
		t.Error(err)
	}
	d := &InMAP{
		InitFuncs: []DomainManipulator{
			cfg.RegularGrid(ctmdata, pop, popIndices, mr, emis),
			cfg.MutateGrid(mutator, ctmdata, pop, mr, emis, nil),
		},
	}
	if err := d.Init(); err != nil {
		t.Error(err)
	}

	for _, c := range d.Cells() {
		val := c.M2u - c.M2d + (*c.above)[0].M2d*(*c.above)[0].Dz/c.Dz
		if absDifferent(val, 0, testTolerance) {
			t.Error(c.Layer, val, c.M2u, c.M2d, (*c.above)[0].M2d)
		}
	}
}
