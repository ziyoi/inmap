/*
Copyright © 2013 the InMAP authors.
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

import (
	"math"
	"testing"
)

const E = 1000000. // emissions

// Tests whether the cells correctly reference each other
func TestCellAlignment(t *testing.T) {
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
	d.testCellAlignment2(t)
}

func (d *InMAP) testCellAlignment2(t *testing.T) {
	const testTolerance = 1.e-8
	for _, cell := range *d.cells {
		var westCoverage, eastCoverage, northCoverage, southCoverage float64
		var aboveCoverage, belowCoverage, groundLevelCoverage float64
		for _, w := range *cell.west {
			westCoverage += w.info.coverFrac
			if !w.boundary {
				pass := false
				for _, e := range *w.east {
					if e.Cell == cell.Cell {
						pass = true
						if different(w.info.diff, e.info.diff, testTolerance) {
							t.Errorf("Kxx doesn't match")
						}
						if different(w.info.centerDistance, e.info.centerDistance, testTolerance) {
							t.Errorf("Dx doesn't match")
							break
						}
					}
				}
				if !pass {
					t.Errorf("Failed for Cell %v West", cell)
				}
			}
		}
		for _, e := range *cell.east {
			eastCoverage += e.info.coverFrac
			if !e.boundary {
				pass := false
				for _, w := range *e.west {
					if w.Cell == cell.Cell {
						pass = true
						if different(e.info.diff, w.info.diff, testTolerance) {
							t.Errorf("Kxx doesn't match")
						}
						if different(e.info.centerDistance, w.info.centerDistance, testTolerance) {
							t.Errorf("Dx doesn't match")
						}
						break
					}
				}
				if !pass {
					t.Errorf("Failed for Cell %v East", cell)
				}
			}
		}
		for _, n := range *cell.north {
			northCoverage += n.info.coverFrac
			if !n.boundary {
				pass := false
				for _, s := range *n.south {
					if s.Cell == cell.Cell {
						pass = true
						if different(n.info.diff, s.info.diff, testTolerance) {
							t.Errorf("Kyy doesn't match")
						}
						if different(n.info.centerDistance, s.info.centerDistance, testTolerance) {
							t.Errorf("Dy doesn't match")
						}
						break
					}
				}
				if !pass {
					t.Errorf("Failed for Cell %v  North", cell)
				}
			}
		}
		for _, s := range *cell.south {
			southCoverage += s.info.coverFrac
			if !s.boundary {
				pass := false
				for _, n := range *s.north {
					if n.Cell == cell.Cell {
						pass = true
						if different(s.info.diff, n.info.diff, testTolerance) {
							t.Errorf("Kyy doesn't match")
						}
						if different(s.info.centerDistance, n.info.centerDistance, testTolerance) {
							t.Errorf("Dy doesn't match")
						}
						break
					}
				}
				if !pass {
					t.Errorf("Failed for Cell %v South", cell)
				}
			}
		}
		for _, a := range *cell.above {
			aboveCoverage += a.info.coverFrac
			if !a.boundary {
				pass := false
				for _, b := range *a.below {
					if b.Cell == cell.Cell {
						pass = true
						if different(a.info.diff, b.info.diff, testTolerance) {
							t.Errorf("Kzz doesn't match above (layer=%v, "+
								"KzzAbove=%v, KzzBelow=%v)", cell.Layer,
								b.info.diff, a.info.diff)
						}
						if different(a.info.centerDistance, b.info.centerDistance, testTolerance) {
							t.Errorf("Dz doesn't match")
						}
						break
					}
				}
				if !pass {
					t.Errorf("Failed for Cell %v Above", cell)
				}
			}
		}
		for _, b := range *cell.below {
			belowCoverage += b.info.coverFrac
			pass := false
			if cell.Layer == 0 && b.Cell == cell.Cell {
				pass = true
			} else {
				for _, a := range *b.above {
					if a.Cell == cell.Cell {
						pass = true
						if different(b.info.diff, a.info.diff, testTolerance) {
							t.Errorf("Kzz doesn't match below")
						}
						if different(b.info.centerDistance, a.info.centerDistance, testTolerance) {
							t.Errorf("Dz doesn't match")
						}
						break
					}
				}
			}
			if !pass {
				t.Errorf("Failed for Cell %v  Below", cell)
			}
		}
		// Assume upper cells are never higher resolution than lower cells
		for _, g := range *cell.groundLevel {
			groundLevelCoverage += g.info.coverFrac
			g2 := g
			pass := false
			for {
				if g2.above.len() == 0 {
					pass = false
					break
				}
				if g2.Cell == (*g2.above)[0].Cell {
					pass = false
					break
				}
				if g2.Cell == cell.Cell {
					pass = true
					break
				}
				g2 = (*g2.above)[0]
			}
			if !pass {
				t.Errorf("Failed for Cell %v GroundLevel", cell)
			}
		}
		const tolerance = 1.0e-10
		if different(westCoverage, 1, tolerance) {
			t.Errorf("cell %v, west coverage %g!=1", cell, westCoverage)
		}
		if different(eastCoverage, 1, tolerance) {
			t.Errorf("cell %v, east coverage %g!=1", cell, eastCoverage)
		}
		if different(southCoverage, 1, tolerance) {
			t.Errorf("cell %v, south coverage %g!=1", cell, southCoverage)
		}
		if different(northCoverage, 1, tolerance) {
			t.Errorf("cell %v, north coverage %g!=1", cell, northCoverage)
		}
		if different(belowCoverage, 1, tolerance) {
			t.Errorf("cell %v, below coverage %g!=1", cell, belowCoverage)
		}
		if different(aboveCoverage, 1, tolerance) {
			t.Errorf("cell %v, above coverage %g!=1", cell, aboveCoverage)
		}
		if different(groundLevelCoverage, 1, tolerance) {
			t.Errorf("cell %v, groundLevel coverage %g!=1", cell, groundLevelCoverage)
		}
	}
}

func different(a, b, tolerance float64) bool {
	if 2*math.Abs(a-b)/math.Abs(a+b) > tolerance || math.IsNaN(a) || math.IsNaN(b) {
		return true
	}
	return false
}

func absDifferent(a, b, tolerance float64) bool {
	if math.Abs(a-b) > tolerance {
		return true
	}
	return false
}
