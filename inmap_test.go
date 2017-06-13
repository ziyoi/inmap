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

package inmap_test

import (
	"math"
	"testing"
	"time"

	"github.com/ctessum/geom"
	"github.com/gonum/floats"
	"github.com/spatialmodel/inmap"
	"github.com/spatialmodel/inmap/science/chem/simplechem"
	"github.com/spatialmodel/inmap/science/drydep/simpledrydep"
	"github.com/spatialmodel/inmap/science/wetdep/emepwetdep"
)

const E = 1000000. // emissions

func TestConverge(t *testing.T) {
	const (
		testTolerance = 1.e-8
		timeout       = 10 * time.Second
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

	convergences := []inmap.DomainManipulator{inmap.SteadyStateConvergenceCheck(2, cfg.PopGridColumn, nil),
		inmap.SteadyStateConvergenceCheck(-1, cfg.PopGridColumn, nil)}
	convergenceNames := []string{"fixed", "criterion"}
	expectedConcentration := []float64{0.348963630316385, 83.8101077058609}

	for i, conv := range convergences {

		iterations := 0

		d := &inmap.InMAP{
			InitFuncs: []inmap.DomainManipulator{
				cfg.RegularGrid(ctmdata, pop, popIndices, mr, emis),
				inmap.SetTimestepCFL(),
			},
			RunFuncs: []inmap.DomainManipulator{
				inmap.Calculations(inmap.AddEmissionsFlux()),
				inmap.Calculations(
					simpledrydep.DryDeposition(simplechem.SimpleDryDepIndices),
					emepwetdep.WetDeposition(simplechem.EMEPWetDepIndices),
				),
				conv,
				func(_ *inmap.InMAP) error {
					iterations++
					return nil
				},
			},
		}
		if err := d.Init(); err != nil {
			t.Error(err)
		}
		timeoutChan := time.After(timeout)
		doneChan := make(chan int)
		go func() {
			if err := d.Run(); err != nil {
				t.Error(err)
			}
			doneChan <- 0
		}()
		select {
		case <-timeoutChan:
			t.Errorf("%s timed out after %d iterations.", convergenceNames[i], iterations)
		case <-doneChan:
			t.Logf("%s completed after %d iterations.", convergenceNames[i], iterations)
		}

		o, err := inmap.NewOutputter("", false, map[string]string{"PrimPM25": "PrimaryPM25"}, nil)
		if err != nil {
			t.Error(err)
		}

		r, err := d.Results(o)
		if err != nil {
			t.Error(err)
		}
		results := r["PrimPM25"]
		total := floats.Sum(results)
		if different(total, expectedConcentration[i], testTolerance) {
			t.Errorf("%s concentration (%v) doesn't equal %v", convergenceNames[i], total, expectedConcentration[i])
		}
	}
}

func BenchmarkRun(b *testing.B) {
	const testTolerance = 1.e-8

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
		b.Error(err)
	}
	d := &inmap.InMAP{
		InitFuncs: []inmap.DomainManipulator{
			cfg.RegularGrid(ctmdata, pop, popIndices, mr, emis),
			cfg.MutateGrid(mutator, ctmdata, pop, mr, emis, nil),
			inmap.SetTimestepCFL(),
		},
		RunFuncs: []inmap.DomainManipulator{
			inmap.Calculations(inmap.AddEmissionsFlux()),
			inmap.Calculations(
				inmap.UpwindAdvection(),
				inmap.Mixing(),
				inmap.MeanderMixing(),
				simpledrydep.DryDeposition(simplechem.SimpleDryDepIndices),
				emepwetdep.WetDeposition(simplechem.EMEPWetDepIndices),
				simplechem.Chemistry(),
			),
			inmap.SteadyStateConvergenceCheck(1000, cfg.PopGridColumn, nil),
		},
	}
	if err = d.Init(); err != nil {
		b.Error(err)
	}
	if err = d.Run(); err != nil {
		b.Error(err)
	}

	o, err := inmap.NewOutputter("", false, map[string]string{"TotalPopDeaths": "coxHazard(loglogRR(TotalPM25), TotalPop, MortalityRate)"}, nil)
	if err != nil {
		b.Error(err)
	}

	r, err := d.Results(o)
	if err != nil {
		b.Error(err)
	}
	results := r["TotalPopDeaths"]
	totald := floats.Sum(results)
	const expectedDeaths = 6.614182415997713e-06

	if different(totald, expectedDeaths, testTolerance) {
		b.Errorf("Deaths (%v) doesn't equal %v", totald, expectedDeaths)
	}
}

// TestBigM2d checks whether the model can run stably with a high rate of
// convective mixing.
func TestBigM2d(t *testing.T) {
	cfg, ctmdata, pop, popIndices, mr := inmap.VarGridTestData()
	ctmdata.Data["M2d"].Data.Scale(100)
	ctmdata.Data["M2u"].Data.Scale(100)

	emis := inmap.NewEmissions()
	emis.Add(&inmap.EmisRecord{
		SOx:  E,
		NOx:  E,
		PM25: E,
		VOC:  E,
		NH3:  E,
		Geom: geom.Point{X: -3999, Y: -3999.},
	}) // ground level emissions

	d := &inmap.InMAP{
		InitFuncs: []inmap.DomainManipulator{
			cfg.RegularGrid(ctmdata, pop, popIndices, mr, emis),
			inmap.SetTimestepCFL(),
		},
		RunFuncs: []inmap.DomainManipulator{
			inmap.Calculations(inmap.AddEmissionsFlux()),
			inmap.Calculations(
				inmap.UpwindAdvection(),
				inmap.Mixing(),
				inmap.MeanderMixing(),
				simpledrydep.DryDeposition(simplechem.SimpleDryDepIndices),
				emepwetdep.WetDeposition(simplechem.EMEPWetDepIndices),
				simplechem.Chemistry(),
			),
			inmap.SteadyStateConvergenceCheck(-1, cfg.PopGridColumn, nil),
		},
	}
	if err := d.Init(); err != nil {
		t.Error(err)
	}
	if err := d.Run(); err != nil {
		t.Error(err)
	}

	o, err := inmap.NewOutputter("", false, map[string]string{"TotalPM25": "TotalPM25"}, nil)
	if err != nil {
		t.Error(err)
	}

	r, err := d.Results(o)
	if err != nil {
		t.Error(err)
	}
	results := r["TotalPM25"]
	sum := floats.Sum(results)
	if math.IsNaN(sum) {
		t.Errorf("concentration sum is NaN")
	}
}
