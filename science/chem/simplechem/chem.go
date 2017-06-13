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

// Package simplechem contains a simplified atmospheric chemistry mechanism.
package simplechem

import (
	"github.com/spatialmodel/inmap"
	"github.com/spatialmodel/inmap/science/drydep/simpledrydep"
	"github.com/spatialmodel/inmap/science/wetdep/emepwetdep"
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

// SimpleDryDepIndices provides array indices for use with package simpledrydep.
func SimpleDryDepIndices() (simpledrydep.SOx, simpledrydep.NH3, simpledrydep.NOx, simpledrydep.VOC, simpledrydep.PM25) {
	return simpledrydep.SOx{igS}, simpledrydep.NH3{igNH}, simpledrydep.NOx{igNO}, simpledrydep.VOC{igOrg}, simpledrydep.PM25{ipOrg, iPM2_5, ipNH, ipS, ipNO}
}

// EMEPWetDepIndices provides array indices for use with package emepwetdep.
func EMEPWetDepIndices() (emepwetdep.SO2, emepwetdep.OtherGas, emepwetdep.PM25) {
	return emepwetdep.SO2{igS}, emepwetdep.OtherGas{igNH, igNO, igOrg}, emepwetdep.PM25{ipOrg, iPM2_5, ipNH, ipS, ipNO}
}

// Chemistry returns a function that calculates the secondary formation of PM2.5.
// It explicitly calculates formation of particulate sulfate
// from gaseous and aqueous SO2.
// It partitions organic matter ("gOrg" and "pOrg"), the
// nitrogen in nitrate ("gNO and pNO"), and the nitrogen in ammonia ("gNH" and
// "pNH) between gaseous and particulate phase
// based on the spatially explicit partioning present in the baseline data.
// The function arguments represent the array indices of each chemical species.
func Chemistry() inmap.CellManipulator {
	return func(c *inmap.Cell, Δt float64) {
		// All SO4 forms particles, so sulfur particle formation is limited by the
		// SO2 -> SO4 reaction.
		ΔS := c.SO2oxidation * c.Cf[igS] * Δt
		c.Cf[ipS] += ΔS
		c.Cf[igS] -= ΔS
		// NH3 / pNH4 partitioning
		totalNH := c.Cf[igNH] + c.Cf[ipNH]
		c.Cf[ipNH] = totalNH * c.NHPartitioning
		c.Cf[igNH] = totalNH * (1 - c.NHPartitioning)

		// NOx / pN0 partitioning
		totalNO := c.Cf[igNO] + c.Cf[ipNO]
		c.Cf[ipNO] = totalNO * c.NOPartitioning
		c.Cf[igNO] = totalNO * (1 - c.NOPartitioning)

		// VOC/SOA partitioning
		totalOrg := c.Cf[igOrg] + c.Cf[ipOrg]
		c.Cf[ipOrg] = totalOrg * c.AOrgPartitioning
		c.Cf[igOrg] = totalOrg * (1 - c.AOrgPartitioning)
	}
}
