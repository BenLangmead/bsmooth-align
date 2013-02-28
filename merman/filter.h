/*
 * Copyright 2011, Ben Langmead <blangmea@jhsph.edu>
 *
 * This file is part of Merman.
 *
 * Merman is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Merman is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Merman.  If not, see <http://www.gnu.org/licenses/>.
 */

enum {
	// Read was too short for the index being used
	FILTER_TOO_SHORT_FOR_INDEX         = 0x01,

	// Read was shorter than the minimum length allowed by alignment parameters
	// directly governing alignment length (e.g. Merman's --minlen)
	FILTER_TOO_SHORT_FOR_MINLEN_PARAMS = 0x02,
	
	// Read was shorter than the minimum length allowed by alignment parameters
	// for scoring thresholds.  E.g. if the minimum score is 10 and the match
	// bonus is 1, a 9-nt read can't possibly align
	FILTER_TOO_SHORT_FOR_SCORE_PARAMs  = 0x04,
	
	// An upstream QC step flagged the read as failing QC
	FILTER_QC                          = 0x08,
	
	// Too many characters in the read are non-A/C/G/T
	FILTER_TOO_MANY_AMBIGS             = 0x10,
	
	// Too many low quality positions and/or homopolymers
	FILTER_BAD_QUALITIES               = 0x20
};
