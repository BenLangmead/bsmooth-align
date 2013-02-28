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

#include <iostream>
#include "assert_helpers.h"
#include "ds.h"
#include "sstring.h"
#include "random_source.h"

/**
 * Simple parameterized random nucleotide generator.
 */
class RandNucGen {

public:

	RandNucGen(
		float atboost,     // over-represent A/Ts by given fraction (can be -ve)
		float aboost,      // over-represent As versus Ts by given fraction (can be -ve)
		float cboost,      // over-represent Cs versus Gs by given fraction (can be -ve)
		float singletonNf, // chance that a single N appears at any given position
		float stretchNf,   // chance that a stretch of Ns will appear at any given position
		int stretchNl,     // min number of Ns in a stretch
		int stretchNh)     // max number of Ns in a stretch
	{
		init(atboost, aboost, cboost, singletonNf, stretchNf, stretchNl, stretchNh);
	}

	/**
	 *
	 */
	void init(
		float atboost,     // over-represent A/Ts by given fraction (can be -ve)
		float aboost,      // over-represent As versus Ts by given fraction (can be -ve)
		float cboost,      // over-represent Cs versus Gs by given fraction (can be -ve)
		float singletonNf, // chance that a single N appears at any given position
		float stretchNf,   // chance that a stretch of Ns will appear at any given position
		int stretchNl,     // min number of Ns in a stretch
		int stretchNh)     // max number of Ns in a stretch
	{
		frac[0] = frac[1] = frac[2] = frac[3] = 0.25f;
		if(frac[0] + atboost <= 1.0f) atboost = 1.0f - frac[0];
		frac[0] += atboost;
		frac[3] += atboost;
		frac[1] -= atboost;
		frac[2] -= atboost;
		if(frac[0] + aboost <= 1.0f) aboost = 1.0f - frac[0];
		frac[0] += aboost;
		frac[3] -= aboost;
		if(frac[1] + cboost <= 1.0f) cboost = 1.0f - frac[1];
		frac[1] += cboost;
		frac[2] -= cboost;
		assert_lt(frac[0] + frac[1] + frac[2] + frac[3], 1.001f);
		assert_gt(frac[0] + frac[1] + frac[2] + frac[3], 0.999f);
		singletonNFrac = singletonNf;
		stretchNFrac   = stretchNf;
		stretchNLo     = stretchNl;
		stretchNHi     = stretchNh;
	}
	
	/**
	 * Install another character or stretch of characters in s.
	 */
	template<typename T>
	void nextStretch(T& s) {
		
	}

protected:

	float frac[4];
	float singletonNFrac;
	float stretchNFrac;
	int   stretchNLo;
	int   stretchNHi;
};

/**
 *
 */
class MermanRef {

public:

	MermanRef() { }
	
	void init() {
	}

protected:

	EList<SStringExpandable<uint8_t> > strs;
};
