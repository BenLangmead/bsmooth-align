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

#ifndef SEEDS_H_
#define SEEDS_H_

#include <stdint.h>
#include "ds.h"

/**
 * Encapsulates a spaced-seed strategy.
 */
class SeedSet {
public:
	virtual ~SeedSet() { }
	/**
	 * Return the number of coincident seed hits that are necessary for
	 * an i-mismatch hit.
	 */
	virtual int thresh(int i) const = 0;

	/**
	 * Given i, return true iff i-mismatch hits are guaranteed to be
	 * found.
	 */
	virtual bool guaranteed(int i) const = 0;

	/**
	 * Populate a list with the appropriate seeds.
	 */
	virtual void getSeeds(EList<uint64_t>& s, int upto, int rounds) const = 0;

	/**
	 * Create a set of seeds that adheres to the user's requested
	 * specificity or, if specified, their N-choose-K preference.
	 */
	static SeedSet* create(int i, int specificty, int chooseN, int chooseK);
};

/**
 * Concrete subclass of SeedSet that calculates thresholds and seeds
 * using a simple n-choose-k spaced-seed approach.
 */
class SeedSetChoose : public SeedSet {

public:

	SeedSetChoose(size_t x, size_t y);
	
	virtual ~SeedSetChoose() { }

	/**
	 * Return the number of coincident seed hits that are necessary for
	 * an i-mismatch hit.
	 */
	virtual int thresh(int i) const {
		if(i >= (int)threshs_.size()) return 1;
		else return threshs_[i];
	}
	
	/**
	 * Which nu
	 */
	virtual bool guaranteed(int i) const {
		return i < (int)threshs_.size();
	}
	
	/**
	 * Populate a list with the appropriate seeds.
	 */
	virtual void getSeeds(EList<uint64_t>& s, int upto, int rounds) const;
	
	/**
	 * Sanity check this seeding scheme.
	 */
	bool repOk() const {
		for(size_t i = 1; i < threshs_.size(); i++) {
			assert_leq(threshs_[i], threshs_[i-1]);
		}
		return true;
	}

protected:

	/**
	 * Calculate n choose k.
	 */
	unsigned long long choose(size_t n, size_t k) {
		if (k > n)   return 0;		
		if (k > n/2) k = n-k; // Take advantage of symmetry
		long double accum = 1;
		for(size_t i = 1; i <= k; i++) {
			accum = accum * (n-k+i) / i;
		}
		return static_cast<unsigned long long>(accum + 0.5); // avoid rounding error
	}

	size_t x_, y_;
	EList<uint64_t> masks_;
	EList<int> threshs_;
};

#endif /*SEEDS_H_*/
