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
#include "seeds.h"

using namespace std;

SeedSet* SeedSet::create(int i, int specificity, int chooseN, int chooseK) {
	if(chooseN > 0 && chooseK > 0) {
		if(chooseN < chooseK) {
			cerr << "Error: When specifying an n-choose-k seed strategy, N must be >= K." << endl;
			cerr << "N was " << chooseN << " and K was " << chooseK << "." << endl;
			throw 1;
		}
		if(chooseN == chooseK) {
			chooseN = chooseK = 1;
		}
		return new SeedSetChoose(chooseN, chooseK);
	}
	if(i == 0) {
		return new SeedSetChoose(1, 1);
	}
	else if(i == 1) {
		return new SeedSetChoose(2, 1);
	}
	else if(i == 2) {
		if(specificity < 1) {
			return new SeedSetChoose(4, 2);
		} else if(specificity == 1) {
			return new SeedSetChoose(5, 3);
		} else if(specificity == 2) {
			return new SeedSetChoose(6, 4);
		}
	}
	else if(i == 3) {
		if(specificity < 1) {
			return new SeedSetChoose(4, 1);
		} else if(specificity == 1) {
			return new SeedSetChoose(5, 2);
		} else if(specificity == 2) {
			return new SeedSetChoose(6, 3);
		}
	}
	else if(i == 4) {
		if(specificity < 1) {
			return new SeedSetChoose(5, 1);
		} else if(specificity == 1) {
			return new SeedSetChoose(6, 2);
		} else if(specificity == 2) {
			return new SeedSetChoose(7, 3);
		}
	}
	else if(i == 5) {
		if(specificity < 1) {
			return new SeedSetChoose(6, 1);
		} else if(specificity == 1) {
			return new SeedSetChoose(7, 2);
		} else if(specificity == 2) {
			return new SeedSetChoose(8, 3);
		}
	}
	else {
		// Not fully sensitive
		return new SeedSetChoose(7, 2);
	}
	return NULL;
}

/**
 * Shuffle columns of the spaced-seed matrix randomly.  Make sure not
 * to involve columns outside of the index window, or we could
 * compromise the sensitivity inside the window.
 */
static void randomizeSeeds(EList<uint64_t>& s, int upto, int rounds) {
	for(int i = 0; i < rounds; i++) {
		int r1 = rand() % upto;
		int r2;
		do { r2 = rand() % upto; } while(r2 == r1);
		for(size_t j = 0; j < s.size(); j++) {
			if(((s[j] >> r1) & 1llu) != ((s[j] >> r2) & 1llu)) {
				s[j] ^= ((1llu << r1) | (1llu << r2));
			}
		}
	}
}

/**
 * x choose y
 */
SeedSetChoose::SeedSetChoose(size_t x, size_t y) : x_(x), y_(y) {
	assert_gt(x, 0);
	assert_gt(y, 0);
	assert_leq(x, 10);
	assert_leq(y, x);
	// Generate all possible masks for x choose y
	for(size_t i1  = 0;    i1  < (y < 1  ? 1    : x); i1++ ) {
	for(size_t i2  = i1+1; i2  < (y < 2  ? i1+2 : x); i2++ ) {
	for(size_t i3  = i2+1; i3  < (y < 3  ? i2+2 : x); i3++ ) {
	for(size_t i4  = i3+1; i4  < (y < 4  ? i3+2 : x); i4++ ) {
	for(size_t i5  = i4+1; i5  < (y < 5  ? i4+2 : x); i5++ ) {
	for(size_t i6  = i5+1; i6  < (y < 6  ? i5+2 : x); i6++ ) {
	for(size_t i7  = i6+1; i7  < (y < 7  ? i6+2 : x); i7++ ) {
	for(size_t i8  = i7+1; i8  < (y < 8  ? i7+2 : x); i8++ ) {
	for(size_t i9  = i8+1; i9  < (y < 9  ? i8+2 : x); i9++ ) {
	for(size_t i10 = i9+1; i10 < (y < 10 ? i9+2 : x); i10++) {
		uint64_t pat = 0;
		if(y >= 10) pat |= 1llu << i10;
		if(y >= 9)  pat |= 1llu << i9;
		if(y >= 8)  pat |= 1llu << i8;
		if(y >= 7)  pat |= 1llu << i7;
		if(y >= 6)  pat |= 1llu << i6;
		if(y >= 5)  pat |= 1llu << i5;
		if(y >= 4)  pat |= 1llu << i4;
		if(y >= 3)  pat |= 1llu << i3;
		if(y >= 2)  pat |= 1llu << i2;
		if(y >= 1)  pat |= 1llu << i1;
		size_t fill = x;
		uint64_t final = pat;
		while(fill < 64) {
			final <<= x;
			final |= pat;
			fill += x;
		}
		masks_.push_back(final);
	}}}}}}}}}}
	// Confirm we made the right number of masks
	assert_eq(choose(x, y), masks_.size());
	assert_gt(masks_.size(), 0);
#ifndef NDEBUG
	for(size_t offset = 0; offset < 64-x; offset += x) {
		for(size_t i = 0; i < masks_.size(); i++) {
			size_t set = 0;
			for(size_t k = 0; k < x; k++) {
				if((((masks_[i] >> offset) >> k) & 1) != 0) set++;
			}
			assert_eq(y, set);		
			for(size_t j = i+1; j < masks_.size(); j++) {
				assert_neq((masks_[i] >> offset), masks_[j]);
			}
		}
	}
#endif
	// Record the 
	while(x >= y) {
		threshs_.push_back((int)choose(x--, y));
		assert_gt(threshs_.back(), 0);
	}
	assert(repOk());
}

/**
 * Populate a list with the appropriate seeds.
 */
void SeedSetChoose::getSeeds(EList<uint64_t>& s, int upto, int rounds) const {
	s = masks_;
	randomizeSeeds(s, upto, rounds);
}
