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

#ifndef MER_ENT_H_
#define MER_ENT_H_

#include <cassert>
#include <iostream>
#include <stdint.h>
#include "assert_helpers.h"

/**
 * An entry in a simple seed-based mer index.  Each entry is 10 bytes,
 * but the struct is declared as packed (this is GNU-specific) to
 * prevent it from being padded.
 */
struct mer_ent {

	mer_ent() {
		key = value = key2 = seed = 0;
	}

	mer_ent(const mer_ent& me) {
		key = me.key;
		value = me.value;
		key2 = me.key2;
		seed = me.seed;
	}

	mer_ent& operator=(const mer_ent& me) {
		key = me.key;
		value = me.value;
		key2 = me.key2;
		seed = me.seed;
		return (*this);
	}

	/**
	 * Return true iff this mer_ent is less than rhs.  Sort first by
	 * seed, then by key2 (most recently extracted), then by key, then
	 * by value
	 */
	int operator< (const mer_ent &rhs) const {
		if(seed < rhs.seed) return 1;
		if(seed > rhs.seed) return 0;
		if(key2 < rhs.key2) return 1;
		if(key2 > rhs.key2) return 0;
		if(key < rhs.key) return 1;
		if(key > rhs.key) return 0;
		return (value < rhs.value)? 1 : 0;
	}

	int operator== (const mer_ent &rhs) const {
		return(key == rhs.key &&
		       key2 == rhs.key2 &&
		       seed == rhs.seed &&
		       value == rhs.value);
	}

	mer_ent(uint64_t k, uint32_t v, int s) {
		init(k, v, s);
	}

	/// Initialize with 64-bit k, v, and s
	void init(uint64_t k, uint32_t v, int s) {
		assert_eq((k >> 40llu), 0llu);
		assert_geq(s, 0);
		key = (uint32_t)k;
		key2 = (k >> 32llu); // most recently extracted bits go in key2
		assert_eq(getKey(), k);
		value = v;
		seed = s;
	}

	uint64_t getKey() const {
		return ((uint64_t)key2 << 32llu) | (uint64_t)key;
	}

	uint32_t key;   // extracted subsequence, first 16 chars
	uint32_t value; // offset into reference
	uint8_t  key2;  // extracted subsequence, next 4 chars
	uint8_t  seed     : 8; // which seed?

} __attribute__((__packed__));

/**
 * An entry in a simple seed-based mer index.  Each entry is 10 bytes,
 * but the struct is declared as packed (this is GNU-specific) to
 * prevent it from being padded.
 */
struct mer_ent_sm {

	mer_ent_sm() {
		key = value = key2 = 0;
	}

	mer_ent_sm& operator=(const mer_ent& me) {
		key = (me.key >> 16);
		value = me.value;
		key2 = me.key2;
		// Ignore seed and bottom 16 bits of key; they're now encoded
		// in the merJump_ index.
		return (*this);
	}

	/**
	 * Return true iff this mer_ent is less than rhs.  Sort first by
	 * seed, then by key2 (most recently extracted), then by key, then
	 * by value
	 */
	int operator< (const mer_ent_sm &rhs) const {
		if(key < rhs.key) return 1;
		if(key > rhs.key) return 0;
		if(key2 < rhs.key2) return 1;
		if(key2 > rhs.key2) return 0;
		return (value < rhs.value)? 1 : 0;
	}

	int operator== (const mer_ent_sm &rhs) const {
		return(key == rhs.key &&
		       key2 == rhs.key2 &&
		       value == rhs.value);
	}

	uint32_t getKey() const {
		return key | ((uint32_t)key2 << 16);
	}

	uint16_t key;   // extracted subsequence, first 16 chars
	uint32_t value; // offset into reference
	uint8_t  key2;  // extracted subsequence, next 4 chars

} __attribute__((__packed__));

#endif /* MER_ENT_H_ */
