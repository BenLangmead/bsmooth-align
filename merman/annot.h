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

#ifndef ANNOT_H_
#define ANNOT_H_

#include <stdint.h>
#include <map>
#include <iostream>
#include <fstream>

/**
 * Encapsulates a sorted list of reference positions that are annotated
 * somehow (e.g. as a SNP).
 */
class AnnotationMap {
public:
	typedef std::pair<uint32_t, uint32_t> U32Pair;
	typedef std::pair<char, char> CharPair;
	typedef std::map<U32Pair, CharPair> AnnotMap;
	typedef std::map<U32Pair, CharPair>::const_iterator Iter;

	AnnotationMap(const char *fname) {
		fname_ = fname;
		parse();
	}

	/**
	 * Get the first entry in the map_ at a reference coordinate
	 * greater than or equal to the given coordinate.
	 */
	Iter lower_bound(const U32Pair& h) const {
		return map_.lower_bound(h);
	}

	/**
	 * Return a beginning iterator.
	 */
	Iter begin() const {
		return map_.begin();
	}

	/**
	 * Return an end iterator.
	 */
	Iter end() const {
		return map_.end();
	}

protected:

	/**
	 * Parse an annotation-map file.
	 */
	void parse();

	/// filename of file containing the annotation map
	const char *fname_;
	/// maps reference positions to character annotations
	AnnotMap map_;
};

#endif /* ANNOT_H_ */
