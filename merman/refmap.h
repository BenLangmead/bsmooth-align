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

#ifndef REFMAP_H_
#define REFMAP_H_

#include <stdint.h>
#include <cassert>
#include <iostream>
#include <fstream>
#include "ds.h"
#include "sstring.h"

typedef SStringExpandable<char, 256> TRefStr;

class ReferenceMap {
	typedef std::pair<uint32_t, uint32_t> U32Pair;

public:
	ReferenceMap(const char *fname, bool parseNames) {
		fname_ = fname;
		parseNames_ = parseNames;
		parse();
	}

	/**
	 * Give a reference coordinate in the index, translate it into a
	 * new reference coordinate via the reference map supplied by the
	 * user.
	 */
	void map(U32Pair& h) const;

	/**
	 * Return true iff we have a name for reference with id 'i'.
	 */
	bool hasName(size_t i) const {
		if(!parseNames_) return false;
		return !names_[i].empty();
	}

	/**
	 * Get the name for reference with id 'i'.
	 */
	const TRefStr& getName(size_t i) const {
		assert(parseNames_);
		assert(hasName(i));
		return names_[i];
	}

protected:

	/**
	 * Parse a reference-map file.
	 */
	void parse();

	const char *fname_;
	EList<U32Pair> map_;
	bool parseNames_;
	EList<TRefStr> names_;
};

#endif /* REFMAP_H_ */
