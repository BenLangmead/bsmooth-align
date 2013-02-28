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

#ifndef TOKENIZE_H_
#define TOKENIZE_H_

#include <string>
#include <sstream>
#include "ds.h"

/**
 * Split string s according to given delimiters.  Mostly borrowed
 * from C++ Programming HOWTO 7.3.
 */
static inline void tokenize(const std::string& s,
                            const char *delims,
                            EList<std::string>& ss,
                            size_t max = 9999)
{
	std::string::size_type lastPos = s.find_first_not_of(delims, 0);
	std::string::size_type pos = s.find_first_of(delims, lastPos);
	while (std::string::npos != pos || std::string::npos != lastPos) {
		ss.push_back(s.substr(lastPos, pos - lastPos));
		lastPos = s.find_first_not_of(delims, pos);
		pos = s.find_first_of(delims, lastPos);
		if(ss.size() == (max - 1)) {
			pos = std::string::npos;
		}
	}
}

#endif /*TOKENIZE_H_*/
