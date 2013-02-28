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

#include "annot.h"

using namespace std;

/**
 * Parse an annotation-map file.
 */
void AnnotationMap::parse() {
	ifstream in(fname_);
	if(!in.good() && in.is_open()) {
		cerr << "Could not open annotation file " << fname_ << endl;
		throw 1;
	}
	while(in.peek() != EOF) {
		U32Pair pos;
		CharPair an;
		in >> pos.first >> pos.second >> an.first >> an.second;
		map_[pos] = an;
		while(isspace(in.peek())) in.get();
	}
	in.close();
}
