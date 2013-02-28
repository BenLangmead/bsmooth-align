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

#include "refmap.h"

using namespace std;

/**
 * Given a refid,offset pair in the index space, transform it into the
 * reference coordinate space according to the reference mappings
 * provided by the user.
 */
void ReferenceMap::map(U32Pair& h) const {
	if(h.first >= map_.size()) {
		cerr << "Could not find a reference-map entry for reference "
		     << h.first << " in map file \"" << fname_ << "\"" << endl;
		throw 1;
	}
	h.second += map_[h.first].second;
	h.first = map_[h.first].first;
}

/**
 * Parse a reference-map file.
 */
void ReferenceMap::parse() {
	ifstream in(fname_);
	if(!in.good() || !in.is_open()) {
		cerr << "Could not open reference map file " << fname_ << endl;
		throw 1;
	}
	while(in.peek() != EOF) {
		if(in.peek() == '>') {
			// This appears to be a name line
			in.get(); // chop off the initial '>'
			uint32_t off;
			in >> off;
			in.get(); // chop off tab
			char buf[1024];
			in.getline(buf, 1023);
			if(parseNames_) {
				if(names_.size() <= off) names_.resize(off+1);
				names_[off] = string(buf);
			}
			continue;
		}
		uint32_t id, off;
		in >> id;
		if(in.eof()) break;
		in >> off;
		map_.resize(map_.size()+1);
		map_.back().first = id;
		map_.back().second = off;
		while(isspace(in.peek())) in.get();
	}
	if(map_.empty()) {
		cerr << "Warnings: no entries in refmap file " << fname_
		     << " (" << names_.size() << " names)" << endl;
	}
	in.close();
}
