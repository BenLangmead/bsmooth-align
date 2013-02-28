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
#include "alphabet.h"
#include "hit_set.h"
#include "annot.h"
#include "refmap.h"
#include "mer_index.h"
#include "threading.h"
#include "ref.h"
#include "orient.h"

using namespace std;

/**
 * Report up to 'khits' hits from this HitSet.
 */
void HitSet::reportUpTo(
	const Read& rd,
	AlignOutput& os,
	int khits,
	const ReferenceSet& refs,
	const ReferenceMap* rmap, // can be NULL
	bool xformCoords,
	const AnnotationMap *amap)
{
	khits = min(khits, (int)size());
	const BTDnaString& seq = rd.seq;
	const BTString& qual = rd.qual;
	BTDnaString seqrc;
	BTString qualr;
	bool seqrcInit = false;
	assert(!seq.empty());
	assert_eq(seq.length(), qual.length());
	assert(seq.isPrintable());
	assert(qual.isPrintable());
	assert(rd.name.isPrintable());
	for(int i = 0; i < khits; i++) {
		HitSetEnt& h = ents[i];
		bool printfw = orientPrintFw(h.orient);
		if(!orientPrintFw(h.orient) && !seqrcInit) {
			// Lazily initialize seqrc and qualr
			size_t seqlen = seq.length();
			seqrc.resize(seqlen);
			qualr.resize(seqlen);
			for(size_t j = 0; j < seqlen; j++) {
				seqrc.set(rd.color ? seq[j] : comp(seq[j]), seqlen - j - 1);
				qualr.set(qual[j], seqlen - j - 1);
			}
			seqrcInit = true;
		}
		// If necessary, transform the reference coordinates into the
		// user-specified coordinate system
		if(rmap != NULL && xformCoords) rmap->map(h.h);
		if(amap != NULL) {
			AnnotationMap::Iter ai = amap->lower_bound(h.h);
			bool added = false;
			const size_t len = seq.length();
			for(; ai != amap->end(); ai++) {
				assert_geq(ai->first.first, h.h.first);
				if(ai->first.first != h.h.first) {
					// Different chromosome
					break;
				}
				if(ai->first.second >= h.h.second + len) {
					// Doesn't fall into alignment
					break;
				}
				if(ai->second.first != 'S') {
					// Not a SNP annotation
					continue;
				}
				size_t off = ai->first.second - h.h.second;
				// snpc <- reference char from SNP annotation
				char snpc = toupper(ai->second.second);
				// c <- query character at given offset
				char c;
				if(!printfw) {
					c = toupper(seqrc[off]);
					off = len - off - 1;
				} else {
					c = toupper(seq[off]);
				}
				bool mm = false;
				for(size_t j = 0; j < h.size(); j++) {
					if(h[j].pos == off) {
						assert_eq(EDIT_TYPE_MM, h[j].type);
						// Force the mismatch's reference base to match
						// the true reference base - the previous
						// Bowtie phase may introduce mistakes here
						h[j].chr = snpc;
						mm = true;
						break;
					}
				}
				if(mm) continue;
				if(isDna(c)) {
					// If c matches the reference char, then this isn't an
					// edit
					if(c != snpc) {
						// Extract original reference character from annotation
						h.expand();
						h.back().pos = (uint32_t)off;
						h.back().type = EDIT_TYPE_SNP;
						h.back().chr = snpc;
						added = true;
					}
				}
			}
		}
		h.sort();
		// Set up refNamePtr
		char refName[1024];
		const char* refNamePtr = refName;
		if(rmap != NULL) {
			refNamePtr = rmap->getName(h.h.first).toZBuf();
		} else {
			refNamePtr = refs.getName(h.h.first).toZBuf();
		}
		// Set up qchr fields of Edits
		for(size_t j = 0; j < h.edits.size(); j++) {
			Edit& e = h.edits[j];
			if(e.qchr == 0)
				e.qchr = (printfw ? seq[e.pos] : comp(seq[e.pos]));
		}
		// TODO: handle 5' clipping
		os.printAlignment(
			rd,
			(h.seq.length() > 0 ? h.seq : (printfw ? seq : seqrc)),
			(printfw ? qual : qualr),
			true,
			h.clip5,
			h.clip3,
			h.h.first,
			refNamePtr,
			h.h.second,
			h.orient,
			h.edits,
			h.aedits,
			h.cedits,
			h.ccedits,
			h.oms,
			h.mapq,
			h.stratum,
			h.cost,
			h.ccost);
	}
}

ostream& operator << (ostream& os, const HitSetEnt& hs) {
	os << "\t" << hs.h.first << ":" << hs.h.second;
	return os;
}

//ostream& operator << (ostream& os, const HitSet& hs) {
//	os << hs.name << ":" << hs.seq << ":" << hs.qual << endl;
//	for(size_t i = 0; i < hs.ents.size(); i++) os << hs.ents[i];
//	return os;
//}
