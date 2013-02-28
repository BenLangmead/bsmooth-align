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
#include "ref.h"
#include "edit.h"
#include "hit_set.h"
#include "check.h"
#include "orient.h"
#include "sstring.h"

using namespace std;

/**
 * Sanity check a hit by seeing if the covered nucleotides agree with
 * the read nucleotides and edits.  The read sequence and the edit
 * characters have already been converted
 *
 * The trimmed, decoded read/quality sequence is in h.seq/h.qual.
 * Trimming has already been applied, and h.h.second has already been
 * updated accordingly.
 */
bool sanityCheckHit(const ReferenceSet& refs, const HitSetEnt& h) {
	int ei = 0;
	size_t nlen = h.seq.length();
	assert_eq(refs[h.h.first].fwWatsonIdx, h.h.first); // must be forward reference strand
	EList<Edit> edits = h.edits;
	if(!orientPrintFw(h.orient)) {
		Edit::invertPoss(edits, nlen);
	}
	// Get the reference sequence aligned to, from the strand aligned to
	SStringExpandable<char, 256> rf;
	size_t rfoff = h.horig.second;
	if(!is5pLeft(h.orient)) {
		assert_geq(rfoff, h.clip3);
	} else {
		assert_geq(rfoff, h.clip5);
	}
	for(size_t i = 0; i < nlen; i++) {
		rf.append(refs[h.horig.first].seq.charAt(rfoff + i, false));
	}
	// If reference sequence was crick, we reverse complement it before
	// comparing it to h.seq.
	if(refs[h.horig.first].crick) {
		rf.reverse();
		// Reverse-complement reference, including reverse
		// complementing the IUPAC codes.
		for(size_t i = 0; i < nlen; i++) {
			rf.set(asc2dnacomp[(int)rf[i]], i);
		}
	}
	// Clipping has already been applied to h.seq/h.qual and the
	// reference offset has been adjusted accordingly.  We shouldn't
	// need to do any adjusting here.
	int clipLeft  = 0; //is5pLeft(h.orient) ? h.clip5 : h.clip3;
	int clipRight = 0; //is5pLeft(h.orient) ? h.clip3 : h.clip5;
	// Move left-to-right along forward reference strand
	for(size_t i = 0; i < nlen; i++) {
		char readc, refc;
		// refoff <- offset of leftmost character on forward reference strand involved in the alignment
		refc = toupper(rf[i]);
		readc = toupper(h.seq.toChar(i+clipLeft));
		bool mm = false;
		bool err = false;
		while(ei < (int)edits.size() && edits[ei].pos < i) ei++;
		while(ei < (int)edits.size() && edits[ei].pos == i) {
			if(edits[ei].isMismatch()) {
				if(readc != (char)edits[ei].qchr) {
					cerr << "edits[" << ei << "].qchr is " << (char)edits[ei].qchr
					<< " but toupper(h.seq.toChar(" << i+clipLeft << ") == " << readc
					<< endl;
					err = true;
				}
				if(refc != (char)edits[ei].chr) {
					cerr << "edits[" << ei << "].chr is " << (char)edits[ei].chr
						 << " but toupper(refs[" << h.h.first << "].seq.charAt("
						 << h.h.second << " + " << i+clipLeft << ", false) == " << refc
						 << endl;
					err = true;
				}
				refc = readc;
				mm = true;
			}
			ei++;
		}
		bool relevant = (i >= (size_t)clipLeft) && (i < (nlen - (size_t)clipRight));
		if(relevant && (err || !ambigCompat(refc, readc, false))) {
			// Print friendly side-by-side comparison before dying
			string rds, rfs, mms;
			ei = 0;
			for(size_t j = 0; j < nlen; j++) {
				size_t refoff = h.h.second;
				int rfc = toupper(refs[h.h.first].seq.charAt(refoff + j, false));
				int rdc = toupper(h.seq.toChar(j+clipLeft));
				mm = false;
				while(ei < (int)edits.size() && edits[ei].pos < j) ei++;
				while(ei < (int)edits.size() && edits[ei].pos == j) {
					if(edits[ei].isMismatch()) mm = true;
					ei++;
				}
				bool relevantIn =
					(j >= (size_t)clipLeft) &&
					(j  < (nlen - (size_t)clipRight));
				rds += rdc;
				rfs += rfc;
				if(relevantIn) {
					if(mm) mms += " ";
					else   mms += "|";
				} else {
					mms += "C";
					continue;
				}
			}
			cerr << "Read: " << rds << endl;
			cerr << "      " << mms << endl;
			cerr << " Ref: " << rfs << endl;
			assert(false);
		}
	}
	return true;
}
