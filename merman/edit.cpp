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
#include "edit.h"
#include "orient.h"

using namespace std;

/**
 * Print a list of edits, separated by commas.
 */
void Edit::print(ostream& os, const EList<Edit>& edits) {
	for(size_t i = 0; i < edits.size(); i++) {
		os << edits[i];
		if(i < edits.size()-1) {
			os << ',';
		}
	}
}

/**
 * Flip all the edits.pos fields so that they're with respect to
 * the other end of the read (of length 'sz').
 */
void Edit::invertPoss(EList<Edit>& edits, size_t sz) {
	// Invert elements
	edits.reverse();
	// Invert all the .pos's
	for(size_t i = 0; i < edits.size(); i++) {
		assert_lt(edits[i].pos, sz);
		edits[i].pos = (uint32_t)(sz - edits[i].pos - (edits[i].type == EDIT_TYPE_INS ? 0 : 1));
	}
}

/**
 * Flip all the edits.pos fields so that they're with respect to
 * the other end of the read (of length 'sz') and complement all
 * DNA characters.
 */
void Edit::invertCompPoss(EList<Edit>& edits, size_t sz, bool color) {
	// Invert elements
	edits.reverse();
	// Invert all the .pos's
	for(size_t i = 0; i < edits.size(); i++) {
		Edit& e = edits[i];
		assert_lt(e.pos, sz);
		e.pos = (uint32_t)(sz - e.pos - (e.type == EDIT_TYPE_INS ? 0 : 1));
		if(!color) {
			e.chr = compDna(e.chr);
			e.qchr = compDna(e.qchr);
		}
	}
}

/**
 * For now, we pretend that the alignment is in the forward orientation
 * and that the Edits are listed from left- to right-hand side.
 */
void Edit::printQAlign(std::ostream& os,
                       const BTDnaString& read,
                       const EList<Edit>& edits)
{
	printQAlign(os, "", read, edits);
}

/**
 * For now, we pretend that the alignment is in the forward orientation
 * and that the Edits are listed from left- to right-hand side.
 */
void Edit::printQAlign(std::ostream& os,
                       const char *prefix,
                       const BTDnaString& read,
                       const EList<Edit>& edits)
{
	size_t eidx = 0;
	os << prefix << "Rd: ";
	// Print read
	for(size_t i = 0; i < read.length(); i++) {
		bool ins = false, del = false, mm = false;
		while(eidx < edits.size() && edits[eidx].pos == i) {
			if(edits[eidx].isInsert()) {
				ins = true;
				os << '-';
			} else if(edits[eidx].isDelete()) {
				del = true;
				//assert_eq(edits[eidx].qchr, read.toChar(i));
				os << read.toChar(i);
			} else {
				mm = true;
				assert(edits[eidx].isMismatch());
				//assert_eq(edits[eidx].qchr, read.toChar(i));
				os << read.toChar(i);
			}
			eidx++;
		}
		if(!del && !mm) os << read.toChar(i);
	}
	os << endl;
	os << prefix << "    ";
	eidx = 0;
	// Print match bars
	for(size_t i = 0; i < read.length(); i++) {
		bool ins = false, del = false, mm = false;
		while(eidx < edits.size() && edits[eidx].pos == i) {
			if(edits[eidx].isInsert()) {
				ins = true;
				os << ' ';
			} else if(edits[eidx].isDelete()) {
				del = true;
				os << ' ';
			} else {
				mm = true;
				assert(edits[eidx].isMismatch());
				os << ' ';
			}
			eidx++;
		}
		if(!del && !mm) os << '|';
	}
	os << endl;
	os << prefix << "Rf: ";
	eidx = 0;
	// Print reference
	for(size_t i = 0; i < read.length(); i++) {
		bool ins = false, del = false, mm = false;
		while(eidx < edits.size() && edits[eidx].pos == i) {
			if(edits[eidx].isInsert()) {
				ins = true;
				os << (char)edits[eidx].chr;
			} else if(edits[eidx].isDelete()) {
				del = true;
				os << '-';
			} else {
				mm = true;
				assert(edits[eidx].isMismatch());
				os << (char)edits[eidx].chr;
			}
			eidx++;
		}
		if(!del && !mm) os << read.toChar(i);
	}
	os << endl;
}
/**
 * For now, we pretend that the alignment is in the forward orientation
 * and that the Edits are listed from left- to right-hand side.
 */
void Edit::printQAlign(std::ostream& os,
                       const char *prefix,
                       const BTDnaString& read,
                       const EList<Edit>& edits,
                       bool fw)
{
	EList<Edit>& eedits = const_cast<EList<Edit>&>(edits);
	if(!fw) Edit::invertPoss(eedits, read.length());
	Edit::printQAlign(os, prefix, read, eedits);
	if(!fw) Edit::invertPoss(eedits, read.length());
}

/**
 * Given a read string and some edits, generate and append the
 * corresponding reference string to 'ref'.
 */
void Edit::toRef(const BTDnaString& read,
                 const EList<Edit>& edits,
                 BTDnaString& ref)
{
	size_t eidx = 0;
	// Print reference
	for(size_t i = 0; i < read.length(); i++) {
		bool ins = false, del = false, mm = false;
		while(eidx < edits.size() && edits[eidx].pos == i) {
			if(edits[eidx].isInsert()) {
				ins = true;
				ref.appendChar((char)edits[eidx].chr);
			} else if(edits[eidx].isDelete()) {
				del = true;
			} else {
				mm = true;
				assert(edits[eidx].isMismatch());
				ref.appendChar((char)edits[eidx].chr);
			}
			eidx++;
		}
		if(!del && !mm) {
			ref.append(read[i]);
		}
	}
}

/**
 * Sanity check a set of edits with respect to a read sequence.
 */
bool Edit::repOk(
	const EList<Edit>& edits,
	const BTDnaString& s,
	size_t clipLeft,
	size_t clipRight)
{
	for(size_t i = 0; i < edits.size(); i++) {
		const Edit& e = edits[i];
		size_t pos = e.pos;
		if(i > 0) {
			// poss must range from low to high
			assert_geq(pos, edits[i-1].pos);
		}
		bool del = false, ins = false, mm = false;
		while(i < edits.size() && edits[i].pos == pos) {
			const Edit& ee = edits[i];
			if(ee.pos >= s.length()-clipRight) {
				printQAlign(std::cerr, "  ", s, edits);
				assert_lt(ee.pos, s.length()-clipRight);
			}
			if(ee.pos < clipLeft) {
				printQAlign(std::cerr, "  ", s, edits);
				assert_geq(ee.pos, clipLeft);
			}
			if(ee.qchr != '-') {
				if(!ee.isDelete() && !ee.isMismatch()) {
					printQAlign(std::cerr, "  ", s, edits);
					assert(ee.isDelete() || !ee.isMismatch());
				}
				if((char)ee.qchr != s.toChar(ee.pos)) {
					printQAlign(std::cerr, "  ", s, edits);
					assert_eq((char)ee.qchr, s.toChar(ee.pos));
				}
			}
			if(ee.isMismatch()) {
				assert(!mm);
				mm = true;
				assert(!del);
			} else if(ee.isInsert()) {
				ins = true;
				assert(!mm);
			} else if(ee.isDelete()) {
				assert(!mm);
				assert(!del);
				del = true;
			}
			i++;
		}
	}
	return true;
}

/**
 * Sanity check a set of edits with respect to a read sequence (which
 * may have aligned as reverse-complement).
 *
 * TODO: This is made much more complicated by the need to take both
 * orientation and clipping into account.  Need to rethink and rewrite
 * to be less clunky.
 */
bool Edit::repOk(
	const EList<Edit>& edits,
	const BTDnaString& s,
	int orient,
	size_t clip5,
	size_t clip3)
{
	// Don't actually know how to handle clipping on the 5' end
	assert_eq(0, clip5);
	bool fivepleft = is5pLeft(orient);
	size_t clipLeft  = fivepleft ? clip5 : clip3;
	size_t clipRight = fivepleft ? clip3 : clip5;
	EList<Edit>& eedits = const_cast<EList<Edit>&>(edits);
	if(!fivepleft) Edit::invertPoss(eedits, s.length());
	Edit::repOk(eedits, s, clipLeft, clipRight);
	if(!fivepleft) Edit::invertPoss(eedits, s.length());
	return true;
}

/**
 * Merge second argument into the first.  Assume both are sorted to
 * begin with.
 */
void Edit::merge(EList<Edit>& dst, const EList<Edit>& src) {
	size_t di = 0, si = 0;
	while(di < dst.size()) {
		if(src[si].pos < dst[di].pos) {
			dst.insert(src[si], di);
			si++; di++;
		} else if(src[si].pos == dst[di].pos) {
			// There can be two inserts at a given position, but we
			// can't merge them because there's no way to know their
			// order
			assert(src[si].isInsert() != dst[di].isInsert());
			if(src[si].isInsert()) {
				dst.insert(src[si], di);
				si++; di++;
			} else if(dst[di].isInsert()) {
				di++;
			}
		}
	}
	while(si < src.size()) dst.push_back(src[si++]);
}

ostream& operator<< (ostream& os, const Edit& e) {
	os << e.pos << ":" << (char)e.chr << ">" << (char)e.qchr;
	return os;
}
