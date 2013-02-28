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

#ifndef EDIT_H_
#define EDIT_H_

#include <iostream>
#include <stdint.h>
#include "assert_helpers.h"
#include "filebuf.h"
#include "sstring.h"
#include "ds.h"

/**
 * 3 types of edits; mismatch (substitution), insertion in the
 * reference, deletion in the reference.
 */
enum {
	EDIT_TYPE_INS = 1,
	EDIT_TYPE_DEL,
	EDIT_TYPE_MM,
	EDIT_TYPE_SNP
};

/**
 * Encapsulates an edit between the read sequence and the reference
 * sequence.
 */
struct Edit {
	
	Edit() : pos(1023), reserved(0) { }
	
	Edit(int po, int ch, int qc, int ty) :
	chr(ch), qchr(qc), type(ty), pos(po), reserved(0)
	{
		assert_in(toupper(chr), "ACGTMRWSYKVHDBN-");
		assert_in(qchr, "ACGTN-");
	}
	
	/**
	 * Write Edit to an OutFileBuf.
	 */
	void serialize(OutFileBuf& fb) const {
		assert_eq(4, sizeof(Edit));
		fb.writeChars((const char*)this, 4);
	}
	
	/**
	 * Read Edit from a FileBuf.
	 */
	void deserialize(FileBuf& fb) {
		fb.get((char*)this, 4);
	}
	
	/**
	 * Edit less-than overload.
	 */
	int operator< (const Edit &rhs) const {
		if(pos < rhs.pos) return 1;
		if(pos > rhs.pos) return 0;
		if(type < rhs.type) return 1;
		if(type > rhs.type) return 0;
		if(chr < rhs.chr) return 1;
		if(chr > rhs.chr) return 0;
		return (qchr < rhs.qchr)? 1 : 0;
	}
	
	/**
	 * Edit equals overload.
	 */
	int operator== (const Edit &rhs) const {
		return(pos  == rhs.pos &&
			   chr  == rhs.chr &&
			   qchr == rhs.qchr &&
			   type == rhs.type);
	}
	
	/**
	 * Return true iff this Edit is initialized.
	 */
	bool initialized() const {
		return pos != 1023;
	}
	
	/**
	 * Return true iff this Edit is an initialized insertion.
	 */
	bool isInsert() const {
		return initialized() && type == EDIT_TYPE_INS;
	}
	
	/**
	 * Return true iff this Edit is an initialized deletion.
	 */
	bool isDelete() const {
		return initialized() && type == EDIT_TYPE_DEL;
	}
	
	/**
	 * Return true if this Edit is either an initialized deletion or an
	 * initialized insertion.
	 */
	bool isGap() const {
		return initialized() && (type == EDIT_TYPE_DEL || type == EDIT_TYPE_INS);
	}
	
	/**
	 * Return true iff this Edit is an initialized mismatch.
	 */
	bool isMismatch() const {
		return initialized() && type == EDIT_TYPE_MM;
	}
	
	/**
	 * Flip all the edits.pos fields so that they're with respect to
	 * the other end of the read (of length 'sz').
	 */
	static void invertPoss(EList<Edit>& edits, size_t sz);

	/**
	 * Flip all the edits.pos fields so that they're with respect to
	 * the other end of the read (of length 'sz') and complement all
	 * DNA characters.
	 */
	static void invertCompPoss(EList<Edit>& edits, size_t sz, bool color);
	
	/**
	 * Given a read string and some edits, generate and append the
	 * corresponding reference string to 'ref'.
	 */
	static void toRef(
		const BTDnaString& read,
		const EList<Edit>& edits,
		BTDnaString& ref);
	
	/**
	 * Given a string and its edits with respect to some other string,
	 * print the alignment between the strings with the strings stacked
	 * vertically, with vertical bars denoting matches.
	 */
	static void printQAlign(
		std::ostream& os,
		const BTDnaString& read,
		const EList<Edit>& edits);
	
	/**
	 * Given a string and its edits with respect to some other string,
	 * print the alignment between the strings with the strings stacked
	 * vertically, with vertical bars denoting matches.  Add 'prefix'
	 * before each line of output.
	 */
	static void printQAlign(
		std::ostream& os,
		const char *prefix,
		const BTDnaString& read,
		const EList<Edit>& edits);
	
	/**
	 * Given a string and its edits with respect to some other string,
	 * print the alignment between the strings with the strings stacked
	 * vertically, with vertical bars denoting matches.  Add 'prefix'
	 * before each line of output.  Invert edits first if fw is false.
	 */
	static void printQAlign(
		std::ostream& os,
		const char *prefix,
		const BTDnaString& read,
		const EList<Edit>& edits,
		bool fw);
	
	/**
	 * Sanity check a set of edits with respect to a read sequence.
	 */
	static bool repOk(
		const EList<Edit>& edits,
		const BTDnaString& s,
		size_t clipLeft,
		size_t clipRight);
	
	/**
	 * Sanity check a set of edits with respect to a read sequence (which
	 * may have aligned as reverse-complement).
	 */
	static bool repOk(
		const EList<Edit>& edits,
		const BTDnaString& s,
		int orient,
		size_t clip5,
		size_t clip3);
	
	uint32_t chr       :  8; // reference character involved (for subst and ins)
	uint32_t qchr      :  8; // read character involved (for subst and del)
	uint32_t type      :  4; // 1 -> mm, 2 -> SNP, 3 -> ins, 4 -> del
	uint32_t pos       : 10; // position w/r/t search root
	uint32_t reserved  :  2; // padding
	
	friend std::ostream& operator<< (std::ostream& os, const Edit& e);
	
	/**
	 * Print a comma-separated list of Edits to given output stream.
	 */
	static void print(std::ostream& os, const EList<Edit>& edits);
	
	/**
	 * Merge second argument into the first.  Assume both are sorted to
	 * begin with.
	 */
	static void merge(EList<Edit>& dst, const EList<Edit>& src);
};

#endif /* EDIT_H_ */
