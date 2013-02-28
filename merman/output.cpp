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

#include <string>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include "edit.h"
#include "hit_set.h"
#include "threading.h"
#include "filebuf.h"
#include "output.h"
#include "orient.h"

using namespace std;

enum {
	SAM_FLAG_PAIRED = 1,
	SAM_FLAG_MAPPED_PAIRED = 2,
	SAM_FLAG_UNMAPPED = 4,
	SAM_FLAG_MATE_UNMAPPED = 8,
	SAM_FLAG_QUERY_STRAND = 16,
	SAM_FLAG_MATE_STRAND = 32,
	SAM_FLAG_FIRST_IN_PAIR = 64,
	SAM_FLAG_SECOND_IN_PAIR = 128,
	SAM_FLAG_NOT_PRIMARY = 256,
	SAM_FLAG_FAILS_CHECKS = 512,
	SAM_FLAG_DUPLICATE = 1024
};

/**
 * Flush all output streams.
 */
void AlignOutput::flush() {
	for(size_t i = 0; i < alOuts_.size(); i++) alOuts_[i].flush();
	unalOut().flush();
	maxOut().flush();
}

/**
 * C++ version char* style "itoa":
 */
template<typename T>
char* toa10(T value, char* result) {
	// Check that base is valid
	char* out = result;
	bool neg = false;
	if(!(value > 0 || value == 0)) {
		neg = true;
		value = -value;
	}
	T quotient = value;
	do {
		*out = "0123456789"[ quotient % 10 ];
		++out;
		quotient /= 10;
	} while ( quotient );

	// Only apply negative sign for base 10
	if(neg) *out++ = '-';
	std::reverse( result, out );
	*out = 0; // terminator
	return out;
}

/**
 * Print a single Edit object.
 */
static void printBowtieEdit(
	OutFileBuf& os,
	const Edit& e,
	bool colorize = false)
{
	assert(e.initialized());
	char buf[100];
	toa10<size_t>(e.pos, buf);
	os.writeChars(buf, strlen(buf));
	if(e.type == EDIT_TYPE_SNP) os.write('S');
	assert(isprint((char)e.chr));
	os.write(':');
	if(colorize) {
		os.writeChars(dna2colstr[e.chr]);
	} else {
		os.write(e.chr);
	}
	os.write('>');
	assert_neq((char)e.qchr, (char)e.chr);
	assert(isprint((char)e.qchr));
	assert_neq(e.qchr, 0);
	if(colorize) {
		os.writeChars(dna2colstr[e.qchr]);
	} else {
		os.write(e.qchr);
	}
}

/**
 * Print an alignment in a format similar to the default Bowtie format.  The
 * edit string is formatted differently to include information about colorspace
 * edits and ambiguous reference nucleotides.  Also, the user can ask that the
 * orientation field (2nd column) contain code flexible enough to describe the
 * bisulfite Watson and Crick strands and their reverse complements.
 */
void BowtieOutput::printAlignment(
	const Read& rd,
	const BTDnaString& seq,
	const BTString& qual,
	bool preClipped,  // sequence is already clipped?
	size_t clip5,
	size_t clip3,
	size_t refIdx,
	const char * refName,
	size_t refOff,
	int orient,
	const EList<Edit>& edits,
	const EList<Edit>& aedits,
	const EList<Edit>& cedits,
	const EList<Edit>& ccedits,
	size_t oms,
	size_t mapq,
	int stratum,
	int cost,
	int ccost)
{
	const BTString& name = rd.name;
	ThreadSafe t(&lock_);
	OutFileBuf& os = alOut(refIdx);
	char buf[100];
	os.writeChars(name.toZBuf(), name.length()); os.write('\t');
	if(orientEx_) {
		os.writeString(orientStr(orient));
	} else {
		os.write(is5pLeft(orient) ? '+' : '-');
	}
	os.write('\t');
	if(refIdx_) {
		toa10<size_t>(refIdx, buf);
		os.writeChars(buf, strlen(buf));
	} else {
		if(fullRef_) os.writeChars(refName);
		else os.writeCharsUptoSpace(refName);
	}
	os.write('\t');
	toa10<size_t>(refOff, buf);
	os.writeChars(buf, strlen(buf)); os.write('\t');
#ifndef NDEBUG
	for(size_t i = 0; i < seq.length(); i++) {
		assert_in(seq.toChar(i), "ACGTN");
	}
#endif
	// Apply clipping
	size_t trimLeft = 0;
	size_t trimRight = 0;
	if(!preClipped) {
		trimLeft  = is5pLeft(orient) ? clip5 : clip3;
		trimRight = is5pLeft(orient) ? clip3 : clip5;
	}
	for(size_t i = trimLeft; i < seq.length()-trimRight; i++) {
		os.write("ACGTN"[(int)seq[i]]);
	}
	os.write('\t');
#ifndef NDEBUG
	for(size_t i = 0; i < qual.length(); i++) {
		assert_geq(qual[i], 33);
	}
#endif
	// Apply clipping; TODO: In colorspace mode, we need a really silly check
	// to account for the fact that (a) we're printing COLOR qualities and
	// NUCLEOTIDE characters, (b) we have already trimmed the sequence but NOT
	// THE COLOR QUALITIES.
	size_t quallen = qual.length()-trimRight;
	if(quallen > seq.length()-trimRight) {
		quallen = seq.length()+1;
	}
	for(size_t i = trimLeft; i < quallen; i++) {
		os.write((int)qual[i]);
	}
	os.write('\t');
	toa10<size_t>(oms, buf);
	os.writeChars(buf, strlen(buf)); os.write('\t');
	int editsPrinted = 0;
	// Print nucleotide edits
	for(size_t i = 0; i < edits.size(); i++) {
		editsPrinted++;
		printBowtieEdit(os, edits[i]);
		if(i < edits.size()-1) os.write(',');
	}
	// Print ambiguous nucleotide resolutions
	os.write(';');
	for(size_t i = 0; i < aedits.size(); i++) {
		editsPrinted++;
		printBowtieEdit(os, aedits[i]);
		if(i < aedits.size()-1) os.write(',');
	}
	// Print color miscalls
	os.write(';');
	for(size_t i = 0; i < cedits.size(); i++) {
		editsPrinted++;
		printBowtieEdit(os, cedits[i], true);
		if(i < cedits.size()-1) os.write(',');
	}
	// Print color edits
	os.write(';');
	for(size_t i = 0; i < ccedits.size(); i++) {
		editsPrinted++;
		printBowtieEdit(os, ccedits[i], true);
		if(i < ccedits.size()-1) os.write(',');
	}
	if(printCost_) {
		os.write('\t');
		toa10<int>(stratum, buf);
		os.writeChars(buf, strlen(buf));
		os.write('\t');
		toa10<unsigned int>(cost, buf);
		os.writeChars(buf, strlen(buf));
	}
	os.write('\n');
}

/**
 * Print an alignment record for a read that failed to align.
 */
void BowtieOutput::printUnalignedRead(
	const Read& rd,
	const BTDnaString& seq,
	const BTString& qual,
	int filtered)
{
	const BTString& name = rd.name;
	ThreadSafe t(&lock_);
	unalOut().writeChars(name.toZBuf()); unalOut().write('\t');
	unalOut().write('*'); unalOut().write('\t');
	unalOut().write('*'); unalOut().write('\t');
	unalOut().write('0'); unalOut().write('\t');
	if(true) {
		maxOut().write('*'); maxOut().write('\t');
		maxOut().write('*'); maxOut().write('\t');
	} else {
		unalOut().writeChars(seq.toZBuf()); unalOut().write('\t');
		unalOut().writeChars(qual.toZBuf()); unalOut().write('\t');
	}
	unalOut().write('0'); unalOut().write('\t');
	unalOut().write('-'); // edits
	if(this->printCost_) {
		unalOut().write('\t');
		unalOut().write('0');
		unalOut().write('\t');
		unalOut().write('0');
	}
	unalOut().write('\n');
}

/**
 * Print an alignment record for a read that aligned too many times.
 */
void BowtieOutput::printMaxedRead(
	const Read& rd,
	const BTDnaString& seq,
	const BTString& qual,
	int stratum,
	size_t oms)
{
	const BTString& name = rd.name;
	ThreadSafe t(&lock_);
	char buf[100];
	maxOut().writeChars(name.toZBuf()); maxOut().write('\t');
	maxOut().write('*'); maxOut().write('\t');
	maxOut().write('*'); maxOut().write('\t');
	maxOut().write('0'); maxOut().write('\t');
	if(true) {
		maxOut().write('*'); maxOut().write('\t');
		maxOut().write('*'); maxOut().write('\t');
	} else {
		maxOut().writeCharsUptoSpace(seq.toZBuf());   maxOut().write('\t');
		maxOut().writeCharsUptoSpace(qual.toZBuf()); maxOut().write('\t');
	}
	toa10<size_t>(oms, buf);
	maxOut().writeChars(buf, strlen(buf));
	maxOut().write('\t');
	unalOut().write('-'); // edits
	if(this->printCost_) {
		maxOut().write('\t');
		assert_lt(stratum, 10);
		assert_geq(stratum, 0);
		toa10<int>(stratum, buf);
		maxOut().writeChars(buf, strlen(buf));
		maxOut().write('\t');
		maxOut().write('0');
	}
	maxOut().write('\n');
}

/**
 * Print a read name in a way that doesn't violate SAM's character constraints.
 * \*|[!-()+-<>-~][!-~]*
 */
template<typename T>
static void printSamReadName(OutFileBuf& o, const T& name) {
	for(size_t i = 0; i < Class_sstr_len<T>::sstr_len(name); i++) {
		if(isspace(name[i])) {
			return;
		}
		o.write(name[i]);
	}
}

/**
 * Print a reference name in a way that doesn't violate SAM's character
 * constraints. \*|[!-()+-<>-~][!-~]*
 */
template<typename T>
static void printSamRefName(OutFileBuf& o, const T& name) {
	for(size_t i = 0; i < Class_sstr_len<T>::sstr_len(name); i++) {
		if(isspace(name[i])) {
			return;
		}
		o.write(name[i]);
	}
}


/**
 * Print an alignment in SAM format.
 *
 * If a bisulfite protocol that yields all 4 strands is used, we must
 * distinguish among 4 (not just +/-) in the SAM output somehow.  We use a SAM
 * optional flag: (XB:Z, which can take values BW, BC, BWR, and BCR).
 */
void SamOutput::printAlignment(
	const Read& rd,
	const BTDnaString& seq,
	const BTString& qual,
	bool preClipped,  // sequence is already clipped?
	size_t clip5,
	size_t clip3,
	size_t refIdx,
	const char * refName,
	size_t refOff,
	int orient,
	const EList<Edit>& edits,
	const EList<Edit>& aedits,
	const EList<Edit>& cedits,
	const EList<Edit>& ccedits,
	size_t oms,
	size_t mapq,
	int stratum,
	int cost,
	int ccost)
{
	const BTString& name = rd.name;
	ThreadSafe t(&lock_);
	char buf[100];
	OutFileBuf& os = alOut(refIdx);
	// QNAME
	printSamReadName(os, name);
	os.write('\t');
	// FLAG
	int flags = 0;
	if(!orientFw(orient)) {
		flags |= SAM_FLAG_QUERY_STRAND;
	}
	toa10<int>(flags, buf);
	os.writeChars(buf, strlen(buf));
	os.write('\t');
	// RNAME
	if(refIdx_) {
		toa10<size_t>(refIdx, buf);
		printSamRefName(os, (const char *)buf);
	} else {
		if(fullRef_) {
			os.writeChars(refName);
		} else {
			printSamRefName(os, refName);
		}
	}
	os.write('\t');
	// POS
	toa10<size_t>(refOff+1, buf);
	os.writeChars(buf, strlen(buf));
	os.write('\t');
	// MAPQ
	toa10<size_t>(mapq, buf);
	os.writeChars(buf, strlen(buf));
	os.write('\t');
	// CIGAR
	toa10<size_t>(seq.length(), buf);
	os.writeChars(buf, strlen(buf));
	os.write('M');
	os.write('\t');
	// MRNM
	os.write('*');
	os.write('\t');
	// MPOS
	os.write('0');
	os.write('\t');
	// ISIZE
	os.write('0');
	os.write('\t');
	// SEQ; TODO: use SAM clipping instead of hard clipping here
	size_t trimLeft = 0;
	size_t trimRight = 0;
	if(!preClipped) {
		trimLeft  = is5pLeft(orient) ? clip5 : clip3;
		trimRight = is5pLeft(orient) ? clip3 : clip5;
	}
	for(size_t i = trimLeft; i < seq.length()-trimRight; i++) {
		os.write("ACGTN"[(int)seq[i]]);
	}
	os.write('\t');
	// QUAL; TODO: use SAM clipping instead of clipping here
	for(size_t i = trimLeft; i < qual.length()-trimRight; i++) {
		os.write((int)qual[i]);
	}

	// Write MD string.  Easy when there's no gaps.
	bool inclAmbig = true;
	os.writeChars("\tMD:Z:", 6);
	unsigned int lastpos = 0;
	int nm = 0;
	if(inclAmbig) {
		int eidx = 0, aidx = 0;
		if(is5pLeft(orient)) {
			while(eidx < (int)edits.size() || aidx < (int)aedits.size()) {
				Edit e;
				if(eidx == (int)edits.size() ||
				   (aidx < (int)aedits.size() &&
				    aedits[aidx].pos < edits[eidx].pos))
				{
					e = aedits[aidx];
					aidx++;
				} else {
					e = edits[eidx];
					eidx++;
					if(e.type != EDIT_TYPE_SNP) nm++;
				}
				unsigned int pos = e.pos;
				unsigned int run = pos - lastpos;
				lastpos = pos+1;
				toa10<unsigned int>(run, buf);
				os.writeChars(buf, strlen(buf));
				os.write(e.chr);
			}
		} else {
			eidx = (int)edits.size()-1;
			aidx = (int)aedits.size()-1;
			while(eidx >= 0 || aidx >= 0) {
				Edit e;
				if(eidx < 0 ||
				   (aidx >= 0 && aedits[aidx].pos > edits[eidx].pos))
				{
					e = aedits[aidx];
					aidx--;
				} else {
					e = edits[eidx];
					eidx--;
					if(e.type != EDIT_TYPE_SNP) nm++;
				}
				unsigned int pos = (unsigned int)seq.length()-e.pos-1;
				unsigned int run = pos - lastpos;
				lastpos = pos+1;
				toa10<unsigned int>(run, buf);
				os.writeChars(buf, strlen(buf));
				os.write(e.chr);
			}
		}
	} else {
		if(is5pLeft(orient)) {
			for(int i = 0; i < (int)edits.size(); i++) {
				const Edit& e = edits[i];
				unsigned int pos = e.pos;
				unsigned int run = pos - lastpos;
				lastpos = pos+1;
				toa10<unsigned int>(run, buf);
				os.writeChars(buf, strlen(buf));
				os.write(e.chr);
				if(e.type != EDIT_TYPE_SNP) nm++;
			}
		} else {
			for(int i = (int)edits.size()-1; i >= 0; i--) {
				const Edit& e = edits[i];
				unsigned int pos = (unsigned int)seq.length()-e.pos-1;
				unsigned int run = pos - lastpos;
				lastpos = pos+1;
				toa10<unsigned int>(run, buf);
				os.writeChars(buf, strlen(buf));
				os.write(e.chr);
				if(e.type != EDIT_TYPE_SNP) nm++;
			}
		}
	}
	toa10<unsigned int>((unsigned int)seq.length()-lastpos, buf);
	os.writeChars(buf, strlen(buf));
	
	// Add optional alignment cost field
	os.writeChars("\tAS:i:");
	toa10<int>(-cost, buf);
	os.writeChars(buf, strlen(buf));
	
	// Add optional edit distance field
	os.writeChars("\tNM:i:");
	toa10<int>(nm, buf);
	os.writeChars(buf, strlen(buf));

	// Print bisulfite strand flag
	if(bisulfite_) {
		os.writeChars("\tXB:Z:");
		os.writeChars(orientStr(orient));
	}

	// Possibly print CS:Z and CQ:Z strings
	if(rd.color && cscq_) {
		os.writeChars("\tCS:Z:");
		rd.writeOrigSeq(os);
		os.writeChars("\tCQ:Z:");
		rd.writeOrigQual(os);
		os.writeChars("\tYC:i:");
		toa10<int>(-ccost, buf);
		os.writeChars(buf, strlen(buf));
	}
	
	os.write('\n');
}

void SamOutput::printUnalignedRead(
	const Read& rd,
	const BTDnaString& seq,
	const BTString& qual,
	int filtered)
{
	ThreadSafe t(&lock_);
	const BTString& name = rd.name;
	char buf[100];
	OutFileBuf& os = alOut(0);
	// QNAME
	printSamReadName(os, name);
	os.write('\t');
	// FLAG
	int flags = SAM_FLAG_UNMAPPED;
	toa10<int>(flags, buf);
	os.writeChars(buf, strlen(buf));
	os.write('\t');
	// RNAME
	os.write('*');
	os.write('\t');
	// POS
	os.write('0');
	os.write('\t');
	// MAPQ
	os.write('0');
	os.write('\t');
	// CIGAR
	os.write('*');
	os.write('\t');
	// MRNM
	os.write('*');
	os.write('\t');
	// MPOS
	os.write('0');
	os.write('\t');
	// ISIZE
	os.write('0');
	os.write('\t');
	// SEQ; TODO: use SAM clipping instead of hard clipping here
	for(size_t i = 0; i < seq.length(); i++) {
		os.write("ACGTN"[(int)seq[i]]);
	}
	os.write('\t');
	// QUAL; TODO: use SAM clipping instead of clipping here
	for(size_t i = 0; i < qual.length(); i++) {
		os.write((int)qual[i]);
	}
	// Add optional alignment cost field
	os.writeChars("\tYM:i:0");
	
	// Add optional field describing why the read was filtered.  See filtered.h
	// Read was too short for the index being used
	// FILTER_TOO_SHORT_FOR_INDEX         = 0x01,
	//
	// Read was shorter than the minimum length allowed by alignment parameters
	// directly governing alignment length (e.g. Merman's --minlen)
	// FILTER_TOO_SHORT_FOR_MINLEN_PARAMS = 0x02,
	//
	// Read was shorter than the minimum length allowed by alignment parameters
	// for scoring thresholds.  E.g. if the minimum score is 10 and the match
	// bonus is 1, a 9-nt read can't possibly align
	// FILTER_TOO_SHORT_FOR_SCORE_PARAMs  = 0x04,
	//
	// An upstream QC step flagged the read as failing QC
	// FILTER_QC                          = 0x08,
	//
	// Too many characters in the read are non-A/C/G/T
	// FILTER_TOO_MANY_AMBIGS             = 0x10,
	//
	// Too many low quality positions and/or homopolymers
	// FILTER_BAD_QUALITIES               = 0x20
	//
	if(filtered > 0) {
		os.writeChars("\tXF:i:");
		toa10<int>(filtered, buf);
		os.writeChars(buf, strlen(buf));
	}
	os.write('\n');
}

void SamOutput::printMaxedRead(
	const Read& rd,
	const BTDnaString& seq,
	const BTString& qual,
	int stratum,
	size_t oms)
{
	ThreadSafe t(&lock_);
	const BTString& name = rd.name;
	char buf[100];
	OutFileBuf& os = alOut(0);
	// QNAME
	printSamReadName(os, name);
	os.write('\t');
	// FLAG
	int flags = SAM_FLAG_UNMAPPED;
	toa10<int>(flags, buf);
	os.writeChars(buf, strlen(buf));
	os.write('\t');
	// RNAME
	os.write('*');
	os.write('\t');
	// POS
	os.write('0');
	os.write('\t');
	// MAPQ
	os.write('0');
	os.write('\t');
	// CIGAR
	os.write('*');
	os.write('\t');
	// MRNM
	os.write('*');
	os.write('\t');
	// MPOS
	os.write('0');
	os.write('\t');
	// ISIZE
	os.write('0');
	os.write('\t');
	// SEQ; TODO: use SAM clipping instead of hard clipping here
	for(size_t i = 0; i < seq.length(); i++) {
		os.write("ACGTN"[(int)seq[i]]);
	}
	os.write('\t');
	// QUAL; TODO: use SAM clipping instead of clipping here
	for(size_t i = 0; i < qual.length(); i++) {
		os.write((int)qual[i]);
	}
	// Add optional alignment cost field
	os.writeChars("\tYM:i:1");
	os.write('\n');
}

extern string cmdline;

static void printSamRefName(ostream& os, const string& name) {
	for(size_t i = 0; i < name.length(); i++) {
		if(isspace(name[i])) {
			break;
		}
		os << name[i];
	}
}

/**
 * Print SAM header to given output stream.
 */
static void printSamHeader(
	OutFileBuf& of,
	const EList<std::string>& names,
	const EList<size_t>& lens)
{
	ostringstream os;
	os << "@HD\tVN:1.0\tSO:unsorted" << endl;
	os << "@PG\tID:Merman\tVN:" << MERMAN_VERSION << "\tCL:\"" << cmdline << "\"" << endl;
	assert_eq(names.size(), lens.size());
	for(size_t i = 0; i < names.size(); i++) {
		os << "@SQ\tSN:"; printSamRefName(os, names[i]);
		os << "\tLN:" << lens[i] << endl;
	}
	of.writeString(os.str());
}

/**
 * Print the SAM header to
 */
void SamOutput::printHeader(
	const EList<std::string>& names,
	const EList<size_t>& lens)
{
	ThreadSafe t(&lock_);
	printSamHeader(alOuts_[0], names, lens);
}
