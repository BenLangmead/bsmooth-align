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

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <iostream>
#include <string>
#include "threading.h"
#include "hit_set.h"
#include "edit.h"
#include "ds.h"
#include "sstring.h"
#include "read.h"

/**
 * Abstract parent for output classes.
 */
class AlignOutput {
public:
	AlignOutput(const std::string& of) {
		MUTEX_INIT(lock_);
		alOuts_.resize(1);
		if(of != "-") alOuts_.back().setFile(of);
	}
	virtual ~AlignOutput() { flush(); }
	/// Print a header
	virtual void printHeader(
		const EList<std::string>& names,
		const EList<size_t>& lens) = 0;
	/// Print an alignment
	virtual void printAlignment(
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
			int ccost) = 0;
	/// Print a read that failed to align at all
	virtual void printUnalignedRead(
			const Read& rd,
			const BTDnaString& seq,
			const BTString& qual,
			int filtered) = 0;
	/// Print a read that aligned to too many places for the -m limit
	virtual void printMaxedRead(
			const Read& rd,
			const BTDnaString& seq,
			const BTString& qual,
			int stratum,
			size_t oms) = 0;
	/// Flush the output buffer
	virtual void flush();
protected:
	/// Get the output buffer for refenence 'idx'
	virtual OutFileBuf& alOut(size_t idx) {
		assert(!alOuts_.empty());
		if(idx < alOuts_.size()) return alOuts_[idx];
		return alOuts_[0];
	}
	/// Get the output buffer for unaligned reads
	virtual OutFileBuf& unalOut() {
		return alOuts_[0];
	}
	/// Get the output buffer for maxed-out reads
	virtual OutFileBuf& maxOut() {
		return alOuts_[0];
	}
	MUTEX_T lock_;
	EList<OutFileBuf> alOuts_;
	OutFileBuf unalOut_;
	OutFileBuf maxOut_;
};

/**
 * Bowtie-like verbose output.
 */
class BowtieOutput : public AlignOutput {
public:
	BowtieOutput(
		const std::string& of,
		bool fullRef,
		bool printCost,
		bool refIdx,
		bool orientEx) :
		AlignOutput(of),
		printCost_(printCost),
		fullRef_(fullRef),
		refIdx_(refIdx),
		orientEx_(orientEx) { }
	virtual ~BowtieOutput() { }
	/// Print a header
	virtual void printHeader(
		const EList<std::string>& names,
		const EList<size_t>& lens) { }
	/// Print an alignment
	virtual void printAlignment(
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
			int ccost);
	/// Print a read that failed to align at all
	virtual void printUnalignedRead(
			const Read& rd,
			const BTDnaString& seq,
			const BTString& qual,
			int filtered);
	/// Print a read that aligned to too many places for the -m limit
	virtual void printMaxedRead(
			const Read& rd,
			const BTDnaString& seq,
			const BTString& qual,
			int stratum,
			size_t oms);
protected:
	bool printCost_;
	bool fullRef_;
	bool refIdx_;
	bool orientEx_;
};

/**
 * No output.
 */
class StubOutput : public AlignOutput {
public:
	StubOutput() : AlignOutput("-") { }
	virtual ~StubOutput() { }
	/// Print a header
	virtual void printHeader(
		const EList<std::string>& names,
		const EList<size_t>& lens) { }
	/// Print an alignment
	virtual void printAlignment(
			const Read&,
			const BTDnaString&,
			const BTString&,
			bool,
			size_t, // clip5
			size_t, // clip3
			size_t,
			const char *,
			size_t,
			int,
			const EList<Edit>&,
			const EList<Edit>&,
			const EList<Edit>&,
			const EList<Edit>&,
			size_t,
			size_t,
			int,
			int,
			int) { }
	/// Print a read that failed to align at all
	virtual void printUnalignedRead(
			const Read& rd,
			const BTDnaString& seq,
			const BTString& qual,
			int filtered) { };
	/// Print a read that aligned to too many places for the -m limit
	virtual void printMaxedRead(
			const Read&,
			const BTDnaString&,
			const BTString&,
			int,
			size_t) { }
};

/**
 * SAM output.
 */
class SamOutput : public AlignOutput {

public:

	SamOutput(
		const std::string& of,
		bool fullRef,
		bool refIdx,
		bool bisulfite,
		bool cscq) :
		AlignOutput(of),
		fullRef_(fullRef),
		refIdx_(refIdx),
		bisulfite_(bisulfite),
		cscq_(cscq) { }
	
	virtual ~SamOutput() { }
	/// Print a header
	virtual void printHeader(
		const EList<std::string>& names,
		const EList<size_t>& lens);
	/// Print an alignment
	virtual void printAlignment(
			const Read& rd,
			const BTDnaString& seq,
			const BTString& qual,
			bool preClipped,  // sequence is already clipped?
			size_t clip5,
			size_t clip3,
			size_t refIdx,
			const char * refName,
			size_t refOff,
			int	orient,
			const EList<Edit>& edits,
			const EList<Edit>& aedits,
			const EList<Edit>& cedits,
			const EList<Edit>& ccedits,
			size_t oms,
			size_t mapq,
			int stratum,
			int cost,
			int ccost);
	/// Print a read that failed to align at all
	virtual void printUnalignedRead(
			const Read& rd,
			const BTDnaString& seq,
			const BTString& qual,
			int filtered);
	/// Print a read that aligned to too many places for the -m limit
	virtual void printMaxedRead(
			const Read& rd,
			const BTDnaString& seq,
			const BTString& qual,
			int stratum,
			size_t oms);

protected:

	bool fullRef_;   // print full reference name?
	bool refIdx_;    // print reference index instead of name?
	bool bisulfite_; // print XB:Z flag?
	bool cscq_;      // print CS:Z and CQ:Z flags for colorspace alignments?
};

#endif /* OUTPUT_H_ */
