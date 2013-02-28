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

#ifndef ALIGNER_H_
#define ALIGNER_H_

#include <stdint.h>
#include "random_source.h"

// Forward declarations
class ReferenceMap;
class AnnotationMap;
class HitSet;
class Read;
class ReferenceSet;
class AlignOutput;
template<typename T> class EList;

/**
 * Structure containing parameters that define which alignments are
 * valid and reportable, and how to decode them.
 */
struct AlignParams {

	void init(
		int smms,
		int emms,
		int *penc,
		int minl,
		int adjmin,
		bool fw,
		bool rc,
		bool reqa,
		bool disi,
		bool maqr,
		bool ignq,
		int slen,
		int islen,
		int kh,
		int mh,
		bool msamp,
		bool stra,
		int snpp,
		int rdopen,
		int rdex,
		int rfopen,
		int rfex,
		int gapb,
		bool exdec)
	{
		seedMms = smms;
		e2eMms = emms;
		for(int i = 0; i < 1024; i++) penceil[i] = penc[i];
		minLen = minl;
		adjMinLen = adjmin;
		alignfw = fw;
		alignrc = rc;
		requireAnnot = reqa;
		disallowIupac = disi;
		maqRound = maqr;
		ignoreQuals = ignq;
		seedLen = slen;
		iSeedLen = islen;
		khits = kh;
		mhits = mh;
		msample = msamp;
		strata = stra;
		snpPen = snpp;
		readOpenPen = rdopen;
		readExtendPen = rdex;
		refOpenPen = rdopen;
		refExtendPen = rfex;
		gapBarrier = gapb;
		exDecEnds = exdec;
	}

	// Which alignments are valid?
	int  seedMms;       // if a seeded strategy is used, this is seed mismatch ceiling
	int  e2eMms;        // if an end-to-end strategy is used, this is the mismatch ceiling
	int  penceil[1024]; // penality ceiling for alignments of various lengths
	int  minLen;        // reads shorter than this should not be aligned
	int  adjMinLen;     // 
	bool alignfw;       // true -> don't align the forward read
	bool alignrc;       // true -> don't align the reverse complement of the
	                    // read.  In 1-strand non-bisulfite mode, this means
						// that only alignments to the forward strand will be
						// considered
	bool requireAnnot;  // Only report alignments overlapping annotations.
	bool disallowIupac; // 
	bool maqRound;      // round qualities to the nearest 10 and cap at 30 a la Maq
	bool ignoreQuals;   // ignore qualities
	int  seedLen;       // period = iSeedLen - seedWidth + 1;
	int  iSeedLen;      // sampling window
	// Which should we report?
	int  khits;   // -k mode
	int  mhits;   // -m/-M mode
	bool msample; // -M mode
	bool strata;
	// How should we decode?
	int  snpPen;
	// Gaps
	int  readOpenPen;
	int  readExtendPen;
	int  refOpenPen;
	int  refExtendPen;
	int  gapBarrier;
	bool exDecEnds;
};

/**
 * The query() method fills this in to summarize for the caller what
 * happened.
 */
struct AlignResult {
	bool maxed; // whether # of alignments exceeded -m limit
	int64_t hits; // # hits reported
	int64_t seedHits; // # seed hits
	bool bail; // whether to bail after this attempt due to an error

	void clear() {
		maxed = bail = false;
		hits = seedHits = 0;
	}
};

/**
 * Abstract parent for classes that align reads via a query() func
 * like the one below.  We build in a notion of threads so that the
 * details of how concrete subclasses manage per-thread state are
 * hidden from users of this class.
 */
class Aligner {
public:

	Aligner(int nthreads) : nthreads_(nthreads) { }

	virtual ~Aligner() { }

	virtual void query(
		Read& rd,
		const ReferenceSet& refs,
		const ReferenceMap* rmap,
		const AnnotationMap* amap,
		HitSet& hits,
		AlignOutput& os,
		AlignResult& res,
		const AlignParams& p,
		bool randomize,
		RandomSource& rnd,
		int threadid) = 0;

	int maxThreads() { return nthreads_; }

protected:
	int nthreads_;
};

#endif /* ALIGNER_H_ */
