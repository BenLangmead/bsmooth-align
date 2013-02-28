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

#ifndef ALIGN_NAIVE_H_
#define ALIGN_NAIVE_H_

#include "aligner.h"
#include "color_dec.h"

/**
 * Per-thread state for naive searching.
 */
struct NaivePerThreadState {
	void clear() {
		edits.clear();
		aedits.clear();
		cedits.clear();
		ccedits.clear();
	}
	EList<Edit> edits;
	EList<Edit> aedits;
	EList<Edit> cedits;
	EList<Edit> ccedits;
	ColorspaceDecoder dec_;
};

class NaiveAligner : public Aligner {

public:

	NaiveAligner(
		int nthreads,
		bool bisulfiteC,
		bool bisulfiteCpG,
		bool maqRound) :
		Aligner(nthreads),
		bisulfiteC_(bisulfiteC),
		bisulfiteCpG_(bisulfiteCpG),
		maqRound_(maqRound)
	{
		threadState_.resize(nthreads);
		for(int i = 0; i < nthreads; i++) threadState_[i].clear();
	}

	virtual ~NaiveAligner() { }

	/**
	 * Query the iupac index with the given read, requiring that all
	 * reported alignments must be in a better stratum than
	 * 'worstStratum'.
	 */
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
		int threadid);

protected:
	EList<Edit> nedits_, aedits_, cedits_, ccedits_;
	EList<NaivePerThreadState> threadState_;
	bool bisulfiteC_;
	bool bisulfiteCpG_;
	bool maqRound_;
};

#endif /* ALIGN_NAIVE_H_ */
