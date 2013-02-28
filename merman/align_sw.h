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

#ifndef ALIGN_SW_H_
#define ALIGN_SW_H_

#include "aligner.h"

class SWAligner : public Aligner {

public:

	SWAligner(int nthreads) : Aligner(nthreads) { }

	virtual ~SWAligner () { }

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
};

#endif /* ALIGN_SW_H_ */
