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

/*
 * Classes and routines for decoding color-to-color alignments into
 * base-to-base alignments (and related information).  Two decoders are
 * implemented: one for ungapped alignments, and one for gapped (or ungapped
 * - less efficiently than the ungapped decoder) alignments.
 *
 * Both decoders take a colorspace read (5' end must be on the left) and a
 * reference sequence (which must be arranged to correspond to read characters
 * from 5' to 3', even if the read aligns to the reverse strand) and outputs:
 *
 * 1. A decoded nucleotide sequence, arranged from read's 5' to 3'
 * 2. A decoded quality-value sequence, arranged from read's 5' to 3'
 * 3. A list of Edits describing the decoded nucleotide edits, indexed from 5'
 *    to 3'
 * 4. A list of Edits describing the color mismatches that were classified as
 *    miscalls during decoding, indexed from 5' to 3'
 * 5. A list of Edits describing all color edits, including both those induced
 *    by sequencer miscalls and those induced by nucleotide variants.  Also
 *    indexed from 5' to 3'.
 */

#ifndef COLOR_DEC_H_
#define COLOR_DEC_H_

#include <stdint.h>
#include <string>
#include <utility>
#include <limits>
#include "alphabet.h"
#include "sstring.h"
#include "ds.h"
#include "edit.h"
#include "qual.h"
#include "random_source.h"

typedef int TScore;
#define MAX_SCORE std::numeric_limits<TScore>::max()

// 5-bit pop count
static int alts[] = {
	-1, 1, 1, 2, 1, 2, 2, 3,
	 1, 2, 2, 3, 2, 3, 3, 4,
	 1, 2, 2, 3, 2, 3, 3, 4,
	 2, 3, 3, 4, 3, 4, 4, 5
};
// Index of lowest set bit
static int firsts[] = {
	-1, 0, 1, 0, 2, 0, 1, 0,
	 3, 0, 1, 0, 2, 0, 1, 0,
	 4, 0, 1, 0, 2, 0, 1, 0,
	 3, 0, 1, 0, 2, 0, 1, 0
};

/**
 * Given a mask with up to 5 bits, return an index corresponding to a
 * set bit in the mask, randomly chosen from among all set bits.
 */
static inline int randFromMask(RandomSource& rand, int mask) {
	assert_gt(mask, 0);
	if(alts[mask] == 1) {
		// only one to pick from, pick it via lookup table
		return firsts[mask];
	}
	assert_gt(mask, 0);
	assert_lt(mask, 16);
	int r = rand.nextU32() % alts[mask];
	assert_geq(r, 0);
	assert_lt(r, alts[mask]);
	// could do the following via lookup table too
	for(int i = 0; i < 5; i++) {
		if((mask & (1 << i)) != 0) {
			if(r == 0) return i;
			r--;
		}
	}
	std::cerr << "Shouldn't get here" << std::endl;
	throw 1;
	return -1;
}

/**
 * All of the incoming information to a gapped-decoder cell.
 */
struct UngappedCell {
	// Incoming diagonal edges; [to][from]
	uint16_t mask[4]; // mask encoding which incoming paths carry the best score
	TScore best[4]; // the best incoming score per [to]
};

/**
 * Helper classes for decoding colorspace alignments into nucleotide
 * space alignments.
 */
class ColorspaceDecoderUngapped {
public:
	void backtrack(const BTDnaString& read,
	               size_t readi,
	               size_t readf,
	               BTString& ref,
	               size_t refi,
	               size_t reff,
	               BTDnaString& dseq,
	               EList<Edit>& nedits,
	               EList<Edit>& aedits,
	               EList<Edit>& cedits,
	               EList<Edit>& ccedits,
	               RandomSource& rand);

	TScore decode(const BTDnaString& read,
	              const BTString& qual,
	              size_t readi,
	              size_t readf,
	              BTString& ref,
	              size_t refi,
	              size_t reff,
	              int snpPen,
	              BTDnaString& dseq,
	              EList<Edit>& nedits,
	              EList<Edit>& aedits,
	              EList<Edit>& cedits,
	              EList<Edit>& ccedits,
	              RandomSource& rand);

private:
	//
	// Dynamic programming table; good for colorspace reads up to 1024
	// colors in length.
	//
	// 0 -> A, 1 -> C, 2 -> G, 3 -> T, 4 -> min(A, C, G, T),
	// 5 -> backtrack mask
	EList<UngappedCell> table_;
};

enum {
	GAPPED_BT_DIAG,
	GAPPED_BT_REF_OPEN,
	GAPPED_BT_REF_EXTEND,
	GAPPED_BT_READ_OPEN,
	GAPPED_BT_READ_EXTEND
};

/**
 * A bitmask encoding which backtracking paths out of a particular cell
 * correspond to optimal subpaths.
 */
struct GappedCellMask {

	/**
	 * Set all flags to 0, indicating there is no way to backtrack from
	 * this cell to an optimal answer.
	 */
	void clear() {
		*((uint16_t*)this) = 0;
	}

	/**
	 * Return true iff there are no backward paths recorded in this
	 * mask.
	 */
	bool empty() const {
		return *((uint16_t*)this) == 0;
	}

	/**
	 * Return true iff it's possible to extend a gap in the reference
	 * in the cell below this one.
	 */
	bool refExtendPossible() const {
		return rfop || rfex;
	}

	/**
	 * Return true iff it's possible to open a gap in the reference
	 * in the cell below this one (false implies that only extension
	 * is possible).
	 */
	bool refOpenPossible() const {
		return diag || rdop || rdex;
	}

	/**
	 * Return true iff it's possible to extend a gap in the read
	 * in the cell to the right of this one.
	 */
	bool readExtendPossible() const {
		return rdop || rdex;
	}

	/**
	 * Return true iff it's possible to open a gap in the read in the
	 * cell to the right of this one (false implies that only extension
	 * is possible).
	 */
	bool readOpenPossible() const {
		return diag || rfop || rfex;
	}

	/**
	 * Select a path for backtracking from this cell.  If there is a
	 * tie among eligible paths, break it randomly.  Return value is
	 * a pair where first = a flag indicating the backtrack type (see
	 * enum defining GAPPED_BT_* above), and second = a selection for
	 * the read character for the next row up.  second should be
	 * ignored if the backtrack type is a gap in the read.
	 */
	std::pair<int, int>
	randBacktrack(RandomSource& rand) const
	{
		std::pair<int, int> ret;
		ret.second = -1;
		int i = ((diag != 0) << 0) |
		        ((rfop != 0) << 1) |
		        ((rfex != 0) << 2) |
		        ((rdop != 0) << 3) |
		        ((rdex != 0) << 4);
		ret.first = randFromMask(rand, i);
		assert_lt(ret.first, 5);
		assert(ret.first == GAPPED_BT_DIAG ||
		       ret.first == GAPPED_BT_REF_OPEN ||
		       ret.first == GAPPED_BT_REF_EXTEND ||
		       ret.first == GAPPED_BT_READ_OPEN ||
		       ret.first == GAPPED_BT_READ_EXTEND);
		if(ret.first < 3) {
			// Must choose character for next row
			if(ret.first == GAPPED_BT_DIAG) {
				assert(diag != 0);
				ret.second = randFromMask(rand, diag);
			} else if(ret.first == GAPPED_BT_REF_OPEN) {
				assert(rfop != 0);
				ret.second = randFromMask(rand, rfop);
			} else if(ret.first == GAPPED_BT_REF_EXTEND) {
				assert(rfex != 0);
				ret.second = randFromMask(rand, rfex);
			}
		}
		return ret;
	}

	uint16_t diag     : 4;
	uint16_t rfop     : 4;
	uint16_t rfex     : 4;
	uint16_t rdop     : 1;
	uint16_t rdex     : 1;
	uint16_t reserved : 2;
};

/**
 * A score
 */
struct GappedScore {

	/**
	 * Gapped scores are invalid until proven valid.
	 */
	GappedScore() {
		gaps = 0;
		score = 0;
	}

	/**
	 * Return an invalid GappedScore.
	 */
	static GappedScore INVALID() {
		GappedScore s;
		s.invalidate();
		return s;
	}

	/**
	 * Return true iff gapped score is valid (i.e., represents one or
	 * more paths leading to a valid partial alignment).
	 */
	bool valid() const {
		return gaps != std::numeric_limits<int32_t>::min();
	}

	/**
	 * Make this score invalid (and therefore <= all other scores).
	 */
	void invalidate() {
		gaps = std::numeric_limits<int32_t>::min();
	}

	/**
	 * Return true iff this score is > score o.
	 * Note: An "invalid" score is <= all other scores.
	 */
	bool operator>(const GappedScore& o) const {
		if(!o.valid()) {
			if(!valid()) {
				// both invalid
				return false;
			} else {
				// I'm valid, other is invalid
				return true;
			}
		} else if(!valid()) {
			// I'm invalid, other is valud
			return false;
		}
		if(gaps > o.gaps) return true;
		if(gaps < o.gaps) return false;
		return score > o.score;
	}

	/**
	 * Scores are equal iff they're bitwise equal.
	 */
	bool operator==(const GappedScore& o) const {
		return score == o.score && gaps == o.gaps;
	}

	/**
	 * Return true iff this score is >= score o.
	 */
	bool operator>=(const GappedScore& o) const {
		return operator==(o) || operator>(o);
	}

	/**
	 * Return true iff this score is < score o.
	 */
	bool operator<(const GappedScore& o) const {
		return !operator>=(o);
	}

	/**
	 * Return true iff this score is <= score o.
	 */
	bool operator<=(const GappedScore& o) const {
		return !operator>(o);
	}

	/**
	 * Calculate difference between two GappedScores.
	 */
	GappedScore operator-(const GappedScore& o) const {
		GappedScore s;
		s.gaps = gaps - o.gaps;
		s.score = score - o.score;
		return s;
	}

	/**
	 * Calculate sum of two GappedScores.
	 */
	GappedScore operator+(const GappedScore& o) const {
		GappedScore s;
		s.gaps = gaps + o.gaps;
		s.score = score + o.score;
		return s;
	}

	/**
	 * Add given GappedScore into this one.
	 */
	GappedScore operator+=(const GappedScore& o) {
		gaps += o.gaps;
		score += o.score;
		return (*this);
	}

	/**
	 * Subtract given GappedScore from this one.
	 */
	GappedScore operator-=(const GappedScore& o) {
		gaps -= o.gaps;
		score -= o.score;
		return (*this);
	}

	/**
	 * Calculate difference between two GappedScores.
	 */
	GappedScore operator-(int o) const {
		GappedScore s;
		s.gaps = gaps;
		s.score = score - o;
		return s;
	}

	// # gaps encountered so far, unless that number exceeds the
	// target, in which case the score becomes invalid and therefore <=
	// all other scores
	int32_t gaps;
	// Score accumulated so far (penalties are subtracted starting at 0)
	int32_t score;
};

std::ostream& operator<<(std::ostream& os, const GappedScore& o);

/**
 * Encapsulates all information needed to encode the optimal subproblem
 * at a gapped-decoder cell.  Also encapsulates a few helper functions
 * that update.
 */
struct GappedCell {

	void updateHoriz(
		const GappedCell& left,
		int totGaps,
		int readOpenPen,
		int readExtendPen);

	void updateDiag(
		const GappedCell& uc,
		int snpPhred,
		int refMask,
		int prevColor,
		int prevQual);

	void updateVert(
		const GappedCell& up,
		int totGaps,
		int prevColor,
		int prevQual,
		int refOpenPen,
		int refExtendPen);

	void updateDiagVert(
		const GappedCell& dc,
		const GappedCell& uc,
		int totGaps,
		int snpPhred,
		int refMask,
		int prevColor, // color b/t this row, one above
		int prevQual,
		int refOpenPen,
		int refExtendPen);

	/**
	 * Clear this cell so that it's ready for updates.
	 */
	void clear() {
		// Initially, best scores are all invalid
		best[0] = best[1] = best[2] = best[3] = GappedScore::INVALID();
		// Initially, there's no way to backtrack from this cell
		mask[0].clear();
		mask[1].clear();
		mask[2].clear();
		mask[3].clear();
	}

	/**
	 * Select a 'to' character for backtracking purposes.  Break ties
	 * randomly.
	 */
	int randTo(RandomSource& rand) const {
		int msk = 0;
		GappedScore bst = GappedScore::INVALID();
		for(int i = 0; i < 4; i++) {
			if(best[i] >= bst) {
				if(best[i] > bst) {
					bst = best[i];
					msk = 0;
				}
				msk |= (1 << i);
			}
		}
		return randFromMask(rand, msk);
	}

	// Best incoming score for each 'to' character
	GappedScore best[4];
	// Mask for tied-for-best incoming paths for each 'to' character
	GappedCellMask mask[4];
};

class ColorspaceDecoderGapped {
public:
	GappedScore backtrack(
		const BTDnaString& read,
		const BTString& qual,
		size_t readi,
		size_t readf,
		BTString& ref,
		size_t refi,
		size_t reff,
		int snpPhred,
		GappedScore decodeScore,
		int readGaps,
		int refGaps,
		int readOpenPen,
		int readExtendPen,
		int refOpenPen,
		int refExtendPen,
		int gapBarrier,
		BTDnaString& decoded,
		EList<Edit>& nedits,
		EList<Edit>& aedits,
		EList<Edit>& cedits,
		EList<Edit>& ccedits,
		int lastC,
		RandomSource& rand);

	TScore decode(
		const BTDnaString& read,
		const BTString& qual,
		size_t readi,
		size_t readf,
		BTString& ref, // this will be overwritten
		size_t refi,
		size_t reff,
		int snpPhred,
		int maxCost,
		int readGaps,
		int refGaps,      // # of reference gaps in
		const EList<Edit>& edits, // colorspace-to-colorspace edits
		int readOpenPen,  // penalty for opening a new gap in the read
		int readExtendPen,// penalty for extending a gap in the read
		int refOpenPen,   // penalty for opening a new gap in the reference
		int refExtendPen, // penalty for extending a gap in the reference
		int gapBarrier,   // # bases on either side of alignment that must be gap-free
		BTDnaString& decoded,
		EList<Edit>& nedits,
		EList<Edit>& aedits,
		EList<Edit>& cedits,
		EList<Edit>& ccedits,
		RandomSource& rand);

private:
	EList<EList<GappedCell> > table_;
	SStringFixed<uint16_t> delAllow_;
	SStringFixed<uint16_t> insAllow_;
};

class ColorspaceDecoder {

public:
	TScore decode(
		const BTDnaString& read,
		const BTString& qual,
		size_t readi, // offset of first character within 'read' to consider
		size_t readf, // offset of last char (exclusive) in 'read' to consider
		BTString& ref,    // reference sequence, as masks
		size_t refi,  // offset of first character within 'ref' to consider
		size_t reff,  // offset of last char (exclusive) in 'ref' to consider
		int snpPen,
		int maxCost,
		int readGaps,
		int refGaps,
		const EList<Edit>& edits,
		int readOpenPen,  // penalty for opening a new gap in the read
		int readExtendPen,// penalty for extending a gap in the read
		int refOpenPen,   // penalty for opening a new gap in the reference
		int refExtendPen, // penalty for extending a gap in the reference
		int gapBarrier,
		bool exEnds,
		bool maqRound,
		BTDnaString& dseq,
		BTString& dqual,
		EList<Edit>& nedits,
		EList<Edit>& aedits,
		EList<Edit>& cedits,
		EList<Edit>& ccedits,
		RandomSource& rand);

private:
	ColorspaceDecoderUngapped ungapped_;
	ColorspaceDecoderGapped   gapped_;
};

#endif /* COLOR_DEC_H_ */
