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

#ifndef REF_H_
#define REF_H_

#include "alphabet.h"
#include "ds.h"
#include "sstring.h"

typedef SStringExpandable<char, 256> TRefStr;

/**
 * Encapsulates a nucleotide reference string and its corresponding
 * colorspace sequence.  Member function typically take a 'bool color'
 * parameter, which the client uses to specify whether the query should
 * be resolved with respect to the colorspace version.
 */
class RefString {

public:
	/**
	 * Return the length of this reference string in either color or
	 * nucleotides.  There's one less color than there are nucleotides.
	 */
	size_t length(bool color) const {
		assert_gt(seq.length(), 0);
		if(color) return seq.length()-1;
		else return seq.length();
	}

	/**
	 * Get nucleotide or color at position i.  If prevChar is specified
	 * and the character at position 'i' is ambiguous, prevChar is the
	 * unambiguous character it resolved to.  resolve the phasing
	 * between an ambiguous color and its successor
	 */
	char charAt(size_t i, bool color, int prevCharI = -1, int prevCharIp = -1) const {
		assert_lt(i, length(color));
		if(color) {
			int m1;
			if(prevCharI != -1) {
				assert(isDna(prevCharI));
				m1 = 1 << asc2dna[prevCharI];
				assert(isAmbigNuc(seq[i]));
			} else {
				m1 = asc2dnamask[(int)seq[i]];
			}
			int m2;
			if(prevCharIp != -1) {
				assert(isDna(prevCharIp));
				m2 = 1 << asc2dna[prevCharIp];
				assert(isAmbigNuc(seq[i+1]));
			} else {
				m2 = asc2dnamask[(int)seq[i+1]];
			}
			assert_gt(m1, 0);
			assert_lt(m1, 16);
			assert_gt(m2, 0);
			assert_lt(m2, 16);
			int col = dnamasks2colormask[m1][m2];
			assert_gt(col, 0);
			assert_lt(col, 16);
#ifndef NDEBUG
			if(mask2popcnt[m1] == 1 && mask2popcnt[m2] == 1) {
				int c1 = asc2dna[(int)(prevCharI  == -1 ? seq[i]   : prevCharI)];
				int	c2 = asc2dna[(int)(prevCharIp == -1 ? seq[i+1] : prevCharIp)];
				assert_leq(c1, 4);
				assert_leq(c2, 4);
				int col2 = 1 << dinuc2color[c1][c2];
				assert_eq(col, col2);
			}
#endif
			return mask2dna[col];
		} else {
			// Nucleotide case is easy
			return seq[i];
		}
	}

	/**
	 * Return true iff the nucleotide or color at position i is
	 * unmatchable.  If color is true, then we consider the color
	 * unmatchable if either of the overlapping nucleotides is
	 * unmatchable.
	 */
	bool isUnmatchableAt(size_t i, bool color) const {
		if(color) {
			return isUnmatchableNuc(seq[i]) || isUnmatchableNuc(seq[i+1]);
		} else {
			return isUnmatchableNuc(seq[i]);
		}
	}

	/**
	 * Return true iff the nucleotide or color at position i is
	 * ambiguous.  If color is true, then we consider the color
	 * ambiguous if either of the overlapping nucleotides is ambiguous.
	 */
	bool isAmbigAt(size_t i, bool color) const {
		if(color) {
			return isAmbigNuc(seq[i]) || isAmbigNuc(seq[i+1]);
		} else {
			return isAmbigNuc(seq[i]);
		}
	}

	/**
	 * Convert a stretch of the reference to a reference "mask".  Abort
	 * and return false if any of the reference characters encountered
	 * were unmatchable.
	 */
	bool toRefMask(
		BTString& rf,  // install result here (clobbers existing content)
		size_t refi,   // offset of first reference nucleotide to extract
		size_t reff,   // offset just past the last reference nucleotide to extract
		bool fw) const // false -> reference stretch is reverse complemented first
	{
		size_t reflen = reff-refi;
		rf.resize(reflen);
		for(size_t i = 0; i < reflen; i++) {
			char c = toupper(charAt(refi+i, false));
			if(isUnmatchableNuc(c)) {
				// Bail
				return false;
			}
			size_t ii = i;
			if(!fw) {
				c = compDna(c);
				ii = reflen - i - 1;
			}
			rf.set(asc2dnamask[(int)c], ii);
			assert_neq(0, asc2dnamask[(int)c]);
		}
		return true;
	}

	/**
	 * Apply or un-apply sequence transformations stored in 'xforms'
	 * field (e.g. CG -> YG for bisulfite alignment)
	 */
	void toggleXforms() {
		for(size_t i = 0; i < xforms_.size(); i++) {
			char tmp = seq[xforms_[i].first];
			seq.set(xforms_[i].second, xforms_[i].first);
			xforms_[i].second = tmp;
		}
	}

	/**
	 * Add a sequence transformation such that it can be undone.  The
	 * caller should pass the original character, not the new
	 * character.
	 */
	void addXform(uint32_t off, char old) {
		xforms_.push_back(std::make_pair(off, old));
	}
	
	/**
	 * Return the list of transformations associated with this
	 * RefString.
	 */
	const EList<std::pair<uint32_t, char> > getXforms() {
		return xforms_;
	}

	TRefStr seq;  // sequence

protected:

	EList<std::pair<uint32_t, char> > xforms_;
};

/**
 * Parameters regarding the reference strings.
 */
struct ReferenceParams {

	void init(
		bool reqa,
		bool bisc,
		bool biscpg,
		bool gencr,
		bool genrc,
		bool watcrirc,
		int ent1,
		int ent2)
	{
		requireAnnot = reqa;
		bisulfiteC = bisc;
		bisulfiteCpG = biscpg;
		genCrick = gencr;
		genRevcomps = genrc;
		watsonCrickRc = watcrirc;
		entThresh.first = ent1;
		entThresh.second = ent2;
	}

	bool requireAnnot;  // Only keep portions of the reference where a read
	                    // could overlap an annotation.
	bool bisulfiteC;    // bisulfite-treat references s.t. every C becomes a Y
	bool bisulfiteCpG;  // bisulfite-treat references s.t. every CpG C becomes
	                    // a Y and every non-CpG C becomes a T
	// See also: AlignParams::alignfw and AlignParams::alignrc
	bool genCrick;      // true -> don't align to the watson strand.  This
	                    // means the same as noreadfw in non-bisulfite
						// alignment modes and in bisulfite 2-strand mode, but
						// in the 4-strand bisulfite mode this suppresses
						// alignments to either the BSW or BSWR strands.
	bool genRevcomps;   // When building and indexing reference, take all the
	                    // input strands (including the generated Crick strand
						// if applicable)
	bool watsonCrickRc; // true iff the transformed Watson and Crick strands
	                    // are still reverse complementary.  This is useful for
						// sanity checking.  E.g. if this is true,
						// ReferenceParams::genCrick and AlignParams::alignrc
						// should not both be true; if they are then we are
						// probably double-reporting alignments.
	std::pair<int, int> entThresh; // Entropy-in-window constraints
};

/**
 * Encapsulates a reference sequence.
 */
struct Reference {
	uint32_t  off;         // Offset into global space of concatenated refs
	uint32_t  idx;         // Index in the global list of References
	RefString seq;         // Sequence
	TRefStr   name;        // Name
	bool      crick;       // Is this a crick strand?  I.e. C, BC or BCR?
	bool      rc;          // Is this a reverse-comped strand?  I.e. WR?
	uint32_t  fwWatsonIdx; // Index of forward Watson version of this reference
	bool      index;       // Whether to index this sequence
};

/**
 * Encapsulates the set of all references being aligned against.
 */
class ReferenceSet {

public:

	ReferenceSet() :
		refs_(),
		offToRef_(),
		totRefLen_(),
		entBuf_()
	{ }

	/**
	 * Return the number of reference sequences indexed.
	 */
	size_t numRefs() const { return refs_.size(); }
	
	/**
	 * Get the most recently added Reference.
	 */
	const Reference& back() const {
		assert(!refs_.empty());
		return refs_.back();
	}

	/**
	 * Get the name of the ith reference sequence.
	 */
	const TRefStr& getName(size_t i) const {
		assert_lt(i, refs_.size());
		return refs_[i].name;
	}

	/**
	 * Read-only accessor for Reference objects.
	 */
	const Reference& operator[](size_t i) const {
		assert_lt(i, refs_.size());
		return refs_[i];
	}

	/**
	 * Clear out all References.
	 */
	void clear() {
		refs_.clear();
	}

	/**
	 * Return number of reference strings contained.
	 */
	size_t size() const {
		return refs_.size();
	}

	/**
	 * Return true iff the set of reference strings is empty.
	 */
	bool empty() const {
		return size() == 0;
	}

	/**
	 * Add the contents of a FASTA file as index sequences.
	 */
	void addOrigReferenceFasta(const char *s, const ReferenceParams& p);

	/**
	 * Add a string as an index sequences.
	 */
	void addOrigReferenceString(const char *s, const ReferenceParams& p);

	/**
	 * Re-add all the index sequences as their reverse complements
	 * either before or after any sequence transformations were applied
	 * (e.g. bisulfite transformations).
	 */
	void addReferenceRevComps(const ReferenceParams& p, bool postXforms, int crick, int rc);

	/**
	 * Remove all the forward references; useful when --nofw and
	 * --rcref are both specified.
	 */
	void removeWatsonOrCrickReferences(bool watson = true);

	/**
	 * Return the index of the reference that contains the given global
	 * offset.
	 */
	const Reference& refWithOff(size_t off) const {
		assert(!offToRef_.empty());
		size_t idx = offToRef_.leqBound(off);
		assert_gt(idx, 0);
		assert_leq(idx, refs_.size());
		const Reference *ref = &refs_[idx-1];
		assert_lt(ref->idx, refs_.size());
		assert_lt(ref->off, totRefLen_);
		assert_leq(ref->off, off);
		assert_gt(ref->off + ref->seq.length(false), off);
		return *ref;
	}
	
	/**
	 * Check that this ReferenceSet is internally consistent.
	 */
	bool repOk() const {
		for(size_t i = 0; i < refs_.size(); i++) {
			assert_lt(refs_[i].idx, refs_.size());
			assert_lt(refs_[i].off, totRefLen_);
		}
		return true;
	}

protected:

	/**
	 * Re-add all the index sequences as their reverse complements.
	 */
	void addReference(const ReferenceParams& p);

	EList<Reference> refs_;
	ELMap<size_t, Reference*> offToRef_;
	size_t totRefLen_;
	EList<int> entBuf_;
};

#endif /* REF_H_ */
