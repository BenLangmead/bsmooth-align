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

#ifndef HIT_SET_H_
#define HIT_SET_H_

#include <algorithm>
#include "assert_helpers.h"
#include "filebuf.h"
#include "edit.h"
#include "alphabet.h"
#include "annot.h"
#include "refmap.h"
#include "ref.h"
#include "read.h"

/**
 * Encapsulates a hit contained within a HitSet that can be
 * (de)serialized to/from FileBufs.  Used for chaining.
 *
 * FIXME: Most of the members are now vestigial and won't chain with
 * anything.
 */
struct HitSetEnt {
	typedef std::pair<uint32_t,uint32_t> U32Pair;

	HitSetEnt() { }

	HitSetEnt& operator= (const HitSetEnt& e) {
		h = e.h;
		horig = e.horig;
		seedoff = e.seedoff;
		orient = e.orient;
		stratum = e.stratum;
		cost = e.cost;
		oms = e.oms;
		mapq = e.mapq;
		clip3 = e.clip3;
		clip5 = e.clip5;
		seq = e.seq;
		cseq = e.cseq;
		qual = e.qual;
		cqual = e.cqual;
		edits = e.edits;
		aedits = e.aedits;
		cedits = e.cedits;
		ccedits = e.ccedits;
		return *this;
	}

	/**
	 * Write binary representation of HitSetEnt to an OutFileBuf.
	 */
	void serialize(OutFileBuf& fb) const {
		fb.writeChars((const char*)&h.first, 4);
		fb.writeChars((const char*)&h.second, 4);
		fb.writeChars((const char*)&seedoff, 2);
		assert_range(0, 3, (int)orient);
		fb.write(orient);
		assert_geq(stratum, 0);
		assert_lt(stratum, 4);
		fb.write(stratum);
		assert_eq(stratum, (cost >> 14));
		fb.writeChars((const char*)&cost, 2);
		fb.writeChars((const char*)&oms, 4);
		fb.writeChars((const char*)&clip3, 2);
		fb.writeChars((const char*)&clip5, 2);
		// Write nucleotide sequence
		assert_eq(seq.length(), qual.length());
		uint16_t seqLen = seq.length();
		fb.writeChars((const char*)&seqLen, 2);
		fb.writeChars((const char*)seq.toZBuf(), seqLen);
		fb.writeChars((const char*)qual.toZBuf(), seqLen);
		// Write color sequence (if present)
		assert_eq(cseq.length(), cqual.length());
		uint16_t cseqLen = cseq.length();
		fb.writeChars((const char*)&cseqLen, 2);
		fb.writeChars((const char*)cseq.toZBuf(), cseqLen);
		fb.writeChars((const char*)cqual.toZBuf(), cseqLen);
		// Write nucleotide edits
		size_t sz = edits.size();
		fb.writeChars((const char*)&sz, 4);
		for(size_t i = 0; i < edits.size(); i++) {
			edits[i].serialize(fb);
		}
		sz = aedits.size();
		fb.writeChars((const char*)&sz, 4);
		for(size_t i = 0; i < aedits.size(); i++) {
			aedits[i].serialize(fb);
		}
		sz = cedits.size();
		fb.writeChars((const char*)&sz, 4);
		for(size_t i = 0; i < cedits.size(); i++) {
			cedits[i].serialize(fb);
		}
		sz = ccedits.size();
		fb.writeChars((const char*)&sz, 4);
		for(size_t i = 0; i < ccedits.size(); i++) {
			ccedits[i].serialize(fb);
		}
	}

	/**
	 * Repopulate a HitSetEnt from its binary representation in FileBuf.
	 */
	void deserialize(FileBuf& fb) {
		fb.get((char*)&h.first, 4);
		fb.get((char*)&h.second, 4);
		fb.get((char*)&seedoff, 2);
		orient = fb.get();
		assert_range(0, 3, (int)orient);
		stratum = fb.get();
		assert_geq(stratum, 0);
		assert_lt(stratum, 4);
		fb.get((char*)&cost, 2);
		assert_eq(stratum, (cost >> 14));
		fb.get((char*)&oms, 4);
		fb.get((char*)&clip3, 2);
		fb.get((char*)&clip5, 2);
		// Read nucleotide and quality sequences
		// sz = length of nucs string
		size_t seqLen;
		fb.get((char*)&seqLen, 2);
		if(seqLen > 0) {
			seq.resize(seqLen);
			fb.get((char*)seq.wbuf(), seqLen);
			qual.resize(seqLen);
			fb.get((char*)qual.wbuf(), seqLen);
		}
		// Read color sequence
		// sz = length of color string
		size_t cseqLen;
		fb.get((char*)&cseqLen, 2);
		if(cseqLen > 0) {
			cseq.resize(cseqLen);
			fb.get((char*)cseq.wbuf(), cseqLen);
			cqual.resize(cseqLen);
			fb.get((char*)cqual.wbuf(), cseqLen);
		}
		// sz = length of edits list
		uint32_t sz = 0;
		fb.get((char*)&sz, 4);
		assert_lt(sz, 1024);
		edits.resize(sz);
		for(uint32_t i = 0; i < sz; i++) {
			edits[i].deserialize(fb);
		}
		// sz = length of aedits list
		fb.get((char*)&sz, 4);
		assert_lt(sz, 1024);
		aedits.resize(sz);
		for(uint32_t i = 0; i < sz; i++) {
			aedits[i].deserialize(fb);
		}
		// sz = length of cedits list
		fb.get((char*)&sz, 4);
		assert_lt(sz, 1024);
		cedits.resize(sz);
		for(uint32_t i = 0; i < sz; i++) {
			cedits[i].deserialize(fb);
		}
		// sz = length of ccedits list
		fb.get((char*)&sz, 4);
		assert_lt(sz, 1024);
		ccedits.resize(sz);
		for(uint32_t i = 0; i < sz; i++) {
			ccedits[i].deserialize(fb);
		}
	}

	/**
	 * Less than operator.  Break HitSetEnt ties by:
	 *  - Stratum, then
	 *  - Cost, then
	 *  - Position, then
	 *  - Orientation
	 */
	int operator< (const HitSetEnt &rhs) const {
		if(stratum < rhs.stratum) return 1;
		if(stratum > rhs.stratum) return 0;
		if(clip3+clip5 < rhs.clip3+rhs.clip5) return 1;
		if(clip3+clip5 > rhs.clip3+rhs.clip5) return 0;
		if(cost < rhs.cost) return 1;
		if(cost > rhs.cost) return 0;
		if(clip3 < rhs.clip3) return 1;
		if(clip3 > rhs.clip3) return 0;
		if(clip5 < rhs.clip5) return 1;
		if(clip5 > rhs.clip5) return 0;
		if(h < rhs.h) return 1;
		if(h > rhs.h) return 0;
		return (orient < rhs.orient)? 1 : 0;
	}

	/**
	 * Greater than operator.
	 */
	int operator> (const HitSetEnt &rhs) const {
		if(stratum < rhs.stratum) return 0;
		if(stratum > rhs.stratum) return 1;
		if(clip3+clip5 < rhs.clip3+rhs.clip5) return 0;
		if(clip3+clip5 > rhs.clip3+rhs.clip5) return 1;
		if(cost < rhs.cost) return 0;
		if(cost > rhs.cost) return 1;
		if(clip3 < rhs.clip3) return 0;
		if(clip3 > rhs.clip3) return 1;
		if(clip5 < rhs.clip5) return 0;
		if(clip5 > rhs.clip5) return 1;
		if(h < rhs.h) return 0;
		if(h > rhs.h) return 1;
		return (orient <= rhs.orient)? 0 : 1;
	}

	/**
	 * Equality comparison operator.
	 */
	bool operator== (const HitSetEnt &rhs) const {
		return(stratum == rhs.stratum &&
		       cost == rhs.cost &&
		       clip3 == rhs.clip3 &&
		       clip5 == rhs.clip5 &&
		       orient == rhs.orient &&
		       h == rhs.h);
	}

	/**
	 * Equality comparison operator.
	 */
	bool operator!= (const HitSetEnt &rhs) const {
		return !(*this == rhs);
	}
	
	/**
	 * Indexing returns edits.
	 */
	Edit& operator[](size_t x) {
		return edits[x];
	}

	/**
	 * Indexing returns edits.
	 */
	const Edit& operator[](size_t x) const {
		return edits[x];
	}

	/**
	 * Another way to get at an edit.
	 */
	Edit& editAt(unsigned i) {
		return edits[i];
	}

	/**
	 * Another way to get at a const edit.
	 */
	const Edit& editAt(unsigned i) const {
		return edits[i];
	}

	/**
	 * Get the ith color edit.
	 */
	Edit& colorEditAt(unsigned i) {
		return cedits[i];
	}

	/**
	 * Another way to get at an edit.
	 */
	const Edit& colorEditAt(unsigned i) const {
		return cedits[i];
	}

	/**
	 * Return the front entry.
	 */
	Edit& front() {
		return edits.front();
	}

	/**
	 * Return the back entry.
	 */
	Edit& back() {
		return edits.back();
	}

	/**
	 * Expand the entry list by one.
	 */
	void expand() {
		edits.resize(edits.size() + 1);
	}

	/**
	 * Sort edits by position
	 */
	void sort() {
		edits.sort(); aedits.sort(); cedits.sort(); ccedits.sort();
	}

	/**
	 * Return number of edits.
	 */
	size_t size() const {
		return edits.size();
	}

	/**
	 * Return true iff there are no elements in the 'edits' array.
	 */
	bool empty() const {
		return edits.empty();
	}
	
	/**
	 * Reset state.
	 */
	void clear() {
		seq.clear();
		cseq.clear();
		qual.clear();
		cqual.clear();
		edits.clear();
		aedits.clear();
		cedits.clear();
		ccedits.clear();
	}
	
	/**
	 * Check that this HitSetEnt is internally consistent.
	 */
	bool repOk() {
		assert_eq(qual.length(), seq.length());
		assert_eq(cseq.length(), cqual.length());
		for(size_t i = 0; i < edits.size(); i++) {
			assert_lt(edits[i].pos, seq.length());
		}
		for(size_t i = 0; i < aedits.size(); i++) {
			assert_lt(aedits[i].pos, seq.length());
		}
		for(size_t i = 0; i < cedits.size(); i++) {
			assert_lt(cedits[i].pos, cseq.length());
		}
		for(size_t i = 0; i < ccedits.size(); i++) {
			assert_lt(ccedits[i].pos, cseq.length());
		}
		return true;
	}

	/**
	 * Write HitSetEnt to an output stream.
	 */
	friend std::ostream& operator << (std::ostream& os, const HitSetEnt& hse);

	U32Pair h;        // reference coordinates
	U32Pair horig;    // reference coordinates for reference where we originally aligned
	uint16_t seedoff; // offset of leftmost seed char from h
	uint8_t orient;   // orientation
	int8_t stratum;   // stratum
	uint16_t cost;    // cost, including stratum (if colorspace, = decoded cost)
	uint16_t ccost;   // color-to-color alignment cost, including stratum
	uint32_t oms;     // # others
	uint32_t mapq;    // # others
	uint16_t clip3;   // bases clipped from 3' end
	uint16_t clip5;   // bases clipped from 5' end
	BTDnaString seq;  // decoded bases for colorspace reads, or original bases for nucleotide reads
	BTDnaString cseq; // original bases for colorspace reads
	BTString qual;    // decoded qualities for colorspace reads, or original qualities for nucleotide reads
	BTString cqual;   // original qualities for colorspace reads
	EList<Edit> edits; // edits to get from reference to subject
	EList<Edit> aedits; // resolutions of ambiguous nucleotides
	EList<Edit> cedits; // color edits to get from reference to subject
	EList<Edit> ccedits; // color edits to get from reference to subject
};

class AlignOutput; // forward decl
class Read; // forward decl

/**
 * Encapsulates a set of hits that can be (de)serialized to/from
 * FileBufs.  Used for chaining.
 */
struct HitSet {

	typedef std::pair<uint32_t,uint32_t> U32Pair;

	HitSet() {
		maxedStratum = -1;
	}

	HitSet(FileBuf& fb) {
		deserialize(fb);
	}

	/**
	 * Copy the contents of h into this HitSet.
	 */
	HitSet& operator= (const HitSet& h) {
		//name = h.name;
		//seq = h.seq;
		//qual = h.qual;
		maxedStratum = h.maxedStratum;
		//color = h.color;
		ents = h.ents;
		return *this;
	}

	/**
	 * Write binary representation of HitSet to an OutFileBuf.
	 */
	void serialize(OutFileBuf& fb) const {
	#if 0
		fb.write(color ? 1 : 0);
		size_t i = name.length();
		assert_gt(i, 0);
		fb.writeChars((const char*)&i, 4);
		fb.writeChars(name.toZBuf(), i);
		i = seq.length();
		assert_gt(i, 0);
		assert_lt(i, 1024);
		fb.writeChars((const char*)&i, 4);
		for(size_t j = 0; j < i; j++) {
			fb.write("ACGTN"[(int)seq[j]]);
		}
		fb.writeChars(qual.toZBuf(), i);
		i = ents.size();
		fb.writeChars((const char*)&i, 4);
		for(size_t j = 0; j < ents.size(); j++) {
			ents[j].serialize(fb);
		}
		fb.write(maxedStratum);
	#endif
	}

	/**
	 * Repopulate a HitSet from its binary representation in FileBuf.
	 */
	void deserialize(FileBuf& fb) {
	#if 0
		color = (fb.get() != 0 ? true : false);
		uint32_t sz = 0;
		if(fb.get((char*)&sz, 4) != 4) {
			name.clear();
			seq.clear();
			return;
		}
		assert_gt(sz, 0);
		assert_lt(sz, 1024);
		name.resize(sz);
		fb.get(name.wbuf(), sz);
		fb.get((char*)&sz, 4);
		assert_gt(sz, 0);
		assert_lt(sz, 1024);
		seq.resize(sz);
		fb.get(seq.wbuf(), sz);
		qual.resize(sz);
		fb.get(qual.wbuf(), sz);
		fb.get((char*)&sz, 4);
		if(sz > 0) {
			ents.resize(sz);
			for(size_t i = 0; i < sz; i++) {
				ents[i].deserialize(fb);
			}
		} else {
			ents.clear();
		}
		maxedStratum = fb.get();
	#endif
	}

	/**
	 * Return true iff this HitSet is initialized with a read (but not
	 * necessarily any alignments).
	 */
	//bool initialized() const {
	//	return seq.length() != 0;
	//}

	/**
	 * Return true iff this HitSet has no hits.
	 */
	bool empty() const {
		return ents.empty();
	}

	/**
	 * Return number of entries in this HitSet.
	 */
	size_t size() const {
		return ents.size();
	}

	/**
	 * Remove all entries from this HitSet.
	 */
	void clear() {
		maxedStratum = -1;
		ents.clear();
	}
	
	void reset() { clear(); }

	/**
	 * Return the front entry.
	 */
	HitSetEnt& front() {
		assert(!ents.empty());
		return ents.front();
	}

	/**
	 * Return the back entry.
	 */
	HitSetEnt& back() {
		assert(!ents.empty());
		return ents.back();
	}

	/**
	 * Remove the backmost entry.
	 */
	void pop_back() {
		assert(!ents.empty());
		ents.pop_back();
	}
	
	/**
	 * Expand the entry list by one.
	 */
	void expand() {
		ents.resize(ents.size() + 1);
		ents.back().clear();
		assert(ents.back().empty());
	}

	/**
	 * Resize entry list
	 */
	void resize(size_t sz) {
		ents.resize(sz);
	}

	/**
	 * Sort hits by stratum/penalty.
	 */
	void sort() { ents.sort(); }

	/**
	 * Return true iff Ents are sorted
	 */
	bool sorted() const {
		if(ents.empty()) return true;
		for(size_t i = 0; i < ents.size()-1; i++) {
			if(!(ents[i] < ents[i+1])) return false;
		}
		return true;
	}

	/**
	 * Remove the ith hit from the HitSet, shifting all subsequent hits
	 * up by one.
	 */
	void remove(size_t i) {
		ents.remove(i);
	}

	/**
	 * Return true if strata are uniform across hits; assert if they're
	 * not.
	 */
	bool uniformStrata() const {
		for(size_t i = 0; i < ents.size()-1; i++) {
			assert(ents[i].stratum == ents[i+1].stratum);
		}
		return true;
	}

	/**
	 * Indexing returns entries.
	 */
	HitSetEnt& operator[](size_t x) {
		return ents[x];
	}

	/**
	 * Indexing returns entries.
	 */
	const HitSetEnt& operator[](size_t x) const {
		return ents[x];
	}

	/**
	 * Apply a reference mappings to all the contained hits.
	 */
	void applyReferenceMap(const ReferenceMap& map) {
		for(size_t i = 0; i < ents.size(); i++) map.map(ents[i].h);
	}

	/**
	 * Clear out all the strings and all the entries.
	 */
	void clearAll() {
		//name.clear();
		//seq.clear();
		//qual.clear();
		ents.clear();
		maxedStratum = -1;
		//color = false;
	}

	/**
	 * Given a new hit, including its position, orientation, cost, etc,
	 * determine whether it is strictly better than an existing hit
	 * and, if so, replace the existing hit with the new one.
	 *
	 * BTL: This is only really used for chaining; Merman itself is not
	 * going to arrive at two different alignments for the same seed
	 * hit.  It's pretty much vestigial, hence the throw 1 at the
	 * beginning.
	 */
	bool tryReplacing(
		U32Pair h,
		bool orient,
		uint16_t cost,
		uint32_t oms,
		uint32_t mapq,
		uint16_t seedoff,
		uint16_t clip3,
		uint16_t clip5,
		size_t& pos)
	{
		throw 1;
#if 0
		uint64_t newSeedId =
			((uint64_t)h.second + seedoff) |
			 (uint64_t)h.first << 32 | (fw ? (1llu << 63) : 0);
		
		for(size_t i = 0; i < ents.size(); i++) {
			uint64_t seedId = ((uint64_t)ents[i].h.second + ents[i].seedoff) |
			                  (uint64_t)ents[i].h.first << 32 |
			                  (ents[i].fw ? (1llu << 63) : 0);
			if(seedId == newSeedId) {
				if(clip3 + clip5 < ents[i].clip3 + ents[i].clip5 ||
				   (clip3 + clip5 == ents[i].clip3 + ents[i].clip5 &&
				    cost < ents[i].cost))
				{
					// New hit at same position is better!  Replace it.
					ents[i].h       = h;
					ents[i].seedoff = seedoff;
					ents[i].fw      = fw;
					ents[i].stratum = cost >> 14;
					ents[i].cost    = cost;
					ents[i].oms     = oms;
					ents[i].mapq    = 40;
					ents[i].clip3   = clip3;
					ents[i].clip5   = clip5;
					pos = i;
					return true;
				} else {
					pos = 0xffffffff;
					return true;
				}
			}
		}
#endif
		return false;
	}

	/**
	 * Report up to 'khits' hits from this HitSet.
	 */
	void reportUpTo(
		const Read& rd,
		AlignOutput& os,
		int khits,
		const ReferenceSet& refs,
		const ReferenceMap *rmap,
		bool xformCoords,
		const AnnotationMap *amap);

	/**
	 * Print this HitSet.
	 */
	friend std::ostream& operator << (std::ostream& os, const HitSet& hs);

	//BTString         name;   // name
	//BTDnaString      seq; // original sequence
	//BTString         qual;   // qualities
	int8_t           maxedStratum;
	EList<HitSetEnt> ents;
	//bool             color; // whether read was orginally in colorspace
};

#endif /* HIT_SET_H_ */
