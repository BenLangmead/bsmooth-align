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

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <limits>
#include <set>
#include "mer_index.h"
#include "qual.h"
#include "filebuf.h"
#include "read.h"
#include "color.h"
#include "parsort.h"
#include "color_dec.h"
#include "ds.h"
#include "globals.h"
#include "aligner.h"
#include "check.h"
#include "orient.h"

using namespace std;
typedef std::pair<uint32_t, uint32_t> U32Pair;

/**
 * Append a DNA character to a 64-bit buffer.  Convert it to two bits
 * first.
 */
static inline void appendChar(char c, uint64_t& mer, bool comp = false) {
	if(!isDna(c)) c = 3;
	mer <<= 2llu;
	int ic = asc2dna[(int)c];
	assert(ic >= 0 && ic < 4);
	mer |= (uint64_t)(comp ? (ic ^ 3) : ic);
}

/**
 * Append a DNA character to a 64-bit buffer.
 */
static inline void append(int c, uint64_t& mer, bool comp = false) {
	assert(c >= 0 && c < 4);
	mer <<= 2llu;
	mer |= comp ? c^3 : c;
}

/**
 * Find the first element of mers_ that is not < the query.  This
 * is hand-coded because it's performance-critical.
 */
const mer_ent_sm* MerIndex::lowerBound(
	uint16_t key1,
	uint16_t key2,
	uint8_t key3,
	int seed,
	mer_ent_sm*& end) const
{
	assert(sorted_);
	if(mers_[seed][key1] == NULL) {
		end = NULL;
		return NULL;
	}
	size_t top = 0;
	size_t bot = mersLen_[seed][key1];
	end = &mers_[seed][key1][bot];
	while(bot > top) {
		size_t mid = top + ((bot - top) >> 1);
		if(mers_[seed][key1][mid].key < key2 ||
		   (mers_[seed][key1][mid].key == key2 &&
		    mers_[seed][key1][mid].key2 < key3))
		{
			// mers_[mid] < query
			if(top == mid) {
				top++; break;
			}
			top = mid;
		} else {
			bot = mid;
		}
	}
	const mer_ent_sm* ret = &mers_[seed][key1][top];
	if(ret == end || ret->key != key2 || ret->key2 != key3) ret = NULL;
	return ret;
}

/**
 * Align the read against the indexed, IUPAC-containing reference.
 */
void MerIndex::queryHelper(
	Read& r,
	const ReferenceSet& refs,
	const ReferenceMap* rmap,
	const AnnotationMap* amap,
	HitSet& hits,
	AlignOutput& os,
	AlignResult& res,
	const AlignParams& pa,
	bool randomize,
	RandomSource& rnd,
	int threadid)
{
	assert_geq(threadid, 0);
	assert_lt((size_t)threadid, threadState_.size());
	MerIndexPerThreadState& st = threadState_[threadid];
	assert(!empty());
	assert(!refs.empty());
	assert(sorted_);
	assert(pa.adjMinLen == -1 || pa.adjMinLen >= pa.seedLen);
	assert(pa.adjMinLen == -1 || pa.adjMinLen >= pa.iSeedLen);
	int e2eMms = pa.e2eMms;
	int seedMms = pa.seedMms;
	int color = r.color ? 1 : 0;
	int adjMinLen = pa.adjMinLen;
	if(adjMinLen == -1) adjMinLen = (int)r.seq.length();
	assert_geq((int)r.seq.length(), seedWidth_);
	// Check whether this read was already completely maxed out by an
	// upstream tool.
	if(hits.maxedStratum == 0) {
		res.maxed = true;
		os.printMaxedRead(r, r.seq, r.qual, hits.maxedStratum, pa.mhits);
		return;
	}
	// Does arrangement of Ns in the read disqualifiy it from aligning?
	if(nFilter(r, adjMinLen, seedMms, e2eMms)) {
		if(verbose)
			cout << "  filtered " << r.name << " due to Ns" << endl;
		os.printUnalignedRead(r, r.seq, r.qual, 0 /* not filtered */);
		return;
	}
	bool e2e = (e2eMms != INT_MAX);
	int bestStrat = 0;
	int worstStrat = seedMms;
	int maxHitCost = INT_MAX, minHitCost = INT_MAX;
	bool atCapacity = (pa.mhits == INT_MAX && (int)hits.size() == pa.khits);
	assert(pa.mhits == INT_MAX || (int)hits.size() <= pa.khits);
	// Establish bounds
	if(!hits.empty()) {
#ifndef NDEBUG
		bool first = true;
		int lastStrat = -1;
		int nhits = 0;
		uint16_t lastCost;
		for(size_t i = 0; i < hits.size(); i++) {
			if(first) {
				first = false;
				lastStrat = hits[i].stratum;
				lastCost = hits[i].cost;
			}
			nhits++;
			assert_leq(nhits, pa.mhits);
			assert_geq(hits[i].cost, lastCost);
			if(pa.strata) assert_eq(lastStrat, hits[i].stratum);
			lastStrat = hits[i].stratum;
		}
#endif
		if(pa.strata || atCapacity) {
			worstStrat = hits.back().stratum;
		}
		minHitCost = hits.front().cost;
		maxHitCost = hits.back().cost;
	}
	if(hits.maxedStratum > 0) {
		worstStrat = min(worstStrat, hits.maxedStratum-1);
	}
	if(e2e) {
		e2eMms = min(e2eMms, worstStrat);
	} else {
		seedMms = min(seedMms, worstStrat);
	}
	int64_t seedHits = 0;
	// Extract relevant subsequences from Read
	extractReadmers(r, st.readmers, st.readmerstrs, !pa.alignfw, !pa.alignrc);
	const size_t readmersSz = st.readmers.size();
	if(readmersSz == 0) {
		if(verbose)
			cout << "  filtered " << r.name << " due to lack of readmers" << endl;
		os.printUnalignedRead(r, r.seq, r.qual, 0 /* not filtered */);
		return;
	}
	// Map from reference offsets to the number of seed hits occurring
	// at that offset
	int stratThreshs[] = {
		sset_->thresh(0), sset_->thresh(1),
		sset_->thresh(2), sset_->thresh(3),
		sset_->thresh(4), sset_->thresh(5),
		sset_->thresh(6), sset_->thresh(7),
		sset_->thresh(8)
	};
#ifndef NDEBUG
	for(int i = 1; i < 9; i++) {
		assert_leq(stratThreshs[i], stratThreshs[i-1]);
	}
#endif
	int seedsHit = 0;
	assert_leq(readmersSz, st.merItsSz);
	for(size_t mi = 0; mi < readmersSz; mi++) {
		const ReadMer& m = st.readmers[mi];
		assert(m.repOk());
		st.merIts[mi] = NULL; // seeds are misses until proven hits
		st.merItEnds[mi] = NULL;
		assert(m.fw  || pa.alignrc);
		assert(!m.fw || pa.alignfw);
		// st.merIts[mi] now points to the first element of the reference mer
		// list mers_[m.s] that is not less than the query mer m.mer.  
		st.merIts[mi] = lowerBound(m.mer, m.s, st.merItEnds[mi]);
		if(st.merIts[mi] == NULL) continue;
		assert(st.merItEnds[mi] != NULL);
		seedsHit++; seedHits++;
	}
	// Can we achieve even the worst permitted alignment?
	if(seedsHit < stratThreshs[worstStrat]) {
		// No; bail
		if(verbose)
			cout << "  filtered " << r.name << " due to lack of seed hits" << endl;
		os.printUnalignedRead(r, r.seq, r.qual, 0 /* not filtered */);
		return;
	}
	// Sort seed hits into full per-stratum hits for all strata under
	// consideration
	for(int i = 0; i < 7; i++) {
		st.stratHits[i].clear();
	}
	ASSERT_ONLY(uint64_t lastMinValue = std::numeric_limits<uint64_t>::min());
	while(seedsHit > 0) {
		int numTied = 0;
		uint64_t minValue = std::numeric_limits<uint64_t>::max();
		int minI = -1;
		ASSERT_ONLY(set<int> seedsTied);
		// Get the next seed hit for each seed shape.  For each hit, calculate
		// where it is.  Across all hits, calculate which is the furthest
		// upstream on the Watson strand.  That is the next place we'll extend.
		for(size_t mi = 0; mi < readmersSz; mi++) {
			if(st.merIts[mi] == NULL) {
				// No hits for seed of this shape
				continue;
			}
			const ReadMer& m = st.readmers[mi];
			// Skip over any entries where the offset into the reference is
			// such that if falls off the beginning of the refernce.  We could
			// filter these out later, but this also avoids
			// (st.merIts[mi]->value - m.i) yielding a negative number, which
			// is troublesome.
			while(st.merIts[mi]->value < (size_t)m.i) {
				uint32_t mmer = (uint32_t)(m.mer >> 16);
				st.merIts[mi]++;
				if(st.merIts[mi]->getKey() != mmer ||
				   st.merIts[mi] == st.merItEnds[mi])
				{
					st.merIts[mi] = st.merItEnds[mi] = NULL;
					seedsHit--;
					break;
				}
			}
			if(st.merIts[mi] == NULL) {
				// See if statement in previous loop for how this happens
				continue;
			}
			st.merItVals[mi] = (st.merIts[mi]->value - m.i) << 1llu;
			if(!m.fw) st.merItVals[mi] |= 1llu;
			assert_lt(st.merItVals[mi], 0x200000000llu);
			
			if(st.merItVals[mi] < minValue) {
				minValue = st.merItVals[mi];
				assert(lastMinValue == 0 || minValue > lastMinValue);
				minI = m.i;
				numTied = 1;
				ASSERT_ONLY(seedsTied.clear());
				ASSERT_ONLY(seedsTied.insert(m.s));
			} else if(st.merItVals[mi] == minValue) {
				assert_eq(m.i, minI);
				numTied++;
				assert(seedsTied.find(m.s) == seedsTied.end());
			}
		}
		// No winner?  We must be done.
		if(minValue == std::numeric_limits<uint64_t>::max()) {
			assert_eq(0, seedsHit);
			break; // Stop trying to find extensions
			// Break from while(seedsHit > 0) loop
		}
		assert_lt(minValue, 0x200000000llu);
		ASSERT_ONLY(lastMinValue = minValue);
		// Advance all the iterators that contributed to this hit.
		// TODO: Can't we avoid this loop in the case where numTied==1 or is
		// otherwise small?
		for(size_t mi = 0; mi < readmersSz; mi++) {
			if(st.merIts[mi] != NULL) {
				assert_geq(st.merItVals[mi], minValue);
				if(st.merItVals[mi] == minValue) {
					// This seed shape contributed to the hit
					const ReadMer& m = st.readmers[mi];
					uint32_t mmer = (uint32_t)(m.mer >> 16);
					st.merIts[mi]++;
					// Now advance down the the st.merIts[mi] list
					if(st.merIts[mi]->getKey() != mmer ||
					   st.merIts[mi] == st.merItEnds[mi])
					{
						st.merIts[mi] = st.merItEnds[mi] = NULL;
						seedsHit--;
						// TODO: bail at this point?
					} else {
#ifndef NDEBUG
						const ReadMer& mm = st.readmers[mi];
						// Check that the next hit for this readmer isn't
						// the same as the last
						uint64_t tmp = (st.merIts[mi]->value - mm.i) << 1llu;
						if(!mm.fw) tmp |= 1llu;
						assert_gt(tmp, minValue);
#endif
						seedHits++;
					}
				}
			}
		}
#ifndef NDEBUG
		bool inSomeStratum = false;
		for(int i = 0; i < 9; i++) {
			if(numTied == stratThreshs[i]) {
				inSomeStratum = true;
				break;
			}
		}
		assert(color || inSomeStratum);
#endif
		// Iterate through strata
		for(int j = bestStrat; j <= worstStrat; j++) {
			if(numTied >= stratThreshs[j]) {
#if 0
				for(size_t k = 0; k < st.stratHits[j].size(); k++) {
					assert_neq((uint64_t)st.stratHits[j][k], minValue);
				}
#endif
				// Commit this hit to the list of hits of the given stratum.
				// 'minValue' contains the reference offset w/r/t the Watson
				// reference strand.
				st.stratHits[j].push_back(minValue);
				break;
			}
		}
#if 0
		std::map<uint64_t, int> vals;
		for(size_t i = 0; i < 7; i++) {
			for(size_t j = 0; j < st.stratHits[i].size(); j++) {
				uint64_t val = st.stratHits[i][j];
				if(vals.find(val) != vals.end()) {
					cerr << "Hit " << val << " appeared multiple times in stratHits" << endl;
					cerr << "First in stratum " << vals[val] << endl;
					cerr << "Then in stratum " << i << endl;
					throw 1;
				}
				vals[val] = (int)i;
			}
		}
#endif
		if(seedsHit < stratThreshs[worstStrat]) {
			// Stop looking for extensions
			break;
		}
	}
	unsigned seedOverhang = max<int>((int)r.seq.length() - iSeedLen_, 0);
	unsigned minSeedOverhang = max<int>(adjMinLen - iSeedLen_, 0);
	// Go from "best" hits (hit by the most shapes) to worst (hit by the fewest
	// shapes)
	for(int strat = bestStrat; strat <= worstStrat; strat++) {
		bool done = false;
		// For a given stratum (# of shapes hit), visit the hits in a random
		// order.
		const size_t stratSz = st.stratHits[strat].size();
		if(stratSz == 0) {
			continue;
		}
		assert_gt(stratSz, 0);
		if(randomize) {
			// Optionally randomize the order in which we visit the hits
			for(size_t i = stratSz; i > 1; i--) {
				size_t dst = rnd.nextU32() % i;
				if(dst != i-1) {
					swap(st.stratHits[strat][i-1], st.stratHits[strat][dst]);
				}
			}
		}
		for(size_t it = 0; it < stratSz; it++) {
			uint64_t val = st.stratHits[strat][it];
			assert_lt(val, 0x200000000llu);
			const bool fw = ((val & 1llu) == 0);
			const uint32_t off = (uint32_t)(val >> 1llu);
			// Look up the reference info
			const Reference& ref = refs.refWithOff(off);
			int orient = (ref.crick ? ORIENT_CRICK : ORIENT_WATSON) |
			             (fw        ? ORIENT_FW    : ORIENT_RC);
			int fpl = is5pLeft(orient);
			// Calculate the seed's leftmost character's offset
			// into the reference sequence
			uint32_t roff = off - ref.off;
			if(roff + iSeedLen_ > ref.seq.length(color)) {
				// Seed falls off end
				continue;
			}
			// roff = offset of leftmost base on forward reference
			// strand involved in this hypothetical alignment
			int adjAllen = (int)r.seq.length();
			if(!fw) {
				if(roff < minSeedOverhang) continue; // fell off LHS
				if(roff < seedOverhang) {
					// Overhang doesn't entirely fit, but minimum fits
					adjAllen -= (seedOverhang - roff);
					assert_geq(adjAllen, iSeedLen_);
					roff = 0;
				} else {
					// Overhang fits
					roff -= seedOverhang;
				}
			} else {
				if(roff + r.seq.length() > ref.seq.length(color)) {
					if(roff + adjMinLen > ref.seq.length(color)) {
						// Minimum overhang doesn't fit
						continue;
					}
					// Overhang doesn't entirely fit, but minimum fits
					adjAllen -= (int)((roff + r.seq.length()) - ref.seq.length(color));
				}
			}
			int allen = adjAllen + (color ? 1 : 0);
			// 'allen' and 'roff' have now been adjusted based on the
			// reference sequence boundaries.  We might adjust them
			// again after alignment.
			st.edits.clear();
			st.aedits.clear();
			st.cedits.clear();
			st.ccedits.clear();
			int qualPen = 0, mms = 0, smms = 0;
			U32Pair p = make_pair(ref.idx, roff /* might change */);
			if(rmap != NULL) rmap->map(p);
			roff = p.second;
			AnnotationMap::Iter ait;
			bool overlapsAnnots = false; // Does it overlap annotations?
			bool compatAnnots = false; // Are there any compatible annots?
			// Copy annotations that overlap candidate alignment
			// into annots[] array.  5' end is in annots[0].
			char annots[1024];
			memset(annots, 0, allen);
			if(pa.requireAnnot && !overlapsAnnots) continue; // Reject this seed placement
			bool good = false;
			
			// A note re: ambiguous nucleotides that lead to ambiguous colors.
			// If each position is treated independantly and the position's
			// color is the colorspace mask encodoing the union of all possible
			// nucleotide pairings, we have the additional problem of "phasing"
			// the decision made at the first affected color and the second
			// affected color.  This is simple if ambiguous nucleotides are
			// never adjacent (and we handle this), but harder if they're
			// adjacent (we print an error and abort).
			
			// Iterate from 5' to 3' end of read
			int j;
			char prevChar = -1;
			for(j = 0; j < adjAllen; j++) {
				char annot = annots[j];
				size_t refj;
				if(fw) refj = p.second + j;
				else   refj = p.second + (adjAllen - j - 1);
				// Get reference character; re-encode color as nuc if necessary
				int refc = toupper(
					ref.seq.charAt(
						refj,
						color,
						fw ? prevChar : -1,
						fw ? -1 : prevChar));
				prevChar = -1;
				int refcfw = (ref.crick && !color) ? compDna(refc) : refc;
				int qc = toupper(r.charAt(j));
				int qcfw = (!fpl && !color) ? compDna(qc) : qc;
				bool refUnambig = isUnambigNuc(refc);
				if(isUnmatchableNuc(refc) ||
				   (isAmbigNuc(refc) && pa.disallowIupac))
				{
					// Overlaps unmatchable character in the reference
					break; // for(j = 0; j < adjAllen; j++)
				}
				assert(qc == 'N' || isUnambigNuc(qc));
				bool match = false;
				if(refUnambig && annot == 0) {
					// Unambiguous, unannotated reference position
					// versus unambiguous read position
					match = ((!fw && !color) ? compDna(refc) : refc) == qc;
				} else if(isUnambigNuc(qc)) {
					// Either ambiguous or annotated (or both)
					// reference position versus unambiguous read
					// position
					match = ambigCompat(refc, ((fw || color) ? qc : compDna(qc)), color);
					if(match && isAmbigNuc(refc) && color) {
						int refa = toupper(ref.seq.charAt(refj + (fw ? 1 : 0), false));
						if(isAmbigNuc(refa)) {
							int refu = toupper(ref.seq.charAt(refj + (fw ? 0 : 1), false));
							assert(isAmbigNuc(refa));
							assert(isUnambigNuc(refu));
							int prevCharN = nuccol2nuc[asc2dna[refu]][asc2dna[qc]];
							assert_lt(prevCharN, 4);
							assert_geq(prevCharN, 0);
							prevChar = "ACGT"[prevCharN];
						}
					}
					if(isAmbigNuc(refc) && !color) {
						if(match) {
							st.aedits.push_back(Edit(j, refcfw, qcfw, EDIT_TYPE_MM));
						}
					}
					if(annot != 0) {
						if(match && qc != annot) {
							// Matched, but didn't match the
							// reference allele
							compatAnnots = true;
						} else if(!match) {
							// Try again with annotation char
							match = (qc == annot);
						}
					}
				} else {
					// Ambiguous read character; match == false
				}
				if(!match) {
					if(!e2e && j < seedLen_) {
						assert_leq(worstStrat, seedMms);
						// Blew our -n budget?
						if(smms+1 > worstStrat) {
							break; // for(j = 0; j < adjAllen; j++)
						}
						smms++;
					} else if(e2e) {
						assert_leq(worstStrat, e2eMms);
						// Blew our -v budget?
						if(mms+1 > worstStrat) {
							break; // for(j = 0; j < adjAllen; j++)
						}
						mms++;
					}
					int q = phredcToPhredq(r.qualAt(j, true));
					if(pa.ignoreQuals) q = 30;
					qualPen += mmPenalty(maqRound_, q);
					// Blew our -e budget?
					if(!pa.ignoreQuals && pa.penceil[j] != -1 && qualPen > pa.penceil[j]) {
						break; // for(j = 0; j < adjAllen; j++)
					}
					if(color) {
						st.ccedits.push_back(Edit(j, refcfw, qcfw, EDIT_TYPE_MM));
					} else {
						st.edits.push_back(Edit(j, refcfw, qcfw, EDIT_TYPE_MM));
					}
				}
				// Add ambiguous base resolutions to st.aedits
				
				if(j+1 >= adjMinLen) {
					// We've reached the minimum length, so we're
					// valid now.  Keep extending.
					good = true;
				}
			} // loop over 5' to 3' read chars
			if(pa.requireAnnot && !overlapsAnnots) continue; // Reject this seed placement
			if(pa.requireAnnot && !compatAnnots) continue; // Reject this seed placement
			if(!good) {
				continue;
			}
			// Calculate stratum and cost
			int stratum = e2e ? mms : smms;
			assert_leq(stratum, worstStrat);
			assert_geq(stratum, bestStrat);
			int cost = stratum << 14 | qualPen;
			size_t replPos = 0;
			bool replaced = false;
			// Any bases we trimmed in the decode step come off now
			// for the ref.rc case.  For the !fw case, they
			// already came off
			if(j < adjAllen) {
				// Any bases we trimmed in the extend step come off now
				adjAllen = j;
				allen = adjAllen + color;
			}
			uint16_t clip3 = r.seq.length() - adjAllen;
			uint16_t seedoff = allen - iSeedLen_;
			assert_geq(allen, iSeedLen_);
			// Before we try replacing an element, see if our stratum
			// means that we can clear all of them
			if(pa.strata && !hits.empty() && stratum < hits.front().stratum) {
				hits.ents.clear();
			}
#if 0
			// BTL: this is only really useful in the phased-out
			// chaining mode.  I may try to rescue this some day but
			// for now it's depricated.
			else if(!hits.empty() && hits.tryReplacing(p, fw, cost, 0, seedoff, clip3, 0, replPos)) {
				if(replPos != 0xffffffff) {
					replaced = true;
					assert_lt(replPos, hits.size());
				} else {
					// Reject: it was no better than an existing
					// alignment at that position
					continue;
				}
			}
#endif
			// If there is an -m limit, check if the current candidate
			// will put us over the limit
			if((int)hits.size() >= pa.mhits) {
				assert(!hits.empty());
				assert(!pa.strata || stratum <= hits.front().stratum);
				if(pa.strata && stratum == hits.front().stratum) {
					assert_geq(hits.front().stratum, bestStrat);
					if(hits.front().stratum == bestStrat) {
						// We exceeded the mhits limit in the best
						// possible stratum; abort
						res.maxed = true;
						hits.maxedStratum = bestStrat;
						if(!pa.msample) {
							hits.ents.clear();
						}
						strat = worstStrat+1;
						done = true;
						break; // for(size_t it = 0; it < stratSz; it++)
					}
					// We exceeded the mhits limit in 'stratum', so
					// that stratum is disqualified
					assert_gt(stratum, 0);
					assert_leq(stratum, worstStrat);
					if(hits.maxedStratum == -1 || stratum < hits.maxedStratum) {
						hits.maxedStratum = stratum;
					}
					res.maxed = true;
					if(!pa.mhits) {
						hits.ents.clear(); // reject all hits so far
					}
					// refrain from any further searching within
					// 'stratum'
					worstStrat = stratum - 1;
					assert_geq(worstStrat, bestStrat);
					if(e2e) e2eMms = worstStrat;
					else seedMms = worstStrat;
				} else if(!pa.strata) {
					// Exceeded unstratified -m limit, reject read
					assert_eq((int)hits.size(), pa.mhits);
					res.maxed = true;
					hits.maxedStratum = 0;
					if(!pa.msample) {
						hits.ents.clear(); // reject all hits so far
					}
					strat = worstStrat+1;
					done = true;
				}
				break; // for(size_t it = 0; it < stratSz; it++)
			}
			// If reporting is stratified and the current candidate is
			// better than any previous candidate, we can raise the bar
			// on future candidates
			if(pa.strata) {
				// Possibly change e2eMms and/or seedMms
				if(e2e) e2eMms = min(e2eMms, mms);
				else    seedMms = min(seedMms, smms);
				worstStrat = min(worstStrat, stratum);
			}
			mms = e2e? mms : smms;
			// Now add the hit
			size_t elt;
			if(replaced) {
				elt = replPos;
			} else {
				elt = hits.size();
				hits.expand();
				if(elt == 0) {
					//hits.seq = r.seq;
					//hits.qual = r.qual;
				}
			}
			if(clip3 > 0 && !fpl && !fw) {
				// Any bases we trimmed in the extend step come off now
				p.second += clip3;
			}
			// Decode if it's in colorspace
			// Determine if this is a valid alignment along its entire
			// length
			HitSetEnt& hit = hits[elt];
			hit.clip3 = clip3;
			hit.clip5 = 0;
			if(color) {
				hit.seq.clear();
				hit.qual.clear();
				hit.cseq = r.seq;
				hit.cqual = r.qual;
				size_t readi = 0;
				size_t readf = allen - 1;
				size_t refi = p.second;
				size_t reff = p.second + allen;
				BTString refMask;
				if(!ref.seq.toRefMask(refMask, refi, reff, fw)) {
					// Couldn't make mask string, probably because
					// there was an unmatchable character
					continue;
				}
				hit.edits.clear();
				hit.aedits.clear();
				hit.cedits.clear();
				hit.ccedits.clear();
				hit.ccost = qualPen;
				hit.cost = st.dec_.decode(
					r.seq,  // ASCII colors, '0', '1', '2', '3', '.'
					r.qual, // ASCII quals, Phred+33 encoded
					readi, // offset of first character within 'read' to consider
					readf, // offset of last char (exclusive) in 'read' to consider
					refMask, // reference sequence, as masks
					0, // offset of first character within 'ref' to consider
					reff-refi, // offset of last char (exclusive) in 'ref' to consider
					pa.snpPen, // penalty incurred by a SNP
					MAX_SCORE, // max cost
					0, // # gaps on read side of colorspace-to-colorspace alignment
					0, // # gaps on ref side of colorspace-to-colorspace alignment
					st.ccedits, // color-to-color edits
					pa.readOpenPen,   // penalty for opening a new gap in the read
					pa.readExtendPen, // penalty for extending a gap in the read
					pa.refOpenPen,    // penalty for opening a new gap in the reference
					pa.refExtendPen,  // penalty for extending a gap in the reference
					pa.gapBarrier,    // # of alignment chars on either side that cannot have gaps
					pa.exDecEnds,     // true -> trim either end of nucleotide string
					pa.maqRound,      // true -> use Maq-like rounding
					hit.seq,    // decoded nucleotides appended here
					hit.qual,   // decoded qualities appended here
					hit.edits,  // destination for decoded nucleotide edits
					hit.aedits, // destination for resolved ambiguous bases
					hit.cedits, // destination for decoded color miscalls
					hit.ccedits,// destination for decoded color edits
					r.rand);
				assert(!e2e || mms <= worstStrat);
				assert(!e2e || hit.ccedits.size() <= (size_t)mms);
				assert_leq(hit.cedits.size(), hit.ccedits.size());
				//assert_leq(hit.edits.size(), hit.ccedits.size());
				assert_eq((int)hit.seq.length(), allen + (pa.exDecEnds ? -2 : 0));
				assert_eq(hit.seq.length(), hit.qual.length());
				assert(hit.repOk());
				if(pa.exDecEnds) {
					// It also shifts the alignment's offset up by 1
					p.second++;
					allen -= 2;
					adjAllen -= 2;
					// Ensure we haven't trimmed to less than the minimum
					assert_geq(hit.seq.length()+1, (size_t)adjMinLen);
				} else {
					// Ensure we haven't trimmed to less than the minimum
					assert_geq(hit.seq.length()-1, (size_t)adjMinLen);
				}
				// hit.seq, hit.qual, hit.edits, hit.cedits and
				// hit.ccedits are all arranged from 5' to 3' already.  The
				// only correction is to complement the query and reference
				// bases for the nucleotide edits.
				if(!fpl) {
					for (size_t i = 0; i < hit.edits.size(); i++) {
						hit.edits[i].qchr  = compDna(hit.edits[i].qchr);
						hit.edits[i].chr   = compDna(hit.edits[i].chr);
					}
					for (size_t i = 0; i < hit.aedits.size(); i++) {
						hit.aedits[i].qchr = compDna(hit.aedits[i].qchr);
						hit.aedits[i].chr  = compDna(hit.aedits[i].chr);
					}
				}
#ifndef NDEBUG
				// All the edits in cedits should also be in ccedits
				size_t ci = 0, cci = 0;
				while(ci < hit.cedits.size() && cci < hit.ccedits.size()) {
					const EList<Edit> cce = hit.ccedits;
					const EList<Edit> ce = hit.cedits;
					while(!cce[cci].isMismatch() ||
						  (cce[cci].isMismatch() && cce[cci].pos < ce[ci].pos))
					{
						cci++;
					}
					assert_lt(ci, hit.cedits.size());
					assert_lt(cci, hit.ccedits.size());
					assert_eq(hit.cedits[ci].pos,  hit.ccedits[cci].pos);
					assert_eq(hit.cedits[ci].qchr, hit.ccedits[cci].qchr);
					ci++; cci++;
				}
				assert_eq(ci, hit.cedits.size());
#endif
			} else {
				hit.seq = r.seq;
				hit.qual = r.qual;
				hit.cseq.clear();
				hit.cqual.clear();
				hit.edits = st.edits;
				hit.aedits = st.aedits;
				hit.cedits.clear();
				hit.ccedits.clear();
				hit.cost = qualPen;
				assert(hit.repOk());
			}
			
			// hit.seq is either the original base string (if !color) or the
			// decoded base string (if color) in the same arrangement as r.seq
			// (i.e. 5' to 3').  hit.edits, .cedits and .ccedits are all
			// arranged from 5' to 3' and hit.edits has already been
			// complemented if necessary.
			
			hit.stratum = stratum;
			hit.seedoff = seedoff;
			hit.oms = 0;
			hit.mapq = 40;
			hit.horig = p;
			if(clip3 > 0 && !fpl && fw) {
				// Any bases we trimmed in the extend step come off now
				p.second += clip3;
			}
			hit.h = p;
			hit.orient = orient;
			if(ref.crick) {
				// If we hit a revcomped reference string, we have to adjust our
				// offset to be w/r/t the forward reference string.  We'll
				// correct seq, qual, and edits fields later
				int coloradj = (color ? 1 : 0);
				int readadj = (int)r.seq.length() + coloradj;
				if(color && pa.exDecEnds) readadj -= 2;
				if(allen < readadj && fw) {
					readadj = allen;
				} else {
					hit.horig.second += hit.clip3;
				}
				// See cases 2b, 3b, 4b, 4d in comment in orient.h:
				assert_leq((int)(roff + readadj + coloradj), (int)(ref.seq.length(false)));
				hit.h.second = (uint32_t)ref.seq.length(false) - (roff + readadj + coloradj);
				assert_lt(hit.h.second, 0xffff0000);
			}
			if(ref.crick || ref.rc) {
				assert_neq(hit.h.first, ref.fwWatsonIdx);
				hit.h.first = ref.fwWatsonIdx;
			}
			// Correct h.seq, h.cseq and h.*edits in the event that the reverse
			// complement aligned 
			assert(hit.repOk());
			assert(!hit.seq.empty());
			// Correct h.seq, h.cseq and h.*edits in the event that
			// the printed read should be the reverse complement.
			// See cases 1b (WR), 2b (C), 3b, 4b (BC), 4c (BWR) in
			// comment in orient.h
			if(!orientPrintFw(hit.orient)) {
				// Whether to invert the positions of the
				// edits before sanity-checking the edit
				// list
				if(!hit.seq.empty()) {
					hit.seq.reverseComp(false);
					hit.qual.reverse();
					if(!color) {
						if(!hit.seq.empty())   hit.seq.trimBegin(hit.clip3);
						if(!hit.qual.empty())  hit.qual.trimBegin(hit.clip3);
					}
					// Make sure read and edits are compatible
					assert(Edit::repOk(hit.edits, hit.seq, hit.orient, 0, 0));
				} else assert(!color);
				if(!hit.cseq.empty()) {
					assert(color);
					hit.cseq.reverseComp(true);
					hit.cqual.reverse();
					if(!color) {
						if(!hit.cseq.empty())  hit.cseq.trimBegin(hit.clip3);
						if(!hit.cqual.empty()) hit.cqual.trimBegin(hit.clip3);
					}
					// Make sure read and edits are compatible
					assert(Edit::repOk(hit.ccedits, hit.cseq, hit.orient, 0, 0));
				} else assert(!color);
			} else if(!color) {
				if(!hit.seq.empty())   hit.seq.trimEnd(hit.clip3);
				if(!hit.qual.empty())  hit.qual.trimEnd(hit.clip3);
				if(!hit.cseq.empty())  hit.cseq.trimEnd(hit.clip3);
				if(!hit.cqual.empty()) hit.cqual.trimEnd(hit.clip3);
			}
			// Make sure read, edits and reference are compatible.  Can't
			// handle alignments to a reference with ref.rc == true in
			// bisulfite mode yet because the mismatches will be legitimately
			// off.
			assert(hit.repOk());
			assert(sanityCheckHit(refs, hit));
			// If there is a -k limit and no -m limit, and more
			// hits than are allowed by the -k limit...
			assert(pa.mhits != INT_MAX || (int)hits.size() <= pa.khits);
			if(pa.mhits == INT_MAX && (int)hits.size() == pa.khits) {
				// Truncate to -k hits
				hits.sort();
				done = true;
				break;
			}
			// Update maxHitCost, minHitCost
			if(hits.size() == 1) {
				maxHitCost = cost;
				minHitCost = cost;
			} else {
				if(cost > maxHitCost) maxHitCost = cost;
				if(cost < minHitCost) minHitCost = cost;
			}
			assert(!pa.strata || hits.uniformStrata());
		} // for(size_t it = 0; it < stratSz; it++)
		if(done) {
			break;
		}
	} // for(int strat = bestStrat; strat <= worstStrat; strat++)
	res.seedHits = seedHits;
	hits.sort();
	// Make sure that an empty hitset for a maxed-out read gets
	// counted as maxed-out, not unaligned
	bool abort = true;
	if(hits.maxedStratum != -1 || res.maxed) {
		assert_leq(hits.maxedStratum, 3);
		res.maxed = true;
		if(pa.msample) {
			assert_eq(pa.mhits, (int)hits.size());
			abort = false;
			for(size_t i = 0; i < hits.size(); i++) {
				hits[i].mapq = 0;
			}
		} else {
			os.printMaxedRead(
				r, r.seq, r.qual, hits.maxedStratum, pa.mhits);
		}
	} else if(hits.empty()) {
		os.printUnalignedRead(r, r.seq, r.qual, 0 /* not filtered */);
	} else {
		abort = false;
	}
	if(abort) {
		return;
	}
	assert_leq((int)hits.size(), pa.mhits);
	if(pa.strata) assert(hits.uniformStrata());
	res.hits = hits.size();
	hits.reportUpTo(r, os, pa.khits, refs, rmap, false, amap);
}

/**
 * Align the read against the indexed, IUPAC-containing reference.
 */
void MerIndex::query(
	Read& r,
	const ReferenceSet& refs,
	const ReferenceMap* rmap,
	const AnnotationMap* amap,
	HitSet& hits,
	AlignOutput& os,
	AlignResult& res,
	const AlignParams& p,
	bool randomize,
	RandomSource& rnd,
	int threadid)
{
	assert_geq(threadid, 0);
	if(naiveCheck_) {
		HitSet naiveHits; // = r.hitset;
		AlignResult nres; // stub
		StubOutput nos; // stub
		assert(r.repOk());
		naive_.query(r, refs, rmap, amap, naiveHits, nos, nres, p, false, rnd, threadid);
		queryHelper (r, refs, rmap, amap,      hits,  os,  res, p, false, rnd, threadid);
		naiveHits.sort();
		// Check results against naive results
		size_t off = 0, noff = 0;
		ostringstream diffout;
		HitSet hitscp; // = r.hitset;
		while(off < hitscp.size() && noff < naiveHits.size()) {
			if(hitscp[off] == naiveHits[noff]) {
				hitscp.remove(off);
				naiveHits.remove(noff);
			} else if(hitscp[off] < naiveHits[noff]) {
				off++;
			} else {
				noff++;
			}
		}
		if(hitscp.empty() && naiveHits.empty()) return;
		res.bail = true;
		os.flush();
		if(!hitscp.empty()) {
			cerr << "Merman but not naive:" << endl;
			hitscp.reportUpTo(r, os, INT_MAX, refs, NULL, false, NULL);
			for (size_t i = 0; i < hitscp[0].seq.length(); i++) {
				cerr << refs[hitscp[0].h.first].seq.charAt(hitscp[0].h.second + i, false);
			}
			cerr << endl;
			os.flush();
			cerr << "--------------------" << endl;
		}
		if(!naiveHits.empty()) {
			cerr << "Naive but not merman:" << endl;
			naiveHits.reportUpTo(r, os, INT_MAX, refs, NULL, false, NULL);
			for (size_t i = 0; i < naiveHits[0].seq.length(); i++) {
				cerr << refs[naiveHits[0].h.first].seq.charAt(naiveHits[0].h.second + i, false);
			}
			cerr << endl;
			os.flush();
			cerr << "--------------------" << endl;
		}
		assert(false);
	} else {
		queryHelper(r, refs, rmap, amap, hits, os, res, p, true, rnd, threadid);
	}
}

/**
 * Extract read mer subsequences from the read and store them in a
 * list, sorted such that ReadMers with greater entropy among the 4
 * bases are first.
 */
void MerIndex::extractReadmers(const Read& r,
                               EList<ReadMer>& readmers,
                               EList<string>& merstrs,
                               bool nofw, bool norc) const
{
	readmers.clear();
	for(int s = 0; s < (int)seeds_.size(); s++) {
		// For both forward and reverse-complement orientations
		for(int o = (nofw ? 1 : 0); o < (norc ? 1 : 2); o++) {
			bool fw = (o == 0);
			if(merverbose && (norc || !fw)) {
				uint64_t seed = seeds_[s];
				for(int j = 0; j < seedWidth_; j++) {
					if((seed & 1llu) != 0) {
						cout << "1";
					} else {
						cout << "0";
					}
					seed >>= 1;
				}
				cout << ", seed: " << s << endl;
			}
			// For each placement of the seed width within the seed portion
			// of the read
			assert(iSeedLen_ >= seedWidth_);
			assert((int)r.seq.length() >= iSeedLen_);
			for(int i = 0; i < (iSeedLen_ - seedWidth_ + 1); i++) {
				uint64_t seed = seeds_[s];
				readmers.expand();
				readmers.back().mer = 0;
				if(merverbose) {
					merstrs.push_back("");
				}
				int cnts[] = {0, 0, 0, 0};
				bool aborted = false;
				int cnt = 0;
				// Construct the seed string from the read
				for(int j = 0; j < seedWidth_; j++) {
					if((seed & 1llu) != 0) {
						cnt++;
						int c = fw ? r.seq[i+j] :
						             r.seq[iSeedLen_ - (i+j) - 1];
						assert_range(0, 4, c);
						if(c < 4) {
							cnts[(int)c]++;
							append(c, readmers.back().mer, !fw && !r.color);
							if(merverbose) merstrs.back().push_back("0123"[(int)c]);
						} else {
							// N in the read guarantees that this seed
							// never hits
							readmers.pop_back();
							aborted = true;
							break;
						}
					}
					else if(merverbose) {
						if(merverbose) merstrs.back().push_back(' ');
					}
					if(cnt == 20) break;
					seed >>= 1llu;
				}
				// Add range information
				if(!aborted) {
					int min = 9999;
					int max = 0;
					for(int j = 0; j < 4; j++) {
						if(cnts[j] > max) max = cnts[j];
						if(cnts[j] < min) min = cnts[j];
					}
					assert(max >= min);
					readmers.back().range = max - min;
					readmers.back().fw = fw;
					readmers.back().s = s;
					readmers.back().i = i;
				}
				else if(merverbose) {
					merstrs.pop_back();
				}
			}
		}
	}
	if(readmers.size() > 1) {
		// Move the highest-information seeds to the beginning
		readmers.sort();
	}
	// Print them
	if(merverbose) {
		if(r.color) {
			for(size_t i = 0; i < merstrs.size(); i++) {
				assert_leq((int)merstrs[i].length(), seedWidth_);
				int ons = 0;
				for(size_t j = 0; j < merstrs[i].length(); j++) {
					if(merstrs[i][j] != ' ') ons++;
					printColor(merstrs[i][j]);
				}
				assert_leq(ons, 20);
				cout << endl;
			}
		}
	}
}

/**
 * Return true iff the placement of Ns in the read dictates that it
 * won't align anywhere.
 */
bool MerIndex::nFilter(const Read& r,
                       int adjMinLen,
                       int seedMms,
                       int e2eMms) const
{
	// First thing's first: check if the arrangement of Ns in the read
	// disqualifies it from aligning
	int mms = 0;
	int smms = 0;
	for(int i = 0; i < adjMinLen; i++) {
		assert_range(0, 4, (int)r.seq[i]);
		bool mm = r.seq[i] == 4;
		if(mm) {
			// For the read, we assume any character that isn't
			// unambiguous DNA incurs a mismatch
			if(++smms > seedMms && (int)i < seedLen_) {
				// No alignments possible
				return true;
			}
			if(++mms > e2eMms) {
				// No alignments possible
				return true;
			}
		}
	}
	return false;
}

/**
 * Extract all the subsequences from the reference that will make up
 * the index.
 *
 * tid = thread id, 0-based
 * nt = number of threads
 */
pair<size_t, size_t>
MerIndex::extractMers(
	const ReferenceSet& refs, int tid, int nt, bool install, bool color)
{
	assert_gt(nt, 0);
	assert_geq(tid, 0);
	assert_lt(tid, nt);
	// If there's just one thread, we extract a subsequence every
	// 'period_' positions in the reference.  If there are multiple
	// threads, each thread handles a staggered subset of the positions
	// and the period becomes 'period_ * nt'.
	//
	// Note that period_ = iSeedLen_ - seedWidth_ + 1
	int period = period_ * nt;
	size_t tot = 0, extra = 0;
	mer_ent lMers[1024]; // temporary, large mer_ents
	size_t nMers = 0;
	size_t lOcc = 0;
	size_t runlen = 0;
	EList<uint64_t> refmers;
	// This list (nucmasks) is useful in the case of a colorspace
	// reference with ambiguity.
	EList<char> nucmasks;
	EList<string> refmerstrs;
	for(size_t ri = 0; ri < refs.size();
	    runlen += refs[ri].seq.length(false), ri++)
	{
		if(!refs[ri].index) {
			// This reference is flagged to indicate we shouldn't
			// index it.
			continue;
		}
		const RefString& s = refs[ri].seq;
		size_t len = s.length(color);
		// Initialize i to 'tid * period_' so that multiple threads are
		// staggered.
		for(size_t i = (tid * period_); i < len; i += period) {
			// Index the mer with window anchored at i
			// For each seed...
			if(i + seedWidth_ > len) break;
			for(int j = 0; j < (int)seeds_.size(); j++) {
				uint64_t seed = seeds_[j];
				// Clear list where we temporarily store the
				// subsequences extracted from the reference.
				refmers.clear();
				if(merverbose) refmerstrs.clear();
				// We set this to false as soon as we know that there
				// will be more than once subsequence copy generated at
				// this position (i.e. because the mask overlaps an
				// IUPAC code).
				refmers.push_back(0llu);
				if(merverbose) refmerstrs.push_back("");
				int cnt = 0;
				size_t subtot = 1, subextra = 1;
				// true iff the previous position in the mask was set to 1
				bool lastSet = false;
				// Extract the Dna characters
				for(int k = 0; k < seedWidth_; k++) {
					if((seed & 1llu) != 0) {
						seed >>= 1llu;
						cnt++; // count of set mask bits
						char c = 4;
						char cik = toupper(s.charAt(i+k, color));
						if(!color) {
							// We're in nucleotide space
							if(s.isUnmatchableAt(i+k, color)) {
								// The nucleotide is unmatchable; this
								// invalidates all of the subsequences we were
								// extracting
								refmers.clear();
								subtot = subextra = 0;
								if(merverbose) refmerstrs.clear();
								break;
							} else if(s.isAmbigAt(i+k, color)) {
								// The nucleotide is ambiguous
								int num = 0;
								int nucalts[] = {4, 4, 4, 4};
								decodeAmbigNuc(cik, num, nucalts);
								assert_gt(num, 1);
								assert_neq(nucalts[0], nucalts[1]);
								subextra *= num;
								// Create new sets of refmers with all but the
								// first alternative
								const size_t oldsz = refmers.size();
								for(int l = 1; l < num; l++) {
									for(size_t m = 0; m < oldsz; m++) {
										refmers.push_back(refmers[m]);
										append(nucalts[l], refmers.back());
										if(merverbose) {
											refmerstrs.push_back(refmerstrs[m]);
											refmerstrs.back().push_back("ACGT"[nucalts[l]]);
											assert_leq((int)refmerstrs.back().length(), seedWidth_);
										}
									}
								}
								// Now add the first alternative to the
								// original refmers
								for(size_t m = 0; m < oldsz; m++) {
									append(nucalts[0], refmers[m]);
									if(merverbose) {
										refmerstrs[m].push_back("ACGT"[nucalts[0]]);
										assert_leq((int)refmerstrs[m].length(), seedWidth_);
									}
								}
							} else {
								assert(isDna(cik));
								c = asc2dna[(int)cik];
								for(size_t it = 0; it < refmers.size(); it++) {
									append(c, refmers[it]);
									if(merverbose) {
										refmerstrs[it].push_back("ACGT"[(int)c]);
										assert_leq((int)refmerstrs[it].length(), seedWidth_);
									}
								}
							}
						} else {
							// We're in colorspace
							if(s.isUnmatchableAt(i+k, color)) {
								// The color is unmatchable; this invalidates
								// all subsequences coming from this reference
								// position
								refmers.clear();
								subtot = subextra = 0;
								if(merverbose) refmerstrs.clear();
								break;
							} else if(s.isAmbigAt(i+k, color)) {
#if 0
								// Ambiguous color
								int num = 0;
								int calts1[] = { 4, 4, 4, 4 };
								// Get the upstream nucleotide
								char n1 = toupper(s.charAt(i+k, false));
								// nextToo <- true iff the next character is also in
								// the subsequence to be extracted
								bool nextToo = ((seed & 1llu) != 0);
								if(!isAmbigNuc(n1) &&        // not ambiguous
								   nextToo &&                // next one's also in the mask
								   k < seedWidth_-1 &&       // we're not at the edge of the mask
								   i+k+1 < s.length(true) && // we're not at the edge of the reference
								   cnt != 20)                // we haven't already filled our extracted mer
								{
									// This color is the upstream member of the
									// pair of colors that correspond to an
									// ambiguous nucleotide, and the next color is
									// going into the same seed
									char am = toupper(s.charAt(i+k+1, false));
									char n2 = toupper(s.charAt(i+k+2, false));
									assert(isAmbigNuc(am));
									assert(!isAmbigNuc(n2));
									int calts2[] = { 4, 4, 4, 4 };
									decodeAmbigNucTrio(n1, am, n2, num, calts1, calts2);
									k++;
									seed >>= 1;
									cnt++;
									// Create new sets of refmers with all but the
									// first alternative
									const size_t oldsz = refmers.size();
									subextra *= num;
									for(int l = 1; l < num; l++) {
										for(size_t m = 0; m < oldsz; m++) {
											refmers.push_back(refmers[m]);
											append(calts1[l], refmers.back());
											append(calts2[l], refmers.back());
											if(merverbose) {
												refmerstrs.push_back(refmerstrs[m]);
												refmerstrs.back().push_back("ACGT"[calts1[l]]);
												refmerstrs.back().push_back("ACGT"[calts2[l]]);
												assert_leq((int)refmerstrs[m].length(), seedWidth_);
											}
										}
									}
									// Now add the first alternative to the
									// original refmers
									for(size_t m = 0; m < oldsz; m++) {
										append(calts1[0], refmers[m]);
										append(calts2[0], refmers[m]);
										if(merverbose) {
											refmerstrs[m].push_back("ACGT"[calts1[0]]);
											refmerstrs[m].push_back("ACGT"[calts2[0]]);
											assert_leq((int)refmerstrs[m].length(), seedWidth_);
										}
									}
								} else {
									// Get the dinucleotide
									char n2 = toupper(s.charAt(i+k+1, false));
									decodeAmbigNucPair(n1, n2, num, calts1);
									// Create new sets of refmers with all but the
									// first alternative
									const size_t oldsz = refmers.size();
									subextra *= num;
									for(int l = 1; l < num; l++) {
										for(size_t m = 0; m < oldsz; m++) {
											refmers.push_back(refmers[m]);
											append(calts1[l], refmers.back());
											if(merverbose) {
												refmerstrs.push_back(refmerstrs[m]);
												refmerstrs.back().push_back("ACGT"[calts1[l]]);
												assert_leq((int)refmerstrs.back().length(), seedWidth_);
											}
										}
									}
									// Now add the first alternative to the
									// original refmers
									for(size_t m = 0; m < oldsz; m++) {
										append(calts1[0], refmers[m]);
										if(merverbose) {
											refmerstrs[m].push_back("ACGT"[calts1[0]]);
											assert_leq((int)refmerstrs[m].length(), seedWidth_);
										}
									}
								}
#else
								bool upambig = s.isAmbigAt(i+k,   false);
								bool dnambig = s.isAmbigAt(i+k+1, false);
								assert(upambig || dnambig);
								// nextSet <- true iff the next color is also
								// in the subsequence to be extracted
								bool nextSet = ((seed & 1llu) != 0);
								// If the downstream nucleotide is ambiguous
								// and the next mask bit is set, we need to
								// propagate information about the downstream
								// nucleotide associated with each color we
								// generate to the next iteration.
								bool recordDnNucs = nextSet && dnambig;
								// If the upstream nucleotide is ambiguous and
								// the *previous* mask bit was set, we need to
								// take downstream-nucleotide information from
								// the last iteration into account when
								// deciding which refmers to extend with which
								// colors.
								bool extendUpNucs = lastSet && upambig;
								// Nucleotide upstream of color
								int upnuc = s.charAt(i+k,   false);
								// Nucleotide downstream of color
								int dnnuc = s.charAt(i+k+1, false);
								int mks1[4], mks2[4];
								bool cs[4];
								decodeAmbigColorNucPair(
									upnuc, // nucleotide 1
									dnnuc, // nucleotide 2
									mks1,  // upstream nucleotide alternatives
									mks2,  // downstream nucleotide alternatives
									cs);   // cs[i] == true iff color i is possibile
								// We may want to clear nucmasks so that we can
								// install new masks at some point, but we
								// shouldn't yet because we may still need it
								size_t oldsz = refmers.size();
								if(extendUpNucs) {
									// Extend mers in refmers list, with the
									// complication that the added color's
									// preceding nucleotide must be compatible
									// with the last-nucleotide mask for the
									// refmer
									assert_eq(nucmasks.size(), oldsz);
									// For each refmer in need of extending
									for(size_t mi = 0; mi < oldsz; mi++) {
										bool first = true;
										assert_range(1, 15, (int)nucmasks[mi]);
										int firstci = -1;
										// For each color we might extend with
										for(int ci = 0; ci < 4; ci++) {
											if((nucmasks[mi] & mks1[ci]) != 0) {
												// Extend this refmer with this
												// color
												if(first) {
													// Replace mask at this
													// position with a new mask
													firstci = ci;
													first = false;
												} else {
													// We're extending the same
													// refmer in more than one
													// way and this is not the
													// first, so we need to
													// create a new refmer
													subextra++;
													refmers.push_back(refmers[mi]);
													nucmasks.push_back(mks2[ci]);
													append(ci, refmers.back());
													if(merverbose) {
														refmerstrs.push_back(refmerstrs[mi]);
														refmerstrs.back().push_back("ACGT"[ci]);
														assert_leq((int)refmerstrs.back().length(), seedWidth_);
													}
												}
											}
										}
										assert(!first); // must have extended
										assert_range(0, 3, firstci);
										nucmasks[mi] = mks2[firstci];
										append(firstci, refmers[mi]);
										if(merverbose) {
											refmerstrs[mi].push_back("ACGT"[(int)firstci]);
											assert_leq((int)refmerstrs[mi].length(), seedWidth_);
										}
									}
								} else {
									// Extend mers in refmers list; ignore nucmasks
									// For each refmer in need of extending
									if(recordDnNucs) {
										nucmasks.resize(oldsz);
									}
									for(size_t mi = 0; mi < oldsz; mi++) {
										bool first = true;
										int firstci = -1;
										for(int ci = 0; ci < 4; ci++) {
											if(mks1[ci] != 0) {
												if(first) {
													// Replace mask at this
													// position with a new mask
													firstci = ci;
													first = false;
												} else {
													// We're extending the same
													// refmer in more than one
													// way and this is not the
													// first, so we need to
													// create a new refmer
													subextra++;
													refmers.push_back(refmers[mi]);
													if(recordDnNucs) {
														nucmasks.push_back(mks2[ci]);
													}
													append(ci, refmers.back());
													if(merverbose) {
														refmerstrs.push_back(refmerstrs[mi]);
														refmerstrs.back().push_back("ACGT"[ci]);
														assert_leq((int)refmerstrs.back().length(), seedWidth_);
													}
												}
											}
										}
										assert(!first); // must have extended
										assert_range(0, 3, firstci);
										if(recordDnNucs) {
											nucmasks[mi] = mks2[firstci];
										}
										append(firstci, refmers[mi]);
										if(merverbose) {
											refmerstrs[mi].push_back("ACGT"[(int)firstci]);
											assert_leq((int)refmerstrs[mi].length(), seedWidth_);
										}
									}
								}
#endif
							} else {
								// Neither the nucleotide at i+k nor the
								// nucleotide at i+k+1 is ambiguous.
								assert(isDna(cik));
								c = asc2dna[(int)cik];
								for(size_t it = 0; it < refmers.size(); it++) {
									append(c, refmers[it]);
									if(merverbose) {
										refmerstrs[it].push_back("ACGT"[(int)c]);
										assert_leq((int)refmerstrs[it].length(), seedWidth_);
									}
								}
							}
						}
#if 0
						// Case 3: Nucleotide/color is ambiguous
						else if(!color) {
						} else {
#endif
						lastSet = true;
					} else {
						seed >>= 1llu;
						if(merverbose) {
							for(size_t it = 0; it < refmerstrs.size(); it++) {
								refmerstrs[it].push_back(' ');
								assert_leq((int)refmerstrs[it].length(), seedWidth_);
							}
						}
						lastSet = false;
					}
					assert_leq(cnt, 20);
					if(cnt == 20) break;
				}
				tot += subtot;
				extra += subextra;
				// refmers is now populated with all the appropriate seed
				// keys for this window
				for(size_t it = 0; it < refmers.size(); it++) {
					assert_lt(lOcc, 1024);
					lMers[lOcc].init(refmers[it], (uint32_t)(runlen + i), j);
					//dump << refmers[it] << '\t' << (runlen+i) << '\t' << j << endl;
					lOcc++;
					assert_leq(lOcc, 1024);
					if(lOcc == 1024) {
						if(!install) countMers(lMers, 1024);
						else {
							installMers(lMers, 1024);
							nMers += 1024;
							assert_leq(nMers, mersSz_);
						}
						lOcc = 0;
					}
				}
				if(merverbose && !install) {
					string spaces;
					spaces.resize(i, ' ');
					for(size_t it = 0; it < refmerstrs.size(); it++) {
						assert_leq((int)refmerstrs[it].length(), seedWidth_);
						cout << spaces;
						if(color) {
							int ons = 0;
							for(size_t k = 0; k < refmerstrs[it].length(); k++) {
								if(refmerstrs[it][k] != ' ') ons++;
								printColor(refmerstrs[it][k]);
							}
							assert_leq(ons, 20);
							cout << endl;
						} else {
							cout << refmerstrs[it] << endl;
						}
					}
				}
			} // loop over seeds
		} // loop over positions
		if(merverbose) {
			if(color) {
				for(size_t j = 0; j < s.length(true); j++) {
					printColor(s.charAt(j, true));
				}
				cout << endl;
			}
			cout << s.seq << endl;
		}
	}
	// Flush remaining merss
	if(lOcc > 0) {
		assert_lt(lOcc, 1024);
		if(!install) countMers(lMers, lOcc);
		else {
			installMers(lMers, lOcc);
			nMers += lOcc;
			assert(nt > 1 || nMers == mersSz_);
		}
	}
	assert(nt > 1 || nMers == mersSz_);
	//dump.close();
	return make_pair(extra, tot);
}

/**
 * Threadsafe function for installing a list of mers into the global
 * mers_ array.
 */
void MerIndex::installMers(const mer_ent* mers, size_t sz) {
	ThreadSafe t(&lock_);
	for(size_t i = 0; i < sz; i++) {
		uint8_t s = mers[i].seed;
		assert_lt(s, seeds_.size());
		uint16_t key = (uint16_t)mers[i].key;
		uint32_t occ = mersOcc_[s][key];
		assert_lt(occ, mersLen_[s][key]);
		mers_[s][key][occ] = mers[i];
		mersOcc_[s][key]++;
	}
}

/**
 * Threadsafe function for incrementing counters according to the
 * elements in a list of mers.
 */
void MerIndex::countMers(const mer_ent* mers, size_t sz) {
	ThreadSafe t(&lock_);
	for(size_t i = 0; i < sz; i++) {
		uint8_t s = mers[i].seed;
		assert_lt(s, seeds_.size());
		uint16_t key = (uint16_t)mers[i].key;
		mersLen_[s][key]++;
	}
}

/**
 * Threadsafe function for incrementing counters according to the
 * elements in a list of mers.
 */
void MerIndex::allocateMers() {
	for(size_t i = 0; i < seeds_.size(); i++) {
		for(size_t j = 0; j < 256*256; j++) {
			if(mersLen_[i][j] > 0) {
				uint32_t len = mersLen_[i][j];
				mers_[i][j] = new mer_ent_sm[len];
				mersSz_ += len;
				assert(mers_[i][j] != NULL);
			} else {
				mers_[i][j] = NULL;
			}
		}
	}
}

/**
 * Checks if the mersOcc_, mersLen_, and mers_ arrays seem to be
 * consistent with each other.
 */
bool MerIndex::sanityCheckMers() const {
	for(size_t i = 0; i < seeds_.size(); i++) {
		for(size_t j = 0; j < 256*256; j++) {
			assert_eq(mersOcc_[i][j], mersLen_[i][j]);
			if(mersLen_[i][j] > 0) {
				assert(mers_[i][j] != NULL);
			} else {
				assert(mers_[i][j] == NULL);
			}
		}
	}
	return true;
}

/**
 * Sort the mer_ent list by key and
 */
void MerIndex::sort(int nt) {
	sanityCheckMers();
	EList<pair<mer_ent_sm*,mer_ent_sm*> > ps;
	for(size_t i = 0; i < seeds_.size(); i++) {
		for(size_t j = 0; j < 256*256; j++) {
			if(mers_[i][j] == NULL) {
				assert_eq(0, mersLen_[i][j]);
				continue;
			}
			assert_gt(mersLen_[i][j], 0);
			if(nt == 1) {
				// Go ahead and sort it here in the master thread
				std::sort<mer_ent_sm*>(
					mers_[i][j],                    // begin
					mers_[i][j] + mersLen_[i][j]);  // end
				assert(reallySorted(i, j));
				assert(reallySorted(0, 0)); // corruption?
			} else {
				// Add an element to the worklist for one of the sort
				// threads to work on
				ps.push_back(make_pair(
					mers_[i][j],                    // begin
					mers_[i][j] + mersLen_[i][j])); // end
			}
		}
	}
	if(nt > 1) workingListParallelSort(ps, nt);
	assert(reallySorted());
	sorted_ = true;
}

/**
 * Sort the mer_ent list by key and
 */
bool MerIndex::reallySorted(size_t i, size_t j) const {
	if(mersLen_[i][j] < 2) return true;
	for(size_t k = 0; k < mersLen_[i][j]-1; k++) {
		assert(!(mers_[i][j][k+1] < mers_[i][j][k]));
	}
	return true;
}

/**
 * Sort the mer_ent list by key and
 */
bool MerIndex::reallySorted() const {
	for(size_t i = 0; i < seeds_.size(); i++) {
		for(size_t j = 0; j < 256*256; j++) {
			assert(reallySorted(i, j));
		}
	}
	return true;
}
