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

#include <utility>
#include "align_naive.h"
#include "read.h"
#include "ref.h"
#include "ds.h"
#include "refmap.h"
#include "annot.h"
#include "qual.h"
#include "check.h"
#include "orient.h"

using namespace std;

typedef pair<uint32_t, uint32_t> U32Pair;

void NaiveAligner::query(
	Read& rd,
	const ReferenceSet& refs,
	const ReferenceMap* rmap,
	const AnnotationMap* amap,
	HitSet& hits,
	AlignOutput& /*os*/,
	AlignResult& /*res*/,
	const AlignParams& pa,
	bool randomize,
	RandomSource& rnd,
	int threadid)
{
	assert_geq(threadid, 0);
	assert_lt((size_t)threadid, threadState_.size());
	NaivePerThreadState& st = threadState_[threadid];
	bool e2e = (pa.e2eMms != INT_MAX);
	// For each reference sequence...
	for(size_t r = 0; r < refs.size(); r++) {
		const Reference& ref = refs[r];
		const RefString& refseq = ref.seq;
		int color = rd.color ? 1 : 0;
		// For each strand
		for(int o = (pa.alignfw ? 0 : 1); o < (pa.alignrc ? 2 : 1); o++) {
			bool fw = (o == 0);
			int orient = (ref.crick ? ORIENT_CRICK : ORIENT_WATSON) |
			             (fw        ? ORIENT_FW    : ORIENT_RC);
			int fpl = is5pLeft(orient);
			// minLen = minimum allowed length of the aligned portion
			// of the read in nucleotides.
			int adjMinLen = (pa.minLen == -1 ? (int)rd.seq.length() : pa.adjMinLen);
			assert_leq(adjMinLen, (int)rd.seq.length());
			// Need overhang because trimming might allow us to align a
			// read that, without trimming, would have gone off the end
			int minOverhang = (int)rd.seq.length() - adjMinLen;
			assert_geq(minOverhang, 0);
			// Scan reference
			for(int j = -minOverhang; j <= ((int)refseq.length(false) - pa.minLen); j++) {
				st.edits.clear();
				st.aedits.clear();
				st.ccedits.clear();
				int qualPen = 0;
				if(fw && j < 0) {
					// 5' end of read can't hang off end of the reference
					continue;
				}
				if(!fw && j + rd.seq.length() > refseq.length(color)) {
					// 5' end of read can't hang off end of the reference
					continue;
				}
				int adjAllen = (int)rd.seq.length();
				// Decoded nucleotide alignment is one character longer
				// than colorspace alignment
				int allen = (int)rd.seq.length() + color;
				int refi = j;
				U32Pair p = make_pair(r, refi);
				int smms = 0; // # mismatches in seed
				int mms = 0; // # mismatches overall
				bool good = false; // current extenison is valid
				bool overlapsAnnots = false; // Does it overlap annotations?
				bool compatAnnots = false; // Are there any compatible annots?
				AnnotationMap::Iter ait;
				// Translate into target coordinate system before
				// looking up annotations in annotation map
				if(rmap != NULL) rmap->map(p);
				// Copy annotations that overlap candidate alignment
				// into annots[] array.  5' end is in annots[0].
				char annots[1024];
				memset(annots, 0, allen);
				if(amap != NULL) {
					ait = amap->lower_bound(p);
					while(ait != amap->end() &&
					      ait->first.first == p.first &&
					      ait->first.second < p.second + allen)
					{
						overlapsAnnots = true;
						size_t off = ait->first.second - p.second;
						if(!fw) off = allen - off - 1;
						// annots[0] is annotation for 5'-most base
						annots[off] = ait->second.second;
						ait++;
					}
				}
				int k;
				int prevChar = -1;
				// For each posision involved in the alignment, going
				// from 5' to 3' on the read
				for(k = 0; k < adjAllen; k++) {
					// See if we can extend through this position
					int loff = k; // offset w/r/t reference
					if(!fw) {
						// Going right-to-left on the reference
						loff = (int)rd.seq.length() - k - 1;
					}
					int fpoff = k; // offset w/r/t 5' end of read
					bool match = false;
					if((size_t)(j+loff) >= refseq.length(color)) {
						// Not caught by the j loop bounds?
						break;
					}
					int refc = toupper(refseq.charAt(j+loff, color, fw ? prevChar : -1, fw ? -1 : prevChar));
					assert(prevChar == -1 || isUnambigNuc(refc));
					prevChar = -1;
					if(isUnmatchableNuc(refc) ||
					   (isAmbigNuc(refc) && pa.disallowIupac))
					{
						// Covers an unmatchable character
						break;
					}
					// Get next character (going from read 5' to 3')
					char readc = toupper(rd.charAt5p(fpoff, fw));
					if(isUnambigNuc(readc)) {
						if(isUnambigNuc(refc) && annots[fpoff] == 0) {
							// Unambiguous, unannotated reference position
							// versus unambiguous read position
							match = (asc2dna[(int)readc] == asc2dna[(int)refc]);
						} else if(isUnambigNuc(readc)) {
							// Either ambiguous or annotated (or both)
							// reference position versus unambiguous read
							// position
							match = ambigCompatNuc(refc, readc);
							if(match && color) {
								int refa = toupper(refseq.charAt(j+loff + (fw ? 1 : 0), false));
								if(isAmbigNuc(refa)) {
									int refu = toupper(refseq.charAt(j+loff + (fw ? 0 : 1), false));
									if(isAmbigNuc(refu)) {
										cerr << "Error: Encountered two ambiguous reference nucleotides in a row" << endl;
										throw 1;
									}
									assert(isAmbigNuc(refa));
									assert(isUnambigNuc(refu));
									int prevCharN = nuccol2nuc[asc2dna[refu]][asc2dna[(int)readc]];
									assert_lt(prevCharN, 4);
									assert_geq(prevCharN, 0);
									prevChar = "ACGT"[prevCharN];
								}
							}
							if(annots[fpoff] != 0) {
								if(match && readc != annots[fpoff] ) {
									// Matched, but didn't match the
									// reference allele
									compatAnnots = true;
								}
								else if(!match) {
									match = (readc == annots[fpoff] );
								}
							}
						}
					}
					if(!match) {
						assert_neq(refc, readc);
						if(fpoff < pa.seedLen) {
							// Blew the -n budget?
							if(smms+1 > pa.seedMms) break;
							smms++;
						}
						// Blew the -v budget?
						if(mms+1 > pa.e2eMms) break;
						mms++;
						int q = phredcToPhredq(rd.qualAt(k, fw));
						if(pa.ignoreQuals) q = 30;
						qualPen += mmPenalty(maqRound_, q);
						// Blew the -e budget?
						if(!pa.ignoreQuals && pa.penceil[k] != -1 && qualPen > pa.penceil[k]) break;
						// Accept the edit
						if(!color) {
							st.edits.push_back(Edit(fpoff, refc, readc, EDIT_TYPE_MM));
						} else {
							st.ccedits.push_back(Edit(fpoff, refc, readc, EDIT_TYPE_MM));
						}
					}
					if(adjMinLen >= 0 ? (k+1 >= adjMinLen) : (k == adjAllen-1)) {
						// We've reached the minimum length, so we're
						// valid now.  Keep extending.
						good = true;
					}
				} // loop over alignment chars
				if(pa.requireAnnot && (!overlapsAnnots || !compatAnnots)) {
					continue;
				}
				assert_gt(k, -1);
				if(good) {
					// Length of actual (possbly trimmed) color alignment
					adjAllen = k;
					// Length of actual (possbly trimmed) nucleotide alignment
					allen = adjAllen + color;
					assert_geq(allen, pa.minLen);
					assert_geq(adjAllen, pa.adjMinLen);
					assert_geq(adjAllen, pa.iSeedLen);
					// Calculate # chars trimmed
					uint16_t clip3 = rd.seq.length() - adjAllen;
					// Leave edits alone
					hits.expand();
					HitSetEnt& hit = hits.back();
					hit.clip5 = 0;
					hit.clip3 = clip3;
					if(clip3 > 0 && !fpl && !fw) {
						// Any bases we trimmed in the extend step come off now
						p.second += clip3;
					}
					if(color) {
						hit.seq.clear();
						hit.qual.clear();
						hit.cseq = rd.seq;
						hit.cqual = rd.qual;
						// Decode colorspace alignment
						size_t readi = 0;
						size_t readf = allen-1;
						size_t reff = p.second + allen;
						int64_t rfi = p.second;
						if(reff > refseq.length(false)) {
							int diff = (int)reff - (int)refseq.length(true);
							readf -= diff;
							reff -= diff;
							allen -= diff;
						}
						if(rfi < 0) {
							int diff = -(int)rfi;
							allen -= diff;
							readi += diff;
							rfi += diff;
						}
						if(pa.minLen != -1 && reff-rfi < (size_t)pa.minLen) {
							continue;
						}
						BTString refMask;
						if(!refseq.toRefMask(refMask, rfi, reff, fw)) {
							// Couldn't make mask string, probably because
							// there was an unmatchable character
							continue;
						}
#ifndef NDEBUG
						for(size_t i = 0; i < (reff-rfi); i++) {
							assert_gt(refMask[i], 0);
							assert_lt(refMask[i], 16);
						}
#endif
						hit.edits.clear();
						hit.aedits.clear();
						hit.cedits.clear();
						hit.ccedits.clear();
						qualPen = (int)st.dec_.decode(
							rd.seq,  // ASCII colors, '0', '1', '2', '3', '.'
							rd.qual, // ASCII quals, Phred+33 encoded
							readi, // offset of first character within 'read' to consider
							readf, // offset of last char (exclusive) in 'read' to consider
							refMask, // reference sequence, as masks
							0, // offset of first character within 'ref' to consider
							reff-rfi, // offset of last char (exclusive) in 'ref' to consider
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
							pa.exDecEnds, // true -> trim either end of nucleotide string
							pa.maqRound, // true -> use Maq-like rounding
							hits.back().seq,    // decoded nucleotides appended here
							hits.back().qual,   // decoded qualities appended here
							hit.edits, // destination for decoded nucleotide edits
							hit.aedits, // destination for resolved ambiguous bases
							hit.cedits, // destination for decoded color miscalls
							hit.ccedits, // destination for decoded color edits
							rd.rand);
						hit.ccost = qualPen;
						//assert_geq(hit.ccedits.size(), st.ccedits.size());
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
						}
						// hit.seq, hit.qual, hit.edits, hit.cedits and
						// hit.ccedits are all arranged from 5' to 3' already.  The
						// only correction is to complement the query and reference
						// bases for the nucleotide edits.
						if((!fw) != refs[r].crick) {
							for (size_t i = 0; i < hit.edits.size(); i++) {
								hit.edits[i].qchr = compDna(hit.edits[i].qchr);
								hit.edits[i].chr = compDna(hit.edits[i].chr);
							}
							for (size_t i = 0; i < hit.aedits.size(); i++) {
								hit.aedits[i].qchr = compDna(hit.aedits[i].qchr);
								hit.aedits[i].chr = compDna(hit.aedits[i].chr);
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
						// !color
						hit.seq = rd.seq;
						hit.qual = rd.qual;
						// hit.seq and hit.qual will be trimmed below,
						// after they're possibly reverse-complemented
						hit.cseq.clear();
						hit.cqual.clear();
						hit.edits = st.edits;
						hit.aedits = st.aedits;
						hit.cedits.clear();
						hit.ccedits.clear();
						hit.cost = qualPen;
						assert(hit.repOk());
					}
					// Clipping should already have been applied by now
					if(good) {
						hit.cost = qualPen;
						hit.stratum = e2e ? mms : smms;
						hit.oms = 0;
						hit.mapq = 40;
						hit.horig = p;
						if(clip3 > 0 && !fpl && fw) {
							// Any bases we trimmed in the extend step come off now
							p.second += clip3;
						}
						hit.h = p;
						hit.orient = (ref.crick ? ORIENT_CRICK : ORIENT_WATSON) |
						             (fw        ? ORIENT_FW    : ORIENT_RC);
						if(ref.crick) {
							// If we hit a revcomped reference string, we have to adjust our
							// offset to be w/r/t the forward reference string.  We'll
							// correct seq, qual, and edits fields later
							int coloradj = (color ? 1 : 0);
							int readadj = (int)rd.seq.length() + coloradj;
							if(color && pa.exDecEnds) readadj -= 2;
							if(allen < readadj && fw) {
								readadj = allen;
							} else {
								hit.horig.second += hit.clip3;
							}
							// See cases 2b, 3b, 4b, 4d in comment in orient.h:
							assert_leq(j + readadj + coloradj, (int)(ref.seq.length(false)));
							hit.h.second = (uint32_t)refseq.length(false) - (j + readadj + coloradj);
							assert_lt(hit.h.second, 0xffff0000);
						}
						if(ref.crick || ref.rc) {
							assert_neq(hit.h.first, ref.fwWatsonIdx);
							hit.h.first = ref.fwWatsonIdx;
							assert_eq(refs[hit.h.first].fwWatsonIdx, hit.h.first);
						}
						assert_eq(refs[hit.h.first].fwWatsonIdx, hit.h.first);
						assert(hit.repOk());
						assert(!hit.seq.empty());
						//
						// TODO: make the following block a member function of hit
						//
						// Correct h.seq, h.cseq and h.*edits in the event that
						// the printed read should be the reverse complement.
						// See cases 1b (WR), 2b (C), 3b, 4b (BC), 4c (BWR) in
						// comment in orient.h:
						if(!orientPrintFw(hit.orient)) {
							// Whether to invert the positions of the
							// edits before sanity-checking the edit
							// list
							if(!hit.seq.empty()) {
								hit.seq.reverseComp(false);
								hit.qual.reverse();
								if(!color) {
									if(!hit.seq.empty())  hit.seq.trimBegin(hit.clip3);
									if(!hit.qual.empty()) hit.qual.trimBegin(hit.clip3);
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
						assert_eq((int)hit.seq.length(), allen);
						assert_eq((int)hit.qual.length(), allen);
						// Make sure read, edits and reference are compatible.  Can't
						// handle alignments to a reference with ref.rc == true in
						// bisulfite mode yet because the mismatches will be legitimately
						// off.
						assert(hit.repOk());
						// TODO: This can't currently work in both
						// nucleotide-space and colorspace because we
						// truncate the colorspace reads before we get
						// here, but not so with nucleotide space.
						assert(sanityCheckHit(refs, hit));
#ifndef NDEBUG
						for(size_t i = 0; i < hits.size(); i++) {
							for(size_t l = l+1; l < hits.size(); l++) {
								assert_neq(hits[i], hits[l]);
							}
						}
#endif
					} else {
						// Never mind, we rejected the alignment based on the
						// decoding.
						hits.pop_back();
					}
				}
			}
		}
	}
	if(pa.strata) {
		int bestStrat = 999;
		for(size_t i = 0; i < hits.size(); i++) {
			if(hits[i].stratum < bestStrat) {
				bestStrat = hits[i].stratum;
			}
		}
		for(size_t i = 0; i < hits.size(); i++) {
			if(hits[i].stratum != bestStrat) {
				hits.remove(i);
				i--;
			}
		}
	}
}
