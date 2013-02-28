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

#include <iostream>
#include <sstream>
#include <string>
#include "alphabet.h"
#include "filebuf.h"
#include "globals.h"
#include "ref.h"

using namespace std;

/**
 * Add the sequences in a FASTA file as references.  We assume these sequences
 * represent the forward Watson strand(s).
 */
void ReferenceSet::addOrigReferenceFasta(
	const char *ss,           // name of the FASTA file
	const ReferenceParams& p) // config variables related to the reference
{
	// Add a new Reference object
	TRefStr *seq = NULL, *name = NULL;
	// Open the file
	FILE *f = fopen(ss, "r");
	if(f == NULL) {
		cerr << "Could not open reference fasta file " << ss << " for reading" << endl;
		throw 1;
	}
	// Nest in a FileBuf object for easier manipulation and buffering
	FileBuf fb(f);
	size_t origRefsSz = numRefs();
	while(!fb.eof()) {
		// Parse name line
		if(fb.peek() == '>') {
			if(seq == NULL || !seq->empty()) {
				// Add a new Reference object
				if(seq != NULL) {
					// Finish adding the Reference we just finished
					// parsing
					refs_.back().off = (uint32_t)totRefLen_;
					addReference(p);
				}
				refs_.expand();
				refs_.back().idx = (uint32_t)refs_.size()-1;
				refs_.back().crick = false;
				refs_.back().rc = false;
				refs_.back().fwWatsonIdx = (uint32_t)refs_.size()-1;
				refs_.back().index = true;
			}
			seq = &refs_.back().seq.seq;
			name = &refs_.back().name;
			if(!name->empty()) {
				name->clear();
			}
			fb.get();
			// Install the name in refs_.back().name
			while(fb.peek() != '\n' && fb.peek() != '\r') {
				assert(isprint(fb.peek()));
				name->append(fb.peek());
				fb.get();
			}
			while(fb.peek() == '\n' || fb.peek() == '\r') fb.get();
		}
		while(fb.peek() != '>' && !fb.eof()) {
			// Install the sequence (taken to be all the alphabetical
			// characters) in refs_.back().seq.
			int c = fb.peek();
			if(c == '-') c = 'N';
			if(isalpha(c)) {
				assert(isprint(c));
				if(seq != NULL) {
					seq->append(c);
				} else {
					// Oops, saw sequence characters before we saw a
					// '>' line.
				}
			}
			fb.get();
		}
		// Parse sequence
	}
	if(seq != NULL) {
		refs_.back().off = (uint32_t)totRefLen_;
		addReference(p);
	}
	if(origRefsSz == numRefs()) {
		cerr << "Warning: fasta input file had no sequences" << endl;
	}
	fb.close();
}

/**
 * Add a string to the list of indexed references.  We assume that these
 * sequences represent the forward Watson strand(s).
 */
void ReferenceSet::addOrigReferenceString(
	const char *ss,           // nucleotide string for reference
	const ReferenceParams& p) // 
{
	// Add a new Reference object
	if(ss[0] == '\n') return;
	refs_.expand();
	refs_.back().seq.seq.install(ss);
	refs_.back().idx = (uint32_t)refs_.size()-1;    // this reference's position in the list
	refs_.back().fwWatsonIdx = (uint32_t)refs_.size()-1;  // index of the forward version of this sequence
	refs_.back().index = true;            // yes, index me
	refs_.back().off = (uint32_t)totRefLen_;        // this reference's offset into the concatenated
	                                      // string containing all references in order
	refs_.back().crick = false;           // no, this does not represent the crick strand
	refs_.back().rc = false;              // no, this is not a reverse-comp sequence
	refs_.back().name = tostring(refs_.size()); // copy gets same name as original
	addReference(p);
}

/**
 * Mark all forward references as non-indexable.
 */
void ReferenceSet::removeWatsonOrCrickReferences(bool watson) {
	for(size_t i = 0; i < refs_.size(); i++) {
		assert(refs_[i].index);
		if(refs_[i].crick == (!watson)) refs_[i].index = false;
	}
}

/**
 * Take all of the reference sequences currently in the reference list,
 * add a new copy of it onto the end of the list then reverse
 * complement the copy.
 *
 * The most common reason to do this is if we're aligning bisulfite
 * treated reads.  Recall that the "bisulfite reference" typically
 * comprises either two strands (BiSulfite Watson = BSW and BiSulfite
 * Crick = BSC) or four strands (BSW, BSC, reverse complement of BSW =
 * BSWR, and reverse complement of BSC = BSCR).  For the two-strand
 * case, the two strands are not reverse complements of each other.
 * For this reason, we can't simply align the read and its reverse
 * complement to a single reference strand and call it a day - rather,
 * we have to construct and index both strands and align the read
 * against them separately.  This is also true of the four-strand case,
 * though the additional two strands *can* be queried by taking the
 * reverse complement of the read, so there is no need to construct and
 * index all four strands.
 */
void ReferenceSet::addReferenceRevComps(
	const ReferenceParams& p,
	bool postXforms,
	int crick, // whether resulting sequence should have crick=true, -1 means "inherit"
	int rc)    // whether resulting sequence should have rc=true, -1 means "inherit"
{
	size_t oldsz = refs_.size();
	for(size_t i = 0; i < oldsz; i++) {
		// Add a new slot at the end
		refs_.expand();
		// Let copy's sequence be the reverse complement of the
		// *original* sequence.  This means we have to undo any
		// modifications applied to the original first, before creating
		// the copy.
		if(!postXforms) refs_[i].seq.toggleXforms(); 
		revcompDna(refs_[i].seq.seq, refs_.back().seq.seq, false);
		// Length should be same as the original
		assert_eq(refs_[i].seq.length(true), refs_.back().seq.length(true));
		assert_eq(refs_[i].seq.length(false), refs_.back().seq.length(false));
#ifndef NDEBUG
		// Sequence should be same as the original
		for(size_t j = 0; j < refs_[i].seq.length(false); j++) {
			assert(isprint(refs_.back().seq.seq[j]));
		}
#endif
		// Copy is now created so we can re-apply the modifications
		// that were applied to the original reference
		if(!postXforms) refs_[i].seq.toggleXforms();
		// TODO: can we free the xform memory at this point?
		refs_.back().idx = (uint32_t)refs_.size()-1; // this reference's position in the list
		refs_.back().off = (uint32_t)totRefLen_;     // this reference's offset into the concatenated
		                                   // string containing all references in order
		if(crick != -1) {
			refs_.back().crick = (crick == 1);
		} else {
			refs_.back().crick = refs_[i].crick;
		}
		if(rc != -1) {
			refs_.back().rc = (rc == 1);
		} else {
			refs_.back().rc = refs_[i].rc;
		}
		refs_.back().name = refs_[i].name; // copy gets same name as original
		refs_.back().fwWatsonIdx = (uint32_t)i;      // index of the forward version of this sequence
		refs_.back().index = true;         // yes, index me
		addReference(p);                   // add the new reverse-complemented reference
	}
}

/**
 * Extract mers from the reference at the back of the refs_ list.
 */
void ReferenceSet::addReference(const ReferenceParams& p) {
	// Copy sequence into a new Reference object
	TRefStr& s = refs_.back().seq.seq;
	const size_t len = s.length();
	if(verbose) cout << "  reading ref " << refs_.back().name << endl;
	assert_gt(len, 0);
	if(p.bisulfiteC) {
		if(verbose) cout << "    bisulfite-treating Cs in reference" << endl;
		// Turn all the Cs in the reference into ys
		int num = 0;
		for(size_t i = 0; i < len; i++) {
			int c = s[i];
			assert(isprint(c));
			int C = toupper(c);
			if(C == 'C') {
				// Found a C in the reference; remember the old
				// character
				refs_.back().seq.addXform((uint32_t)i, c);
				s.set('y', i); // C or T
				num++;
			}
		}
		if(verbose) cout << "      treated " << num << " Cs" << endl;
	} else if(p.bisulfiteCpG) {
		if(verbose) cout << "    bisulfite-treating CpGs in reference" << endl;
		int cgNum = 0, cNum = 0;
		// Look for CpG dinucleotides and non-CpG C nucleotides
		for(size_t i = 0; i < len - 1; i++) {
			int c1 = s[i];
			int c2 = s[i+1];
			assert(isprint(c1));
			assert(isprint(c2));
			int C1 = toupper(c1);
			int C2 = toupper(c2);
			if(C1 == 'C' && C2 == 'G') {
				// Found a CpG
				refs_.back().seq.addXform((uint32_t)i, c1);
				s.set('y', i); // C or T
				cgNum++;
			} else if(C1 == 'C') {
				refs_.back().seq.addXform((uint32_t)i, c1);
				s.set('t', i); // just t
				cNum++;
			}
		}
		// Don't forget the last base
		int cl = s[len-1];
		int Cl = toupper(cl);
		if(Cl == 'C') {
			refs_.back().seq.addXform((uint32_t)len-1, cl);
			s.set('t', len-1);
			cNum++;
		}
		if(verbose) {
			cout << "      treated " << cgNum << " CpGs and " << cNum
			     << " non-CpG Cs"<< endl;
		}
	}
	// Mask out regions that fail to reach the entropy threshold
	bool masked = false;
	if(p.entThresh.first > 0) {
		if(verbose) cout << "    entropy-masking reference" << endl;
		entBuf_.resizeAndFill(p.entThresh.first, 0);
		int cur = 0;
		bool full = false; // whether all elts in the buffer
		int sum = 0;
		string msg;
		int ranges = 0;
		long bases = 0;
		for(size_t i = 0; i < len; i++) {
			// Were we below the threshold in the last iteration?
			bool wasAbove = sum > p.entThresh.second;
			assert_lt((size_t)cur, entBuf_.size());
			assert(full || entBuf_[cur] == 0);
			// Adjust sum to remove oldest contribution
			sum -= entBuf_[cur];
			// Add new contribution
			assert(isprint(s[i]));
			entBuf_[cur] = (isUnmatchableNuc(s[i]) ? -9999 : iupacAlleles(s[i]));
			sum += entBuf_[cur];
			// full <- true iff rotating buffer has a full window worth
			// of data
			if(cur == p.entThresh.first-1) full = true;
			if(full) {
				bool isBelow = sum > p.entThresh.second;
				if(isBelow) masked = true;
				if(isBelow && !wasAbove) {
					// Mask entire window
					assert(i >= (size_t)p.entThresh.first);
					for(size_t j = i - p.entThresh.first + 1; j <= i; j++) {
						if(verbose) msg += s[j];
						s.set('-', j);
						bases++;
					}
					ranges++;
				} else if (isBelow) {
					// Just mask most recent character
					if(verbose) msg += s[i];
					s.set('-', i);
					bases++;
				} else if (!isBelow && wasAbove && verbose) {
					cout << "Entropy-masked: " << msg << endl;
					msg.clear();
				}
			}
			cur = (cur + 1) % p.entThresh.first;
		}
		if(verbose) {
			cout << "      masked " << ranges << " ranges, " << bases
			     << " bases" << endl;
		}
		if(false && verbose && !msg.empty()) {
			cout << "Entropy-masked: " << msg << endl;
		}
	}
	// Associate this reference with its global offset, for when we'd
	// like to take a hit in the global space and resolve it to a
	// particular reference sequence.
	ASSERT_ONLY(size_t osz = offToRef_.size());
	offToRef_.insert(refs_.back().off, &refs_.back());
	assert_eq(osz, offToRef_.ltBound(refs_.back().off));
	assert_eq(offToRef_.get(osz), &refs_.back());
	totRefLen_ += len;
}
