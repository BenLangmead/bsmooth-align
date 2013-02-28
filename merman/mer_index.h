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

#ifndef MER_INDEX_H_
#define MER_INDEX_H_

#include <string>
#include <algorithm>
#include <utility>
#include "alphabet.h"
#include "mer_ent.h"
#include "read.h"
#include "edit.h"
#include "hit_set.h"
#include "annot.h"
#include "refmap.h"
#include "threading.h"
#include "seeds.h"
#include "output.h"
#include "ds.h"
#include "ref.h"
#include "align_naive.h"
#include "random_source.h"

/**
 * Comparator for finding a desired Reference in a list of
 * Reference's sorted by offset.
 */
struct Reference_Compare {
	bool operator() (const Reference& a, const unsigned int& b) {
		return a.off <= b;
	}
};

/**
 * Encapsulates a subsequence extracted from somewhere in the indexable
 * region of the read.
 */
struct ReadMer {
	ReadMer() {
		mer = 0llu;
		range = i = s = 999;
		fw = true;
	}
	int operator< (const ReadMer &rhs) const {
		return range < rhs.range;
	}
	
	/**
	 * Check that ReadMer is internally consistent.
	 */
	bool repOk() const {
		return true;
	}

	uint64_t mer; // sequence extracted
	// difference in occurrence count between most- and least-often
	// occurring characters; vaguely approximates entropy
	int range;
	bool fw; // true = extracted from forward read, false = rev comp
	int i; // offset into indexable region
	int s; // index of seed
};

/**
 * Per-thread state for threads querying a MerIndex.  We could put
 * these things in the query() function itself, but that would incur
 * fresh memory allocations for each call to query().
 */
struct MerIndexPerThreadState {

public:

	MerIndexPerThreadState() {
		merIts = NULL;
		merItEnds = NULL;
		merItVals = NULL;
		merItsSz = 0;
	}

	MerIndexPerThreadState(size_t seeds, int isl, int sw) {
		init(seeds, isl, sw);
	}

	void init(size_t seeds, int isl, int sw) {
		assert_geq(isl, sw);
		merItsSz = seeds * 2 * (isl - sw + 1);
		merIts = new const mer_ent_sm*[merItsSz];
		merItEnds = new mer_ent_sm*[merItsSz];
		merItVals = new uint64_t[merItsSz];
	}

	~MerIndexPerThreadState() {
		delete[] merIts;
		delete[] merItEnds;
		delete[] merItVals;
	}

	EList<ReadMer>     readmers;
	EList<std::string> readmerstrs;
	EList<int64_t>     stratHits[7];
	const mer_ent_sm** merIts;
	mer_ent_sm**       merItEnds;
	uint64_t*          merItVals;
	size_t             merItsSz;

	EList<Edit> edits;
	EList<Edit> aedits;
	EList<Edit> cedits;
	EList<Edit> ccedits;
	
	ColorspaceDecoder dec_;
};


/**
 * Reads in the reference string (which may include IUPAC codes),
 * extracts and installs mers in an array of mer_index'es, and then
 * sorts the array.
 */
class MerIndex : public Aligner {
public:

	MerIndex(
		const AlignParams& p,
		const ReferenceParams& rp,
		int readLen,
		int width,
		int chooseN,
		int chooseK,
		int specificity,
		int begin,
		bool naiveCheck,
		int nthreads) :
			Aligner(nthreads),
			sset_(SeedSet::create(p.seedMms, specificity, chooseN, chooseK)),
			naive_(nthreads, rp.bisulfiteC, rp.bisulfiteCpG, p.maqRound)
	{
		sset_->getSeeds(seeds_, width, 100);
		readLen_ = readLen;
		seedWidth_ = width;
		seedLen_ = p.seedLen;
		iSeedLen_ = p.iSeedLen;
		assert_geq(iSeedLen_, seedWidth_);
		period_ = iSeedLen_ - seedWidth_ + 1;
		assert_gt(period_, 0);
		sorted_ = false;
		maqRound_ = p.maqRound;
		begin_ = begin;
		bisulfiteC_ = rp.bisulfiteC;
		bisulfiteCpG_ = rp.bisulfiteCpG;
		naiveCheck_ = naiveCheck;
		mersLen_ = new uint32_t*[seeds_.size()];
		mersOcc_ = new uint32_t*[seeds_.size()];
		mers_ = new mer_ent_sm**[seeds_.size()];
		for(int i = 0; i < (int)seeds_.size(); i++) {
			mersLen_[i] = new uint32_t[256*256];
			mersOcc_[i] = new uint32_t[256*256];
			mers_[i] = new mer_ent_sm*[256*256];
			memset(mersLen_[i], 0, (256*256)*sizeof(uint32_t));
			memset(mersOcc_[i], 0, (256*256)*sizeof(uint32_t));
			memset(mers_[i], 0, (256*256)*sizeof(mer_ent_sm*));
		}
		mersSz_ = 0llu;
		threadState_.resize(nthreads);
		for(int i = 0; i < nthreads; i++) {
			threadState_[i].init(seeds_.size(), iSeedLen_, seedWidth_);
		}
		MUTEX_INIT(lock_);
	}

	~MerIndex() {
		for(int i = 0; i < (int)seeds_.size(); i++) {
			for(int j = 0; j < (256*256); j++) {
				if(mers_[i][j] != NULL) {
					delete[] mers_[i][j];
					mers_[i][j] = 0;
				}
			}
			delete[] mersLen_[i];
			delete[] mersOcc_[i];
			delete[] mers_[i];
		}
		delete[] mersLen_;
		delete[] mersOcc_;
		delete[] mers_;
		delete sset_;
	}

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
		int nthreads);

	/**
	 * Return true iff the mer_end list is sorted.
	 */
	bool sorted() const { return sorted_; }

	/**
	 * Sort the mer_ent list by key and
	 */
	void sort(int nt = 1);

	/**
	 * Return true iff the index contains no mers.
	 */
	bool empty() const { return mersSz_ == 0llu; }

	/**
	 * Return the number of spaced seeds.
	 */
	size_t numSeeds() const { return seeds_.size(); }

	/**
	 * Extract the subsequences that will constitute the index and call
	 * installMers to install them into the index.
	 */
	std::pair<size_t, size_t> extractMers(
		const ReferenceSet& refs,
		int tid,
		int nt,
		bool install,
		bool color);

	/**
	 * Just count (do not extract) subsequences and return total.
	 */
	std::pair<size_t, size_t> countMers(int tid, int nt, bool color) const;

	/**
	 * Having counted the number of mers in each category, allocate
	 * memory for them.
	 */
	void allocateMers();

	/**
	 * Return number of subsequences in the index.
	 */
	size_t size() const { return mersSz_; }

	/**
	 * Checks if the mersOcc_, mersLen_, and mers_ arrays seem to be
	 * consistent with each other.
	 */
	bool sanityCheckMers() const;

	int getISeedLen() const  { return iSeedLen_; }
	int getSeedLen() const   { return seedLen_; }
	int getSeedWidth() const { return seedWidth_; }

protected:

	/**
	 * Find the first element of mers_ that is not < the query.  This
	 * is performance-critical.
	 */
	const mer_ent_sm* lowerBound(
		uint16_t key1,
		uint16_t key2,
		uint8_t key3,
		int seed,
		mer_ent_sm*& end) const;

	/**
	 * Find the first element of mers_ that is not < the query.
	 */
	const mer_ent_sm* lowerBound(uint64_t key,
	                             int seed,
	                             mer_ent_sm*& end) const
	{
		assert(sorted_);
		return lowerBound(key, (key >> 16), (key >> 32), seed, end);
	}

	/**
	 * Query the iupac index with the given read, requiring that all
	 * reported alignments must be in a better stratum than
	 * 'worstStratum'.
	 */
	void queryHelper(
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
		int nthreads);

	/**
	 * Assert that mers_ lists are really sorted.
	 */
	bool reallySorted() const;

	/**
	 * Assert that the i,jth element of the mers_ list is really sorted.
	 */
	bool reallySorted(size_t i, size_t j) const;

	/**
	 * Install a batch of extracted mers into the mers_ array using
	 * synchronization.
	 */
	void installMers(const mer_ent* mers, size_t sz);

	/**
	 * Threadsafe function for incrementing counters according to the
	 * elements in a list of mers.
	 */
	void countMers(const mer_ent* mers, size_t sz);

	/**
	 * Return true iff the placement of Ns in the read dictates that it
	 * won't align.
	 */
	bool nFilter(const Read& r, int minLen, int seedMms, int e2eMms) const;

	/**
	 * Extract read mer subsequences from the read and store them in a
	 * list.  Store "more promising" (less homogenous) subsequences
	 * earlier in the list, with the thought that they're more likely
	 * to lead to an alignment if one exists.
	 */
	void extractReadmers(const Read& r,
	                     EList<ReadMer>& readmers,
	                     EList<std::string>& merstrs,
	                     bool nofw, bool norc) const;

	int readLen_;
	int seedWidth_;
	EList<uint64_t> seeds_;
	SeedSet *sset_;
	int seedLen_;
	int iSeedLen_;
	int period_;
	// 4 MB
	uint32_t **mersLen_;
	// 4 MB
	uint32_t **mersOcc_;
	// 4/8 MB (depending on 32/64-bit pointers)
	mer_ent_sm* **mers_;
	uint64_t mersSz_;
	bool maqRound_;
	bool sorted_;
	int begin_;
	bool bisulfiteC_;
	bool bisulfiteCpG_;
	bool naiveCheck_;
	EList<MerIndexPerThreadState> threadState_;
	NaiveAligner naive_;
	MUTEX_T lock_;
};

/**
 * Thread that performs MerIndex indexing on a subset of the reference
 * mers where the subset is defined by the thread id tid and total
 * number of threads nt.
 *
 * Note that indexing is divided into a counting phase followed by a
 * mer extraction phase.  Threads split and join separately for the
 * phases.
 */
class MerIndexThread {
public:

	void runIndex(const ReferenceSet* refs,
	              MerIndex *ind,
	              int tid,
	              int nt,
	              bool color)
	{
		refs_ = refs;
		ind_ = ind; tid_ = tid; nt_ = nt;
		color_ = color;
		THREAD_CREATE(thread_, MerIndexThread::startIndexThread, this);
	}

	void runCount(const ReferenceSet* refs,
	              MerIndex *ind,
	              int tid,
	              int nt,
	              bool color)
	{
		refs_ = refs;
		ind_ = ind; tid_ = tid; nt_ = nt;
		cnt_ = std::make_pair(0, 0);
		color_ = color;
		THREAD_CREATE(thread_, MerIndexThread::countSubseqsThread, this);
	}

	/**
	 * Wait until this thread is finished before returning.
	 */
	std::pair<size_t, size_t> join() {
		THREAD_JOIN(thread_);
		return cnt_;
	}

private:

	/**
	 * Do a stride's worth of indexing.
	 */
	void indexWork() {
#ifndef NDEBUG
		std::pair<size_t, size_t> cnt;
		cnt =
#endif
		ind_->extractMers(*refs_, tid_, nt_, true, color_);
		assert_eq(cnt.first, cnt_.first);
		assert_eq(cnt.second, cnt_.second);
	}

	/**
	 * Do a stride's worth of subsequence counting.
	 */
	void countWork() {
		cnt_ = ind_->extractMers(*refs_, tid_, nt_, false, color_);
	}

	/**
	 * Start the work of a single indexing thread.
	 */
	static void* startIndexThread(void *obj) {
		reinterpret_cast<MerIndexThread*>(obj)->indexWork();
		return NULL;
	}

	/**
	 * Start the work of a single counting thread.
	 */
	static void* countSubseqsThread(void *obj) {
		reinterpret_cast<MerIndexThread*>(obj)->countWork();
		return NULL;
	}

	const ReferenceSet *refs_;
	int tid_;
	int nt_;
	MerIndex *ind_;
	THREAD_T thread_;
	std::pair<size_t, size_t> cnt_;
	bool color_;
};

#endif /* MER_INDEX_H_ */
