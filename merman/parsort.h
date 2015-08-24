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

#ifndef PARSORT_H_
#define PARSORT_H_

#include <utility>
#include <algorithm>
#include <memory>
#include "threading.h"
#include "ds.h"

/**
 * Create a thread that sorts a partition of a longer list.
 */
template<typename T>
class DivideMergeSortThread {
public:

	/**
	 * Sort a partition of the list, depending on the thread id.
	 */
	void run(T* begin, T* end, int tid, int nt) {
		tid_ = tid;
		nt_ = nt;
		begin_ = begin;
		end_ = end;
		a_ = begin + (((end - begin)/nt) * tid);
		b_ = (tid == nt-1) ? end : (begin + (((end - begin)/nt) * (tid+1)));
		THREAD_CREATE(thread_, DivideMergeSortThread::startThread, this);
	}

	/**
	 * Wait until this thread is finished before returning.
	 */
	void join() { THREAD_JOIN(thread_); }

private:

	/**
	 * Do a stride's worth of sorting.
	 */
	void work() { std::sort(a_, b_); }

	/**
	 * Start the work of a single thread.
	 */
	static void* startThread(void *obj) {
		reinterpret_cast<DivideMergeSortThread*>(obj)->work();
		return NULL;
	}

	int tid_;
	int nt_;
	T* begin_, *end_;
	T* a_, *b_;
	THREAD_T thread_;
};

/**
 * Create a thread that performs an in-place merge of adjacent, sorted
 * partitions of a longer list.
 */
template<typename T>
class DivideMergeMergeThread {
public:

	/**
	 * Sort a partition of the list, depending on the thread id.
	 */
	void run(T* begin, T* end, int tid, int nt) {
		tid_ = tid;
		nt_ = nt;
		begin_ = begin;
		end_ = end;
		a1_ = begin + (((end - begin)/nt) * tid);
		b1_ = begin + (((end - begin)/nt) * (tid+1));
		a2_ = b1_;
		assert_leq(tid, nt-2);
		b2_ = (tid == nt-2) ? end : (begin + (((end - begin)/nt) * (tid+2)));
		THREAD_CREATE(thread_, DivideMergeMergeThread::startThread, this);
	}

	/**
	 * Wait until this thread is finished before returning.
	 */
	void join() { THREAD_JOIN(thread_); }

private:

	/**
	 * Do a stride's worth of sorting.
	 */
	void work() { std::inplace_merge(a1_, b1_, b2_); }

	/**
	 * Start the work of a single thread.
	 */
	static void* startThread(void *obj) {
		reinterpret_cast<DivideMergeMergeThread*>(obj)->work();
		return NULL;
	}

	int tid_;
	int nt_;
	T* begin_, *end_;
	T* a1_, *b1_, *a2_, *b2_;
	THREAD_T thread_;
};

/**
 * Parallel sort where all threads are working on one large chunk of
 * memory and proceed in waves of subsequence sorts and merges.
 */
template <typename T>
void divideAndMergeParallelSort(T* begin, T* end, int nt = 1) {
	using namespace std;
	typedef EList<DivideMergeSortThread<T> > SortThreadVec;
	typedef EList<DivideMergeMergeThread<T> > MergeThreadVec;
	typedef typename EList<DivideMergeSortThread<T> >::iterator SortThreadVecIter;
	typedef typename EList<DivideMergeMergeThread<T> >::iterator MergeThreadVecIter;
	assert_geq(nt, 1);
	// Round nt down to nearest power of 2
	for(int i = 30; i >= 1; i--) {
		if(nt > (1 << i)) {
			nt = (1 << i);
			break;
		}
	}
	SortThreadVec sthreads;
	sthreads.resize(nt);
	SortThreadVecIter sit;
	int tid = 0;
	for(sit = sthreads.begin(); sit != sthreads.end(); sit++) {
		sit->run(begin, end, tid++, nt);
	}
	for(sit = sthreads.begin(); sit != sthreads.end(); sit++) sit->join();
	while(nt > 1) {
		nt >>= 1;
		MergeThreadVec mthreads;
		mthreads.resize(nt);
		MergeThreadVecIter mit;
		tid = 0;
		for(mit = mthreads.begin(); mit != mthreads.end(); mit++) {
			mit->run(begin, end, tid, nt<<1);
			tid += 2;
		}
		for(mit = mthreads.begin(); mit != mthreads.end(); mit++) mit->join();
	}
}

/**
 * Create a thread that performs an in-place merge of adjacent, sorted
 * partitions of a longer list.
 */
template<typename T>
class WorkingListSortThread {
public:

	/**
	 * Sort a partition of the list, depending on the thread id.
	 */
	void run(const EList<std::pair<T*,T*> >* ps,
	         MUTEX_T* lock,
	         size_t *cur, int tid, int nt)
	{
		ps_ = ps; lock_ = lock;
		cur_ = cur; tid_ = tid; nt_ = nt;
		THREAD_CREATE(thread_, WorkingListSortThread::startThread, this);
	}

	/**
	 * Wait until this thread is finished before returning.
	 */
	void join() { THREAD_JOIN(thread_); }

private:

	/**
	 * Do a stride's worth of sorting.
	 */
	void work() {
		T *begin = NULL, *end = NULL;
#ifndef NDEBUG
		{
			ThreadSafe t(lock_);
		}
#endif
		while(true) {
			{
				ThreadSafe t(lock_);
				if((*cur_) < ps_->size()) {
					begin = (*ps_)[*cur_].first;
					assert(begin != NULL);
					end = (*ps_)[*cur_].second;
					assert(end != NULL);
					(*cur_)++;
				} else {
					break;
				}
			}
			ASSERT_ONLY(MUTEX_T* oldLock = lock_);
			assert(begin != NULL);
			assert(end > begin);
			std::sort(begin, end);
			assert(lock_ == oldLock);
		}
	}

	/**
	 * Start the work of a single thread.
	 */
	static void* startThread(void *obj) {
		reinterpret_cast<WorkingListSortThread*>(obj)->work();
		return NULL;
	}

	int tid_;
	int nt_;
	T* begin_, *end_;
	size_t *cur_;
	THREAD_T thread_;
	MUTEX_T* lock_;
	const EList<std::pair<T*,T*> >* ps_;
};

/**
 * Parallel sort where threads grab regions of memory to be sorted on a
 * first-come-first-served basis.
 */
template <typename T>
void workingListParallelSort(const EList<std::pair<T*,T*> >& ps,
                             int nt = 1)
{
	using namespace std;
	typedef EList<WorkingListSortThread<T> > SortThreadVec;
	size_t cur = 0;
	std::auto_ptr<MUTEX_T> wlpsLock(new MUTEX_T);
	MUTEX_INIT(*wlpsLock.get());
	assert_geq(nt, 1);
	auto_ptr<SortThreadVec> wlpsThreads(new SortThreadVec());
	wlpsThreads->resize(nt);
	int tid = 0;
	for(size_t i = 0; i < (*wlpsThreads).size(); i++) {
		(*wlpsThreads)[i].run(&ps, wlpsLock.get(), &cur, tid++, nt);
	}
	for(size_t i = 0; i < (*wlpsThreads).size(); i++) (*wlpsThreads)[i].join();
}

#endif /* PARSORT_H_ */
