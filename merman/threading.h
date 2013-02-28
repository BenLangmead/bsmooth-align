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

#ifndef THREADING_H_
#define THREADING_H_

#include <algorithm>
#include <stdint.h>

#ifdef USE_PTHREADS
#include <sched.h>
#    include <pthread.h>
#    define THREAD_T               pthread_t
#    define THREAD_CREATE(t, f, v) createThread(t, f, v)
#    define THREAD_JOIN(t)         joinThread(t)
#    define THREAD_YIELD()         sched_yield()
#    define MUTEX_T                pthread_mutex_t
#    define MUTEX_INIT(l)          pthread_mutex_init(&l, NULL)
#    define MUTEX_LOCK(l)          pthread_mutex_lock(&l)
#    define MUTEX_UNLOCK(l)        pthread_mutex_unlock(&l)
#else
#    define THREAD_T        int
#    define THREAD_CREATE(t, f, v) (*f)(v)
#    define THREAD_JOIN(thr)
#    define THREAD_YIELD()	0
#    define MUTEX_T         int
#    define MUTEX_INIT(l)   l = 0
#    define MUTEX_LOCK(l)   l = 1
#    define MUTEX_UNLOCK(l) l = 0
#endif /* USE_SPINLOCK */

/**
 * Wrap a lock; obtain lock upon construction, release upon destruction.
 */
class ThreadSafe {
public:
	ThreadSafe(MUTEX_T* lock) {
		lock_ = lock;
		MUTEX_LOCK(*lock_);
	}
	~ThreadSafe() { MUTEX_UNLOCK(*lock_); }
private:
	MUTEX_T *lock_;
};

/**
 * Synchronized counter, capable of thread-safe inc, dec, add, sub, max
 * and min operations.
 */
class SyncCounter {
public:
	SyncCounter() : cnt_(0) {
		MUTEX_INIT(lock_);
	}
	SyncCounter& operator++() {
		ThreadSafe l(&lock_);
		cnt_++;
		return *this;
	}
	SyncCounter& operator++(int /*i*/) {
		ThreadSafe l(&lock_);
		cnt_++;
		return *this;
	}
	SyncCounter& operator--() {
		ThreadSafe l(&lock_);
		cnt_--;
		return *this;
	}
	SyncCounter& operator--(int /*i*/) {
		ThreadSafe l(&lock_);
		cnt_--;
		return *this;
	}
	SyncCounter& operator+=(int64_t c) {
		ThreadSafe l(&lock_);
		cnt_ += c;
		return *this;
	}
	SyncCounter& operator-=(int64_t c) {
		ThreadSafe l(&lock_);
		cnt_ -= c;
		return *this;
	}
	void max(int64_t c) {
		ThreadSafe l(&lock_);
		cnt_ = std::max<int64_t>(c, cnt_);
	}
	void min(int64_t c) {
		ThreadSafe l(&lock_);
		cnt_ = std::min<int64_t>(c, cnt_);
	}
	int64_t value() const {
		return cnt_;
	}
private:
	MUTEX_T lock_;
	int64_t cnt_;
};

#ifdef USE_PTHREADS
void createThread(pthread_t& thr, void *(*func)(void*), void *vp);
void joinThread(const pthread_t& thr);
#endif

#endif
