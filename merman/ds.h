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

#ifndef DS_H_
#define DS_H_

#include <stdexcept>
#include "assert_helpers.h"

/**
 * A simple fixed-length array of type T, automatically freed in the
 * destructor.
 */
template<typename T>
class AutoArray {
public:
	AutoArray(size_t sz) {
		t_ = NULL;
		t_ = new T[sz];
		memset(t_, 0, sz*sizeof(T));
		sz_ = sz;
	}
	~AutoArray() { if(t_ != NULL) delete[] t_; }
	T& operator[](size_t sz) {
		return t_[sz];
	}
	const T& operator[](size_t sz) const {
		return t_[sz];
	}
	size_t size() const { return sz_; }
private:
	T *t_;
	size_t sz_;
};

/**
 * Expandable list using only heap storage.  Useful for data structures
 * that can be recycled from read to read; the first several reads will
 * tend to force expansions of the memory as some will require
 * as-yet-unneeded amounts of memory, but eventually a steady state
 * will be reached and no more heap operations will be needed.
 *
 * For efficiency reasons, ELists should not be declared on the stack
 * in often-called worker functions; rather, they should be re-used
 * across invocations of such functions.
 */
template <typename T>
class EList {
public:

	/**
	 * Allocate initial default of 128 elements.
	 */
	EList() : list_(NULL), sz_(128), cur_(0) {
		list_ = new T[sz_];
	}

	/**
	 * Initially allocate given number of elements; should be > 0.
	 */
	EList(size_t isz) : list_(NULL), sz_(isz), cur_(0) {
		assert_gt(isz, 0);
		list_ = new T[sz_];
	}

	/**
	 * Copy from another EList.
	 */
	EList(const EList<T>& o) : list_(NULL), sz_(0), cur_(0) {
		*this = o;
	}

	/**
	 * Destructor.
	 */
	~EList() { if(list_ != NULL) delete[] list_; }

	/**
	 * Make this object into a copy of o by allocat
	 */
	EList& operator=(const EList<T>& o) {
		if(o.cur_ == 0) {
			cur_ = 0;
			return *this;
		}
		if(sz_ < o.cur_) expandNoCopy(o.cur_ + 1);
		assert_geq(sz_, o.cur_);
		cur_ = o.cur_;
		for(size_t i = 0; i < cur_; i++) {
			list_[i] = o.list_[i];
		}
		return *this;
	}

	/**
	 * Return number of elements.
	 */
	size_t size() const { return cur_; }

	/**
	 * Return true iff there are no elements.
	 */
	bool empty() const { return cur_ == 0; }

	/**
	 * Add an element to the back and immediately initialize it via
	 * operator=.
	 */
	void push_back(const T& el) {
		assert_leq(cur_, sz_);
		if(cur_ == sz_) {
			T el2 = el;
			expandCopy(sz_+1);
			list_[cur_++] = el2;
		} else {
			list_[cur_++] = el;
		}
	}

	/**
	 * Add an element to the back.  No intialization is done.
	 */
	void expand() {
		assert_leq(cur_, sz_);
		if(cur_ == sz_) expandCopy(sz_+1);
		cur_++;
	}

	/**
	 * Add an element to the back.  No intialization is done.
	 */
	void fill(size_t begin, size_t end, const T& v) {
		assert_leq(begin, end);
		assert_leq(end, cur_);
		for(size_t i = begin; i < end; i++) {
			list_[i] = v;
		}
	}

	/**
	 * If size is less than requested size, resize up to at least sz
	 * and set cur_ to requested sz.
	 */
	void resize(size_t sz) {
		if(sz <= cur_) {
			cur_ = sz;
			return;
		}
		if(sz_ < sz) expandCopy(sz);
		cur_ = sz;
	}

	/**
	 * Resize buffer to 'len' elements and set all to 'el'
	 */
	void resizeAndFill(size_t len, const T& el) {
		resize(len);
		fill(0, len, el);
	}

	/**
	 * Erase element at offset idx.
	 */
	void erase(size_t idx) {
		assert_lt(idx, cur_);
		for(size_t i = idx; i < cur_-1; i++) {
			list_[i] = list_[i+1];
		}
		cur_--;
	}

	/**
	 * Insert value 'el' at offset 'idx'
	 */
	void insert(const T& el, size_t idx) {
		assert_lt(idx, cur_);
		if(cur_ == sz_) {
			T el2 = el;
			expandCopy(sz_+1);
			for(size_t i = cur_; i > idx; i--) {
				list_[i] = list_[i-1];
			}
			list_[idx] = el2;
		} else {
			for(size_t i = cur_; i > idx; i--) {
				list_[i] = list_[i-1];
			}
			list_[idx] = el;
		}
		cur_++;
	}

	/**
	 * Insert contents of list 'l' at offset 'idx'
	 */
	void insert(const EList<T>& l, size_t idx) {
		assert_lt(idx, cur_);
		if(l.cur_ == 0) return;
		if(cur_ + l.cur_ > sz_) expandCopy(cur_ + l.cur_);
		for(size_t i = cur_ + l.cur_ - 1; i > idx + (l.cur_ - 1); i--) {
			list_[i] = list_[i - l.cur_];
		}
		for(size_t i = 0; i < l.cur_; i++) {
			list_[i+idx] = l.list_[i];
		}
		cur_ += l.cur_;
	}

	/**
	 * Remove an element from the top of the stack.
	 */
	void pop_back() {
		assert_gt(cur_, 0);
		cur_--;
	}

	/**
	 * Make the stack empty.
	 */
	void clear() {
		cur_ = 0; // re-use stack memory
		// Don't clear heap; re-use it
	}

	/**
	 * Get the element on the top of the stack.
	 */
	T& back() {
		assert_gt(cur_, 0);
		return list_[cur_-1];
	}

	/**
	 * Reverse list elements.
	 */
	void reverse() {
		if(cur_ > 1) std::reverse(list_, list_ + cur_);
	}

	/**
	 * Get the element on the top of the stack, const version.
	 */
	const T& back() const { return back(); }

	/**
	 * Get the frontmost element (bottom of stack).
	 */
	T& front() {
		assert_gt(cur_, 0);
		return list_[0];
	}

	/**
	 * Get the element on the bottom of the stack, const version.
	 */
	const T& front() const { return front(); }

	/**
	 * Return a reference to the ith element.
	 */
	T& operator[](size_t i) {
		assert_lt(i, cur_);
		return list_[i];
	}

	/**
	 * Return a reference to the ith element.
	 */
	const T& operator[](size_t i) const {
		assert_lt(i, cur_);
		return list_[i];
	}

	/**
	 * Sort contents
	 */
	void sort() {
		if(cur_ > 1) std::sort(list_, list_ + cur_);
	}

	/**
	 * Delete element at position 'idx'; slide subsequent chars up.
	 */
	void remove(size_t idx) {
		assert_lt(idx, cur_);
		assert_gt(cur_, 0);
		for(size_t i = idx; i < cur_-1; i++) {
			list_[i] = list_[i+1];
		}
		cur_--;
	}

private:

	/**
	 * Expand the list_ buffer until it has at least 'thresh' elements.
	 * Expansions are quadratic.  Copy old contents into new buffer
	 * using operator=.
	 */
	void expandCopy(size_t thresh) {
		if(thresh <= sz_) return;
		size_t newsz = (sz_ * 2)+1;
		while(newsz < thresh) newsz *= 2;
		T* tmp = new T[newsz];
		if(list_ != NULL) {
			for(size_t i = 0; i < cur_; i++) {
				// Note: operator= is used
				tmp[i] = list_[i];
			}
			delete[] list_;
		}
		list_ = tmp;
		sz_ = newsz;
	}

	/**
	 * Expand the list_ buffer until it has at least 'thresh' elements.
	 * Expansions are quadratic.  Don't copy old contents over.
	 */
	void expandNoCopy(size_t thresh) {
		if(thresh <= sz_) return;
		size_t newsz = (sz_ * 2)+1;
		while(newsz < thresh) newsz *= 2;
		T* tmp = new T[newsz];
		if(list_ != NULL) {
			delete[] list_;
		}
		list_ = tmp;
		sz_ = newsz;
		assert_gt(sz_, 0);
	}

	T *list_;
	size_t sz_, cur_;
};

/**
 * Expandable set using a heap-allocated sorted array.
 *
 * Note that the copy constructor and operator= routines perform
 * shallow copies (w/ memcpy).
 */
template <typename T>
class ELSet {
public:

	/**
	 * Allocate initial default of 128 elements.
	 */
	ELSet() : list_(NULL), sz_(128), cur_(0) {
		list_ = new T[sz_];
	}

	/**
	 * Initially allocate given number of elements; should be > 0.
	 */
	ELSet(size_t isz) : list_(NULL), sz_(isz), cur_(0) {
		assert_gt(isz, 0);
		list_ = new T[sz_];
	}

	/**
	 * Copy from another ELSet.
	 */
	ELSet(const ELSet<T>& o) : list_(NULL) {
		*this = o;
	}

	/**
	 * Destructor.
	 */
	~ELSet() {
		if(list_ != NULL) delete[] list_;
	}

	/**
	 * Copy contents of given ELSet into this ELSet.
	 */
	ELSet& operator=(const ELSet<T>& o) {
		sz_ = o.sz_;
		cur_ = o.cur_;
		if(list_ != NULL) {
			delete[] list_;
			list_ = NULL;
		}
		list_ = new T[sz_];
		memcpy(list_, o.list_, cur_ * sizeof(T));
		return *this;
	}

	/**
	 * Return number of elements.
	 */
	size_t size() const { return cur_; }

	/**
	 * Return true iff there are no elements.
	 */
	bool empty() const { return cur_ == 0; }

	/**
	 * Insert a new element into the set in sorted order.
	 */
	bool insert(const T& el) {
		size_t i = 0;
		if(cur_ == 0) {
			insert(el, 0);
			return true;
		}
		if(cur_ < 16) {
			// Linear scan
			i = scanLoBound(el);
		} else {
			// Binary search
			i = bsearchLoBound(el);
		}
		if(list_[i] == el) return false;
		insert(el, i);
		return true;
	}

	/**
	 * Calculate
	 */
	size_t lowerBound(const T& el) {
		if(cur_ < 16) {
			// Linear scan
			return scanLoBound(el);
		} else {
			// Binary search
			return bsearchLoBound(el);
		}
	}

	/**
	 * Return true iff this set contains 'el'.
	 */
	bool contains(const T& el) const {
		if(cur_ == 0) return false;
		else if(cur_ == 1) return el == list_[0];
		size_t i;
		if(cur_ < 16) {
			// Linear scan
			i = scanLoBound(el);
		} else {
			// Binary search
			i = bsearchLoBound(el);
		}
		return i != cur_ && list_[i] == el;
	}

	/**
	 * Remove element from set.
	 */
	void remove(const T& el) {
		size_t i;
		if(cur_ < 16) {
			// Linear scan
			i = scanLoBound(el);
		} else {
			// Binary search
			i = bsearchLoBound(el);
		}
		assert(i != cur_ && list_[i] == el);
		erase(i);
	}

	/**
	 * If size is less than requested size, resize up to at least sz
	 * and set cur_ to requested sz.
	 */
	void resize(size_t sz) {
		if(sz <= cur_) return;
		if(sz_ < sz) expandCopy(sz);
	}

	/**
	 * Merge the sorted lists.
	 */
//	void insert(const EList<T>& l) {
//		// TODO
//		throw 1;
//		assert_lt(idx, cur_);
//		if(l.cur_ == 0) return;
//		if(cur_ + l.cur_ > sz_) expandCopy(cur_ + l.cur_);
//		for(size_t i = cur_ + l.cur_ - 1; i > idx + (l.cur_ - 1); i--) {
//			list_[i] = list_[i - l.cur_];
//		}
//		for(size_t i = 0; i < l.cur_; i++) {
//			list_[i+idx] = l.list_[i];
//		}
//		cur_ += l.cur_;
//	}

	/**
	 * Clear set without deallocating (or setting) anything.
	 */
	void clear() { cur_ = 0; }

private:

	/**
	 * Simple linear scan that returns the index of the first element
	 * of list_ that is not less than el, or cur_ if all elements are
	 * less than el.
	 */
	size_t scanLoBound(const T& el) const {
		for(size_t i = 0; i < cur_; i++) {
			if(!(list_[i] < el)) {
				// Shouldn't be equal
				return i;
			}
		}
		return cur_;
	}

	/**
	 * Perform a binary search for the first element that is not less
	 * than 'el'.  Return cur_ if all elements are less than el.
	 */
	size_t bsearchLoBound(const T& el) const {
		size_t hi = cur_;
		size_t lo = 0;
		while(true) {
			if(lo == hi) {
				assert_eq(lo, scanLoBound(el));
				return lo;
			}
			size_t mid = lo + ((hi-lo)>>1);
			assert_neq(mid, hi);
			if(list_[mid] < el) {
				if(lo == mid) {
					assert_eq(hi, scanLoBound(el));
					return hi;
				}
				lo = mid;
			} else {
				hi = mid;
			}
		}
	}

	/**
	 * Return true if sorted, assert otherwise.
	 */
	bool sorted() const {
		if(cur_ <= 1) return true;
		for(size_t i = 0; i < cur_-1; i++) {
			assert_lt(list_[i], list_[i+1]);
		}
		return true;
	}

	/**
	 * Insert value 'el' at offset 'idx'.  It's OK to insert at cur_,
	 * which is equivalent to appending.
	 */
	void insert(const T& el, size_t idx) {
		assert_leq(idx, cur_);
		if(cur_ == sz_) expandCopy(sz_+1);
		for(size_t i = cur_; i > idx; i--) {
			list_[i] = list_[i-1];
		}
		list_[idx] = el;
		cur_++;
		assert(sorted());
	}

	/**
	 * Erase element at offset idx.
	 */
	void erase(size_t idx) {
		assert_lt(idx, cur_);
		for(size_t i = idx; i < cur_-1; i++) {
			list_[i] = list_[i+1];
		}
		cur_--;
		assert(sorted);
	}

	/**
	 * Expand the list_ buffer until it has at least 'thresh' elements.
	 * Expansions are quadratic.
	 */
	void expandCopy(size_t thresh) {
		if(thresh <= sz_) return;
		size_t newsz = sz_ * 2;
		while(newsz < thresh) newsz *= 2;
		T* tmp = new T[newsz];
		if(list_ != NULL) {
			memcpy(tmp, list_, cur_ * sizeof(T));
			delete[] list_;
		}
		list_ = tmp;
		sz_ = newsz;
	}

	T *list_;
	size_t sz_, cur_;
};

/**
 * Expandable set using a heap-allocated sorted array.
 *
 * Note that the copy constructor and operator= routines perform
 * shallow copies (w/ memcpy).
 */
template <typename K, typename V>
class ELMap {
public:

	/**
	 * Allocate initial default of 128 elements.
	 */
	ELMap() : list_(NULL), vals_(NULL), sz_(128), cur_(0) {
		list_ = new K[sz_];
		vals_ = new V[sz_];
	}

	/**
	 * Initially allocate given number of elements; should be > 0.
	 */
	ELMap(size_t isz) : list_(NULL), vals_(NULL), sz_(isz), cur_(0) {
		assert_gt(isz, 0);
		list_ = new K[sz_];
		vals_ = new V[sz_];
	}

	/**
	 * Copy from another ELSet.
	 */
	ELMap(const ELMap<K, V>& o) : list_(NULL), vals_(NULL) { *this = o; }

	/**
	 * Destructor.
	 */
	~ELMap() {
		if(list_ != NULL) delete[] list_;
		if(vals_ != NULL) delete[] vals_;
	}

	/**
	 * Copy contents of given ELSet into this ELSet.
	 */
	ELMap& operator=(const ELMap<K, V>& o) {
		sz_ = o.sz_;
		cur_ = o.cur_;
		if(list_ != NULL) {
			delete[] list_;
			list_ = NULL;
		}
		if(vals_ != NULL) {
			delete[] vals_;
			vals_ = NULL;
		}
		list_ = new K[sz_];
		vals_ = new V[sz_];
		memcpy(list_, o.list_, cur_ * sizeof(K));
		memcpy(vals_, o.vals_, cur_ * sizeof(V));
		return *this;
	}

	/**
	 * Return number of elements.
	 */
	size_t size() const { return cur_; }

	/**
	 * Return true iff there are no elements.
	 */
	bool empty() const { return cur_ == 0; }

	/**
	 * Insert a new element into the set in sorted order.
	 */
	bool insert(const K& el, const V& val) {
		if(cur_ == 0) {
			insert(el, val, 0);
			return true;
		}
		size_t i = ltBound(el);
		if(list_[i] == el) {
			// No duplicate keys
			return false;
		}
		insert(el, val, i);
		return true;
	}

	/**
	 * Calculate index of lowest-index element not less than el.
	 */
	size_t ltBound(const K& el) const {
		if(cur_ == 0) return 0;
		if(cur_ < 16) {
			// Linear scan
			return scanLoBound(el);
		} else {
			// Binary search
			return bsearchLoBound(el);
		}
	}

	/**
	 * Calculate index of lowest-index element not less or equal to el.
	 */
	size_t leqBound(const K& el) const {
		if(cur_ == 0) return 0;
		size_t i;
		if(cur_ < 16) {
			// Linear scan
			i = scanLoBound(el);
		} else {
			// Binary search
			i = bsearchLoBound(el);
		}
		if(list_[i] == el) return i+1;
		return i;
	}

	/**
	 * Calculate index of lowest-index element not less than el.
	 */
	const V& get(size_t i) const { return vals_[i]; }

	/**
	 * Return true iff this set contains 'el'.
	 */
	bool containsKey(const K& el) const {
		size_t i = ltBound(el);
		return i != cur_ && list_[i] == el;
	}

	/**
	 * Remove element from set.
	 */
	void remove(const K& el) {
		size_t i = ltBound(el);
		assert(i != cur_ && list_[i] == el);
		erase(i);
	}

	/**
	 * If size is less than requested size, resize up to at least sz
	 * and set cur_ to requested sz.
	 */
	void resize(size_t sz) {
		if(sz <= cur_) return;
		if(sz_ < sz) expandCopy(sz);
	}

	/**
	 * Clear set without deallocating (or setting) anything.
	 */
	void clear() { cur_ = 0; }

private:

	/**
	 * Simple linear scan that returns the index of the first element
	 * of list_ that is not less than el, or cur_ if all elements are
	 * less than el.
	 */
	size_t scanLoBound(const K& el) const {
		for(size_t i = 0; i < cur_; i++) {
			if(!(list_[i] < el)) {
				// Shouldn't be equal
				return i;
			}
		}
		return cur_;
	}

	/**
	 * Perform a binary search for the first element that is not less
	 * than 'el'.  Return cur_ if all elements are less than el.
	 */
	size_t bsearchLoBound(const K& el) const {
		size_t hi = cur_;
		size_t lo = 0;
		while(true) {
			if(lo == hi) {
				assert_eq(lo, scanLoBound(el));
				return lo;
			}
			size_t mid = lo + ((hi-lo)>>1);
			assert_neq(mid, hi);
			if(list_[mid] < el) {
				if(lo == mid) {
					assert_eq(hi, scanLoBound(el));
					return hi;
				}
				lo = mid;
			} else {
				hi = mid;
			}
		}
	}

	/**
	 * Return true if sorted, assert otherwise.
	 */
	bool sorted() const {
		if(cur_ <= 1) return true;
		for(size_t i = 0; i < cur_-1; i++) {
			assert_lt(list_[i], list_[i+1]);
		}
		return true;
	}

	/**
	 * Insert value 'el' at offset 'idx'.  It's OK to insert at cur_,
	 * which is equivalent to appending.
	 */
	void insert(const K& el, const V& val, size_t idx) {
		assert_leq(idx, cur_);
		if(cur_ == sz_) {
			expandCopy(sz_+1);
			K el2 = el;
			V val2 = val;
			for(size_t i = cur_; i > idx; i--) {
				list_[i] = list_[i-1];
				vals_[i] = vals_[i-1];
			}
			list_[idx] = el2;
			vals_[idx] = val2;
		} else {
			for(size_t i = cur_; i > idx; i--) {
				list_[i] = list_[i-1];
				vals_[i] = vals_[i-1];
			}
			list_[idx] = el;
			vals_[idx] = val;
		}
		cur_++;
		assert(sorted());
	}

	/**
	 * Erase element at offset idx.
	 */
	void erase(size_t idx) {
		assert_lt(idx, cur_);
		for(size_t i = idx; i < cur_-1; i++) {
			list_[i] = list_[i+1];
			vals_[i] = vals_[i+1];
		}
		cur_--;
		assert(sorted);
	}

	/**
	 * Expand the list_ buffer until it has at least 'thresh' elements.
	 * Expansions are quadratic.
	 */
	void expandCopy(size_t thresh) {
		if(thresh <= sz_) return;
		size_t newsz = sz_ * 2;
		while(newsz < thresh) newsz *= 2;
		K* ktmp = new K[newsz];
		V* vtmp = new V[newsz];
		if(list_ != NULL) {
			memcpy(ktmp, list_, cur_ * sizeof(K));
			delete[] list_;
		}
		if(vals_ != NULL) {
			memcpy(vtmp, vals_, cur_ * sizeof(V));
			delete[] vals_;
		}
		list_ = ktmp;
		vals_ = vtmp;
		sz_ = newsz;
	}

	K *list_;
	V *vals_;
	size_t sz_, cur_;
};

#endif /* DS_H_ */
