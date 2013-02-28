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

#ifndef SSTRING_H_
#define SSTRING_H_

#include <string.h>
#include <iostream>
#include <algorithm>
#include <ctype.h>
#include "alphabet.h"
#include "assert_helpers.h"

/**
 * Four kinds of strings defined here:
 *
 * SString:
 *   A fixed-length string using heap memory with size set at construction time
 *   or when install() member is called.
 *
 * S2bDnaString:
 *   Like SString, but stores a list uint32_t words where each word is divided
 *   into 16 2-bit slots interpreted as holding one A/C/G/T nucleotide each.
 *
 * TODO: S3bDnaString allowing N.  S4bDnaString allowing nucleotide masks.
 *
 * SStringExpandable:
 *   A string using heap memory where the size of the backing store is
 *   automatically resized as needed.  Supports operations like append, insert,
 *   erase, etc.
 *
 * SStringFixed:
 *   A fixed-length string using stack memory where size is set at compile
 *   time.
 *
 * All string classes have some extra facilities that make it easy to print the
 * string, including when the string uses an encoded alphabet.  See toZBuf()
 * and toZBufXForm().
 *
 * Global lt, eq, and gt template functions are supplied.  They are capable of
 * doing lexicographical comparisons between any of the three categories of
 * strings defined here.
 */

template<typename T>
class Class_sstr_len {
public:
	static inline size_t sstr_len(const T& s) {
		return s.length();
	}
};

template<unsigned N>
class Class_sstr_len<const char[N]> {
public:
	static inline size_t sstr_len(const char s[N]) {
		return std::strlen(s);
	}
};

template<>
class Class_sstr_len<const char *> {
public:
	static inline size_t sstr_len(const char *s) {
		return std::strlen(s);
	}
};

template<>
class Class_sstr_len<const unsigned char *> {
public:
	static inline size_t sstr_len(const unsigned char *s) {
		return std::strlen((const char *)s);
	}
};

template<typename T1, typename T2>
static inline bool sstr_eq(const T1& s1, const T2& s2) {
	size_t len1 = Class_sstr_len<T1>::sstr_len(s1);
	size_t len2 = Class_sstr_len<T2>::sstr_len(s2);
	if(len1 != len2) return false;
	for(size_t i = 0; i < len1; i++) {
		if(s1[i] != s2[i]) return false;
	}
	return true;
}

template<typename T1, typename T2>
static inline bool sstr_neq(const T1& s1, const T2& s2) {
	return !sstr_eq(s1, s2);
}

/**
 * Return true iff the given suffix of s1 is equal to the given suffix of s2 up
 * to upto characters.
 */
template<typename T1, typename T2>
static inline bool sstr_suf_upto_eq(
	const T1& s1, size_t suf1,
	const T2& s2, size_t suf2,
	size_t upto,
	bool endlt = true)
{
	assert_leq(suf1, Class_sstr_len<T1>::sstr_len(s1));
	assert_leq(suf2, Class_sstr_len<T2>::sstr_len(s2));
	size_t len1 = Class_sstr_len<T1>::sstr_len(s1) - suf1;
	size_t len2 = Class_sstr_len<T2>::sstr_len(s2) - suf2;
	if(len1 > upto) len1 = upto;
	if(len2 > upto) len2 = upto;
	if(len1 != len2) return false;
	for(size_t i = 0; i < len1; i++) {
		if(s1[suf1+i] != s2[suf2+i]) {
			return false;
		}
	}
	return true;
}

/**
 * Return true iff the given suffix of s1 is equal to the given suffix of s2 up
 * to upto characters.
 */
template<typename T1, typename T2>
static inline bool sstr_suf_upto_neq(
	const T1& s1, size_t suf1,
	const T2& s2, size_t suf2,
	size_t upto,
	bool endlt = true)
{
	return !sstr_suf_upto_eq(s1, suf1, s2, suf2, upto, endlt);
}

/**
 * Return true iff s1 is less than s2.
 */
template<typename T1, typename T2>
static inline bool sstr_lt(const T1& s1, const T2& s2, bool endlt = true) {
	size_t len1 = Class_sstr_len<T1>::sstr_len(s1);
	size_t len2 = Class_sstr_len<T2>::sstr_len(s2);
	size_t minlen = (len1 < len2 ? len1 : len2);
	for(size_t i = 0; i < minlen; i++) {
		if(s1[i] < s2[i]) {
			return true;
		} else if(s1[i] > s2[i]) {
			return false;
		}
	}
	if(len1 == len2) return false;
	return (len1 < len2) == endlt;
}

/**
 * Return true iff the given suffix of s1 is less than the given suffix of s2.
 */
template<typename T1, typename T2>
static inline bool sstr_suf_lt(
	const T1& s1, size_t suf1,
	const T2& s2, size_t suf2,
	bool endlt = true)
{
	assert_leq(suf1, Class_sstr_len<T1>::sstr_len(s1));
	assert_leq(suf2, Class_sstr_len<T2>::sstr_len(s2));
	size_t len1 = Class_sstr_len<T1>::sstr_len(s1) - suf1;
	size_t len2 = Class_sstr_len<T2>::sstr_len(s2) - suf2;
	size_t minlen = (len1 < len2 ? len1 : len2);
	for(size_t i = 0; i < minlen; i++) {
		if(s1[suf1+i] < s2[suf2+i]) {
			return true;
		} else if(s1[suf1+i] > s2[suf2+i]) {
			return false;
		}
	}
	if(len1 == len2) return false;
	return (len1 < len2) == endlt;
}

/**
 * Return true iff the given suffix of s1 is less than the given suffix of s2.
 * Treat s1 and s2 as though they have lengths len1/len2.
 */
template<typename T1, typename T2>
static inline bool sstr_suf_lt(
	const T1& s1, size_t suf1, size_t len1,
	const T2& s2, size_t suf2, size_t len2,
	bool endlt = true)
{
	assert_leq(suf1, len1);
	assert_leq(suf2, len2);
	size_t left1 = len1 - suf1;
	size_t left2 = len2 - suf2;
	size_t minleft = (left1 < left2 ? left1 : left2);
	for(size_t i = 0; i < minleft; i++) {
		if(s1[suf1+i] < s2[suf2+i]) {
			return true;
		} else if(s1[suf1+i] > s2[suf2+i]) {
			return false;
		}
	}
	if(left1 == left2) return false;
	return (left1 < left2) == endlt;
}

/**
 * Return true iff the given suffix of s1 is less than the given suffix of s2
 * up to upto characters.
 */
template<typename T1, typename T2>
static inline bool sstr_suf_upto_lt(
	const T1& s1, size_t suf1,
	const T2& s2, size_t suf2,
	size_t upto,
	bool endlt = true)
{
	assert_leq(suf1, Class_sstr_len<T1>::sstr_len(s1));
	assert_leq(suf2, Class_sstr_len<T2>::sstr_len(s2));
	size_t len1 = Class_sstr_len<T1>::sstr_len(s1) - suf1;
	size_t len2 = Class_sstr_len<T2>::sstr_len(s2) - suf2;
	if(len1 > upto) len1 = upto;
	if(len2 > upto) len2 = upto;
	size_t minlen = (len1 < len2 ? len1 : len2);
	for(size_t i = 0; i < minlen; i++) {
		if(s1[suf1+i] < s2[suf2+i]) {
			return true;
		} else if(s1[suf1+i] > s2[suf2+i]) {
			return false;
		}
	}
	if(len1 == len2) return false;
	return (len1 < len2) == endlt;
}

/**
 * Return true iff the given prefix of s1 is less than the given prefix of s2.
 */
template<typename T1, typename T2>
static inline bool sstr_pre_lt(
	const T1& s1, size_t pre1,
	const T2& s2, size_t pre2,
	bool endlt = true)
{
	assert_leq(pre1, Class_sstr_len<T1>::sstr_len(s1));
	assert_leq(pre2, Class_sstr_len<T2>::sstr_len(s2));
	size_t len1 = pre1;
	size_t len2 = pre2;
	size_t minlen = (len1 < len2 ? len1 : len2);
	for(size_t i = 0; i < minlen; i++) {
		if(s1[i] < s2[i]) {
			return true;
		} else if(s1[i] > s2[i]) {
			return false;
		}
	}
	if(len1 == len2) return false;
	return (len1 < len2) == endlt;
}

/**
 * Return true iff s1 is less than or equal to s2.
 */
template<typename T1, typename T2>
static inline bool sstr_leq(const T1& s1, const T2& s2, bool endlt = true) {
	size_t len1 = Class_sstr_len<T1>::sstr_len(s1);
	size_t len2 = Class_sstr_len<T2>::sstr_len(s2);
	size_t minlen = (len1 < len2 ? len1 : len2);
	for(size_t i = 0; i < minlen; i++) {
		if(s1[i] < s2[i]) {
			return true;
		} else if(s1[i] > s2[i]) {
			return false;
		}
	}
	if(len1 == len2) return true;
	return (len1 < len2) == endlt;
}

/**
 * Return true iff the given suffix of s1 is less than or equal to the given
 * suffix of s2.
 */
template<typename T1, typename T2>
static inline bool sstr_suf_leq(
	const T1& s1, size_t suf1,
	const T2& s2, size_t suf2,
	bool endlt = true)
{
	assert_leq(suf1, Class_sstr_len<T1>::sstr_len(s1));
	assert_leq(suf2, Class_sstr_len<T2>::sstr_len(s2));
	size_t len1 = Class_sstr_len<T1>::sstr_len(s1) - suf1;
	size_t len2 = Class_sstr_len<T2>::sstr_len(s2) - suf2;
	size_t minlen = (len1 < len2 ? len1 : len2);
	for(size_t i = 0; i < minlen; i++) {
		if(s1[suf1+i] < s2[suf2+i]) {
			return true;
		} else if(s1[suf1+i] > s2[suf2+i]) {
			return false;
		}
	}
	if(len1 == len2) return true;
	return (len1 < len2) == endlt;
}

/**
 * Return true iff the given prefix of s1 is less than or equal to the given
 * prefix of s2.
 */
template<typename T1, typename T2>
static inline bool sstr_pre_leq(
	const T1& s1, size_t pre1,
	const T2& s2, size_t pre2,
	bool endlt = true)
{
	assert_leq(pre1, Class_sstr_len<T1>::sstr_len(s1));
	assert_leq(pre2, Class_sstr_len<T2>::sstr_len(s2));
	size_t len1 = pre1;
	size_t len2 = pre2;
	size_t minlen = (len1 < len2 ? len1 : len2);
	for(size_t i = 0; i < minlen; i++) {
		if(s1[i] < s2[i]) {
			return true;
		} else if(s1[i] > s2[i]) {
			return false;
		}
	}
	if(len1 == len2) return true;
	return (len1 < len2) == endlt;
}

/**
 * Return true iff s1 is greater than s2.
 */
template<typename T1, typename T2>
static inline bool sstr_gt(const T1& s1, const T2& s2, bool endlt = true) {
	size_t len1 = Class_sstr_len<T1>::sstr_len(s1);
	size_t len2 = Class_sstr_len<T2>::sstr_len(s2);
	size_t minlen = (len1 < len2 ? len1 : len2);
	for(size_t i = 0; i < minlen; i++) {
		if(s1[i] > s2[i]) {
			return true;
		} else if(s1[i] < s2[i]) {
			return false;
		}
	}
	if(len1 == len2) return false;
	return (len1 > len2) == endlt;
}

/**
 * Return true iff the given suffix of s1 is greater than the given suffix of
 * s2.
 */
template<typename T1, typename T2>
static inline bool sstr_suf_gt(
	const T1& s1, size_t suf1,
	const T2& s2, size_t suf2,
	bool endlt = true)
{
	assert_leq(suf1, Class_sstr_len<T1>::sstr_len(s1));
	assert_leq(suf2, Class_sstr_len<T2>::sstr_len(s2));
	size_t len1 = Class_sstr_len<T1>::sstr_len(s1) - suf1;
	size_t len2 = Class_sstr_len<T2>::sstr_len(s2) - suf2;
	size_t minlen = (len1 < len2 ? len1 : len2);
	for(size_t i = 0; i < minlen; i++) {
		if(s1[suf1+i] > s2[suf2+i]) {
			return true;
		} else if(s1[suf1+i] < s2[suf2+i]) {
			return false;
		}
	}
	if(len1 == len2) return false;
	return (len1 > len2) == endlt;
}

/**
 * Return true iff the given prefix of s1 is greater than the given prefix of
 * s2.
 */
template<typename T1, typename T2>
static inline bool sstr_pre_gt(
	const T1& s1, size_t pre1,
	const T2& s2, size_t pre2,
	bool endlt = true)
{
	assert_leq(pre1, Class_sstr_len<T1>::sstr_len(s1));
	assert_leq(pre2, Class_sstr_len<T2>::sstr_len(s2));
	size_t len1 = pre1;
	size_t len2 = pre2;
	size_t minlen = (len1 < len2 ? len1 : len2);
	for(size_t i = 0; i < minlen; i++) {
		if(s1[i] > s2[i]) {
			return true;
		} else if(s1[i] < s2[i]) {
			return false;
		}
	}
	if(len1 == len2) return false;
	return (len1 > len2) == endlt;
}

/**
 * Return true iff s1 is greater than or equal to s2.
 */
template<typename T1, typename T2>
static inline bool sstr_geq(const T1& s1, const T2& s2, bool endlt = true) {
	size_t len1 = Class_sstr_len<T1>::sstr_len(s1);
	size_t len2 = Class_sstr_len<T2>::sstr_len(s2);
	size_t minlen = (len1 < len2 ? len1 : len2);
	for(size_t i = 0; i < minlen; i++) {
		if(s1[i] > s2[i]) {
			return true;
		} else if(s1[i] < s2[i]) {
			return false;
		}
	}
	if(len1 == len2) return true;
	return (len1 > len2) == endlt;
}

/**
 * Return true iff the given suffix of s1 is greater than or equal to the given
 * suffix of s2.
 */
template<typename T1, typename T2>
static inline bool sstr_suf_geq(
	const T1& s1, size_t suf1,
	const T2& s2, size_t suf2,
	bool endlt = true)
{
	assert_leq(suf1, Class_sstr_len<T1>::sstr_len(s1));
	assert_leq(suf2, Class_sstr_len<T2>::sstr_len(s2));
	size_t len1 = Class_sstr_len<T1>::sstr_len(s1) - suf1;
	size_t len2 = Class_sstr_len<T2>::sstr_len(s2) - suf2;
	size_t minlen = (len1 < len2 ? len1 : len2);
	for(size_t i = 0; i < minlen; i++) {
		if(s1[suf1+i] > s2[suf2+i]) {
			return true;
		} else if(s1[suf1+i] < s2[suf2+i]) {
			return false;
		}
	}
	if(len1 == len2) return true;
	return (len1 > len2) == endlt;
}

/**
 * Return true iff the given prefix of s1 is greater than or equal to the given
 * prefix of s2.
 */
template<typename T1, typename T2>
static inline bool sstr_pre_geq(
	const T1& s1, size_t pre1,
	const T2& s2, size_t pre2,
	bool endlt = true)
{
	assert_leq(pre1, Class_sstr_len<T1>::sstr_len(s1));
	assert_leq(pre2, Class_sstr_len<T2>::sstr_len(s2));
	size_t len1 = pre1;
	size_t len2 = pre2;
	size_t minlen = (len1 < len2 ? len1 : len2);
	for(size_t i = 0; i < minlen; i++) {
		if(s1[i] > s2[i]) {
			return true;
		} else if(s1[i] < s2[i]) {
			return false;
		}
	}
	if(len1 == len2) return true;
	return (len1 > len2) == endlt;
}

template<typename T>
static inline const char * sstr_to_cstr(const T& s) {
	return s.toZBuf();
}

template<>
static inline const char * sstr_to_cstr<std::basic_string<char> >(
	const std::basic_string<char>& s)
{
	return s.c_str();
}

const static int BTString_len = 1024;

/**
 * Simple string class.
 */
template<typename T, int S = 1024, int M = 2>
class SStringExpandable {
public:
	SStringExpandable() : cs_(NULL), len_(), sz_() { }

	SStringExpandable(size_t sz) : cs_(NULL), len_(), sz_(sz) {
		expandNoCopy(sz);
	}

	/**
	 * Create an SStringExpandable from another SStringExpandable.
	 */
	SStringExpandable(const SStringExpandable<T, S>& o) {
		*this = o;
	}

	/**
	 * Create an SStringExpandable from a std::basic_string of the
	 * appropriate type.
	 */
	SStringExpandable(const std::basic_string<T>& str) {
		install(str.c_str(), str.length());
	}

	/**
	 * Create an SStringExpandable from an array and size.
	 */
	SStringExpandable(const T* b, size_t sz) {
		install(b, sz);
	}

	/**
	 * Create an SStringExpandable from a zero-terminated array.
	 */
	SStringExpandable(const T* b) {
		install(b, strlen(b));
	}

	virtual ~SStringExpandable() { } // C++ needs this

	/**
	 * Assignment to other SStringFixed.
	 */
	SStringExpandable<T,S>& operator=(const SStringExpandable<T,S>& o) {
		install(o.cs_, o.len_);
		return *this;
	}

	/**
	 * Insert char c before position 'idx'; slide subsequent chars down.
	 */
	void insert(const T& c, size_t idx) {
		assert_lt(idx, len_);
		if(sz_ < len_ + 1) expandCopy((len_ + 1 + S) * M);
		// Move everyone down by 1
		for(int i = len_; i > idx; i--) {
			cs_[i] = cs_[i-1];
		}
		cs_[idx] = c;
		len_++;
	}

	/**
	 * Set character at index 'idx' to 'c'.
	 */
	void set(int c, size_t idx) {
		assert_lt(idx, len_);
		cs_[idx] = c;
	}

	/**
	 * Append char c.
	 */
	void append(const T& c) {
		if(sz_ < len_ + 1) expandCopy((len_ + 1 + S) * M);
		cs_[len_++] = c;
	}

	/**
	 * Append char c.  For compaitibility with std::string.
	 */
	void push_back(const T& c) { append(c); }

	/**
	 * Delete char at position 'idx'; slide subsequent chars up.
	 */
	void remove(size_t idx) {
		assert_lt(idx, len_);
		assert_gt(len_, 0);
		for(size_t i = idx; i < len_-1; i++) {
			cs_[i] = cs_[i+1];
		}
		len_--;
	}

	/**
	 * Retrieve constant version of element i.
	 */
	const T& operator[](size_t i) const {
		assert_lt(i, len_);
		return cs_[i];
	}

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string.
	 */
	virtual void install(const T* b, size_t sz) {
		if(sz_ < sz) {
			expandNoCopy((sz + S) * M);
			memcpy(cs_, b, sz * sizeof(T));
		} else {
			memcpy(cs_, b, sz * sizeof(T));
		}
		len_ = sz;
	}


	/**
	 * Copy all bytes from zero-terminated buffer 'b' into this string.
	 */
	void install(const T* b) { install(b, strlen(b)); }

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string, reversing them
	 * in the process.
	 */
	void installReverse(const char* b, size_t sz) {
		if(sz_ < sz) expandNoCopy((sz + S) * M);
		for(size_t i = 0; i < sz; i++) {
			cs_[i] = b[sz-i-1];
		}
		len_ = sz;
	}

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string, reversing them
	 * in the process.
	 */
	void installReverse(const SStringExpandable<T, S>& b) {
		if(sz_ < b.len_) expandNoCopy((b.len_ + S) * M);
		for(size_t i = 0; i < b.len_; i++) {
			cs_[i] = b.cs_[b.len_ - i - 1];
		}
		len_ = b.len_;
	}

	/**
	 * Reverse the buffer in place.
	 */
	void reverse() {
		for(size_t i = 0; i < (len_ >> 1); i++) {
			T tmp = cs_[i];
			cs_[i] = cs_[len_-i-1];
			cs_[len_-i-1] = tmp;
		}
	}

	/**
	 * Simply resize the buffer.  If the buffer is resized to be
	 * longer, the newly-added elements will contain garbage and should
	 * be initialized immediately.
	 */
	void resize(size_t len) {
		if(sz_ < len) expandCopy((len + S) * M);
		len_ = len;
	}

	/**
	 * Simply resize the buffer.  If the buffer is resized to be
	 * longer, new elements will be initialized with 'el'.
	 */
	void resize(size_t len, const T& el) {
		if(sz_ < len) expandCopy((len + S) * M);
		if(len > len_) {
			for(size_t i = len_; i < len; i++) {
				cs_[i] = el;
			}
		}
		len_ = len;
	}

	/**
	 * Set the first len elements of the buffer to el.
	 */
	void fill(size_t len, const T& el) {
		assert_leq(len, len_);
		for(size_t i = 0; i < len; i++) {
			cs_[i] = el;
		}
	}

	/**
	 * Set all elements of the buffer to el.
	 */
	void fill(const T& el) {
		fill(len_, el);
	}

	/**
	 * Resize buffer to 'len' elements and set all to 'el'
	 */
	void resizeAndFill(size_t len, const T& el) {
		resize(len);
		fill(el);
	}

	/**
	 * Trim len characters from the beginning of the string.
	 */
	void trimBegin(size_t len) {
		assert_leq(len, len_);
		if(len == len_) {
			len_ = 0; return;
		}
		for(size_t i = 0; i < len_-len; i++) {
			cs_[i] = cs_[i+len];
		}
		len_ -= len;
	}

	/**
	 * Trim len characters from the end of the string.
	 */
	void trimEnd(size_t len) {
		if(len >= len_) len_ = 0;
		else len_ -= len;
	}

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string.
	 */
	void append(const T* b, size_t sz) {
		if(sz_ < len_ + sz) expandCopy((len_ + sz + S) * M);
		memcpy(cs_ + len_, b, sz * sizeof(T));
		len_ += sz;
	}

	/**
	 * Copy elements from zero-terminated buffer 'b' into this string.
	 */
	void append(const T* b) {
		append(b, strlen(b));
	}

	/**
	 * Return true iff all characters in the string are printable
	 * according to isprint().
	 */
	bool isPrintable() const {
		for(size_t i = 0; i < len_; i++) {
			if(!isprint(cs_[i])) return false;
		}
		return true;
	}

	/**
	 * Return the length of the string.
	 */
	size_t length() const { return len_; }

	/**
	 * Return the length of the string (for compatibility with std::string)
	 */
	size_t size() const { return len_; }

	/**
	 * Clear the buffer.
	 */
	void clear() { len_ = 0; }

	/**
	 * Return true iff the buffer is empty.
	 */
	bool empty() const { return len_ == 0; }

	/**
	 * Put a terminator in the 'len_'th element and then return a
	 * pointer to the buffer.  Useful for printing.
	 */
	const T* toZBuf() const {
		assert_lt(len_, sz_);
		const_cast<T*>(cs_)[len_] = 0;
		return cs_;
	}

protected:
	/**
	 * Allocate new, bigger buffer and copy old contents into it.  If
	 * requested size can be accommodated by current buffer, do nothing.
	 */
	void expandCopy(size_t sz) {
		if(sz_ >= sz) return; // done!
		T *tmp = new T[sz];
		if(cs_ != NULL) {
			memcpy(tmp, cs_, sizeof(T)*len_);
			delete[] cs_;
		}
		cs_ = tmp;
		sz_ = sz;
	}

	/**
	 * Allocate new, bigger buffer.  If requested size can be
	 * accommodated by current buffer, do nothing.
	 */
	void expandNoCopy(size_t sz) {
		if(sz_ >= sz) return; // done!
		if(cs_ != NULL) delete[] cs_;
		cs_ = new T[sz];
		sz_ = sz;
	}

	T *cs_;      // heap-allocated buffer
	size_t len_; // # filled-in elements
	size_t sz_;  // size capacity of cs_
};

/**
 * Simple string class with in-object storage.
 *
 * All copies induced by, e.g., operator=, the copy constructor,
 * install() and append(), are shallow (using memcpy/sizeof).  If deep
 * copies are needed, use a different class.
 *
 * Reading from an uninitialized element results in an assert as long
 * as NDEBUG is not defined.  If NDEBUG is defined, the result is
 * undefined.
 */
template<typename T, int S = BTString_len>
class SStringFixed {
public:
	SStringFixed() : len_(0) { }

	/**
	 * Create an SStringFixed from another SStringFixed.
	 */
	SStringFixed(const SStringFixed<T, S>& o) {
		*this = o;
	}

	/**
	 * Create an SStringFixed from another SStringFixed.
	 */
	SStringFixed(const std::basic_string<T>& str) {
		install(str.c_str(), str.length());
	}

	/**
	 * Create an SStringFixed from an array and size.
	 */
	SStringFixed(const T* b, size_t sz) {
		install(b, sz);
	}

	/**
	 * Create an SStringFixed from a zero-terminated string.
	 */
	SStringFixed(const T* b) {
		install(b, strlen(b));
	}

	virtual ~SStringFixed() { } // C++ needs this

	/**
	 * Retrieve constant version of element i.
	 */
	const T& operator[](size_t i) const {
		return get(i);
	}

	/**
	 * Retrieve constant version of element i.
	 */
	const T& get(size_t i) const {
		assert_lt(i, len_);
		return cs_[i];
	}

	/**
	 * Assignment to other SStringFixed.
	 */
	SStringFixed<T,S>& operator=(const SStringFixed<T,S>& o) {
		install(o.cs_, o.len_);
		return *this;
	}

	/**
	 * Return true iff all corresponding elements between this string
	 * and the given string are equal (i.e. operator== returns true).
	 * Lengths must also be equal.
	 */
	bool operator==(const SStringFixed<T,S>& o) const {
		if(len_ != o.len_) return false;
		for(size_t i = 0; i < len_; i++) {
			if(cs_[i] != o.cs_[i]) return false;
		}
		return true;
	}

	/**
	 * Return the inverse of operator==.
	 */
	bool operator!=(const SStringFixed<T,S>& o) const {
		return !operator==(o);
	}

	/**
	 * Insert char c before position 'idx'; slide subsequent chars down.
	 */
	void insert(const T& c, size_t idx) {
		assert_lt(len_, S);
		assert_lt(idx, len_);
		// Move everyone down by 1
		for(int i = len_; i > idx; i--) {
			cs_[i] = cs_[i-1];
		}
		cs_[idx] = c;
		len_++;
	}

	/**
	 * Set character at index 'idx' to 'c'.
	 */
	void set(int c, size_t idx) {
		assert_lt(idx, len_);
		cs_[idx] = c;
	}

	/**
	 * Append char c.
	 */
	void append(const T& c) {
		assert_lt(len_, S);
		cs_[len_++] = c;
	}

	/**
	 * Delete char at position 'idx'; slide subsequent chars up.
	 */
	void remove(size_t idx) {
		assert_lt(idx, len_);
		assert_gt(len_, 0);
		for(size_t i = idx; i < len_-1; i++) {
			cs_[i] = cs_[i+1];
		}
		len_--;
	}

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string.
	 */
	virtual void install(const T* b, size_t sz) {
		assert_leq(sz, S);
		memcpy(cs_, b, sz * sizeof(T));
		len_ = sz;
	}

	/**
	 * Copy all bytes from zero-terminated buffer 'b' into this string.
	 */
	void install(const T* b) { install(b, strlen(b)); }

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string, reversing them
	 * in the process.
	 */
	void installReverse(const char* b, size_t sz) {
		assert_leq(sz, S);
		for(size_t i = 0; i < sz; i++) {
			cs_[i] = b[sz-i-1];
		}
		len_ = sz;
	}

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string, reversing them
	 * in the process.
	 */
	void installReverse(const SStringFixed<T, S>& b) {
		assert_leq(b.len_, S);
		for(size_t i = 0; i < b.len_; i++) {
			cs_[i] = b.cs_[b.len_ - i - 1];
		}
		len_ = b.len_;
	}

	/**
	 * Reverse the buffer in place.
	 */
	void reverse() {
		for(size_t i = 0; i < (len_ >> 1); i++) {
			T tmp = cs_[i];
			cs_[i] = cs_[len_-i-1];
			cs_[len_-i-1] = tmp;
		}
	}

	/**
	 * Simply resize the buffer.  If the buffer is resized to be
	 * longer, the newly-added elements will contain garbage and should
	 * be initialized immediately.
	 */
	void resize(size_t len) {
		assert_lt(len, S);
		len_ = len;
	}

	/**
	 * Simply resize the buffer.  If the buffer is resized to be
	 * longer, new elements will be initialized with 'el'.
	 */
	void resize(size_t len, const T& el) {
		assert_lt(len, S);
		if(len > len_) {
			for(size_t i = len_; i < len; i++) {
				cs_[i] = el;
			}
		}
		len_ = len;
	}

	/**
	 * Set the first len elements of the buffer to el.
	 */
	void fill(size_t len, const T& el) {
		assert_leq(len, len_);
		for(size_t i = 0; i < len; i++) {
			cs_[i] = el;
		}
	}

	/**
	 * Set all elements of the buffer to el.
	 */
	void fill(const T& el) {
		fill(len_, el);
	}

	/**
	 * Resize buffer to 'len' elements and set all to 'el'
	 */
	void resizeAndFill(size_t len, const T& el) {
		resize(len);
		fill(el);
	}

	/**
	 * Trim len characters from the beginning of the string.
	 */
	void trimBegin(size_t len) {
		assert_leq(len, len_);
		if(len == len_) {
			len_ = 0; return;
		}
		for(size_t i = 0; i < len_-len; i++) {
			cs_[i] = cs_[i+len];
		}
		len_ -= len;
	}

	/**
	 * Trim len characters from the end of the string.
	 */
	void trimEnd(size_t len) {
		if(len >= len_) len_ = 0;
		else len_ -= len;
	}

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string.
	 */
	void append(const T* b, size_t sz) {
		assert_leq(sz + len_, S);
		memcpy(cs_ + len_, b, sz * sizeof(T));
		len_ += sz;
	}

	/**
	 * Copy bytes from zero-terminated buffer 'b' into this string.
	 */
	void append(const T* b) {
		append(b, strlen(b));
	}

	/**
	 * Return true iff all characters in the string are printable
	 * according to isprint().
	 */
	bool isPrintable() const {
		for(size_t i = 0; i < len_; i++) {
			if(!isprint(cs_[i])) return false;
		}
		return true;
	}

	/**
	 * Return the length of the string.
	 */
	size_t length() const { return len_; }

	/**
	 * Clear the buffer.
	 */
	void clear() { len_ = 0; }

	/**
	 * Return true iff the buffer is empty.
	 */
	bool empty() const { return len_ == 0; }

	/**
	 * Put a terminator in the 'len_'th element and then return a
	 * pointer to the buffer.  Useful for printing.
	 */
	virtual const T* toZBuf() const {
		const_cast<T*>(cs_)[len_] = 0;
		return cs_;
	}

	/**
	 * Put a terminator in the 'len_'th element and then return a
	 * pointer to the buffer.  Useful for printing.
	 */
	const T* toZBufXForm(const char *xform) const {
		ASSERT_ONLY(size_t xformElts = strlen(xform));
		T* printcs = const_cast<T*>(printcs_);
		for(size_t i = 0; i < len_; i++) {
			assert_lt(cs_[i], (int)xformElts);
			printcs[i] = xform[cs_[i]];
		}
		printcs[len_] = 0;
		return printcs_;
	}

	/**
	 * Return a const version of the raw buffer.
	 */
	const T* buf() const { return cs_; }

	/**
	 * Return a writeable version of the raw buffer.
	 */
	T* wbuf() { return cs_; }

protected:
	T cs_[S+1]; // +1 so that we have the option of dropping in a terminating "\0"
	T printcs_[S+1]; // +1 so that we have the option of dropping in a terminating "\0"
	size_t len_;
};

template <typename T, int S>
std::ostream& operator<< (std::ostream& os, const SStringFixed<T, S>& str) {
	os << str.toZBuf();
	return os;
}

template <typename T, int S, int M>
std::ostream& operator<< (std::ostream& os, const SStringExpandable<T, S, M>& str) {
	os << str.toZBuf();
	return os;
}

template<int S = BTString_len>
class SDnaStringFixed : public SStringFixed<char, S> {
public:

	SDnaStringFixed() : SStringFixed<char, S>() { }

	/**
	 * Create an SStringFixed from another SStringFixed.
	 */
	SDnaStringFixed(const SDnaStringFixed<S>& o) :
		SStringFixed<char, S>(o) { }

	/**
	 * Create an SStringFixed from a C++ basic_string.
	 */
	SDnaStringFixed(const std::basic_string<char>& str) :
		SStringFixed<char, S>(str) { }

	/**
	 * Create an SStringFixed from an array and size.
	 */
	SDnaStringFixed(const char* b, size_t sz) :
		SStringFixed<char, S>(b, sz) { }

	/**
	 * Create an SStringFixed from a zero-terminated string.
	 */
	SDnaStringFixed(const char* b) :
		SStringFixed<char, S>(b) { }

	virtual ~SDnaStringFixed() { } // C++ needs this

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string, reverse-
	 * complementing them in the process, assuming an encoding where
	 * 0=A, 1=C, 2=G, 3=T, 4=N.
	 */
	void installReverseComp(const char* b, size_t sz) {
		assert_leq(sz, S);
		for(size_t i = 0; i < sz; i++) {
			this->cs_[i] = (b[sz-i-1] == 4 ? 4 : b[sz-i-1] ^ 3);
		}
		this->len_ = sz;
	}

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string, reverse-
	 * complementing them in the process, assuming an encoding where
	 * 0=A, 1=C, 2=G, 3=T, 4=N.
	 */
	void installReverseComp(const SDnaStringFixed<S>& b) {
		assert_leq(b.len_, S);
		for(size_t i = 0; i < b.len_; i++) {
			this->cs_[i] = (b.cs_[b.len_-i-1] == 4 ? 4 : b.cs_[b.len_-i-1] ^ 3);
		}
		this->len_ = b.len_;
	}

	/**
	 * Either reverse or reverse-complement (depending on "color") this
	 * DNA buffer in-place.
	 */
	void reverseComp(bool color = false) {
		if(color) {
			this->reverse();
		} else {
			// Swap elements on extreme ends, working from outside in
			for(size_t i = 0; i < (this->len_ >> 1); i++) {
				char tmp1 = (this->cs_[i] == 4 ? 4 : this->cs_[i] ^ 3);
				char tmp2 = (this->cs_[this->len_-i-1] == 4 ? 4 : this->cs_[this->len_-i-1] ^ 3);
				this->cs_[i] = tmp2;
				this->cs_[this->len_-i-1] = tmp1;
			}
			// If there's a middle element
			if((this->len_ & 1) != 0) {
				// Do middle element
				char tmp = this->cs_[this->len_ >> 1];
				tmp	= (tmp == 4 ? 4 : tmp ^ 3);
				this->cs_[this->len_ >> 1] = tmp;
			}
		}
	}

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string.
	 */
	virtual void install(const char* b, size_t sz) {
		assert_leq(sz, S);
		memcpy(this->cs_, b, sz);
#ifndef NDEBUG
		for(size_t i = 0; i < sz; i++) {
			assert_leq(this->cs_[i], 4);
			assert_geq(this->cs_[i], 0);
		}
#endif
		this->len_ = sz;
	}

	/**
	 * Copy buffer 'b' of ASCII DNA characters into normal DNA
	 * characters.
	 */
	virtual void installChars(const char* b, size_t sz) {
		assert_leq(sz, S);
		for(size_t i = 0; i < sz; i++) {
			assert_in(toupper(b[i]), "ACGTN-");
			this->cs_[i] = asc2dna[(int)b[i]];
			assert_geq(this->cs_[i], 0);
			assert_leq(this->cs_[i], 4);
		}
		this->len_ = sz;
	}

	/**
	 * Copy buffer 'b' of ASCII color characters into normal DNA
	 * characters.
	 */
	virtual void installColors(const char* b, size_t sz) {
		assert_leq(sz, S);
		for(size_t i = 0; i < sz; i++) {
			assert_in(b[i], "0123.");
			this->cs_[i] = asc2col[(int)b[i]];
			assert_geq(this->cs_[i], 0);
			assert_leq(this->cs_[i], 4);
		}
		this->len_ = sz;
	}

	/**
	 * Copy buffer 'b' of ASCII color or DNA characters into normal DNA
	 * characters.
	 */
	virtual void installCharsAndColors(const char* b, size_t sz) {
		assert_leq(sz, S);
		for(size_t i = 0; i < sz; i++) {
			assert_in(b[i], "ACGTN0123.");
			this->cs_[i] = asc2dnaOrCol[(int)b[i]];
			assert_geq(this->cs_[i], 0);
			assert_leq(this->cs_[i], 4);
		}
		this->len_ = sz;
	}
	
	/**
	 * Convert current contents of cs_ into 2-bit encoded characters.
	 */
	virtual void charOrColorIze() {
		for(size_t i = 0; i < this->len_; i++) {
			assert_in(this->cs_[i], "ACGTN0123.");
			this->cs_[i] = asc2dnaOrCol[(int)this->cs_[i]];
			assert_geq(this->cs_[i], 0);
			assert_leq(this->cs_[i], 4);
		}
	}

	/**
	 * Copy C++ string of ASCII DNA characters into normal DNA
	 * characters.
	 */
	virtual void installChars(const std::basic_string<char>& str) {
		installChars(str.c_str(), str.length());
	}

	/**
	 * Copy C++ string of ASCII color characters into normal DNA
	 * characters.
	 */
	virtual void installColors(const std::basic_string<char>& str) {
		installColors(str.c_str(), str.length());
	}

	/**
	 * Copy C++ string of ASCII DNA or color characters into normal DNA
	 * characters.
	 */
	virtual void installCharsAndColors(const std::basic_string<char>& str) {
		installCharsAndColors(str.c_str(), str.length());
	}

	/**
	 * Set DNA character at index 'idx' to 'c'.
	 */
	void set(int c, size_t idx) {
		assert_lt(idx, this->len_);
		assert_leq(c, 4);
		assert_geq(c, 0);
		this->cs_[idx] = c;
	}

	/**
	 * Return true iff all characters in the string are printable
	 * according to isprint().
	 */
	bool isPrintable() const {
		return true; // because we're going to translate into dna
	}

	/**
	 * Append DNA char c.
	 */
	void append(const char& c) {
		assert_lt(this->len_, S);
		assert_leq(c, 4);
		assert_geq(c, 0);
		this->cs_[this->len_++] = c;
	}

	/**
	 * Set DNA character at index 'idx' to 'c'.
	 */
	void setChar(char c, size_t idx) {
		assert_lt(idx, this->len_);
		assert_in(toupper(c), "ACGTN");
		this->cs_[idx] = asc2dna[(int)c];
	}

	/**
	 * Append DNA character.
	 */
	void appendChar(char c) {
		assert_lt(this->len_, S);
		assert_in(toupper(c), "ACGTN");
		this->cs_[this->len_++] = asc2dna[(int)c];
	}

	/**
	 * Return DNA character corresponding to element 'idx'.
	 */
	char toChar(size_t idx) const {
		assert_geq((int)this->cs_[idx], 0);
		assert_leq((int)this->cs_[idx], 4);
		return "ACGTN"[(int)this->cs_[idx]];
	}

	/**
	 * Retrieve constant version of element i.
	 */
	const char& operator[](size_t i) const {
		return this->get(i);
	}

	/**
	 * Retrieve constant version of element i.
	 */
	const char& get(size_t i) const {
		assert_lt(i, this->len_);
		assert_leq(this->cs_[i], 4);
		assert_geq(this->cs_[i], 0);
		return this->cs_[i];
	}

	/**
	 * Retrieve constant version of element i.  Don't check that it's in the
	 * range of legal DNA chars.
	 */
	const char& getNoCheck(size_t i) const {
		assert_lt(i, this->len_);
		return this->cs_[i];
	}

	/**
	 * Put a terminator in the 'len_'th element and then return a
	 * pointer to the buffer.  Useful for printing.
	 */
	virtual const char* toZBuf() const { return this->toZBufXForm("ACGTN"); }
};

typedef SStringFixed<char, BTString_len> BTString;
typedef SDnaStringFixed<BTString_len> BTDnaString;

#endif /* SSTRING_H_ */
