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

#ifndef ALPHABETS_H_
#define ALPHABETS_H_

#include <iostream>
#include <cassert>
#include <string>
#include <sstream>
#include <stdexcept>
#include <stdint.h>
#include "assert_helpers.h"

/// Convert an ascii char to a DNA category.  Categories are:
/// 0 -> invalid
/// 1 -> unambiguous a, c, g or t
/// 2 -> ambiguous
/// 3 -> unmatchable
extern uint8_t asc2dnacat[];
/// Convert masks to ambiguous nucleotides
extern char mask2dna[];
/// Convert ambiguous ASCII nuceleotide to mask
extern uint8_t asc2dnamask[];
/// Convert mask to # of alternative in the mask
extern int mask2popcnt[];
/// Convert an ascii char to a 2-bit base: 0=A, 1=C, 2=G, 3=T, 4=N
extern uint8_t asc2dna[];
/// Convert an ascii char representing a base or a color to a 2-bit
/// code: 0=A,0; 1=C,1; 2=G,2; 3=T,3; 4=N,.
extern uint8_t asc2dnaOrCol[];
/// Convert a pair of DNA masks to a color mask
extern uint8_t dnamasks2colormask[16][16];

/// Convert an ascii char to a color category.  Categories are:
/// 0 -> invalid
/// 1 -> unambiguous 0, 1, 2 or 3
/// 2 -> ambiguous (not applicable for colors)
/// 3 -> unmatchable
extern uint8_t asc2colcat[];
/// Convert an ascii char to a 2-bit base: 0=A, 1=C, 2=G, 3=T, 4=N
extern uint8_t asc2col[];
/// Convert an ascii char to its DNA complement, including IUPACs
extern char asc2dnacomp[];

/// Convert a pair of 2-bit (and 4=N) encoded DNA bases to a color
extern uint8_t dinuc2color[5][5];
/// Convert a 2-bit nucleotide (and 4=N) and a color to the
/// corresponding 2-bit nucleotide
extern uint8_t nuccol2nuc[5][5];
/// Convert a 4-bit mask into an IUPAC code
extern char mask2iupac[16];

/// Convert an ascii color to an ascii dna char
extern char col2dna[];
/// Convert an ascii dna to a color char
extern char dna2col[];
/// Convert an ascii dna to a color char
extern const char* dna2colstr[];

extern void setIupacsCat(uint8_t cat);

/**
 * Return true iff c is an unambiguous Dna character.
 */
static inline bool isDna(char c) {
	return asc2dnacat[(int)c] == 1;
}

/**
 * Return true iff c is an unambiguous color character (0,1,2,3).
 */
static inline bool isColor(char c) {
	return asc2colcat[(int)c] == 1;
}

/**
 * Return true iff c is an ambiguous Dna character.
 */
static inline bool isAmbigNuc(char c) {
	return asc2dnacat[(int)c] == 2;
}

/**
 * Return true iff c is an ambiguous color character.
 */
static inline bool isAmbigColor(char c) {
	return asc2colcat[(int)c] == 2;
}

/**
 * Return true iff c is an ambiguous character.
 */
static inline bool isAmbig(char c, bool color) {
	return (color ? asc2colcat[(int)c] : asc2dnacat[(int)c]) == 2;
}

/**
 * Return true iff c is an unambiguous DNA character.
 */
static inline bool isUnambigNuc(char c) {
	return asc2dnacat[(int)c] == 1;
}

/**
 * Return true iff c is an unambiguous color character.
 */
static inline bool isUnambigCol(char c) {
	return asc2colcat[(int)c] == 1;
}

/**
 * Return true iff c is an unambiguous character.
 */
static inline bool isUnambig(char c, bool color) {
	return (color ? asc2colcat[(int)c] : asc2dnacat[(int)c]) == 1;
}

/**
 * Return true iff c is an unmatchable (e.g. gap) Dna character.
 */
static inline bool isUnmatchableNuc(char c) {
	return asc2dnacat[(int)c] == 3;
}

/**
 * Return true iff c is an unmatchable (e.g. gap) Dna character.
 */
static inline bool isUnmatchableColor(char c) {
	return asc2colcat[(int)c] == 3;
}

/**
 * Return true iff c is an unmatchable (e.g. gap) Dna character.
 */
static inline bool isUnmatchable(char c, bool color) {
	return (color ? asc2dnacat[(int)c] : asc2colcat[(int)c]) == 3;
}

/**
 * M=A/C
 * R=A/G
 * W=A/T
 * S=C/G
 * Y=C/T
 * K=G/T
 * V=A/C/G
 * H=A/C/T
 * D=A/G/T
 * B=C/G/T
 * N=A/C/G/T
 */
static inline void decodeAmbigNuc(char c , int& num, int *alts) {
	switch(toupper(c)) {
	case 'M': alts[0] = 0; alts[1] = 1; num = 2; break;
	case 'R': alts[0] = 0; alts[1] = 2; num = 2; break;
	case 'W': alts[0] = 0; alts[1] = 3; num = 2; break;
	case 'S': alts[0] = 1; alts[1] = 2; num = 2; break;
	case 'Y': alts[0] = 1; alts[1] = 3; num = 2; break;
	case 'K': alts[0] = 2; alts[1] = 3; num = 2; break;
	case 'V': alts[0] = 0; alts[1] = 1; alts[2] = 2; num = 3; break;
	case 'H': alts[0] = 0; alts[1] = 1; alts[2] = 3; num = 3; break;
	case 'D': alts[0] = 0; alts[1] = 2; alts[2] = 3; num = 3; break;
	case 'B': alts[0] = 1; alts[1] = 2; alts[2] = 3; num = 3; break;
	case 'N': alts[0] = 0; alts[1] = 1; alts[2] = 2; alts[3] = 3; num = 4; break;
	default: {
		std::cerr << "Bad IUPAC code: " << c << ", (int: " << (int)c << ")" << std::endl;
		throw std::runtime_error("");
	}
	}
}

/**
 * Decode a not-necessarily-ambiguous nucleotide.
 */
static inline void decodeNuc(char c , int& num, int *alts) {
	switch(toupper(c)) {
	case 'A': alts[0] = 0; num = 1; break;
	case 'C': alts[0] = 1; num = 1; break;
	case 'G': alts[0] = 2; num = 1; break;
	case 'T': alts[0] = 3; num = 1; break;
	case 'M': alts[0] = 0; alts[1] = 1; num = 2; break;
	case 'R': alts[0] = 0; alts[1] = 2; num = 2; break;
	case 'W': alts[0] = 0; alts[1] = 3; num = 2; break;
	case 'S': alts[0] = 1; alts[1] = 2; num = 2; break;
	case 'Y': alts[0] = 1; alts[1] = 3; num = 2; break;
	case 'K': alts[0] = 2; alts[1] = 3; num = 2; break;
	case 'V': alts[0] = 0; alts[1] = 1; alts[2] = 2; num = 3; break;
	case 'H': alts[0] = 0; alts[1] = 1; alts[2] = 3; num = 3; break;
	case 'D': alts[0] = 0; alts[1] = 2; alts[2] = 3; num = 3; break;
	case 'B': alts[0] = 1; alts[1] = 2; alts[2] = 3; num = 3; break;
	case 'N': alts[0] = 0; alts[1] = 1; alts[2] = 2; alts[3] = 3; num = 4; break;
	default: {
		std::cerr << "Bad IUPAC code: " << c << ", (int: " << (int)c << ")" << std::endl;
		throw std::runtime_error("");
	}
	}
}


/**
 * Given a pair of consecutive nucleotides, one downstream (n1) one upstream
 * (n2), generate all possible compatible triplets of (a) an unambiguous
 * upstream nucleotide compatible with n1, (b) an unambiguous downstream
 * nucleotide compatible with n2, (c) the unambiguous color corresponding to
 * (a) and (b).  Then collapse those triplets down so that each color is
 * represented at most once in the cs array and all of the compatible
 * nucleotides are collapsed into a single mask
 *
 * Input arrays n1s, n2s, and cs must have enough room for 4 elements.
 */
static inline void decodeAmbigColorNucPair(
	char n1,  // nucleotide 1
	char n2,  // nucleotide 2
	int *n1s, // upstream nucleotide alternatives (parallel with n2s, cs)
	int *n2s, // downstream nucleotide alternatives (parallel with n1s, cs)
	bool *cs) // cs[i] == true iff color i is a possibility
{
	int num1, num2;
	int nucalts1[4], nucalts2[4];
	n1s[0] = n1s[1] = n1s[2] = n1s[3] = 0;
	n2s[0] = n2s[1] = n2s[2] = n2s[3] = 0;
	decodeNuc(n1, num1, nucalts1);
	decodeNuc(n2, num2, nucalts2);
	for(int i1 = 0; i1 < num1; i1++) {
		for(int i2 = 0; i2 < num2; i2++) {
			int col = dinuc2color[nucalts1[i1]][nucalts2[i2]];
			assert_range(0, 3, col);
			n1s[col] |= (1 << nucalts1[i1]);
			n2s[col] |= (1 << nucalts2[i2]);
			cs[col] = true;
		}
	}
}

/**
 * Given a series of three nucleotides, where the middle nucleotide is
 * ambiguous, calculate and return the list of corresponding color
 * pairs.
 */
static inline void decodeAmbigNucTrio(char n1, char am, char n2,
                                      int& num,
                                      int *calts1, int *calts2)
{
	int alts[] = {-1, -1, -1, -1};
	switch(am) {
		case 'M': alts[0] = 0; alts[1] = 1; num = 2; break;
		case 'R': alts[0] = 0; alts[1] = 2; num = 2; break;
		case 'W': alts[0] = 0; alts[1] = 3; num = 2; break;
		case 'S': alts[0] = 1; alts[1] = 2; num = 2; break;
		case 'Y': alts[0] = 1; alts[1] = 3; num = 2; break;
		case 'K': alts[0] = 2; alts[1] = 3; num = 2; break;
		case 'V': alts[0] = 0; alts[1] = 1; alts[2] = 2; num = 3; break;
		case 'H': alts[0] = 0; alts[1] = 1; alts[2] = 3; num = 3; break;
		case 'D': alts[0] = 0; alts[1] = 2; alts[2] = 3; num = 3; break;
		case 'B': alts[0] = 1; alts[1] = 2; alts[2] = 3; num = 3; break;
		case 'N': alts[0] = 0; alts[1] = 1; alts[2] = 2; alts[3] = 3; num = 4; break;
		default: {
			std::cerr << "Bad ambiguous nucleotide code: " << am
			          << ", (int: " << (int)am << ")" << std::endl;
			throw std::runtime_error("");
		}
	}
	for(int i = 0; i < num; i++) {
		calts1[i] = dinuc2color[asc2dna[(int)n1]][alts[i]];
		calts2[i] = dinuc2color[asc2dna[(int)n2]][alts[i]];
		assert_lt(calts1[i], 4);
		assert_lt(calts2[i], 4);
	}
}

/**
 * Given a series of three nucleotides, where the middle nucleotide is
 * ambiguous, calculate and return the list of corresponding color
 * pairs.
 */
static inline void decodeAmbigNucPair(char n1, char n2, int& num,
                                      int *calts)
{
	assert(!(isDna(n1) && isDna(n2)));
	int am = isDna(n1) ? n2 : n1;
	int notAm = asc2dna[(int)(isDna(n1) ? n1 : n2)];
	int alts[] = {-1, -1, -1, -1};
	switch(am) {
		case 'M': alts[0] = 0; alts[1] = 1; num = 2; break;
		case 'R': alts[0] = 0; alts[1] = 2; num = 2; break;
		case 'W': alts[0] = 0; alts[1] = 3; num = 2; break;
		case 'S': alts[0] = 1; alts[1] = 2; num = 2; break;
		case 'Y': alts[0] = 1; alts[1] = 3; num = 2; break;
		case 'K': alts[0] = 2; alts[1] = 3; num = 2; break;
		case 'V': alts[0] = 0; alts[1] = 1; alts[2] = 2; num = 3; break;
		case 'H': alts[0] = 0; alts[1] = 1; alts[2] = 3; num = 3; break;
		case 'D': alts[0] = 0; alts[1] = 2; alts[2] = 3; num = 3; break;
		case 'B': alts[0] = 1; alts[1] = 2; alts[2] = 3; num = 3; break;
		case 'N': alts[0] = 0; alts[1] = 1; alts[2] = 2; alts[3] = 3; num = 4; break;
		default: {
			std::cerr << "Bad ambiguous nucleotide code: " << am
			          << ", (int: " << (int)am << ")" << std::endl;
			throw std::runtime_error("");
		}
	}
	assert_lt(notAm, 4);
	for(int i = 0; i < num; i++) {
		calts[i] = dinuc2color[notAm][alts[i]];
		assert_lt(calts[i], 4);
	}
}

/**
 * Return true iff the DNA character 'c' is compatible with the
 * ambiguous DNA character 'iuc'.
 */
bool ambigCompatNuc(char iuc, char c);

/**
 * Return true iff the color character 'c' is compatible with the
 * ambiguous color character 'iuc'.
 */
bool ambigCompatColor(char iuc, char c);

/**
 * Return true iff the DNA character 'c' is compatible with the IUPAC
 * character 'iuc'.
 */
bool ambigCompat(char iuc, char c, bool color);

/**
 * OK, so it's entropy times 2.
 */
static inline int iupacAlleles(char iu) {
	switch(toupper(iu)) {
	case 'A':
	case 'C':
	case 'G':
	case 'T': return 1;
	case 'M':
	case 'R':
	case 'W':
	case 'S':
	case 'Y':
	case 'K': return 2;
	case 'V':
	case 'H':
	case 'D':
	case 'B': return 3;
	case '-':
	case 'N': return 4;
	default: {
		std::cerr << "Bad DNA char: " << iu << std::endl;
		throw std::runtime_error("");
	}
	}
}

/**
 * Return the DNA complement of the given ASCII char.
 */
static inline char compDna(int c) {
	return asc2dnacomp[c];
}

/**
 * Return the DNA complement of the given ASCII char.
 */
static inline char comp(int c) {
	assert_range(0, 4, c);
	return c == 4 ? 4 : c ^ 3;
}

/**
 * Calculate the reverse complement of 'src' and store it in 'dst',
 * assuming dst is empty to begin with.
 */
template<typename T>
void revcomp(const T& src, T& dst, bool color) {
	dst.clear();
	if(color) {
		for(int i = src.length()-1; i >= 0; i--) {
			dst.push_back(src[i]);
		}
	} else {
		for(int i = src.length()-1; i >= 0; i--) {
			dst.push_back(comp(src[i]));
			assert_range(0, 4, dst[dst.length()-1]);
		}
	}
}

/**
 * Calculate the reverse complement of 'src' and store it in 'dst',
 * assuming dst is empty to begin with.
 */
template<typename T>
void revcompDna(const T& src, T& dst, bool color) {
	dst.clear();
	if(color) {
		for(int i = (int)src.length()-1; i >= 0; i--) {
			dst.push_back(src[i]);
		}
	} else {
		for(int i = (int)src.length()-1; i >= 0; i--) {
			assert_neq(0, compDna(src[i]));
			dst.push_back(compDna(src[i]));
			assert(isprint(dst[dst.length()-1]));
		}
	}
}

/**
 * Parse some numeric data via stringstream.
 */
template<typename T>
T parse(const char *s) {
	T tmp;
	std::string str = s;
	std::istringstream ss(str);
	assert(ss.good());
	try {
		ss >> tmp;
		char ch;
		if(!ss) {
			std::cerr << "Error: number parser couldn't handle: \"" << s << "\"" << std::endl;
		} else if(ss >> ch) {
			std::cerr << "Error: number parser couldn't handle: \"" << s << "\"; excess characters" << std::endl;
		}
	} catch(std::exception& e) {
		std::cerr << "Warning: bad argument value: \"" << s << "\"" << std::endl;
	}
	return tmp;
}

/**
 * Parse some numeric data via stringstream.
 */
template<typename T>
std::string tostring(const T& i) {
	std::ostringstream ss;
	ss << i;
	return ss.str();
}

/**
 * Calculate the reverse complement of 'src' and return it.
 */
template<typename T>
T revcomp(const T& src, bool color) {
	T dst;
	revcomp(src, dst, color);
	return dst;
}

#endif /*ALPHABETS_H_*/
