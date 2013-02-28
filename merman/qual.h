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

#ifndef QUAL_H_
#define QUAL_H_

#include <stdint.h>

extern unsigned char qualRounds[];
extern unsigned char solToPhred[];

/// Translate a Phred-encoded ASCII character into a Phred quality
static inline uint8_t phredcToPhredq(char c) {
	return ((uint8_t)c >= 33 ? ((uint8_t)c - 33) : 0);
}

/**
 * Convert a Solexa-scaled quality value into a Phred-scale quality
 * value.
 *
 * p = probability that base is miscalled
 * Qphred = -10 * log10 (p)
 * Qsolexa = -10 * log10 (p / (1 - p))
 * See: http://en.wikipedia.org/wiki/FASTQ_format
 *
 */
static inline uint8_t solexaToPhred(int sol) {
	assert(sol < 256);
	if(sol < -10) return 0;
	return solToPhred[sol+10];
}

class SimplePhredPenalty {
public:
	static uint8_t mmPenalty (uint8_t qual) {
		return qual;
	}
	static uint8_t delPenalty(uint8_t qual) {
		return qual;
	}
	static uint8_t insPenalty(uint8_t qual_left, uint8_t qual_right) {
		return std::max(qual_left, qual_right);
	}
};

class MaqPhredPenalty {
public:
	static uint8_t mmPenalty (uint8_t qual) {
		return qualRounds[qual];
	}
	static uint8_t delPenalty(uint8_t qual) {
		return qualRounds[qual];
	}
	static uint8_t insPenalty(uint8_t qual_left, uint8_t qual_right) {
		return qualRounds[std::max(qual_left, qual_right)];
	}
};

static inline uint8_t mmPenalty(bool maq, uint8_t qual) {
	if(maq) {
		return MaqPhredPenalty::mmPenalty(qual);
	} else {
		return SimplePhredPenalty::mmPenalty(qual);
	}
}

#endif /*QUAL_H_*/
