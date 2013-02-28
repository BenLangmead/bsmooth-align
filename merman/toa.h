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

#ifndef TOA_H_
#define TOA_H_

/**
 * C++ version char* style "itoa":
 */
template<typename T>
static char* toa10(T value, char* result) {
	// Check that base is valid
	char* out = result;
	T quotient = value;
	do {
		*out = "0123456789"[ quotient % 10 ];
		++out;
		quotient /= 10;
	} while ( quotient );

	// Only apply negative sign for base 10
	if (value < 0) *out++ = '-';
	std::reverse( result, out );

	*out = 0; // terminator
	return out;
}

/**
 * Same, but without negative check.
 */
template<typename T>
static char* toa10pos(T value, char* result) {
	// Check that base is valid
	char* out = result;
	T quotient = value;
	do {
		*out = "0123456789"[ quotient % 10 ];
		++out;
		quotient /= 10;
	} while ( quotient );

	// Only apply negative sign for base 10
	std::reverse( result, out );

	*out = 0; // terminator
	return out;
}

#endif /* TOA_H_ */
