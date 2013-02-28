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

#ifndef FILEBUF_H_
#define FILEBUF_H_

#include <cassert>
#include <stdexcept>
#include <stdio.h>
#include <stdint.h>
#include "assert_helpers.h"

/**
 * Simple wrapper for a FILE*, istream or ifstream that reads it in
 * chunks (with fread) and keeps those chunks in a buffer.  It also
 * services calls to get(), peek() and gets() from the buffer, reading
 * in additional chunks when necessary.
 */
class FileBuf {
public:
	FileBuf(size_t sz = 256 * 1024) : sz_(sz) {
		buf_ = new uint8_t[sz_];
		init();
		assert_gt(sz_, 0);
	}

	FileBuf(const FileBuf& /*o*/) {
		throw std::runtime_error("unsupported");
		assert_gt(sz_, 0);
	}
	
	FileBuf& operator=(const FileBuf& /*o*/) {
		throw std::runtime_error("unsupported");
		return *this;
		assert_gt(sz_, 0);
	}
	
	~FileBuf() {
		delete[] buf_;
		close();
	}

	FileBuf(FILE *in, size_t sz = 256 * 1024) : sz_(sz) {
		buf_ = new uint8_t[sz_];
		init();
		in_ = in;
		assert(in_ != NULL);
		assert_gt(sz_, 0);
	}

	void setFile(FILE *in) {
		init();
		in_ = in;
		assert(in_ != NULL);
	}

	bool isOpen() {
		return in_ != NULL;
	}

	void skipWhitespace() {
		while(isspace(peek())) get();
	}

	/**
	 * Close the input stream (if that's possible)
	 */
	void close() {
		if(in_ != stdin && !closed_) fclose(in_);
		closed_ = true;
	}

	/**
	 * Get the next character of input and advance.
	 */
	int get() {
		assert(in_ != NULL);
		int c = peek();
		if(c != -1) cur_++;
		return c;
	}

	/**
	 * Return true iff all input is exhausted.
	 */
	bool eof() {
		return (cur_ == bufSz_) && done_;
	}

	/**
	 * Restore state as though we just started reading the input
	 * stream.
	 */
	void reset() {
		if(in_ != NULL) rewind(in_);
		cur_ = sz_;
		bufSz_ = sz_;
		done_ = false;
	}

	/**
	 * Peek at the next character of the input stream without
	 * advancing.  Typically we can simple read it from the buffer.
	 * Occasionally we'll need to read in a new buffer's worth of data.
	 */
	int peek() {
		assert(in_ != NULL);
		assert(cur_ <= bufSz_);
		if(cur_ == bufSz_) {
			if(done_) {
				// We already exhausted the input stream
				return -1;
			}
			// Read a new buffer's worth of data
			else {
				// Get the next chunk
				assert(in_ != NULL);
				bufSz_ = fread(buf_, 1, sz_, in_);
				cur_ = 0;
				if(bufSz_ == 0) {
					// Exhausted, and we have nothing to return to the
					// caller
					done_ = true;
					return -1;
				} else if(bufSz_ < sz_) {
					// Exhausted
					done_ = true;
				}
			}
		}
		return (int)buf_[cur_];
	}

	/**
	 * Store a string of characters from the input file into 'buf',
	 * until we see a newline, EOF, or until 'len' characters have been
	 * read.
	 */
	size_t gets(char *buf, size_t len) {
		size_t stored = 0;
		while(true) {
			int c = get();
			if(c == -1) {
				// End-of-file
				buf[stored] = '\0';
				assert(done_);
				return stored;
			}
			if(stored == len-1 || c == '\n' || c == '\r') {
				// End of string
				buf[stored] = '\0';
				// Skip over all end-of-line characters
				int pc = peek();
				while(pc == '\n' || pc == '\r') {
					get(); // discard
					pc = peek();
				}
				// Next get() will be after all newline characters
				return stored;
			}
			buf[stored++] = (char)c;
		}
	}

	/**
	 * Store a string of characters from the input file into 'buf',
	 * until we see a newline, EOF, or until 'len' characters have been
	 * read.
	 */
	size_t gets(std::string& buf, size_t len = 0xffffffff) {
		size_t stored = 0;
		while(true) {
			int c = get();
			if(c == -1) {
				// End-of-file
				assert(done_);
				return stored;
			}
			if(stored == len-1 || c == '\n' || c == '\r') {
				// Skip over all end-of-line characters
				int pc = peek();
				while(pc == '\n' || pc == '\r') {
					get(); // discard
					pc = peek();
				}
				// Next get() will be after all newline characters
				return stored;
			}
			buf.push_back((char)c);
		}
	}

	/**
	 * Store a string of characters from the input file into 'buf',
	 * until 'len' characters have been read.
	 */
	size_t get(char *buf, size_t len) {
		size_t stored = 0;
		for(size_t i = 0; i < len; i++) {
			int c = get();
			if(c == -1) {
				assert(done_);
				return i;
			}
			buf[stored++] = (char)c;
		}
		return len;
	}

	/**
	 * Store a string of characters from the input file into 'buf',
	 * until 'len' characters have been read.
	 */
	size_t get(std::string& buf, size_t len) {
		for(size_t i = 0; i < len; i++) {
			int c = get();
			if(c == -1) {
				assert(done_);
				return i;
			}
			buf.push_back((char)c);
		}
		return len;
	}

	/**
	 * Gets characters until a \r or \n is reached, then skips over the
	 * run of \r and \n's.
	 */
	void skipLine() {
		while(true) {
			int c = get();
			if(c == -1) return;
			if(c == '\n' || c == '\r') {
				// Skip over all end-of-line characters
				int pc = peek();
				while(pc == '\n' || pc == '\r') {
					get(); // discard
					pc = peek();
				}
				return;
			}
		}
	}


private:

	void init() {
		in_ = NULL;
		cur_ = bufSz_ = sz_;
		done_ = false;
		closed_ = false;
		// no need to clear buf_[]
	}

	size_t    sz_;
	FILE     *in_;
	size_t    cur_;
	size_t    bufSz_;
	bool      done_;
	uint8_t  *buf_; // (large) input buffer
	bool      closed_;
};

/**
 * Wrapper for a buffered output stream that writes characters and
 * other data types.  This class is *not* synchronized; the caller is
 * responsible for synchronization.
 */
class OutFileBuf {

public:

	/**
	 * Open a new output stream to a file with given name.
	 */
	OutFileBuf(const char *out, bool binary = false) :
		name_(out), cur_(0), closed_(false)
	{
		assert(out != NULL);
		out_ = fopen(out, binary ? "wb" : "w");
		if(out_ == NULL) {
			std::cerr << "Error: Could not open alignment output file " << out << std::endl;
			throw 1;
		}
	}

	/**
	 * Open a new output stream to standard out.
	 */
	OutFileBuf() : name_("cout"), cur_(0), closed_(false) {
		out_ = stdout;
	}

	/**
	 * Open a new output stream to a file with given name.
	 */
	void setFile(const char *out, bool binary = false) {
		assert(out != NULL);
		close();
		out_ = fopen(out, binary ? "wb" : "w");
		if(out_ == NULL) {
			std::cerr << "Error: Could not open alignment output file " << out << std::endl;
			throw 1;
		}
		reset();
	}

	/**
	 * Open a new output stream to a file with given name.
	 */
	void setFile(const std::string& out, bool binary = false) {
		setFile(out.c_str(), binary);
	}

	/**
	 * Write a single character into the write buffer and, if
	 * necessary, flush.
	 */
	void write(char c) {
		assert(!closed_);
		if(cur_ == BUF_SZ) flush();
		buf_[cur_++] = c;
	}

	/**
	 * Write a c++ string up to and not including the first whitespace
	 * character.  If necessary, flush.
	 */
	void writeStringUptoSpace(const std::string& s) {
		assert(!closed_);
		for(size_t i = 0; i < s.length() && !isspace(s[i]); i++) {
			write(s[i]);
		}
		assert_leq(cur_, BUF_SZ);
	}

	/**
	 * Write a c++ string to the write buffer and, if necessary, flush.
	 * If string to be written is larger than the buffer, flush the
	 * buffer then write the string directly to the output stream.
	 */
	void writeString(const std::string& s) {
		assert(!closed_);
		size_t slen = s.length();
		if(cur_ + slen > BUF_SZ) {
			if(cur_ > 0) flush();
			if(slen >= BUF_SZ) {
				fwrite(s.c_str(), slen, 1, out_);
			} else {
				memcpy(&buf_[cur_], s.data(), slen);
				assert_eq(0, cur_);
				cur_ = slen;
			}
		} else {
			memcpy(&buf_[cur_], s.data(), slen);
			cur_ += slen;
		}
		assert_leq(cur_, BUF_SZ);
	}

	/**
	 * Write a c++ string up to and not including the first whitespace
	 * character.  If necessary, flush.
	 */
	void writeCharsUptoSpace(const char *s) {
		assert(!closed_);
		size_t i = 0;
		while(!isspace(s[i]) && s[i] != '\0') write(s[i++]);
		assert_leq(cur_, BUF_SZ);
	}

	/**
	 * Write a char * string to the write buffer and, if necessary, flush.
	 */
	void writeChars(const char * s) {
		assert(!closed_);
		size_t i = 0;
		while(s[i] != '\0') write(s[i++]);
		assert_leq(cur_, BUF_SZ);
	}

	/**
	 * Write a c++ string to the write buffer and, if necessary, flush.
	 */
	void writeChars(const char * s, size_t len) {
		assert(!closed_);
		if(cur_ + len > BUF_SZ) {
			if(cur_ > 0) flush();
			if(len >= BUF_SZ) {
				fwrite(s, len, 1, out_);
			} else {
				memcpy(&buf_[cur_], s, len);
				assert_eq(0, cur_);
				cur_ = len;
			}
		} else {
			memcpy(&buf_[cur_], s, len);
			cur_ += len;
		}
		assert_leq(cur_, BUF_SZ);
	}

	/**
	 * Write any remaining bitpairs and then close the input
	 */
	void close() {
		if(closed_) return;
		if(cur_ > 0) flush();
		closed_ = true;
		if(out_ != stdout) {
			fclose(out_);
		}
	}

	/**
	 * Reset so that the next write is as though it's the first.
	 */
	void reset() {
		cur_ = 0;
		closed_ = false;
	}

	/**
	 * Send all buffered data to the output stream.
	 */
	void flush() {
		if(cur_ == 0) return;
		if(!fwrite((const void *)buf_, cur_, 1, out_)) {
			std::cerr << "Error while flushing and closing output" << std::endl;
			throw 1;
		}
		cur_ = 0;
	}

	/**
	 * Return true iff this stream is closed.
	 */
	bool closed() const {
		return closed_;
	}

	/**
	 * Return the filename.
	 */
	const char *name() {
		return name_;
	}

private:

	static const size_t BUF_SZ = 16 * 1024;

	const char *name_;
	FILE       *out_;
	size_t      cur_;
	char        buf_[BUF_SZ]; // (large) input buffer
	bool        closed_;
};

#endif /*ndef FILEBUF_H_*/
