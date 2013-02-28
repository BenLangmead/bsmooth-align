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

#include <string>
#include "alphabet.h"
#include "read.h"
#include "toa.h"
#include "qual.h"

/* One question for colorspace reads is: do we trim a base from the 5' end by
 * default?  Or do we leave it on by default and require the user to specify
 * -5 N with N > 0?  
 *
 *
 *
 */

using namespace std;

/**
 * Given a quality string that may be on the solexa scale (solexaScale
 * = true) and may use a +64 encoding scheme (sixty4off) , convert to
 * a Phred+33 encoding scheme.
 */
static void convertQualsToPhred33(
	BTString& qual,
	bool solexaScale,
	bool sixty4off)
{
	if(!solexaScale && !sixty4off) return;
	for(size_t i = 0; i < qual.length(); i++) {
		int c = qual[i];
		c -= (sixty4off ? 64 : 33);
		if(solexaScale) c = solexaToPhred(c);
		qual.set((char)(c+33), i);
	}
}

//
// Read
//

/**
 * Clear all fields of a read.
 */
void Read::clearAll() {
	name.clear();
	seq.clear();
	qual.clear();
	origSeq.clear();
	origQual.clear();
	for(int i = 0; i < 3; i++) {
		altSeq[i].clear();
		altQual[i].clear();
	}
	color = false;
	primer = primerColor = primerQual = -1;
	assert(empty());
}

/**
 * Clear just the seq field of the read, rendering it "empty".
 */
void Read::clear() {
	seq.clear();
	assert(empty());
}

/**
 * Write the original read sequence to the given buffer.  Useful for
 * printing colorspace reads under CS:Z.
 */
void Read::writeOrigSeq(OutFileBuf& os) const {
	if(color) {
		if(primer != -1) {
			os.write(primer);
		}
		if(primerColor != -1) {
			os.write(primerColor);
		}
		for(size_t i = 0; i < origSeq.length(); i++) {
			assert_range(0, 4, (int)origSeq[i]);
			os.write("0123."[(int)origSeq[i]]);
		}
	} else {
		os.writeChars(origSeq.toZBuf());
	}
}

/**
 * Write the original read sequence to the given buffer.  Useful for
 * printing colorspace reads under CQ:Z.
 */
void Read::writeOrigQual(OutFileBuf& os) const {
	if(color && primerQual != -1) {
		os.write(primerQual);
	}
	os.writeChars(origQual.toZBuf());
}

ostream& operator << (ostream& os, const Read& e) {
	os << e.name << ":" << e.seq << ":" << e.qual;
	return os;
}

//
// StringReads
//

/**
 * Read in one additional FASTQ read from the input string.  Return
 * false iff reading is finished.
 */
bool StringReads::nextImpl(Read& r) {
	ThreadSafe t(&this->lock_);
	assert(cur_ <= reads_.size());
	if(cur_ == reads_.size()) {
		r.clear();
		return false;
	}
	readfs_.clear();
	tokenize(reads_[cur_], ":", readfs_);
	if(readfs_.size() == 2) {
		r.name = readfs_[0];
		r.origSeq.installCharsAndColors(readfs_[1]);
	} else if(readfs_.size() > 2) {
		r.name = readfs_[0];
		r.origSeq.installCharsAndColors(readfs_[1]);
		r.origQual = readfs_[2];
	} else {
		r.name.clear();
		r.origSeq.installCharsAndColors(readfs_[0]);
	}
	if(r.name.empty()) r.name = tostring(cur_ + begin_);
	if(r.origQual.empty()) r.origQual.resizeAndFill(r.origSeq.length(), 'I');
	cur_++;
	r.id = this->id++;
	r.seq = r.origSeq;
	r.qual = r.origQual;
	return true;
}

//
// FastaReads
//

/**
 * Read in one additional FASTA read from the FileBuf fb_.  Return
 * false iff reading is finished.
 */
bool FastaReads::nextImpl(Read& r) {
	ThreadSafe t(&this->lock_);
	fb_.skipWhitespace();
	if(fb_.eof()) return false;
	r.color = false;
	while(fb_.peek() == '#') fb_.skipLine();
	if(fb_.peek() != '>') {
		cerr << "Could not parse FASTA record; started with character "
		     << ((char)fb_.peek()) << endl;
		throw 1;
	}
	fb_.get();
	r.name.resize(fb_.gets(r.name.wbuf(), 1024)); // parse name
	assert(!fb_.eof());
	r.origSeq.resize(fb_.gets(r.origSeq.wbuf(), 1024)); // parse sequence
	r.origSeq.charOrColorIze();
	r.origQual.resizeAndFill(r.origSeq.length(), 'I');
	r.id = this->id++;
	r.seq = r.origSeq;
	r.qual = r.origQual;
	return true;
}

//
// FastqReads
//

/**
 * Read in one additional FASTQ read from the FileBuf fb_.  Return
 * false iff reading is finished.
 */
bool FastqReads::nextImpl(Read& r) {
	ThreadSafe t(&this->lock_);
	fb_.skipWhitespace();
	if(fb_.eof()) return false;
	if(fb_.peek() != '@') {
		cerr << "Could not parse FASTQ record; started with character "
		     << ((char)fb_.peek()) << endl;
		throw 1;
	}
	fb_.get();
	r.name.resize(fb_.gets(r.name.wbuf(), 1024)); // parse name
	assert(!fb_.eof());
	r.origSeq.resize(fb_.gets(r.origSeq.wbuf(), 1024)); // parse sequence
	assert(!fb_.eof());
	if(fb_.peek() != '+') {
		cerr << "Couldn't parse FASTQ record; second name line started "
		     << "with character " << ((char)fb_.peek()) << endl;
		throw 1;
	}
	fb_.skipLine(); // skip the second name line
	assert(!fb_.eof());
	r.origQual.resize(fb_.gets(r.origQual.wbuf(), 1024));
	convertQualsToPhred33(r.origQual, solexaScale_, sixty4off_);
	r.origSeq.charOrColorIze();
	if(r.origQual.length() != r.origSeq.length()) {
		cerr << "Length of the quality field (" << r.origQual.length()
		     << ") didn't match length of the sequence field (" << r.origSeq.length()
		     << ") for read " << r.name << endl;
		throw 1;
	}
	r.id = this->id++;
	r.seq = r.origSeq;
	r.qual = r.origQual;
	return true;
}

//
// CSFastqReads
//

/**
 * Read in one additional FASTQ read from the FileBuf fb_.  Return
 * false iff reading is finished.
 */
bool CSFastqReads::nextImpl(Read& r) {
	ThreadSafe t(&this->lock_);
	fb_.skipWhitespace();
	if(fb_.eof()) return false;
	r.color = true;
	if(fb_.peek() != '@') {
		cerr << "Could not parse CSFASTQ record; started with character "
		     << ((char)fb_.peek()) << endl;
		throw 1;
	}
	fb_.get();
	r.name.resize(fb_.gets(r.name.wbuf(), 1024)); // parse name
	assert(!fb_.eof());
	r.origSeq.resize(fb_.gets(r.origSeq.wbuf(), 1024)); // parse sequence
	bool foundPrimer = false;
	if(isalpha(r.origSeq.getNoCheck(0)) && isdigit(r.origSeq.getNoCheck(1))) {
		// Remove primer base and initial color
		foundPrimer = true;
		r.primer = r.origSeq.getNoCheck(0);
		r.primerColor = r.origSeq.getNoCheck(1);
		r.origSeq.trimBegin(2);
	} else {
		r.primer = r.primerColor = -1;
	}
	r.origSeq.charOrColorIze();
	assert(!fb_.eof());
	if(fb_.peek() != '+') {
		cerr << "Couldn't parse FASTQ record; second name line started "
		     << "with character " << ((char)fb_.peek()) << endl;
		throw 1;
	}
	fb_.skipLine(); // skip the second name line
	assert(!fb_.eof());
	r.origQual.resize(fb_.gets(r.origQual.wbuf(), 1024));
	convertQualsToPhred33(r.origQual, solexaScale_, sixty4off_);
	if(foundPrimer && r.origQual.length() == r.origSeq.length()+1) {
		r.primerQual = r.origQual[0];
		r.origQual.trimBegin(1);
	} else {
		r.primerQual = -1;
	}
	if(r.origQual.length() != r.origSeq.length()) {
		cerr << "Length of the quality field (" << r.origQual.length()
		     << ") didn't match length of the sequence field (" << r.origSeq.length()
		     << ") for read " << r.name << endl;
		throw 1;
	}
	r.id = this->id++;
	r.seq = r.origSeq;
	r.qual = r.origQual;
	return true;
}

//
// CSFastaReads
//

/**
 * Read in one additional CSFASTA read from the FileBuf fb_.  Return
 * false iff reading is finished.
 */
bool CSFastaReads::nextImpl(Read& r) {
	ThreadSafe t(&this->lock_);
	fb_.skipWhitespace();
	if(fb_.eof()) return false;
	r.clearAll();
	r.color = true;
	// Skip comment lines
	while(fb_.peek() == '#') fb_.skipLine();
	if(fb_.peek() != '>') {
		cerr << "Could not parse CSFASTA record; started with character "
		     << ((char)fb_.peek()) << endl;
		throw 1;
	}
	fb_.get();
	r.name.resize(fb_.gets(r.name.wbuf(), 1024)); // parse name
	assert(!fb_.eof());
	if(!isDna(fb_.peek())) {
		string s;
		fb_.gets(s);
		cerr << "CSFASTA read didn't start with a primer nucleotide:"
		     << endl << s;
		throw 1;
	}
	r.primer = fb_.get(); // skip primer base
	r.origSeq.resize(fb_.gets(r.origSeq.wbuf(), 1024)); // parse name
	if(isalpha(r.origSeq.getNoCheck(0)) && isdigit(r.origSeq.getNoCheck(1))) {
		// Remove primer base and initial color
		r.primer = r.origSeq.getNoCheck(0);
		r.primerColor = r.origSeq.getNoCheck(1);
		r.origSeq.trimBegin(2);
	} else {
		r.primer = r.primerColor = -1;
	}
	r.origSeq.charOrColorIze();
	r.origQual.resizeAndFill(r.origSeq.length(), 'I');
	r.id = this->id++;
	r.seq = r.origSeq;
	r.qual = r.origQual;
	return true;
}

//
// CSFastaAndQVReads
//


/**
 * Read in one additional CSFASTA read and quality string from the
 * FileBufs ffb_ and qfb_.  Return false iff reading is finished.
 */
bool CSFastaAndQVReads::nextImpl(Read& r) {
	ThreadSafe t(&this->lock_);
	ffb_.skipWhitespace();
	qfb_.skipWhitespace();
	if(ffb_.eof()) return false;
	if(qfb_.eof()) return false;
	r.color = true;
	// Skip comment lines
	while(ffb_.peek() == '#') ffb_.skipLine();
	while(qfb_.peek() == '#') qfb_.skipLine();
	if(ffb_.peek() != '>') {
		cerr << "Could not parse CSFASTA record; started with character "
		     << ((char)ffb_.peek()) << endl;
		throw 1;
	}
	if(qfb_.peek() != '>') {
		cerr << "Could not parse QV.qual record; started with character "
		     << ((char)qfb_.peek()) << endl;
		throw 1;
	}
	ffb_.get(); // skip >
	r.name.resize(ffb_.gets(r.name.wbuf(), 1024));
#ifndef NDEBUG
	{
		string tmp;
		qfb_.get(); // skip >
		qfb_.gets(tmp); // parse name
		assert(strcmp(tmp.c_str(), r.name.toZBuf()) == 0);
	}
#else
	qfb_.skipLine();
#endif
	assert(!ffb_.eof());
	assert(!qfb_.eof());
	if(!isDna(ffb_.peek())) {
		string s;
		ffb_.gets(s);
		cerr << "CSFASTA read didn't start with a primer nucleotide:"
		     << endl << s;
		throw 1;
	}
	r.primer = ffb_.get(); // skip primer base
	r.origSeq.resize(ffb_.gets(r.origSeq.wbuf(), 1024)); // parse sequence
	bool foundPrimer = false;
	if(isalpha(r.origSeq.getNoCheck(0)) && isdigit(r.origSeq.getNoCheck(1))) {
		// Remove primer base and initial color
		foundPrimer = true;
		r.primer = r.origSeq.getNoCheck(0);
		r.primerColor = r.origSeq.getNoCheck(1);
		r.origSeq.trimBegin(2);
	} else {
		r.primer = r.primerColor = -1;
	}
	r.origSeq.charOrColorIze();
	int c = qfb_.peek();
	while(c != '\n' && c != '\r') {
		bool neg = false;
		int num = 0;
		while(!isspace(c = qfb_.get()) && !qfb_.eof()) {
			if(c == '-') {
				neg = true;
				assert_eq(num, 0);
			} else {
				assert(isdigit(c));
				num *= 10;
				num += (c - '0');
			}
		}
		if(neg) num = 0;
		// Phred-33 ASCII encode it and add it to the back of the
		// quality string
		r.origQual.append('!' + num);
		if(c == '\t' || c == '\n') break;
		// Skip over next stretch of whitespace
		c = qfb_.peek();
		while(c != '\r' && c != '\n' && isspace(c) && !qfb_.eof()) {
			qfb_.get();
			c = qfb_.peek();
		}
	}
	if(foundPrimer && r.origQual.length() == r.origSeq.length()+1) {
		r.primerQual = r.origQual[0];
		r.origQual.trimBegin(1);
	} else {
		r.primerQual = -1;
	}
	if(r.origQual.length() != r.origSeq.length()) {
		cerr << "Length of the quality field (" << r.origQual.length()
		     << ") didn't match length of the sequence field (" << r.origSeq.length()
		     << ") for read " << r.name << endl;
		throw 1;
	}
	r.id = this->id++;
	r.seq = r.origSeq;
	r.qual = r.origQual;
	return true;
}

//
// ChainReads
//


/**
 * Read in one additional read and chained HitSet from the FileBuf fb_.
 * Return false iff reading is finished.
 */
bool ChainReads::nextImpl(Read& r) {
//	ThreadSafe t(&this->lock_);
//	fb_.peek();
//	if(fb_.eof()) return false;
//	do {
//		r.hitset.deserialize(fb_);
//	} while(!r.hitset.initialized() && !fb_.eof());
//	if(!r.hitset.initialized()) return false;
//	// Copy the name/sequence/quals into r.name/r.patFw/r.qualFw
//	r.id = this->id++;
	return true;
}

/**
 * Construct FastaContinuousReads given filename, beginning index,
 * window length, and window frequency.
 */
FastaContinuousReads::FastaContinuousReads(const string& f,
                                           int begin,
                                           size_t length,
                                           size_t freq,
                                           bool bisulfite,
                                           bool rc,
                                           bool colorize)
{
	fname_ = f;
	begin_ = begin;
	FILE *fd;
	if(f == "-") {
		fd = stdin;
	} else {
		fd = fopen(fname_.c_str(), "r");
	}
	if(fd == NULL) {
		cerr << "Couldn't open fasta file " << f << " for reading" << endl;
		throw 1;
	}
	fb_.setFile(fd);
	resetForNextRecord();
	length_ = length;
	freq_ = freq;
	eat_ = length_ - 1;
	beginning_ = true;
	nameChars_ = bufCur_ = 0llu;
	bisulfite_ = bisulfite;
	rc_ = rc;
	colorize_ = colorize;
}

/**
 * Get the next read window.
 */
bool FastaContinuousReads::nextImpl(Read& r) {
	ThreadSafe t(&this->lock_);
	while(true) {
		r.clearAll();
		fb_.peek();
		if(fb_.eof()) return false;
		int c = fb_.get();
		if(c == '>') {
			resetForNextRecord();
			c = fb_.peek();
			bool sawSpace = false;
			while(c != '\n' && c != '\r') {
				if(!sawSpace) {
					sawSpace = isspace(c);
				}
				if(!sawSpace) {
					r.name.append(c);
				}
				fb_.get();
				c = fb_.peek();
			}
			while(c == '\n' || c == '\r') {
				fb_.get();
				c = fb_.peek();
			}
			r.name.append(':');
		} else {
			c = toupper(c);
			int cat = asc2dnacat[c];
			if(cat == 2) c = 'N';
			if(cat == 0) {
				// Encountered non-DNA, non-IUPAC char; skip it
				continue;
			} else {
				// DNA char
				buf_[bufCur_++] = c;
				if(bufCur_ == 1024) bufCur_ = 0;
				if(eat_ > 0) {
					eat_--;
					// Try to keep readCnt_ aligned with the offset
					// into the reference; that let's us see where
					// the sampling gaps are by looking at the read
					// name
					if(!beginning_) readCnt_++;
					continue;
				}
				for(size_t i = 0; i < length_; i++) {
					if(length_ - i <= bufCur_) {
						c = buf_[bufCur_ - (length_ - i)];
					} else {
						// Rotate
						c = buf_[bufCur_ - (length_ - i) + 1024];
					}
					r.seq.append(c);
					r.qual.append('I');
				}
				char buf[30];
				toa10pos<uint64_t>(readCnt_, buf);
				char *b = buf;
				while(*b != '\0') r.name.append(*b);
				eat_ = freq_-1;
				readCnt_++;
				beginning_ = false;
				break;
			}
		}
	}
	if(rc_) {
		// Reverse complement now
		for(size_t i = 0; i < (r.seq.length() >> 1); i++) {
			char frontc = r.seq[i];
			char frontq = r.qual[i];
			r.seq.set(comp(r.seq[r.seq.length()-i-1]), i);
			r.qual.set(r.qual[r.seq.length()-i-1], i);
			r.seq.set(comp(frontc), r.seq.length()-i-1);
			r.qual.set(frontq, r.seq.length()-i-1);
			if((r.seq.length() & 1) != 0) {
				r.seq.set(comp(r.seq[r.seq.length() >> 1]), r.seq.length() >> 1);
			}
		}
	}
	if(bisulfite_) {
		for(size_t i = 0; i < r.seq.length(); i++) {
			if(toupper(r.seq[i]) == 'C') {
				r.seq.set('t', i);
			}
		}
	}
	if(colorize_) {
		// colorize now
		r.color = true;
		// Reverse complement now
		for(size_t i = 0; i < r.seq.length()-1; i++) {
			int b1 = asc2dna[(int)r.seq[i]];
			int b2 = asc2dna[(int)r.seq[i+1]];
			int c = dinuc2color[b1][b2];
			assert_geq(c, 0);
			assert_leq(c, 4);
			r.seq.setChar("01234"[c], i);
		}
		r.seq.resize(r.seq.length()-1);
	}
	r.origSeq = r.seq;
	r.origQual = r.qual;
	return true;
}

/**
 * Reset state in anticipation for next reference sequence.
 */
void FastaContinuousReads::resetForNextRecord() {
	eat_ = length_-1;
	beginning_ = true;
	bufCur_ = nameChars_ = readCnt_ = 0;
}
