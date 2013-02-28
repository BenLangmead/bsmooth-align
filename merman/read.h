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

#ifndef READ_H_
#define READ_H_

#include <iostream>
#include <string>
#include "alphabet.h"
#include "tokenize.h"
#include "filebuf.h"
#include "hit_set.h"
#include "threading.h"
#include "globals.h"
#include "random_source.h"

/**
 * Encapsulates a read, including name, sequence, qualities,
 * alternative sequences and qualities (for fuzzy reads), etc.
 */
struct Read {
	BTString    name;     // read name
	BTDnaString seq;      // read sequence: A,C,G,T,N
	BTString    qual;     // read qualities: converted to Phred-33
	BTDnaString origSeq;  // original read sequence
	BTString    origQual; // original read qualities

	// Alternative sequences and qualities, for fuzzy alignment
	BTDnaString altSeq[3];
	BTString    altQual[3];

	//HitSet  hitset; // Chained hits
	int64_t id;     // Numeric id
	bool    color;  // true -> color space, false -> nucleotide space
	int    primer;      // primer base; for SOLiD reads
	int    primerColor; // color adjacent to primer base; for SOLiD reads
	int    primerQual;  // quality adjacent to primer base; for SOLiD reads

	RandomSource rand; // Per-read pseudo-random generator

	/**
	 * Return the read sequence character at the given position,
	 * converted to DNA; if rc is true, then return the character from
	 * the reverse complement.
	 */
	char charAt(size_t i, bool fw = true) const {
		assert_lt(i, seq.length());
		char ret;
		if(color) {
			ret = (fw? seq.toChar(i) : seq.toChar(seq.length()-i-1));
		} else {
			ret = (fw? seq.toChar(i) : compDna(seq.toChar(seq.length()-i-1)));
		}
		assert_in(ret, "ACGTN");
		return ret;
	}

	/**
	 * Return the read sequence character at the given position; if rc
	 * is true, then return the character from the reverse complement.
	 */
	char at(size_t i, bool fw = true) const {
		assert_lt(i, seq.length());
		assert_range(0, 4, (int)seq[i]);
		assert_range(0, 4, (int)seq[seq.length()-i-1]);
		char ret;
		if(color) {
			ret = (fw? seq[i] : seq[seq.length()-i-1]);
		} else {
			ret = (fw? seq[i] : comp(seq[seq.length()-i-1]));
		}
		assert_range(0, 4, (int)ret);
		return ret;
	}

	/**
	 * Return the read sequence character at the given position; if rc
	 * is true, then return the character from the reverse complement.
	 */
	char charAt5p(size_t i, bool fw = true) const {
		assert_lt(i, seq.length());
		char ret;
		if(color) {
			ret = seq.toChar(i);
		} else {
			ret = (fw? seq.toChar(i) : compDna(seq.toChar(i)));
		}
		assert_in(ret, "ACGTN");
		return ret;
	}

	/**
	 * Return the quality-value character at the given position; if rc
	 * is true, then return the quality-value character from the
	 * reverse.
	 */
	char qualAt(size_t i, bool fw = true) const {
		assert_lt(i, seq.length());
		return (fw? qual[i] : qual[seq.length()-i-1]);
	}

	/**
	 * Once the read's fields have been set, this is called to
	 * "finalize" the read, which essentially initializes its
	 * pseudo-random generator.
	 */
	void finalize() {
		// Calculate a per-read random seed based on a combination of
		// the read data (incl. sequence, name, quals) and the global
		// seed
		uint32_t rseed = (randseed + 101) * 59 * 61 * 67 * 71 * 73 * 79 * 83;
		size_t qlen = seq.length();
		// Throw all the characters of the read into the random seed
		for(size_t i = 0; i < qlen; i++) {
			int p = (int)seq[i];
			assert_leq(p, 4);
			size_t off = ((i & 15) << 1);
			rseed ^= (p << off);
		}
		// Throw all the quality values for the read into the random
		// seed
		for(size_t i = 0; i < qlen; i++) {
			int p = (int)qual[i];
			assert_leq(p, 255);
			size_t off = ((i & 3) << 3);
			rseed ^= (p << off);
		}
		// Throw all the characters in the read name into the random
		// seed
		size_t namelen = name.length();
		for(size_t i = 0; i < namelen; i++) {
			int p = (int)name[i];
			assert_leq(p, 255);
			size_t off = ((i & 3) << 3);
			rseed ^= (p << off);
		}
		rand.init(rseed);
	}

	/**
	 * Make the Read "empty".
	 */
	void clear();

	/**
	 * Clear all fields in the read.
	 */
	void clearAll();

	/**
	 * Return true iff Read is empty.
	 */
	bool empty() const {
		return seq.length() == 0;
	}

	/**
	 * Trim 'amt' bases from the 3' end of the sequence.
	 */
	void trim3(size_t amt) {
		if(amt == 0) return;
		seq.trimEnd(amt);
		qual.trimEnd(amt);
	}

	/**
	 * Trim 'amt' bases from the 5' end of the sequence.
	 */
	void trim5(size_t amt) {
		if(amt == 0) return;
		seq.trimBegin(amt);
		qual.trimBegin(amt);
	}
	
	/**
	 * Check that this Read is internally consistent.  Return true if it is,
	 * Assert if it isn't.
	 */
	bool repOk() {
		assert_eq(seq.length(), qual.length());
		return true;
	}
	
	/**
	 * Write the original read sequence to the given buffer.  Useful for
	 * printing colorspace reads under CS:Z.
	 */
	void writeOrigSeq(OutFileBuf& os) const;

	/**
	 * Write the original read sequence to the given buffer.  Useful for
	 * printing colorspace reads under CQ:Z.
	 */
	void writeOrigQual(OutFileBuf& os) const;

	friend std::ostream& operator << (std::ostream& os, const Read& e);
};

/**
 * Abstract parent for classes that dish out reads.
 */
class Reads {
public:
	Reads() { MUTEX_INIT(lock_); id = 0; }
	virtual ~Reads() { }
	bool next(Read& r) {
		r.clearAll();
		//r.hitset.clearAll();
		bool ret = nextImpl(r);
		r.finalize();
		if(ret) assert(r.repOk());
		return ret;
	}
	virtual bool nextImpl(Read& r) = 0;
protected:
	int64_t id;
	MUTEX_T lock_;
};

/**
 * Dish out reads from a comma-delimited string.
 */
class StringReads : public Reads {
public:
	StringReads(const std::string& rs, int begin) {
		assert(!rs.empty());
		tokenize(rs, ",", reads_);
		cur_ = 0;
		begin_ = begin;
	}
	virtual ~StringReads() { }
	virtual bool nextImpl(Read& r);
protected:
	int begin_;
	EList<std::string> reads_; // reads
	EList<std::string> readfs_; // fields per read
	size_t cur_;
};

/**
 * Dish out reads from a FASTA file.
 */
class FastaReads : public Reads {
public:
	FastaReads(const std::string& f, int begin, size_t bufsz) :
		begin_(begin), fb_(bufsz), fname_(f)
	{
		using namespace std;
		FILE *fd;
		if(f == "-") {
			fd = stdin;
		} else {
			fd = fopen(fname_.c_str(), "r");
 		}
		if(fd == NULL) {
			cerr << "Couldn't open FASTA file " << f << " for reading" << endl;
			throw 1;
		}
		fb_.setFile(fd);
	}
	virtual ~FastaReads() { }
	virtual bool nextImpl(Read& r);
protected:
	int begin_;
	FileBuf fb_;
	std::string fname_;
};

/**
 * Dish out reads from overlapping windows of a FASTA file.
 */
class FastaContinuousReads : public Reads {
public:
	FastaContinuousReads(const std::string& f, int begin,
	                     size_t length, size_t freq, bool bisulfiteFw,
	                     bool rc, bool colorize);
	virtual ~FastaContinuousReads() { }
	virtual bool nextImpl(Read& r);
protected:

	void resetForNextRecord();

	int begin_;
	FileBuf fb_;
	std::string fname_;
	size_t length_;     /// length of reads to generate
	size_t freq_;       /// frequency to sample reads
	int policy_;        /// policy for handling Ns
	size_t eat_;        /// number of characters we need to skip before
	                    /// we have flushed all of the ambiguous or
	                    /// non-existent characters out of our read
	                    /// window
	bool beginning_;    /// skipping over the first read length?
	char buf_[1024];    /// read buffer
	char nameBuf_[1024];/// read buffer for name of fasta record being
	                    /// split into mers
	size_t nameChars_;  /// number of name characters copied into buf
	size_t bufCur_;     /// buffer cursor; points to where we should
	                    /// insert the next character
	uint64_t readCnt_;  /// number to subtract from readCnt_ to get
	                    /// the pat id to output (so it resets to 0 for
	                    /// each new sequence)
	bool bisulfite_;
	bool rc_;
	bool colorize_;
};

/**
 * Dish out reads from a colorspace-encoded FASTA file.
 */
class CSFastaReads : public Reads {
public:
	CSFastaReads(const std::string& f, int begin, size_t bufsz) :
		begin_(begin), fb_(bufsz), fname_(f)
	{
		using namespace std;
		FILE *fd;
		if(f == "-") {
			fd = stdin;
		} else {
			fd = fopen(fname_.c_str(), "r");
 		}
		if(fd == NULL) {
			cerr << "Couldn't open CSFASTA file " << f << " for reading" << endl;
			throw 1;
		}
		fb_.setFile(fd);
	}

	virtual ~CSFastaReads() { }
	virtual bool nextImpl(Read& r);
protected:
	int begin_;
	FileBuf fb_;
	std::string fname_;
};

/**
 * Dish out reads from a colorspace-encoded FASTA file.
 */
class CSFastaAndQVReads : public Reads {
public:
	CSFastaAndQVReads(const std::string& f, const std::string& q,
	                  int begin, size_t bufsz) :
		begin_(begin), ffb_(bufsz), qfb_(bufsz), fname_(f), qname_(q)
	{
		using namespace std;
		FILE *fd;
		if(f == "-") {
			fd = stdin;
		} else {
			fd = fopen(fname_.c_str(), "r");
 		}
		if(fd == NULL) {
			cerr << "Couldn't open CSFASTA file " << f << " for reading" << endl;
			throw 1;
		}
		FILE *qfd;
		if(q == "-") {
			qfd = stdin;
		} else {
			qfd = fopen(qname_.c_str(), "r");
 		}
		if(qfd == NULL) {
			cerr << "Couldn't open QV file " << q << " for reading" << endl;
			throw 1;
		}
		ffb_.setFile(fd);
		qfb_.setFile(qfd);
	}
	virtual ~CSFastaAndQVReads() { }
	virtual bool nextImpl(Read& r);
protected:
	int begin_;
	FileBuf ffb_, qfb_;
	std::string fname_, qname_;
};

/**
 * Dish out reads from a FASTQ file.
 */
class FastqReads : public Reads {
public:
	FastqReads(
		const std::string& f,
		bool solexaScale,
		bool sixty4off,
		int begin,
		size_t bufsz) :
		begin_(begin),
		fb_(bufsz),
		fname_(f),
		solexaScale_(solexaScale),
		sixty4off_(sixty4off)
	{
		using namespace std;
		FILE *fd;
		if(f == "-") {
			fd = stdin;
		} else {
			fd = fopen(fname_.c_str(), "r");
 		}
		if(fd == NULL) {
			cerr << "Couldn't open FASTQ file " << f << " for reading" << endl;
			throw 1;
		}
		fb_.setFile(fd);
	}
	virtual ~FastqReads() { }
	virtual bool nextImpl(Read& r);
protected:
	int begin_;
	FileBuf fb_;
	std::string fname_;
	bool solexaScale_;
	bool sixty4off_;
};

/**
 * Dish out reads from a colorspace-encoded FASTQ file.
 */
class CSFastqReads : public Reads {
public:
	CSFastqReads(
		const std::string& f,
		bool solexaScale,
		bool sixty4off,
		int begin,
		size_t bufsz) :
		begin_(begin),
		fb_(bufsz),
		fname_(f),
		solexaScale_(solexaScale),
		sixty4off_(sixty4off)
	{
		using namespace std;
		FILE *fd;
		if(f == "-") {
			fd = stdin;
		} else {
			fd = fopen(fname_.c_str(), "r");
 		}
		if(fd == NULL) {
			cerr << "Couldn't open CSFASTQ file " << f << " for reading" << endl;
			throw 1;
		}
		fb_.setFile(fd);
	}
	virtual ~CSFastqReads() { }
	virtual bool nextImpl(Read& r);
protected:
	int begin_;
	FileBuf fb_;
	std::string fname_;
	bool solexaScale_;
	bool sixty4off_;
};

/**
 * Dish out reads from a chaining file.
 */
class ChainReads : public Reads {
public:
	ChainReads(const std::string& f, int begin, size_t bufsz) :
		begin_(begin), fb_(bufsz), fname_(f)
	{
		using namespace std;
		FILE *fd;
		if(f == "-") {
			fd = stdin;
		} else {
			fd = fopen(fname_.c_str(), "r");
 		}
		if(fd == NULL) {
			cerr << "Couldn't open chain file " << f << " for reading" << endl;
			throw 1;
		}
		fb_.setFile(fd);
	}
	virtual ~ChainReads() { }
	virtual bool nextImpl(Read& r);
protected:
	int begin_;
	FileBuf fb_;
	std::string fname_;
};

#endif /*READ_H_*/
