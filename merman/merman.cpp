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

/*
 * merman.cpp
 *
 *  Created on: Jul 12, 2009
 *      Author: Ben Langmead
 *
 * Main driver for the Merman alignment program.  Merman is a spaced
 * seed aligner designed to handle certain problems where there are
 * mismatching query/reference character combinations that should NOT
 * be penalized as a mismatch.  This feature is useful for alignment
 * of bisulfite-treated reads, and the feature could also be for
 * SNP-aware alignment, though that is not currently well supported.
 *
 * Merman is also designed to work equally well with both
 * nucleotide-space and colorspace alignments.
 *
 * The basic spaced-seed approach taken by Merman is very similar to
 * that proposed by the authors of BSMAP [1], but extended to deal with
 * a wider range of spaced-seed types, and also extended to deal with
 * colorspace reads.
 *
 * Because Merman is typically used for alignment of bisulfite-treated
 * reads, it has facilities for generating and aligning to the various
 * bisulfite-modified strands that might be relevant depending on the
 * protocol used: BSW (bisulfite Watson), BSWR (reverse complement of
 * BSW), BSC (bisulfite Crick), and BSCR (reverse complement of BSC).
 * But of course, Merman can also be used to align non-bisulfite
 * treated reads to the usual Watson and Crick strands.
 *
 * A summary of how the various strand modes are handled:
 *
 * In non-bisulfite modes:
 
 *
 * nofw norc nowatson nocrick rcref alignfw alignrc B
 * BISC2   = 
 * BISCPG2 = 
 * BISC4   = 
 * BISCPG4 = 
 *
 * [1]: Xi Y, Li W. BSMAP: whole genome bisulfite sequence MAPping
 * program. BMC Bioinformatics. 2009 Jul 27;10:232
 */

#include <iostream>
#include <sstream>
#include <memory>
#include <limits>
#include <ctime>
#include <stdint.h>
#include <getopt.h>
#include "mer_index.h"
#include "tokenize.h"
#include "timer.h"
#include "read.h"
#include "refmap.h"
#include "annot.h"
#include "threading.h"
#include "random_source.h"
#include "aligner.h"
#include "filter.h"

using namespace std;

static void printUsage(ostream& out) {
	out << "Usage: merman [options]* <ref> <reads> [<outfile>]" << endl
	    << " <ref>     comma-separated list of DNA reference strings or FASTA files containing reference" << endl
	    << "           sequence." << endl
	    << " <reads>   comma-separated list of reads or files containing reads." << endl
	    << " <outfile> file to write alignments to; if not specified, stdout is used." << endl
	    << "Options (defaults in parentheses):" << endl
	    << " Input:" << endl
	    << "  -q/--fastq            reads are FASTQ (default)" << endl
	    << "  -f/--fasta            reads are FASTA (all qualities are set to 30)" << endl
	    << "  -F                    reads are windows in a long fasta file" << endl
	    << "  -c                    reads are specified on the command line" << endl
	   // << "  --chainin             reads are specified in a chain file" << endl
	    << "  -C                    reads are colorspace" << endl
	    << "  -s <int>              skip the first <int> reads" << endl
	    << "  -u <int>              stop aligning after first <int> reads" << endl
	    << "  -3/--trim3 <int>      trim <int> bases from the 3' end of all reads" << endl
	    << "  -5/--trim5 <int>      trim <int> bases from the 5' end of all reads" << endl
		<< "  --phred33-quals       input quals are Phred+33 (default)" << endl
		<< "  --phred64-quals       input quals are Phred+64 (same as --solexa1.3-quals)" << endl
		<< "  --solexa-quals        input quals are from GA Pipeline ver. < 1.3" << endl
		<< "  --solexa1.3-quals     input quals are from GA Pipeline ver. >= 1.3" << endl
	    << " Alignment:" << endl
	    << "  -v <int>              set # mismatches allowed end-to-end" << endl
		<< " OR:" << endl
	    << "  -n/--seedmms <int>    set # mismatches allowed in seed region" << endl
	    << "  -l/--seedlen <int>    set policy seed length" << endl
	    << "  -e/--maxerr <int>     whole-read quality ceiling" << endl
		<< "" << endl
	    << "  -L/--iseedlen <int>   set indexable seed length; should be >= min read len" << endl
	    << "  -w/--seedwidth <int>  set seed width" << endl
	   // << "  --estep <int>         amt to decrease -e ceiling per base for trimmed alignments" << endl
	    << "  --minlen <int>        minimum allowed alignment length (default: entire read)" << endl
	    << "  --nofw/--norc         do not align to forward/reverse-complement reference strand" << endl
	   // << "  -i/--requireiu        only alignments overlapping an IUPAC code are valid" << endl
	    << "  -R/--snprate          SNP rate (Phred); used in colorspace decode (def: 30)" << endl
	    << "  --bisC-2              ref = BSW + BSC; all Cs treated as Ys" << endl
	    << "  --bisCpG-2            ref = BSW + BSC; CpG Cs treated as Ys, other Cs as Ts" << endl
	    << "  --bisC-4              like --bisC-2, but includes BSWR and BSCR strands" << endl
	    << "  --bisCpG-4            like --bisCpG-2, but includes BSWR and BSCR strands" << endl
	    << "  --ignore-quals        read quality values but don't use them for alignment" << endl
		<< "  --n-inclusive         N matches everything (default: N is unmatchable)" << endl
	    << " Reporting:" << endl
	    << "  -a                    report all valid alignments" << endl
	    << "  -k <int>              report best k valid alignments (1)" << endl
	    << "  -m <int>              if read has > <int> valid alignemnts, report none" << endl
	    << "  -M <int>              like -m, but reports 1 random hit (MAPQ=0)" << endl
	    << "  --strata              stratified alignment" << endl
	   // << "  --refmap <file>       use given file to map hit offsets to other coordinates" << endl
	   // << "  --annot <file>        use given file to determine locations of SNP annotations" << endl
	   // << "  --refidx              print reference index, rather than index name" << endl
	   // << "  --fullref             when printing ref name, don't stop at first whitespace" << endl
		<< "  --no-cs-cq            don't print CS:Z and CQ:Z optiotnal flags for colorspace" << endl
	   // << "  --cost                print cost and stratum information in extra fields" << endl
	    << " Performance:" << endl
	    << "  -p/--threads <int>    number of parallel threads to use for indexing & search" << endl
	    << "  --bufsz <int>         size in KB of input buffer for read file" << endl
	    << "  --nk <N>,<K>          specify N and K for seeds; overrides --specificity" << endl
	    << "  --specificity <int>   seed specificity; 1=normal, greater=more speed, more mem" << endl
	    << " Misc:" << endl
	    << "  --progress            show periodic progress information" << endl
	    << "  --srand <int>         integer seed for pseudo-random number generation" << endl
	    << "  --sample <float>      sample reads at rate <float> (0.2=randomly skip 4/5ths)" << endl
	    << "  --rcref               for opposite strand, reverse-comp ref rather than read" << endl
	   // << "  -B/--beginoff <int>   first read/ref should be numbered <int> (0)" << endl
	    << "  -E <int1>,<int2>      mask int1-bp window if > int2 possible alleles (30,50)" << endl
	    << "  --print-color         use color in printed output" << endl
	    << "  --verbose             verbose output (for debugging)" << endl
	   // << "  --merverbose          be verbose about mer extraction (for debugging)" << endl
	    << "  -h/--help             print detailed description of tool and its options" << endl
	    ;
}

// Input formats
enum {
	INPUT_CMDLINE = 1,
	INPUT_FASTQ,
	INPUT_FASTA,
	INPUT_FASTA_CONT,
	INPUT_CHAININ,
	INPUT_CSFASTQ,
	INPUT_CSFASTA,
	INPUT_CSFASTA_AND_QV
};

// Output formats
enum {
	OUTPUT_SAM = 1,
	OUTPUT_BOWTIE
};

static AlignParams ap;
static ReferenceParams rp;

static int readLen;
static int seedWidth;
static int qualCeil;
static int begin;
static bool naiveCheck;
static bool justNaive;
static bool timing;
static bool refidx;
static bool fullref;
static bool samNoCsCq;
static bool refIsStr;
static bool color;
static bool verboseIndex;
static float estep;
static int nthreads;
static int iformat;
static int oformat;
static int specificity;
static bool fcontBis;
static bool fcontRc;
static size_t fastaContLen;
static size_t fastaContFreq;
static size_t trim3;
static size_t trim5;
static bool solexaScale;
static bool sixty4off;
static string ofile;
static string qualFile;
static const char *refmapFile;
static const char *annotFile;
static int64_t readmax;
static int64_t readskip;
static bool justBlowup;
static bool printCost;
static float sampleRate;
static size_t bufsz;
static time_t progressInterval;
static size_t updateEvery;
static pair<int,int> nk;

bool verbose;
bool merverbose;
bool quiet;
bool useColor;
bool progress;
uint32_t randseed;

static void reset() {

	int penceil[1024];
	for(int i = 0; i < 1024; i++) {
		penceil[i] = -1;
	}

	ap.init(
		2,       // seedMms
		numeric_limits<int>::max(), // e2eMms
		penceil, // penceil
		-1,      // minLen
		-1,      // adjMinLen
		true,    // alignfw
		true,    // alignrc
		false,   // requireAnnot
		false,   // disallowIupac
		true,    // maqRound
		false,   // ignoreQuals
		28,      // seedLen
		28,      // iSeedLen
		1,       // khits
		numeric_limits<int>::max(), // mhits
		false,   // msample
		false,   // strata
		30,      // snpPen
		40,      // readOpenPen
		15,      // readExPen
		40,      // refOpenPen
		15,      // refExPen
		5,       // gapBarrier
		true);   // exDecEnds

	rp.init(
		false,   // requireAnnot
		false,   // bisulfiteC
		false,   // bisulfiteCpG
		false,   // genCrick
		false,   // genRevcomps
		true,    // watsonCrickRc
		30,      // entThresh.first
		50);     // entThresh.second

	readLen = 100;
	seedWidth = 24;
	qualCeil = 70;
	begin = 0;
	verbose = false;
	merverbose = false;
	quiet = false;
	useColor = false;
	naiveCheck = false;
	justNaive = false;
	timing = false;
	refidx = false;
	fullref = false;
	samNoCsCq = false;
	refIsStr = false;
	color = false;
	verboseIndex = false;
	estep = 0.0f;
	nthreads = 1;
	iformat = INPUT_FASTQ;
	oformat = OUTPUT_SAM;
	specificity = 1;
	fcontBis = false;
	fcontRc = false;
	fastaContLen = 50;
	fastaContFreq = 1;
	trim3 = 0;
	trim5 = 0;
	solexaScale = false;
	sixty4off = false;
	ofile = "-";
	qualFile = "";
	refmapFile = NULL;
	annotFile = NULL;
	readskip = 0;
	readmax = numeric_limits<int64_t>::max();
	justBlowup = false; // quit after printing blowup?
	printCost = false; // quit after printing blowup?
	randseed = 33;
	sampleRate = 2.0f; // take all reads; if < 1.0, sample
	bufsz = 64;
	progress = false;
	progressInterval = 5;
	updateEvery = 100;
	nk = make_pair(0, 0);
}

// Counters
static SyncCounter alReads;
static SyncCounter unalReads;
static SyncCounter maxReads;
static SyncCounter totAls;
static SyncCounter totUnsampled;
static SyncCounter maxSeedHits;
static SyncCounter totSeedHits;
static SyncCounter nreads;

// Command-line argument constants
enum {
	ARG_STRATA = 256,
	ARG_VERBOSE,
	ARG_MERVERBOSE,
	ARG_PRINTCOLOR,
	ARG_QUIET,
	ARG_NOALIGNFW,
	ARG_NOALIGNRC,
	ARG_GENCRICK,
	ARG_GENREVCOMPS,
	ARG_NOMAQROUND,
	ARG_NAIVE,
	ARG_JUST_NAIVE,
	ARG_CHAININ,
	ARG_CHAINOUT,
	ARG_REFMAP,
	ARG_ANNOT,
	ARG_REFIDX,
	ARG_FULLREF,
	ARG_BISC2,   // BSW, BSC; all Cs to Ys
	ARG_BISCPG2, // BSW, BSC; CpG Cs turned to Ys, other Cs to Ts
	ARG_BISC4,   // BSW, BSC, BSWR, BSCR; all Cs to Ys
	ARG_BISCPG4, // BSW, BSC, BSWR, BSCR; CpG Cs turned to Ys, other Cs to Ts
	ARG_PHRED33,
	ARG_PHRED64,
	ARG_SOLEXA64,
	ARG_VERBOSE_INDEX,
	ARG_RCREF,
	ARG_JUSTBLOWUP,
	ARG_COST,
	ARG_ESTEP,
	ARG_MINLEN,
	ARG_SPECIFICITY,
	ARG_FRC,
	ARG_FBIS,
	ARG_SRAND,
	ARG_SAMPLE,
	ARG_PROGRESS,
	ARG_UPDATEVERY,
	ARG_PROGRESSINTERVAL,
	ARG_BUFSZ,
	ARG_IGNOREQUALS,
	ARG_NK,
	ARG_COLOR_KEEP_ENDS,
	ARG_REF_IS_STR,
	ARG_SAM_NO_CS_CQ,
	ARG_N_IUPAC,
	ARG_OLDFORMAT         // --old-output
};

static const char *short_opts = "hk:m:M:an:l:e:v:w:tiUB:L:qu:s:E:cp:yCR:G:fQ:S3:5:F:";
static struct option long_opts[] = {
	{(char*)"verbose", no_argument, 0, ARG_VERBOSE},
	{(char*)"merverbose", no_argument, 0, ARG_MERVERBOSE},
	{(char*)"help", no_argument, 0, 'h'},
	{(char*)"fasta", no_argument, 0, 'f'},
	{(char*)"seedmms", required_argument, 0, 'n'},
	{(char*)"seedlen", required_argument, 0, 'l'},
	{(char*)"iseedlen", required_argument, 0, 'L'},
	{(char*)"seedwidth", required_argument, 0, 'w'},
	{(char*)"timing", no_argument, 0, 't'},
	{(char*)"tryhard", no_argument, 0, 'y'},
	{(char*)"requireiu", no_argument, 0, 'i'},
	{(char*)"noiumatch", no_argument, 0, 'U'},
	{(char*)"iumatch", no_argument, 0, 'M'},
	{(char*)"beginoff", required_argument, 0, 'B'},
	{(char*)"sam", no_argument, 0, 'S'},
	{(char*)"color", no_argument, 0, 'C'},
	{(char*)"fastq", no_argument, 0, 'q'},
	{(char*)"fasta", no_argument, 0, 'f'},
	{(char*)"qual", required_argument, 0, 'Q'},
	{(char*)"print-color", no_argument, 0, ARG_PRINTCOLOR},
	{(char*)"nomaqround", no_argument, 0, ARG_NOMAQROUND},
	{(char*)"nofw", no_argument, 0, ARG_NOALIGNFW},
	{(char*)"norc", no_argument, 0, ARG_NOALIGNRC},
	{(char*)"gencrick", no_argument, 0, ARG_GENCRICK},
	{(char*)"genrevcomps", no_argument, 0, ARG_GENREVCOMPS},
	{(char*)"strata", no_argument, 0, ARG_STRATA},
	{(char*)"quiet", no_argument, 0, ARG_QUIET},
	{(char*)"naive", no_argument, 0, ARG_NAIVE},
	{(char*)"justnaive", no_argument, 0, ARG_JUST_NAIVE},
	{(char*)"chainin", no_argument, 0, ARG_CHAININ},
	{(char*)"refmap", required_argument, 0, ARG_REFMAP},
	{(char*)"annot", required_argument, 0, ARG_ANNOT},
	{(char*)"refidx", no_argument, 0, ARG_REFIDX},
	{(char*)"fullref", no_argument, 0, ARG_FULLREF},
	{(char*)"rcref", no_argument, 0, ARG_RCREF},
	{(char*)"bisC-2", no_argument, 0, ARG_BISC2},
	{(char*)"bisCpG-2", no_argument, 0, ARG_BISCPG2},
	{(char*)"bisC-4", no_argument, 0, ARG_BISC4},
	{(char*)"bisCpG-4", no_argument, 0, ARG_BISCPG4},
	{(char*)"gap", required_argument, 0, 'G'},
	{(char*)"verboseindex", no_argument, 0, ARG_VERBOSE_INDEX},
	{(char*)"predict-memory", no_argument, 0, ARG_JUSTBLOWUP},
	{(char*)"snprate", required_argument, 0, 'R'},
	{(char*)"cost", no_argument, 0, ARG_COST},
	{(char*)"progress", no_argument, 0, ARG_PROGRESS},
	{(char*)"update-every", required_argument, 0, ARG_UPDATEVERY},
	{(char*)"progress-interval", required_argument, 0, ARG_PROGRESSINTERVAL},
	{(char*)"trim3", required_argument, 0, '3'},
	{(char*)"trim5", required_argument, 0, '5'},
	{(char*)"phred33-quals", no_argument, 0, ARG_PHRED33},
	{(char*)"phred64-quals", no_argument, 0, ARG_PHRED64},
	{(char*)"solexa-quals", no_argument, 0, ARG_SOLEXA64},
	{(char*)"solexa1.3-quals", no_argument, 0, ARG_PHRED64},
	{(char*)"trim5", required_argument, 0, '5'},
	{(char*)"estep", required_argument, 0, ARG_ESTEP},
	{(char*)"minlen", required_argument, 0, ARG_MINLEN},
	{(char*)"specificity", required_argument, 0, ARG_SPECIFICITY},
	{(char*)"Fbis", no_argument, 0, ARG_FBIS},
	{(char*)"Frc", no_argument, 0, ARG_FRC},
	{(char*)"srand", required_argument, 0, ARG_SRAND},
	{(char*)"sample", required_argument, 0, ARG_SAMPLE},
	{(char*)"bufsz", required_argument, 0, ARG_BUFSZ},
	{(char*)"nk", required_argument, 0, ARG_NK},
	{(char*)"ignore-quals", no_argument, 0, ARG_IGNOREQUALS},
	{(char*)"col-keepends", no_argument, 0, ARG_COLOR_KEEP_ENDS},
	{(char*)"ref-str", no_argument, 0, ARG_REF_IS_STR},
	{(char*)"no-cs-cq", no_argument, 0, ARG_SAM_NO_CS_CQ},
	{(char*)"n-inclusive", no_argument, 0, ARG_N_IUPAC},
	{(char*)"old-output", no_argument, 0, ARG_OLDFORMAT},
	{(char*)0, 0, 0, 0} // terminator
};

/**
 * Given a string, parse a pair of number separated by the given
 * delimiter.
 */
template<typename T>
pair<T, T> parsePair(const char *str, const char* delim) {
	string s(str);
	EList<string> ss;
	tokenize(s, delim, ss);
	pair<T, T> ret;
	ret.first = parse<T>(ss[0].c_str());
	ret.second = parse<T>(ss[1].c_str());
	return ret;
}

/**
 * Pretty-print a number that represents a number of bytes.
 */
static void printBytes(size_t bytes, ostream& out) {
	const char * allUnits[] = {"B", "KB", "MB", "GB", "TB"};
	int unit = 0;
	double by = (double)bytes;
	while((by / 1024) > 10) {
		by /= 1024;
		unit++;
		if(unit > 4) break;
	}
	out.setf(ios::fixed);
	out << setprecision(unit > 0 ? 1 : 0) << by << ' ' << allUnits[unit];
}

/**
 * Thread that periodically reports data about alignment progress so far.
 */
class ProgressThread {

public:

	ProgressThread() : dead_(false) {
		t0_ = tlast_ = time(0);
	}
	
	/**
	 * Create the thread and return.
	 */
	void run() {
		THREAD_CREATE(thread_, ProgressThread::startThread, this);
	}
	
	/**
	 * Tell this thread that it's dead.  It will die next time wround
	 * the loop in work().
	 */
	void kill() { dead_ = true; }
	
	/**
	 * Do not return until thread terminates.
	 */
	void join() { THREAD_JOIN(thread_); }
	
	static void reportStats(time_t elapsed = 0) {
		cerr << "# reads processed: " << nreads.value() << endl;
		cerr << "# reads with at least one reported alignment: "
		     << alReads.value() << fixed << setprecision(2)
		     << " (" << (nreads.value() > 0 ? (100.0 * (double)alReads.value() / (double)nreads.value()) : 0)
		     << "%)" << endl;
		cerr << "# reads that failed to align: "
		     << unalReads.value() << fixed << setprecision(2)
		     << " (" << (nreads.value() > 0 ? (100.0 * (double)unalReads.value() / (double)nreads.value()) : 0)
		     << "%)" << endl;
		if(maxReads.value() > 0) {
			cerr << "# reads with alignments suppressed due to -m: "
			     << maxReads.value() << fixed << setprecision(2)
			     << " (" << (nreads.value() > 0 ? (100.0 * (double)maxReads.value() / (double)nreads.value()) : 0)
			     << "%)" << endl;
		}
		if(totUnsampled.value() > 0) {
			cerr << "# reads not sampled due to --sample: "
			     << totUnsampled.value()  << endl;
		}
		cerr << "Max seed hits for a read: " << maxSeedHits.value() << endl;
		cerr << "Average seed hits per read: " << (nreads.value() > 0 ? (totSeedHits.value()/nreads.value()) : 0) << endl;
		if(elapsed > 0) {
			cerr << "Per second: reads: "
			     << fixed << setprecision(2)
			     << ((double)nreads.value()/elapsed)
			     << ", seed hits: "
			     << fixed << setprecision(2)
			     << ((double)totSeedHits.value()/elapsed)
			     << endl;
		}
		cerr << "Reported " << totAls.value()
		     << " alignments to 1 output stream(s)" << endl;
	}

protected:

	/**
	 * Do the periodic progress update.
	 */
	void work() {
		while(!dead_) {
			int ret;
			if((ret = THREAD_YIELD()) != 0) {
				cerr << "Return value " << ret << " from calling pthread_yield!" << endl;
				throw 1;
			}
			time_t curtime = time(0);
			if(curtime - tlast_ >= progressInterval) {
				time_t elapsed = curtime - t0_;
				reportStats(elapsed);
				cerr << "-----" << endl;
				tlast_ = curtime;
			}
		}
	}
	
	/**
	 * Start the work of a single thread.
	 */
	static void* startThread(void *obj) {
		reinterpret_cast<ProgressThread*>(obj)->work();
		return NULL;
	}
	
	bool dead_;
	time_t t0_, tlast_;
	THREAD_T thread_;
};

/**
 * One thread that aligns reads in parallel with other threads.
 */
class SearchThread {

public:

	SearchThread() { reset(); }

	SearchThread(
		int tid,
		int nt,
		MerIndex* ind,
		Reads* rs,
		const ReferenceSet* refs,
		AlignOutput* os,
		const ReferenceMap* rmap,
		const AnnotationMap* amap)
	{
		init(tid, nt, ind, rs, refs, os, rmap, amap);
	}
	
	/**
	 * Set all search thread fields to an uninitialized state.
	 */
	void reset() {
		ind_ = NULL;
		rs_ = NULL;
		refs_ = NULL;
		os_ = NULL;
		rmap_ = NULL;
		amap_ = NULL;
	}

	/**
	 * Initialize the SearchThread with the given field values.
	 */
	void init(
		int tid,
		int nt,
		MerIndex* ind,
		Reads* rs,
		const ReferenceSet* refs,
		AlignOutput* os,
		const ReferenceMap* rmap,
		const AnnotationMap* amap)
	{
		tid_ = tid;   // thread id
		nt_ = nt;     // tot # threads
		ind_ = ind;   // Merman index
		rs_ = rs;     // input reads
		refs_ = refs; // set of references
		os_ = os;     // output sink
		rmap_ = rmap; // reference map for translating b/t coordinate systems
		amap_ = amap; // annotation map
	}
	
	/**
	 * Return true iff the SearchThread fields have been initialized.
	 */
	bool inited() {
		if(ind_ != NULL) {
			assert(rs_ != NULL);
			assert(refs_ != NULL);
			return true;
		}
		return false;
	}

	/**
	 * Create a new thread; run the startThread member function in the
	 * new thread.
	 */
	void run() {
		assert(inited());
		THREAD_CREATE(thread_, SearchThread::startThread, this);
	}

	/**
	 * Wait until this thread is finished before returning.
	 */
	void join() { THREAD_JOIN(thread_); }

private:

	/**
	 * Called by startThread when a new search thread is initialized.  Actually do alignment.
	 */
	void work() {
		Read r;
		HitSet hitset;
		AlignResult res;
		RandomSource rnd;
		int64_t ltotHits = 0, lalReads = 0, lunalReads = 0, lmaxReads = 0;
		int64_t ltotSeedHits = 0, lmaxSeedHits = 0;
		int64_t lunsampled = 0;
		int64_t skipped = 0;
		while(rs_->next(r)) {
			assert(r.repOk());
			if(skipped < readskip) {
				skipped++;
				continue;
			}
			rnd.init(randseed ^ r.rand.nextU32());
			// If sampleRate < 1.0f, apply sampling odds to this read;
			// if it isn't chosen, skip it.
			if(sampleRate < 1.0f && rnd.nextDouble() >= sampleRate) {
				// Not chosen
				lunsampled++;
				continue;
			}
			r.color = /*r.hitset.color =*/ color;
			// If the number of unskipped reads exceeds the readmax
			// ceiling, break out of the loop.  Note that there's a
			// minor race condition here.
			if(nreads.value() + 1 > readmax) {
				if(verbose) {
					cout << "Stopping because readmax " << readmax
					     << " was exceeded" << endl;
				}
				break;
			}
			// Trim as requested by the user.  Could do something more
			// sophisticated here.
			r.trim3(trim3);
			r.trim5(trim5);
			//r.initHitset();
			hitset.reset();
			nreads++;
			// The read must be at least as long as the mer length that
			// was used when building the index.
			if(r.seq.length() < (size_t)ap.iSeedLen) {
				if(!quiet) {
					cerr << "Warning: Skipping read " << r.name
					     << " because length " << r.seq.length()
					     << " was less than indexable seed length "
					     << ap.iSeedLen << endl;
				}
				os_->printUnalignedRead(r, r.seq, r.qual, FILTER_TOO_SHORT_FOR_INDEX);
				continue;
			}
			// The read must be at least as long as the mer length that
			// was used when building the index.
			size_t clen = r.seq.length();
			if(ap.minLen != -1 && clen < (size_t)ap.adjMinLen) {
				if(!quiet) {
					cerr << "Warning: Skipping read " << r.name
					     << " because length " << clen
					     << " was such that alignment length would be less "
					     << "than --minlen: " << ap.minLen << endl;
				}
				os_->printUnalignedRead(r, r.seq, r.qual, FILTER_TOO_SHORT_FOR_MINLEN_PARAMS);
				continue;
			}
			// The read is trimmed and has passed all filters.  Next we
			// align it.
			if(verbose) cout << "  aligning read: " << r << endl;
			if(!ind_->empty()) {
				res.clear(); // clear the alignment results structure
				//assert(iformat == INPUT_CHAININ || r.hitset.maxedStratum == -1);
				ind_->query(r, *refs_, rmap_, amap_, hitset, *os_, res, ap, true, rnd, tid_);
				// Update per-thread counters
				if(res.hits > 0) {
					lalReads++;
					ltotHits += res.hits;
				} else if(res.maxed) lmaxReads++;
				else lunalReads++;
				ltotSeedHits += res.seedHits;
				lmaxSeedHits = max<int64_t>(lmaxSeedHits, res.seedHits);
				if(res.bail) {
					// The aligner signaled that we should bail at this
					// point.
					throw 1;
				}
			} else {
				// If the index is empty, there can't possibly be any
				// hits.  TODO: this seems like something that should
				// cause an error early on.
				hitset.reportUpTo(r, *os_, ap.khits, *refs_, rmap_, false, amap_);
			}
			if((nreads.value()+1) % updateEvery == 0) {
				// Fold per-thread counters into global counters
				unalReads    += lunalReads;
				alReads      += lalReads;
				maxReads     += lmaxReads;
				maxSeedHits.max(lmaxSeedHits);
				totSeedHits  += ltotSeedHits;
				totAls       += ltotHits;
				totUnsampled += lunsampled;
				lunalReads   = 0;
				lalReads     = 0;
				lmaxReads    = 0;
				ltotSeedHits = 0;
				ltotHits     = 0;
				lunsampled   = 0;
			}
		}
		// Update global counters in synchronized fashion
		unalReads    += lunalReads;
		alReads      += lalReads;
		maxReads     += lmaxReads;
		maxSeedHits.max(lmaxSeedHits);
		totSeedHits  += ltotSeedHits;
		totAls       += ltotHits;
		totUnsampled += lunsampled;
	}

	/**
	 * Start the work of a single thread.
	 */
	static void* startThread(void *obj) {
		reinterpret_cast<SearchThread*>(obj)->work();
		return NULL;
	}

	int                  tid_;    // thread id
	int                  nt_;     // tot # threads
	MerIndex            *ind_;    // Merman index
	const ReferenceMap  *rmap_;   // reference map for translating b/t coordinate systems
	const ReferenceSet  *refs_;   // set of references
	const AnnotationMap *amap_;   // annotation map
	Reads               *rs_;     // input reads
	AlignOutput         *os_;     // output sink
	THREAD_T             thread_; // thread object for this thread
};

string cmdline;

/**
 * Parse command-line arguments and options and set static and global
 * variables accordingly.
 */
static void parseCommandLine(int argc, char **argv) {
	int option_index = 0;
	int next_option;
	for(int i = 0; i < argc; i++) {
		cmdline += argv[i];
		if(i < argc-1) cmdline.push_back(' ');
	}
	do {
		next_option = getopt_long(argc, argv, short_opts, long_opts, &option_index);
		switch (next_option) {
			case 'h': {
				printUsage(cout);
				throw 0;
				break;
			}
			case ARG_VERBOSE: verbose = true; break;
			case ARG_MERVERBOSE: merverbose = true; break;
			case ARG_QUIET: quiet = true; break;
			case ARG_COLOR_KEEP_ENDS: ap.exDecEnds = false; break;
			case ARG_REF_IS_STR: refIsStr = true; break;
			case 'q': iformat = INPUT_FASTQ; break;
			case 'f': iformat = INPUT_FASTA; break;
			case 'c': iformat = INPUT_CMDLINE; break;
			case 'Q': qualFile = optarg; break;
			case ARG_CHAININ: iformat = INPUT_CHAININ; break;
			case ARG_REFMAP: refmapFile = optarg; break;
			case ARG_ANNOT: annotFile = optarg; break;
			case ARG_REFIDX: refidx = true; break;
			case ARG_FULLREF: fullref = true; break;
			case ARG_NOALIGNFW: {
				// Don't align the forward representation of the read
				ap.alignfw = false;
				break;
			}
			case ARG_NOALIGNRC: {
				// Don't align the revcomp representation of the read
				ap.alignrc = false;
				break;
			}
			case ARG_GENCRICK: {
				// Omit Watson strand and derived strands
				rp.genCrick = true;
				break;
			}
			case ARG_GENREVCOMPS: {
				// Omit Crick strand and derived strands
				rp.genRevcomps = true;
				break;
			}
			case ARG_FBIS: fcontBis = true; break;
			case ARG_FRC: fcontRc = true; break;
			case 'k': ap.khits = atoi(optarg); break;
			case 'p': nthreads = atoi(optarg); break;
			case 'a': ap.khits = numeric_limits<int>::max(); break;
			case 'm': ap.mhits = atoi(optarg); ap.msample = false; break;
			case 'M': ap.mhits = atoi(optarg); ap.msample = true; break;
			case 'i': ap.requireAnnot = rp.requireAnnot = true; break;
			case 'U': setIupacsCat(3); break; // IUPAC codes unmatchable
			case 'R': ap.snpPen = parse<int>(optarg); break;
			case 'E': rp.entThresh = parsePair<int>(optarg, ","); break;
			case 'B': begin = parse<int>(optarg); break;
			case 'S': oformat = OUTPUT_SAM; break;
			case ARG_OLDFORMAT: oformat = OUTPUT_BOWTIE; break;
			case 'F': {
				iformat = INPUT_FASTA_CONT;
				pair<size_t, size_t> p = parsePair<size_t>(optarg, ",");
				fastaContLen = p.first;
				fastaContFreq = p.second;
				break;
			}
			case ARG_SRAND: randseed = parse<unsigned>(optarg); break;
			case ARG_SAMPLE: sampleRate = parse<float>(optarg); break;
			case ARG_PRINTCOLOR: useColor = true; break;
			case ARG_NOMAQROUND: ap.maqRound = false; break;
			case ARG_NAIVE: naiveCheck = true; break;
			case ARG_JUST_NAIVE: justNaive = true; break;
			case ARG_IGNOREQUALS: ap.ignoreQuals = true; break;
			case ARG_VERBOSE_INDEX: verboseIndex = true; break;
			case ARG_BISC2: {
				rp.bisulfiteC = true;     // All Cs become Ys
				rp.bisulfiteCpG = false;
				rp.genCrick = true;       // BS mode implies rcref mode
				rp.genRevcomps = false;   // Need need for BSWR or BSCR
				rp.watsonCrickRc = false; // Watson and Crick not revcomps
				ap.alignfw = true;        // Only align fw representation
				ap.alignrc = false;
				break;
			}
			case ARG_BISCPG2: {
				rp.bisulfiteC = false;
				rp.bisulfiteCpG = true;   // CpG Cs become Ys, other Cs become Ts
				rp.genCrick = true;       // BS mode implies rcref mode
				rp.genRevcomps = false;   // Need need for BSWR or BSCR
				rp.watsonCrickRc = false; // Watson and Crick not revcomps
				ap.alignfw = true;        // Only align fw representation
				ap.alignrc = false;
				break;
			}
			case ARG_BISC4: {
				rp.bisulfiteC = true;     // All Cs become Ys
				rp.bisulfiteCpG = false;
				rp.genCrick = true;       // BS mode implies rcref mode
				rp.genRevcomps = false;   // Need need for BSWR or BSCR
				rp.watsonCrickRc = false; // Watson and Crick not revcomps
				ap.alignfw = true;        // Align both fw and rc representations
				ap.alignrc = true;
				break;
			}
			case ARG_BISCPG4: {
				rp.bisulfiteC = false;
				rp.bisulfiteCpG = true;   // CpG Cs become Ys, other Cs become Ts
				rp.genCrick = true;       // BS mode implies rcref mode
				rp.genRevcomps = false;   // Need need for BSWR or BSCR
				rp.watsonCrickRc = false; // Watson and Crick not revcomps
				ap.alignfw = true;        // Align both fw and rc representations
				ap.alignrc = true;
				break;
			}
			case ARG_JUSTBLOWUP: justBlowup = true; break;
			case ARG_COST: printCost = true; break;
			case ARG_ESTEP: estep = parse<float>(optarg); break;
			case ARG_PROGRESS: {
#ifndef USE_PTHREADS
				cerr << "Error: --progress requires that merman be compiled with -DUSE_PTHREADS" << endl;
				throw 1;
#endif
				progress = true; break;
			}
			case ARG_UPDATEVERY: updateEvery = parse<size_t>(optarg); break;
			case ARG_PROGRESSINTERVAL: progressInterval = parse<time_t>(optarg); break;
			case ARG_MINLEN: ap.minLen = parse<int>(optarg); break;
			case ARG_SPECIFICITY: specificity = parse<int>(optarg); break;
			case ARG_NK: {
				nk = parsePair<int>(optarg, ","); break;
				if(nk.first == 0 || nk.second == 0) {
					// Print an error message
					cerr << "--nk values must be >= 1" << endl;
					throw 1;
				}
				if(nk.first < nk.second) {
					// Print an error message
					cerr << "First --nk value (" << nk.first << ") must be >= second (" << nk.second << ")" << endl;
					throw 1;
				}
			}
			case ARG_PHRED33:  solexaScale = false; sixty4off = false; break;
			case ARG_PHRED64:  solexaScale = false; sixty4off = true;  break;
			case ARG_SOLEXA64: solexaScale = true;  sixty4off = true;  break;
			case ARG_SAM_NO_CS_CQ: samNoCsCq = true; break;
			case ARG_BUFSZ: bufsz = parse<size_t>(optarg); bufsz *= 1024; break;
			case 't': timing = true; break;
			case 'G': {
				// Set a character to be unmatchable
				asc2dnacat[tolower(optarg[0])] = 3;
				asc2dnacat[toupper(optarg[0])] = 3;
				break;
			}
			case ARG_N_IUPAC: {
				// Make 'N' be an IUPAC character rather than an unmatchable
				asc2dnacat[(int)'N'] = asc2dnacat[(int)'n'] = 2;
				break;
			}
			case 'u': readmax = parse<int64_t>(optarg); break;
			case 's': readskip = parse<int64_t>(optarg); break;
			case 'n': {
				ap.seedMms = atoi(optarg);
				ap.e2eMms = numeric_limits<int>::max();
				break;
			}
			case 'l': ap.seedLen = parse<int>(optarg); break;
			case 'L': ap.iSeedLen = parse<int>(optarg); break;
			case 'w': seedWidth = parse<int>(optarg); break;
			case 'e': qualCeil = parse<int>(optarg); break;
			case '3': trim3 = parse<size_t>(optarg); break;
			case '5': trim5 = parse<size_t>(optarg); break;
			case 'y': break;
			case 'C': color = true; break;
			case 'v': {
				ap.e2eMms = atoi(optarg);
				ap.seedMms = ap.e2eMms;
				break;
			}
			case ARG_STRATA: ap.strata = true; break;
			case -1: break;
			default: {
				cerr << "Unknown option: " << (char)next_option << endl;
				printUsage(cerr);
				throw 1;
			}
		}
	} while(next_option != -1);
	if(ap.e2eMms < numeric_limits<int>::max()) {
		qualCeil = numeric_limits<int>::max() - 10;
		estep = 0.0f;
	}
	ap.adjMinLen = ap.minLen;
	if(color && ap.adjMinLen != -1) {
		// Set ap.adjMinLen accordingly.  The user expects --minlen to
		// govern the minimum length of the final output, so here we
		// calculate the difference between the final output length and
		// the length in colors.
		if(ap.exDecEnds) {
			ap.adjMinLen++;
		} else {
			ap.adjMinLen--;
		}
	}
	// Set up penceil
	if(ap.minLen == -1) {
		if(estep != 0.0f) {
			if(!quiet)
				cerr << "Warning: ignoring --estep because --minLen is not specified." << endl;
		}
		for(int i = 0; i < 1024; i++) ap.penceil[i] = qualCeil;
	} else {
		for(int i = 0; i < 1024; i++) {
			if(i <= ap.minLen) ap.penceil[i] = qualCeil;
			else {
				int diff = i - ap.minLen;
				// Add all the esteps and round to nearest integer
				ap.penceil[i] = (int)(qualCeil + (diff * estep) + 0.5);
			}
		}
		assert_gt(ap.penceil[ap.minLen], -1);
	}
	if(seedWidth < 1) {
		cerr << "Error: -w/--seedwidth must be > 0, was: " << seedWidth << endl;
		throw 1;
	}
	if(iformat == INPUT_FASTQ && color) {
		iformat = INPUT_CSFASTQ;
	}
	if(iformat == INPUT_FASTA && color) {
		iformat = INPUT_CSFASTA;
	}
	if(!qualFile.empty() && iformat != INPUT_CSFASTA && !quiet) {
		cerr << "Warning: Ignoring -Q " << qualFile
			 << " because input is not CSFASTA (-f -C)" << endl;
	}
	if(!qualFile.empty() && iformat == INPUT_CSFASTA) {
		iformat = INPUT_CSFASTA_AND_QV;
	}
	if(optind+1 >= argc) {
		cerr << "Must supply reference and query as arguments" << endl;
		printUsage(cerr);
		throw 1;
	}
	bool e2e = ap.e2eMms < numeric_limits<int>::max();
	if(!e2e && ap.iSeedLen > ap.seedLen) {
		if(!quiet) {
			cerr << "Warning, indexable seed length (-L) was automatically "
				 << "decreased to " << ap.seedLen << " to match policy "
				 << "seed length (-l)" << endl;
		}
		ap.iSeedLen = ap.seedLen;
	}
	if(seedWidth > ap.iSeedLen) {
		if(!quiet) {
			cerr << "Warning, seed width (-w) was automatically decreased to "
				 << ap.iSeedLen << " to match indexable seed length (-L)" << endl;
		}
		seedWidth = ap.iSeedLen;
	}
	if(ap.minLen != -1 && ap.adjMinLen < ap.seedLen) {
		if(!quiet) {
			cerr << "Warning, minimum alignment length (--minlen) was automatically increased to" << endl
				 << ap.seedLen << " to match policy seed length (-l)." << endl;
			if(color) {
				cerr << "Note that --minlen refers to length of the alignment in nucleotides, which is" << endl
					 << "always 1 greater than the length in colors." << endl;
			}
		}
		int adj = ap.adjMinLen - ap.minLen;
		ap.minLen = ap.seedLen + adj;
		ap.adjMinLen = ap.seedLen;
	}
	if(ap.e2eMms != numeric_limits<int>::max() && ap.e2eMms > 10) {
		cerr << "Invalid -v: " << ap.e2eMms << endl;
		throw 1;
	}
	if(ap.seedMms > 10) {
		cerr << "Invalid -n: " << ap.seedMms << endl;
		throw 1;
	}
	/*
	// If the user did not specify --nk, we choose a default.								
	if(nk.first == 0 || nk.second == 0) {
		if(ap.e2eMms == numeric_limits<int>::max()) {
			// Seeded search
			if(ap.seedMms == 0) {
				nk.first = nk.second = 1;
			} else if(ap.seedMms == 1) {
				
			} else if(ap.seedMms == 2) {
			}
		} else {
			// End-to-end search
		}
		cerr << "N-choose-K scheme defaulted to " << nk.first << ", " << nk.second << endl;
	}
	*/
}

// C++ name mangling is disabled for the merman() function to make it
// easier to use Merman as a library.
extern "C" {

int merman(int argc, char **argv);

/**
 * Merman main driver function.  Does the following:
 *
 * 1. Parses command-line options
 */
int merman(int argc, char **argv) {
	reset();
	try {
		parseCommandLine(argc, argv);
		Timer tov(cerr, "Overall time: ", timing);
		EList<string> refstrs;
		ReferenceSet refs;
		EList<string> refnames;
		EList<size_t> reflens;
		string refstr = argv[optind++];
		tokenize(refstr, ",", refstrs);
		auto_ptr<MerIndex> ind(
			new MerIndex(ap, rp, readLen, seedWidth, nk.first, nk.second,
			             specificity, begin, naiveCheck, nthreads));
		{
			Timer t(cerr, "... ", timing);
			if(timing) cerr << "Reading reference sequences..." << endl;
			for(size_t i = 0; i < refstrs.size(); i++) {
				if(timing) {
					cerr << "  Sequence " << (i+1) << " of " << refstrs.size() << endl;
				}
				if(refIsStr) {
					refs.addOrigReferenceString(refstrs[i].c_str(), rp);
				} else {
					refs.addOrigReferenceFasta(refstrs[i].c_str(), rp);
				}
			}
			for(size_t i = 0; i < refs.numRefs(); i++) {
				refnames.push_back(string(refs[i].name.toZBuf()));
				reflens.push_back(refs[i].seq.length(color));
			}
			if(refs.numRefs() == 0) {
				cerr << "Warning: No references were found" << endl;
			}
			if(rp.genCrick) {
				if(timing) {
					cerr << "  Crickizing" << endl;
				}
				// Add the crick strand.  If there were bisulfite
				// transformations to the Watson strand, they are
				// removed from the Watson strand before the Crick copy
				// is made.  Transformations are then applied to the
				// new Crick strand.  This has the effect of correctly
				// producing either Watson / Crick in the non-bisulfite
				// case, or BS Watson / BS Crick in the bisulfite case.
				refs.addReferenceRevComps(rp, false, 1, 0);
			}
			if(rp.genRevcomps) {
				if(timing) {
					cerr << "  Adding reverse comps" << endl;
				}
				// Add reverse complements of all existing references
				// (after the transformations have already been
				// applied).
				refs.addReferenceRevComps(rp, true, -1, 1);
			}
			assert(refs.repOk());
		}

		pair<size_t, size_t> mers = make_pair(0, 0);
		EList<MerIndexThread> threads;
		{
			Timer t(cerr, "... ", timing);
			if(timing) cerr << "Preparing to extract sub-sequences..." << endl;
			// Instantiate and run index threads
			assert_gt(nthreads, 0);
			threads.resize(nthreads);
			for(int i = 0; i < nthreads; i++) {
				threads[i].runCount(&refs, ind.get(), i, nthreads, color);
			}
			for(int i = 0; i < nthreads; i++) {
				pair<size_t, size_t> mrs = threads[i].join();
				mers.first += mrs.first;
				mers.second += mrs.second;
			}
			ind->allocateMers();
		}
		if(timing || verbose || justBlowup) {
			cerr << "Expecting index footprint of ";
			printBytes(mers.first * sizeof(mer_ent), cerr);
			cerr << endl;
			if(mers.first > mers.second) {
				cerr.setf(ios::fixed);
				cerr << "  base footprint is ";
				printBytes(mers.second * sizeof(mer_ent), cerr);
				cerr << endl
				     << "  blowup factor: " << setprecision(2) << ((double)mers.first / (double)mers.second) << endl;
			}
			if(justBlowup) throw 0;
		}
		{
			Timer t(cerr, "... ", timing);
			if(timing) cerr << "Extracting index sub-sequences..." << endl;
			// Instantiate and run index threads
			for(int i = 0; i < nthreads; i++) {
				threads[i].runIndex(&refs, ind.get(), i, nthreads, color);
			}
			for(int i = 0; i < nthreads; i++) threads[i].join();
		}
		assert_eq(mers.first, ind->size());
		if(verbose) {
			cout << "  read " << refs.numRefs() << " reference strings" << endl;
		}
		if(refs.empty() && iformat != INPUT_CHAININ) {
			cerr << "Index is empty; not enough reference sequence supplied" << endl;
			throw 1;
		}
		if(refs.numRefs() == 0 && iformat != INPUT_CHAININ) {
			cerr << "No reference strings provided; aborting..." << endl;
			throw 1;
		}
		{
			Timer t(cerr, "Sorting reference mers: ", timing);
			ind->sort(nthreads); // sort mers
		}
		{
			Timer t(cerr, "... ", timing);
			if(timing) cerr << "Aligning reads..." << endl;
			string rstr = argv[optind++];
			// Instantiate reference map, which translates to new reference
			// coordinate system prior to alignment output
			auto_ptr<ReferenceMap> rmap(
				refmapFile == NULL ? NULL : new ReferenceMap(refmapFile, !refidx));
			// Instantiate annotation map, which encodes SNP locations & alleles
			auto_ptr<AnnotationMap> amap(
				annotFile == NULL ? NULL : new AnnotationMap(annotFile));
			// Instantiate the read-input object
			auto_ptr<Reads> rs(
				(iformat == INPUT_CMDLINE) ?
					(Reads*)new StringReads(rstr, begin) :
					((iformat == INPUT_FASTA) ?
						(Reads*)new FastaReads(rstr, begin, bufsz) :
						((iformat == INPUT_FASTA_CONT) ?
							(Reads*)new FastaContinuousReads(
								rstr, begin, fastaContLen,
								fastaContFreq, fcontBis, fcontRc,
								color) :
							((iformat == INPUT_FASTQ) ?
								(Reads*)new FastqReads(rstr, solexaScale, sixty4off, begin, bufsz) :
									((iformat == INPUT_CHAININ) ?
										(Reads*)new ChainReads(rstr, begin, bufsz) :
											((iformat == INPUT_CSFASTA) ?
												(Reads*)new CSFastaReads(rstr, begin, bufsz) :
													((iformat == INPUT_CSFASTA_AND_QV) ?
														(Reads*)new CSFastaAndQVReads(rstr, qualFile, begin, bufsz) :
														(Reads*)new CSFastqReads(rstr, solexaScale, sixty4off, begin, bufsz))))))));
			// Set output stream
			string of = "-";
			if(optind < argc) of = argv[optind++];
			// Instantiate the alignment-output object
			auto_ptr<AlignOutput> outs(
				(oformat == OUTPUT_SAM) ?
					(AlignOutput*)new SamOutput(of, fullref, refidx, rp.bisulfiteC || rp.bisulfiteCpG, !samNoCsCq) :
					(AlignOutput*)new BowtieOutput(of, fullref, printCost, refidx, rp.bisulfiteC || rp.bisulfiteCpG));
			outs->printHeader(refnames, reflens);
			// Run the progress thread, if requested
			ProgressThread proThread;
			if(progress) proThread.run();
			// Instantiate and run search threads
			EList<SearchThread> sthreads;
			sthreads.resize(nthreads);
			for(int i = 0; i < (int)sthreads.size(); i++) {
				sthreads[i].init(
					i, (int)sthreads.size(), ind.get(), rs.get(), &refs,
					outs.get(), rmap.get(), amap.get());
				sthreads[i].run();
			}
			// Wait until search sthreads are finished
			for(size_t i = 0; i < sthreads.size(); i++) {
				sthreads[i].join();
			}
			if(progress) {
				proThread.kill();
				proThread.join();
			}
			outs->flush();
		}
		if(!quiet) ProgressThread::reportStats();
	} catch(exception& e) {
		cerr << "Command: ";
		for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
		cerr << endl;
		return 1;
	} catch(int e) {
		if(e != 0) {
			cerr << "Command: ";
			for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
			cerr << endl;
		}
		return e;
	}
	return 0;
}

} // extern "C"

