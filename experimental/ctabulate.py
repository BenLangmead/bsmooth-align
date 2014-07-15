#!/usr/bin/env python

"""
ctabulate.py

Author: Ben Langmead
Date: 4/23/2012
Contact: langmea@cs.jhu.edu
"""

import os
import argparse
import pysam
import sys
import string
import unittest
import copy
import site
import errno
import reference
import logging
from collections import defaultdict

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
pytrim_path = os.path.join(base_path, "modules", "btl")
site.addsitedir(pytrim_path)

version = "1.0.0"

parser = argparse.ArgumentParser(
    description='''
Tabulate evidence at CpGs and/or non-CpG Cs in a manner similar to
BSmooth's bsev_tabulate.pl script.  Operates on SORTED BSmooth .bam
files.  BSmooth does not sort them by default so you will generally
have to run 'samtools sort' on the .watson.bam and .crick.bam files
before using them with this tool.

Also requires FASTA of the reference be available.  For large
genomes, parsing the FASTA can be slow.  This can be accelerated by
providing the --fasta-pickle option, which will cache a pickled copy
of the FASTA file for faster loading.

The output is tab-separated and consists of these columns.

1. Chromosome name
2. Offset (1-based)
3. Strand (W=Watson, C=Crick)
4. Methylated qualities; string with once character per base quality
   value associated with a methylated piece of evidence at this
   position.  Empty if there is no methylated evidence.
5. Number of distinct sequencing cycles among all the methylated
   pieces of evidence.
6. Same as 4 but for unmethylated evidence.
7. Same as 5 but for unmethylated evidence.
8. Number of pieces of evidence filtered out because the sequencing
   cycle is untrusted.  (TODO: add option to set this threshold)
9. Number of pieces of evidence filtered out because the source read
   is too short.  See --min-read-len.
10. Number of pieces of evidence filtered out because the allele was
    not a C or T.
11. Number of pieces of evidence filtered out because the source
    alignment's mapping quality (MAPQ) is too low.  See --min-mapq.
12. Number of pieces of evidence filtered out because the base
    quality is too low.  (TODO: add option to set this threshold)

By default a header line is output, but this can be suppressed with
--no-header.''', formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('--fasta', metavar='path', dest='fasta', type=str, nargs='+', required=True,
                    help='FASTA file(s) containing reference genome sequences')
parser.add_argument('--bsbam', metavar='path', dest='bsbam', type=str, nargs='+', required=True,
                    help='Sorted BAM file(s) containing Bsmooth alignment output for bisulfite sequencing reads')
parser.add_argument('--locus', metavar='chr:offset1-offset2', dest='loci', type=str, nargs='+', required=False,
                    help='Genome interval to tabulate.  If specified more than once, all specified intervals are '
                         'tabulated.')
parser.add_argument('--fasta-index', dest='fa_idx', type=str, required=False,
                    help='Restore index of FASTA files from this file, or, if file doesn\'t yet exist, put the index '
                         'there so it can be used in future runs')
parser.add_argument('--fasta-pickle', type=str, required=False,
                    help='If not already present, store pickled copy of reference in given file. '
                         'If pickle file is present, load from file, which is very quick.')
parser.add_argument('--min-mapq', action='store', type=int, default=20,
                    help='Read-level measurements with mapping quality (MAPQ) less than this threshold are filtered '
                         'out')
parser.add_argument('--min-read-len', dest='min_rdl', action='store', type=int, default=40,
                    help='Read-level measurements from reads with length less than this threshold are filtered out')
parser.add_argument('--all-C', action='store_const', const=True, default=False,
                    help='Tabulate all Cs.  Default is to tabulate just CpG Cs')
parser.add_argument('--context', dest='npad', action='store', type=int, default=0,
                    help='When printing context for a site, print this many positions on either side')
parser.add_argument('--merge', action='store_const', const=True, default=False,
                    help='Merge read-level measurements overlapping the C and G in a CpG before printing records.  '
                         'Records will have M in the type column, indicating the records are merged.')
parser.add_argument('--sanity', action='store_const', const=True, default=False,
                    help='Do various sanity checks')
parser.add_argument('--test', action='store_const', const=True, default=False, help='Do unit tests')
parser.add_argument('--profile', action='store_const', const=True, default=False,
                    help='Print profiling info')
parser.add_argument('--verbose', action='store_const', const=True, default=False, help='Be talkative')
parser.add_argument('--version', action='store_const', const=True, default=False, help='Print version and quit')

# FIXES FOR SAM ISSUES

# Number to subtract from Crick offsets - used on at least one occasion to fix
# bad output from BSmooth
parser.add_argument('--crick-add', dest='cri_add', action='store', type=int, default=0,
                    help='Amount to subtract from all Crick-strand offsets')

# Whether to reverse the CIGAR for Crick-strand alignments
parser.add_argument('--crick-cigar-rev', dest='cri_cigar_rev', action='store_const', const=True, default=False,
                    help='Reverse CIGAR string for Crick-strand alignments')

# Whether to reverse-comp the YO:Z field for reads with various combinations of
# Watson/Crick forward/reverse-comp
parser.add_argument('--rc-wat-fw', action='store_const', const=True, default=False,
                    help='Reverse complement YO:Z string for W-strand alignments')
parser.add_argument('--rc-cri-fw', action='store_const', const=True, default=False,
                    help='Reverse complement YO:Z string for C-strand alignments')
parser.add_argument('--rc-wat-rc', action='store_const', const=True, default=False,
                    help='Reverse complement YO:Z string for WR-strand alignments')
parser.add_argument('--rc-cri-rc', action='store_const', const=True, default=False,
                    help='Reverse complement YO:Z string for CR-strand alignments')

_revcomp_trans = string.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(_x):
    return _x[::-1].translate(_revcomp_trans)


def phred33c2i(c):
    return ord(c) - 33


def phred33i2c(i):
    return chr(i + 33)


def parse(pos, wat, rc, cigars, seq, qual, i, mul, qual_conv=phred33c2i):
    
    """ Parse the CIGAR string for the alignment in rec and install the read
        alleles in the multiple alignment "mul".
        
        Some question as to how to represent insertions in the multiple
        alignment.  How do we want to query them?  We might want to know if
        they're nearby.
        
        pos: 0-based reference offset
        wat: True iff read aligned to BS-Watson strand
        rc: True iff read aligned to reverse-complement of BSW or BSC
        cigars: CIGAR operators and run lengths
        seq: read sequence
        qual: quality sequence
        i: row in multiple alignment structure
        mul: multiple alignment structure """
    
    assert cigars is not None

    def ensure_mul(refoff, j):
        if refoff not in mul:
            mul[refoff] = []
        if len(mul[refoff]) <= j:
            diff = j - len(mul[refoff]) + 1
            mul[refoff] += [None] * diff
    
    assert pos >= 0
    last_op = -1
    rdlen = len(seq)
    refoff, rdoff = pos, 0
    refoff_min, refoff_max = sys.maxint, 0
    
    left5 = wat == rc
    
    # Operations (in order of op id): MIDNSHP
    for c in cigars:
        op, run = c
        assert op != last_op
        last_op = op
        assert (0 <= op <= 2) or op == 4, "Don't know how to interpret CIGAR op %d" % op
        assert run > 0
        if op == 0:   # M
            for j in xrange(0, run):
                assert rdoff < rdlen
                cy = rdoff if left5 else rdlen - rdoff - 1
                ensure_mul(refoff, i)
                mul[refoff][i] = (seq[rdoff], qual_conv(qual[rdoff]), cy, wat, rc)
                refoff_min, refoff_max = min(refoff_min, refoff), max(refoff_max, refoff)
                rdoff += 1
                refoff += 1
        elif op == 1:  # I
            ensure_mul(-refoff, i)
            for j in xrange(0, run):
                assert rdoff < rdlen
                cy = rdoff if left5 else rdlen - rdoff - 1
                mul[-refoff][i] = (seq[rdoff], qual_conv(qual[rdoff]), cy, wat, rc)
                rdoff += 1
        elif op == 2:  # D
            cy = rdoff if left5 else rdlen - rdoff - 1
            for j in xrange(0, run):
                ensure_mul(refoff, i)
                mul[refoff][i] = ('-', -1, cy, wat, rc)
                refoff_min, refoff_max = min(refoff_min, refoff), max(refoff_max, refoff)
                refoff += 1
        elif op == 4:  # S
            # Don't change refoff, since the POS field of the SAM record is
            # already shifted according to soft clipping
            rdoff += run
    
    return refoff_min, refoff_max


class TestParse(unittest.TestCase):
    
    def test_parse1(self):
        mul = dict()
        pos = 7
        wat = True
        fw = True
        cigars = [(0, 10)]  # 10M
        seq = "GCATGCAACT"
        qual = "IIIIIIIIII"
        #       0123456789
        #       7890123456
        parse(pos, wat, fw, cigars, seq, qual, 0, mul)
        
        self.assertTrue(0 not in mul)
        self.assertTrue(6 not in mul)
        self.assertTrue(7 in mul)
        
        self.assertEqual(1,   len(mul[7]))
        self.assertEqual('G', mul[7][0][0])
        self.assertEqual(phred33c2i('I'), mul[7][0][1])
        self.assertEqual(0,   mul[7][0][2])
        
        self.assertEqual(1,   len(mul[16]))
        self.assertEqual('T', mul[16][0][0])
        self.assertEqual(phred33c2i('I'), mul[16][0][1])
        self.assertEqual(9,   mul[16][0][2])
        
        # Try backwards version of above
        
        fw = False
        parse(pos, wat, fw, cigars, seq, qual, 1, mul)
        
        self.assertEqual(2,   len(mul[7]))
        self.assertEqual('G', mul[7][1][0])
        self.assertEqual(phred33c2i('I'), mul[7][1][1])
        self.assertEqual(9,   mul[7][1][2])
        
        self.assertEqual(2,   len(mul[16]))
        self.assertEqual('T', mul[16][1][0])
        self.assertEqual(phred33c2i('I'), mul[16][1][1])
        self.assertEqual(0,   mul[16][1][2])
        
        # Try sticking a little soft clipping on either end
        cigars = [(4, 2), (0, 7), (4, 1)]  # 10M
        if len(cigars) > 0 and cigars[0][0] == 4:
            pos += cigars[0][1]
        parse(pos, wat, fw, cigars, seq, qual, 2, mul)
        
        self.assertEqual(2,   len(mul[7]))
        self.assertEqual(2,   len(mul[8]))
        self.assertEqual(3,   len(mul[9]))
        self.assertEqual('A', mul[9][2][0])
        self.assertEqual(phred33c2i('I'), mul[9][2][1])
        self.assertEqual(7,   mul[9][2][2])
    
    def test_parse2(self):
        mul = dict()
        pos = 16
        wat = True
        fw = True
        cigars = [(0, 4), (1, 1), (0, 1), (2, 1), (0, 4)]  # 4M1I1M1D4M
        seq = "GCATGCAACN"
        qual = "ABCDEFGHIJ"
        #
        #       GCATGC-AACN
        #       ABCDEF-GHIJ
        #
        #       012345-6789
        #       789012-3456
        #       MMMMIMDMMMM
        parse(pos, wat, fw, cigars, seq, qual, 0, mul)
        
        self.assertTrue(15 not in mul)
        self.assertTrue(16 in mul)
        self.assertTrue(17 in mul)
        self.assertTrue(18 in mul)
        self.assertTrue(19 in mul)
        self.assertTrue(20 in mul)
        self.assertTrue(-20 in mul)  # insertion just before position 20 in ref
        self.assertTrue(-19 not in mul)
        self.assertTrue(-21 not in mul)
        
        self.assertEqual(1,   len(mul[16]))
        self.assertEqual('G', mul[16][0][0])
        self.assertEqual(32,  mul[16][0][1])
        self.assertEqual(0,   mul[16][0][2])
        
        self.assertEqual(1,   len(mul[19]))
        self.assertEqual('T', mul[19][0][0])
        self.assertEqual(35,  mul[19][0][1])
        self.assertEqual(3,   mul[19][0][2])
        
        # Insertion
        self.assertEqual(1,   len(mul[-20]))
        self.assertEqual('G', mul[-20][0][0])
        self.assertEqual(36,  mul[-20][0][1])
        self.assertEqual(4,   mul[-20][0][2])
        
        self.assertEqual(1,   len(mul[20]))
        self.assertEqual('C', mul[20][0][0])
        self.assertEqual(37,  mul[20][0][1])
        self.assertEqual(5,   mul[20][0][2])
        
        # Deletion
        self.assertEqual(1,   len(mul[21]))
        self.assertEqual('-', mul[21][0][0])
        self.assertEqual(-1,  mul[21][0][1])
        self.assertEqual(6,   mul[21][0][2])
        
        self.assertEqual(1,   len(mul[22]))
        self.assertEqual('A', mul[22][0][0])
        self.assertEqual(38,  mul[22][0][1])
        self.assertEqual(6,   mul[22][0][2])


class BsSummary:
    
    """ Class for compiling read-level methylation measurements overlapping a
        position in the reference genome. """
    
    def __init__(self,
                 max_qual,
                 filt_cy=lambda _: False,
                 filt_mapq=lambda _: False,
                 filt_baseq=lambda _: False,
                 filt_rdl=lambda _: False,
                 context=None):
        self.max_qual = max_qual
        self.filt_cy = filt_cy
        self.filt_mapq = filt_mapq
        self.filt_baseq = filt_baseq
        self.filt_rdl = filt_rdl

        def _filt_nuc(_x, wat):
            # Return true iff the allele doesn't match the strand
            xu = _x.upper()
            if wat:
                return xu != 'C' and xu != 'T'
            else:
                return xu != 'G' and xu != 'A'

        self.filt_nuc = _filt_nuc
        self.nfilt_cy, self.nfilt_rdl, self.nfilt_nuc = 0, 0, 0
        self.nfilt_mapq, self.nfilt_baseq = 0, 0
        self.ms = []
        self.context = context
    
    def add(self, seq, qual, wat, rc, cy, rdlen, mapq):
        """ Add a read-level measurement """
        filt = False
        if self.filt_cy is not None and self.filt_cy(cy):
            self.nfilt_cy += 1
            filt = True
        elif self.filt_rdl is not None and self.filt_rdl(rdlen):
            self.nfilt_rdl += 1
            filt = True
        elif self.filt_nuc is not None and self.filt_nuc(seq, wat):
            self.nfilt_nuc += 1
            filt = True
        elif self.filt_mapq is not None and self.filt_mapq(mapq):
            self.nfilt_mapq += 1
            filt = True
        elif self.filt_baseq is not None and self.filt_baseq(qual):
            self.nfilt_baseq += 1
            filt = True
        qu = qual if qual <= self.max_qual else self.max_qual+1
        sq = seq
        if not wat:
            if sq == 'A':
                sq = 'T'
            if sq == 'G':
                sq = 'C'
        if not filt:
            self.ms.append((sq, qu, cy, wat, rc))
    
    def merge(self, o):
        """ Merge the summary already in o into this summary """
        self.nfilt_cy += o.nfilt_cy
        self.nfilt_rdl += o.nfilt_rdl
        self.nfilt_nuc += o.nfilt_nuc
        self.nfilt_mapq += o.nfilt_mapq
        self.nfilt_baseq += o.nfilt_baseq
        self.ms += o.ms

    @staticmethod
    def summary_header():
        return 'mqual\tmcyc\tuqual\tucyc\tfilt_cy\tfilt_readlen\tfilt_non_ct\tfilt_mapq\tfilt_basequal'

    def summary(self, npad, strat_by=4):
        """ Summarize all the read-level measurements collected here """
        r = defaultdict(int)
        q, c = defaultdict(lambda: defaultdict(int)), defaultdict(lambda: defaultdict(int))
        for m in self.ms:
            seq, qual, cy, wat, rc = m
            if wat != rc:
                cy = -cy
            r[seq] += 1
            q[seq][qual] += 1
            c[seq][cy] += 1
        ucyc, mcyc = len(c['T']), len(c['C'])
        uq, mq = q.get('T', dict()), q.get('C', dict())
        context = self.context or '(unknown)'
        filts = [self.nfilt_cy, self.nfilt_rdl, self.nfilt_nuc, self.nfilt_mapq, self.nfilt_baseq]
        if strat_by is not None:
            us, ms = [], []
            for i in xrange(0, self.max_qual+2, strat_by):
                u, m = 0, 0
                for j in xrange(0, strat_by):
                    u += uq.get(i+j, 0)
                    m += mq.get(i+j, 0)
                us.append(u)
                ms.append(m)
            ret = ms + [mcyc] + us + [ucyc] + filts
            ret = map(str, ret)
            if npad > 0:
                ret = [context] + ret
        else:
            ret = [''.join(map(phred33i2c, mq))] + [str(mcyc)] + \
                  [''.join(map(phred33i2c, uq))] + [str(ucyc)] + map(str, filts)
            if npad > 0:
                ret = [context] + ret
        return ret


def print_summaries(oh, ref, summ_fw, summ_rc, npad, merge=True, head=False, strat_by=None):
    
    """ Given a dictionary mapping genomic position to summary of the
        read-level measurements at that position, and given the name of the
        reference sequence, print a summary with one locus per row. """
    
    ignore = dict()
    for off in sorted(summ_fw.keys() + summ_rc.keys()):
        if head:
            oh.write("ref\toff\tstrand\t%s\n" % BsSummary.summary_header())
            head = False
        if off in summ_fw:
            sm = copy.deepcopy(summ_fw[off])
            merged = False
            if off+1 in summ_rc and merge:
                sm.merge(summ_rc[off+1])
                merged = True
                ignore[off+1] = True
            sm = sm.summary(npad, strat_by=strat_by)
            oh.write("%s\t%d\t%s\t%s\n" % (ref, off+1, 'M' if merged else 'W', "\t".join(sm)))
        if off in summ_rc and off not in ignore:
            sm = summ_rc[off].summary(npad, strat_by=strat_by)
            oh.write("%s\t%d\tC\t%s\n" % (ref, off+1, "\t".join(sm)))


def tab_ival(off,
             ln,
             aln_get,
             fa_get,
             npad,
             just_cpg=True,
             filt_cy=lambda _: False,
             filt_mapq=lambda _: False,
             filt_baseq=lambda _: False,
             filt_rdl=lambda _: False,
             sanity=False,
             verbose=False):

    # TODO: is the purpose of the padding just so we don't miss a CpG
    # at the edge of a substring?  Does the padding ever need to be >1?

    """ For each CpG (and possibly each C) locus, tabulate the read-level
        measurements overlapping the locus.
        
        ref:      reference string ID
        off:      reference offset, 0-based
        ln:       length of reference interval to tabulate
        bsbams:   pysam Samfile objects for each BSmooth .bam file
        fa_idx:   indexed FASTA for reference strings
        just_cpg: True -> just tabulate CpGs, False -> all Cs
        
        aln_get retrieves records, where each record is like:
        - 
        """

    # need a clear expectation about what this function should do if
    # our requested bounds fall off the end of the reference sequence.
    logging.info('  Retrieving reference sequence at %d:%d +- %d' % (off, off+ln, npad))
    assert ln > 0
    refstr, lpad, rpad = fa_get(off, ln, npad, npad)
    assert lpad <= npad
    assert rpad <= npad
    assert len(refstr) == ln + lpad + rpad, '%d != %d + %d + %d' % (len(refstr), ln, lpad, rpad)
    if lpad < npad:
        refstr = ('x' * (npad - lpad)) + refstr
    if rpad < npad:
        refstr += ('x' * (npad - rpad))
    assert len(refstr) == ln + 2 * npad

    # Now calculate bisulfite strings
    logging.info('  Bisulfite converting reference sequence')
    refstr_wat = refstr.replace('CG', 'YG')
    refstr_wat = refstr_wat.replace('C', 'T' if just_cpg else 'Y')
    refstr_cri = refstr.replace('CG', 'CR')
    refstr_cri = refstr_cri.replace('G', 'A' if just_cpg else 'R')
    refstr_trimmed, refstr_wat, refstr_cri = \
        refstr[npad:len(refstr)-npad], \
        refstr_wat[npad:len(refstr_wat)-npad], \
        refstr_cri[npad:len(refstr_cri)-npad]
    
    # Multi-alignment is a dictionary mapping reference offsets to lists.
    # Each list is a list of strings, where each string is either empty,
    # indicating that the read does not overlap that position, a character, or
    # a '-' indicating that the position is spanned by the read but no read
    # character lines up with the reference character.
    mul = dict()
    # For each BSmooth .bam file
    nrecs = 0
    mapqs, rdlens, wats = [], [], []
    n, n_ival = 0, 1
    logging.info('  Parsing alignments')
    for rec in aln_get(off, off+ln):
        n += 1
        assert n <= n_ival
        if n == n_ival:
            logging.debug('    Parsed %d alignments')
            n_ival = max(int(n_ival * 1.1 + 0.5), n_ival+1)
        pos, flag, wat, rc, cigar, alsc, mapq, seq, qual = rec
        if (flag & 4) != 0:
            # This can happen if this is an unmapped mate and the opposite
            # mate is mapped
            continue
        assert cigar is not None
        lo, hi = parse(pos, wat, rc, cigar, seq, qual, nrecs, mul)
        assert lo >= pos
        assert hi >= lo, "hi was %d, lo was %d" % (hi, lo)
        assert hi-lo <= ln, "hi-lo was %d, ln was %d" % (hi-lo, ln)
        if sanity and verbose:
            npos, nma = 0, 0
            rdstr, rfstr, rfstr_c = "", "", ""
            strn = ("W" if wat else "C") + ("R" if rc else "")
            for ps in xrange(lo, hi+1):
                if len(mul[ps]) > nrecs and mul[ps][nrecs] is not None:
                    seq, qu, cy, wat, rc = mul[ps][nrecs]
                    assert len(seq) == 1
                    if seq != '-':
                        refc = fa_get(ps, ps+1)
                        rfstr_c += refc
                        if wat:
                            if refc == 'C':
                                refc = 'T'
                            if seq == 'C':
                                seq = 'T'
                        else:
                            if refc == 'G':
                                refc = 'A'
                            if seq == 'G':
                                seq = 'A'
                        npos += 1
                        if refc == seq:
                            nma += 1
                        rdstr += seq
                        rfstr += refc
            frac = 1.0 * nma / npos
            logging.info("Strand=%s, fraction matchup=%0.3f, %s, %d\nrdstr:   %s\nrfstr:   %s\nrfstr_c: %s" % \
                (strn, frac, str(cigar), alsc, rdstr, rfstr, rfstr_c))
        nrecs += 1
        mapqs.append(mapq)
        rdlens.append(len(seq))
        wats.append(wat)

    logging.info('  Parsed all %d alignments' % n)
    max_qual = 40
    # One big overall summary:
    summ = BsSummary(max_qual, filt_mapq=filt_mapq, filt_cy=filt_cy,
                     filt_baseq=filt_baseq, filt_rdl=filt_rdl, context=None)
    summs_wat, summs_cri = dict(), dict()
    assert len(refstr_wat) >= ln, (len(refstr_wat), ln)
    assert len(refstr_cri) >= ln

    logging.info('  Iterating over %d columns' % ln)
    for j in xrange(0, ln):
        logging.debug('    Column %d' % j)
        pos = off + j  # offset into refstr_trimmed / refstr_wat / refstr_cri
        chr_wat = refstr_wat[j].upper()
        chr_cri = refstr_cri[j].upper()
        ev_wat = chr_wat == 'Y'
        ev_cri = chr_cri == 'R'
        context = refstr[j:j+2*npad+1]
        if ev_wat or ev_cri:
            # Go down the column in the multi-alignment and tabulate all the
            # evidence for this position
            if pos in mul:
                logging.debug('      has evidence')
                if ev_wat:
                    summs_wat[pos] = BsSummary(max_qual,
                                               filt_mapq=filt_mapq,
                                               filt_cy=filt_cy,
                                               filt_baseq=filt_baseq,
                                               filt_rdl=filt_rdl,
                                               context=context)
                if ev_cri:
                    summs_cri[pos] = BsSummary(max_qual,
                                               filt_mapq=filt_mapq,
                                               filt_cy=filt_cy,
                                               filt_baseq=filt_baseq,
                                               filt_rdl=filt_rdl,
                                               context=context)
                for i in xrange(0, nrecs):
                    if (wats[i] and ev_wat) or (not wats[i] and ev_cri):
                        if len(mul[pos]) > i and mul[pos][i] is not None:
                            seq, qu, cy, wat, rc = mul[pos][i]
                            if ev_wat and wats[i]:
                                summ.add(seq, qu, wat, rc, cy, rdlens[i], mapqs[i])
                                summs_wat[pos].add(seq, qu, wat, rc, cy, rdlens[i], mapqs[i])
                            if ev_cri and not wats[i]:
                                summ.add(seq, qu, wat, rc, cy, rdlens[i], mapqs[i])
                                summs_cri[pos].add(seq, qu, wat, rc, cy, rdlens[i], mapqs[i])
    return summ, summs_wat, summs_cri


class TestTabIval(unittest.TestCase):
    
    def test_tab_ival_1(self):
        pos = 7
        flag = 0
        fw = True
        wat = True
        flag |= 16 if fw else 0
        cigars = [(4, 1), (0, 10)]  # 1S10M
        mapq = 40
        alsc = 30
        refseq = "AGCATGCAACT"  # 3 Cs, 2 Gs, no CpGs
        rdseq = "AGCATGCAACT"
        rdqual = "IIIIIIIIIII"
        recs = [(pos, flag, wat, not fw, cigars, alsc, mapq, rdseq, rdqual)]
        summ, summs_fw, summs_rc = tab_ival(1, 9, lambda _x, y: recs,
                                            lambda _x, _y, _lpad, _rpad: (refseq[_x:_x+_y], 0, 0), 0, True)
        self.assertEqual(0, len(summs_fw.keys()))
        self.assertEqual(0, len(summs_rc.keys()))
        summ, summs_fw, summs_rc = tab_ival(1, 9, lambda _x, y: recs,
                                            lambda _x, _y, _lpad, _rpad: (refseq[_x:_x+_y], 0, 0), 0, False)
        self.assertEqual(3, len(summs_fw.keys()))
        self.assertEqual(2, len(summs_rc.keys()))

if sys.argv[-1] == "--test":
    unittest.main(argv=[sys.argv[0]])
    sys.exit()

if sys.argv[-1] == "--version":
    print version
    sys.exit()

args = parser.parse_args()
if args.loci is None and args.rand_loci is None:
    raise RuntimeError("Neither --locus nor --random-loci specified")

logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', datefmt='%m/%d/%y-%H:%M:%S',
                    level=logging.DEBUG if args.verbose else logging.INFO)

# Open the FASTA files in 
logging.info('Parsing FASTA files...')
if args.fasta_pickle is not None:
    fa_idx = reference.ReferencePicklable(args.fasta, args.fasta_pickle)
else:
    fa_idx = reference.ReferenceSimple(args.fasta)
fa_sanity = dict()
if args.sanity:
    from Bio import SeqIO
    for fa_fn in args.fasta:
        handle = open(fa_fn, "rU")
        tmp = SeqIO.to_dict(SeqIO.parse(fa_fn, "fasta"))
        fa_sanity.update(tmp)

# Index the BAM files if needed
logging.info('Opening and indexing BSmooth BAM files...')
for bm in args.bsbam:
    if not os.path.exists(bm + ".bai"):
        logging.info('Indexing BAM file "%s" ... ' % bm)
        disp = pysam.SamtoolsDispatcher("index", None)
        disp(bm)

# Open BAM files
bsbams = [pysam.Samfile(x, "rb") for x in args.bsbam]


def aln_get(ref, off1, off2):
    """ Get all alignments from all bsbams overlapping the given reference
        interval """
    for f in bsbams:
        for rec in f.fetch(ref, off1, off2):
            seq = rec.opt('YO')
            assert seq is not None, str(rec)
            ori = rec.opt('XB')
            wat, rc = False, False
            # Parse Watson/Crick and forward/reverse-comp
            if ori is not None and ori[0] == 'W':
                wat = True
            if ori is not None and ori[-1] == 'R':
                rc = True
            if wat and not rc and args.rc_wat_fw:
                seq = revcomp(seq)
            elif not wat and not rc and args.rc_cri_fw:
                seq = revcomp(seq)
            elif wat and rc and args.rc_wat_rc:
                seq = revcomp(seq)
            elif not wat and rc and args.rc_cri_rc:
                seq = revcomp(seq)
            # Get POS and correct if necessary
            pos = rec.pos
            if not wat:
                pos += args.cri_add
            # Get CIGAR and correct if necessary
            cigar = rec.cigar
            if not wat and args.cri_cigar_rev and cigar is not None:
                cigar = cigar[::-1]
            tags = dict()
            for t in rec.tags:
                tags[t[0]] = t[1]
            alsc = tags.get('AS', None)
            yield pos, rec.flag, wat, rc, cigar, alsc, rec.mapq, seq, rec.qual


def go():
    # Handle the intervals specified by the caller
    if args.loci is not None:
        logging.info('Handling intervals specified by caller...')
        for loc in args.loci:
            logging.info('  %s' % loc)
            ch, off = string.split(loc, ':')
            off1, off2 = string.split(off, '-')
            off1, off2 = int(off1), int(off2)
            assert off1 > 0, "--locus range is 1-based; limits must be >= 1"
            assert off2 > 0, "--locus range is 1-based; limits must be >= 1"
            # Make them 0-based
            ln = off2 - off1
            off1 -= 1
            summ, summs_fw, summs_rc = \
                tab_ival(off1,
                         ln,
                         lambda _x, _y: aln_get(ch, _x, _y),
                         lambda _x, _y, _lpad, _rpad: fa_idx.get(ch, _x, _y, _lpad, _rpad),
                         args.npad,
                         filt_mapq=lambda _x: _x < args.min_mapq,
                         filt_rdl=lambda _x: _x < args.min_rdl,
                         just_cpg=not args.all_C)
            print_summaries(sys.stdout, ch, summs_fw, summs_rc, args.npad, merge=args.merge, head=True)

try:
    pr = None
    if args.profile:
        import cProfile
        import pstats
        import StringIO
        pr = cProfile.Profile()
        pr.enable()
    go()
    if args.profile:
        pr.disable()
        s = StringIO.StringIO()
        sortby = 'tottime'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats(30)
        print s.getvalue()
except IOError, e:
    if e.errno == errno.EPIPE:
        pass  # EPIPE error; do nothing (don't print error message)
    else:
        raise e  # Other error; pass it along
