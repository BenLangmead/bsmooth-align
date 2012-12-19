#!/usr/bin/perl -w

#
# Copyright (C) 2012, Ben Langmead <langmea@cs.jhu.edu>
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# The BSmooth Alignment Pipeline is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# The BSmooth Alignment Pipeline.  If not, see <http://www.gnu.org/licenses/>.
#

##
# BtlBio::Align::Bisulfite::BsEvidence
#
# Encapsulates a piece of bisulfite evidence, be it pro- or anti-methylation.
#
# Author: Ben Langmead
# Email: langmea@cs.jhu.edu
#

package BtlBio::Align::Bisulfite::BsEvidence;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(parse_meth_ev_record
             parse_meth_ev_records_filehandle
             parse_meth_ev_records_filename
             parse_meth_ev_records_dirname);

use strict;
use warnings;
use Carp;

use FindBin qw($Bin);
use lib "$Bin/../../..";
use BtlBio::Alphabet::DNA;
use BtlBio::Alignment::Util;
use BtlBio::Read::Read;
use BtlBio::Util::File;

##
# Create a new bisulfite evidence object.  An evidence object is just an array
# ref with fields for:
#
# 1. Reference name
# 2. Reference offset
# 3. Evidence allele
# 4. Evidence strand (watson/crick)
# 5. Evidence orientation (fw read/rc read)
# 6. Mate id (0/1/2)
# 7. Positional quality value 1
# 8. Positional quality value 2 (if applicable)
# 9. Sequencing cycle
# 10. Alignment length
# 11. Alignment score
# 12. Mapping quality (if applicable)
#
sub new {
	my (
		$class,
		$rdid,    # read id
		$tname,   # ref sequence name
		$toff,    # ref sequence offset (1-based)
		$ev,      # A/C/G/T (w/r/t the strand with C)
		$watson,  # true -> aligned to Watson strand
		$fw,      # true -> fw read aligned
		$flags,   # SAM flags
		$qu1,     # quality value 1
		$qu2,     # quality value 2
		$cy,      # sequencing cycle (0-based offset from 5')
		$allen,   # alignment length
		$score,   # alignment score
		$mapq     # alignment mapping quality
	) = @_;
	$rdid = $rdid || "?";
	defined($flags) || confess("Must define SAM flags when making BsEvidence");
	my $watson_ev = $watson ? $ev : revcomp($ev);
	my ($xqu1, $xqu2) = (ord($qu1)-33, ord($qu2)-33);
	$xqu1 = 0 if $xqu1 < 0;
	$xqu2 = $xqu1 if $xqu2 < 0;
	my $qual_summ = ($xqu1 + $xqu2)/2;
	return bless {
		_rdid      => $rdid,
		_tname     => $tname,
		_toff      => $toff,
		_ev        => $ev,
		_watson_ev => $watson_ev,
		_watson    => $watson,
		_fw        => $fw,
		_flags     => $flags,
		_qu1       => $qu1,
		_qu2       => $qu2,
		_qual_summ => $qual_summ,
		_cy        => $cy,
		_allen     => $allen,
		_score     => $score,
		_mapq      => $mapq
	}, $class;
}

sub rdid      { return $_[0]->{_rdid}   }
sub read_id   { return $_[0]->{_rdid}   }
sub ref_name  { return $_[0]->{_tname}  }
sub ref_off   { return $_[0]->{_toff}   }
sub ev        { return $_[0]->{_ev}     }
sub watson    { return $_[0]->{_watson} }
sub fw        { return $_[0]->{_fw}     }
sub flags     { return $_[0]->{_flags}  }
sub cycle     { return $_[0]->{_cy}     }
sub aln_len   { return $_[0]->{_allen}  }
sub aln_score { return $_[0]->{_score}  }
sub aln_mapq  { return $_[0]->{_mapq}   }
sub qual1     { return $_[0]->{_qu1}    }
sub qual2     { return $_[0]->{_qu2}    }
sub quals     { return ($_[0]->{_qu1}, $_[0]->{_qu2}) }

sub to_record {
	return join("\t", (
		$_[0]->{_tname},
		$_[0]->{_toff},
		$_[0]->{_rdid},
		$_[0]->{_ev},
		$_[0]->{_watson},
		$_[0]->{_fw},
		$_[0]->{_flags},
		$_[0]->{_qu1},
		$_[0]->{_qu2},
		$_[0]->{_cy},
		$_[0]->{_allen},
		$_[0]->{_score},
		$_[0]->{_mapq}));
};

##
# Given a filehandle for a methylation evidence file, parse it into a list of
# BsEvidence objects.
#
sub parse_meth_ev_records_filehandle($$) {
	my ($fh, $list) = @_;
	while(readline $fh) {
		push @$list, parse_meth_ev_record($_);
	}
}

##
# Given a filename for a methylation evidence file, parse it into a list of
# BsEvidence objects.
#
sub parse_meth_ev_records_filename($$) {
	my ($fn, $list) = @_;
	return parse_meth_ev_records_filehandle(openex($fn), $list);
}

##
# Given a directory name for a directory of methylation evidence, parse it into
# a list of BsEvidence objects.
#
sub parse_meth_ev_records_dirname($$) {
	my ($dirn, $list) = @_;
	for my $fn (<$dirn/*.tsv>) { parse_meth_ev_records_filename($fn, $list); }
}

##
# Given a line from a methylation evidence file, parse it into a BsEvidence
# object and return it.
#
sub parse_meth_ev_record($) {
	my $line = shift;
	chomp($line);
	my @ts = split(/\t/, $line, -1);
	my ($tname, $toff, $rdid, $ev, $watson, $fw, $flags, $qu1, $qu2,
	    $cy, $allen, $score, $mapq) = @ts;
	defined($mapq) ||
		croak("Not enough tokens in methylation evidence record:\n$line");
	$toff  == int($toff ) || die "Bad toff:\n$line";
	$mapq  == int($mapq ) || die "Bad mapq:\n$line";
	$cy    == int($cy   ) || die "Bad cy:\n$line";
	$allen == int($allen) || die "Bad allen:\n$line";
	$score == int($score) || die "Bad score:\n$line";
	my $methev = BtlBio::Align::Bisulfite::BsEvidence->new(
		$rdid,
		$tname,
		$toff,
		$ev,
		$watson,
		$fw,
		$flags,
		$qu1,
		$qu2,
		$cy,
		$allen,
		$score,
		$mapq);
	return $methev;
}

#
# Simple tests
#

sub _test_1() {
	print STDERR "Testing parse_meth_ev_record 1 ... ";
	my $line = "ref1	12	read1	T	1	1	2	I	-1	8	9	40	255\n";
	my $ev = parse_meth_ev_record($line);
	$ev->{_rdid}   eq "read1" || die;
	$ev->{_tname}  eq "ref1"  || die;
	$ev->{_toff}   == 12      || die;
	$ev->{_ev}     eq "T"     || die;
	$ev->{_watson} == 1       || die;
	$ev->{_fw}     == 1       || die;
	$ev->{_flags}  == 2       || die;
	$ev->{_qu1}    eq "I"     || die;
	$ev->{_qu2}    eq "-1"    || die;
	$ev->{_cy}     == 8       || die;
	$ev->{_allen}  == 9       || die;
	$ev->{_score}  == 40      || die;
	$ev->{_mapq}   == 255     || die;
	print STDERR "PASSED\n";
}

sub _test_2() {
	print STDERR "Testing parse_meth_ev_record 2 ... ";
	my $line = "ref2	15	rd2	A	0	1	1	3	-1	11	30	-20	0";
	my $ev = parse_meth_ev_record($line);
	$ev->{_rdid}   eq "rd2"  || die;
	$ev->{_tname}  eq "ref2" || die;
	$ev->{_toff}   == 15     || die;
	$ev->{_ev}     eq "A"    || die;
	$ev->{_watson} == 0      || die;
	$ev->{_fw}     == 1      || die;
	$ev->{_flags} == 1       || die;
	$ev->{_qu1}    eq "3"    || die;
	$ev->{_qu2}    eq "-1"   || die;
	$ev->{_cy}     == 11     || die;
	$ev->{_allen}  == 30     || die;
	$ev->{_score}  == -20    || die;
	$ev->{_mapq}   == 0      || die;
	print STDERR "PASSED\n";
}

sub _test() {
	_test_1();
	_test_2();
}

# Call _test() if the module is being run directly.
_test() unless caller();

1;
