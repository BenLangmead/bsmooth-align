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
# BtlBio::Read::Qseq
#
# Routines converting between Qseq files and BtlBio::Read objects.
#
# Author: Ben Langmead
# Email: langmea@cs.jhu.edu
#

package BtlBio::Read::Qseq;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(parse_qseq_read);

use strict;
use warnings;
use Carp;
use FindBin qw($Bin);
use lib "$Bin/../..";
use BtlBio::Read::Read;

our $format_str = "Qseq";
our $format_ext = "_qseq.txt";

##
# Given a filehandle, parses a record from it and returns a BtlBio::Read
# object.  The mateid field is always set to 0.
#
sub parse_qseq_read($) {
	my ($fh) = @_;
	my $orig = "";
	# Parse name line
	my $name = undef;
	my ($machname, $runnum, $lanenum, $tilenum, $xcoord, $ycoord,
	    $index, $matenum, $seq, $qual, $filter);
	while(readline $fh) {
		# Skip empty lines
		next if /^\s*$/;
		$orig .= $_;
		chomp;
		(
			$machname,
			$runnum,
			$lanenum,
			$tilenum,
			$xcoord,
			$ycoord,
			$index,
			$matenum,
			$seq,
			$qual,
			$filter) = split(/\t/, $_, -1);
		defined($filter) || croak("Expected 11 fields, got:\n$orig");
		$name = join(
			"_", $machname, $runnum, $lanenum, $tilenum, $xcoord, $ycoord, $index);
		if($matenum == 1 || $matenum == 2) {
			$name .= "/$matenum";
		}
		last;
	}
	return undef unless defined($name);
	return BtlBio::Read::Read->new(
		$name,
		$seq,
		$qual,
		0,
		0,  # color
		$orig);
}

##
# Given a pair of filehandles, parses a pair of FASTQ records from them and
# returns a pair of BtlBio::Read objects.  Reads parsed from the first
# filehandle argument have mateid set to 1, and reads parsed from the second
# have mateid set to 2.
#
sub parse_qseq_pair($$) {
	my ($fh1, $fh2) = @_;
	my ($rd1, $rd2) = (parse_qseq_read($fh1), parse_qseq_read($fh2));
	defined($rd1) == defined($rd2) ||
		croak("Paired input FASTQ files had different lengths");
	# Remove mate designation if necessary
	$rd1->{_name_canonical} eq $rd2->{_name_canonical} ||
		croak("Paired-up reads had different names: ".
		      $rd1->{_orig}."--\n".$rd2->{_orig});
	return ($rd1, $rd2);
}

#
# Simple tests
#

sub _test_1() {
	print STDERR "Testing parse_qseq_read 1 ... ";
	my $ostr = "";
	open(my $ofh, '>', \$ostr);
	print {$ofh} qq!
MACH1	70	5	2101	4059	2174	A	1	ACGT	acbd	1
!;
	open(my $ifh, '<', \$ostr);
	my $rd1 = parse_qseq_read($ifh);
	my $rd2 = parse_qseq_read($ifh);
	defined($rd2) && die "Expected third call to parse_qseq_read to fail";
	
	$rd1->name eq "MACH1_70_5_2101_4059_2174_A/1" || die "bad name: ".$rd1->name;
	$rd1->seq  eq "ACGT" || die "bad seq";
	$rd1->qual eq "acbd" || die "bad qual";
	print STDERR "PASSED\n";
}

sub _test_2() {
	print STDERR "Testing parse_qseq_read 2 ... ";
	my $ostr = "";
	open(my $ofh, '>', \$ostr);
	my $qseq =
qq!HWI-ST369	70	5	2101	4059	2174	TCGAAA	1	ACCCCTTTCAAAAACCTCTACCTTCTTCTCTCTAGAAATTATCATAACTAAAACGCTTCCTAAAACATCAACCCCTAACTCTAAAACCAAAACTCCTCTA	fffffffffffffeeffde`fffffffefeecff^WcSccZb_^aaaddd_c`VMcb\`c__a_^EGIZT[ZGXRXUSU\_JZW^c\`ccceeebR_TSR	1
HWI-ST369	70	5	2101	4906	2084	TCGAAA	2	TCCCTCAAAAACTAAAATACATCCATATCCTATATCCCTTTCCCCTTACATCGTATAATCTAAAAACTACAAAAACCATTTAAACGTACTTCTAAAAATA	ccc_cbd_aaeeeeeedcdYdebcdacc``c`cccdc]ddddScceebV`bc^`Zd``d^__]TR_`a^BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	1
!;
	print {$ofh} $qseq;
	open(my $ifh, '<', \$ostr);
	my $rd1 = parse_qseq_read($ifh);
	my $rd2 = parse_qseq_read($ifh);
	my $rd3 = parse_qseq_read($ifh);
	defined($rd3) && die "Expected third call to parse_qseq_read to fail";
	
	$rd1->name eq "HWI-ST369_70_5_2101_4059_2174_TCGAAA/1" || die "bad name: ".$rd1->name;
	$rd1->seq  eq "ACCCCTTTCAAAAACCTCTACCTTCTTCTCTCTAGAAATTATCATAACTAAAACGCTTCCTAAAACATCAACCCCTAACTCTAAAACCAAAACTCCTCTA" || die "bad seq";
	$rd1->qual eq "fffffffffffffeeffde`fffffffefeecff^WcSccZb_^aaaddd_c`VMcb\`c__a_^EGIZT[ZGXRXUSU\_JZW^c\`ccceeebR_TSR" || die "bad qual";

	$rd2->name eq "HWI-ST369_70_5_2101_4906_2084_TCGAAA/2" || die "bad name: ".$rd2->name;
	$rd2->seq  eq "TCCCTCAAAAACTAAAATACATCCATATCCTATATCCCTTTCCCCTTACATCGTATAATCTAAAAACTACAAAAACCATTTAAACGTACTTCTAAAAATA" || die "bad seq";
	$rd2->qual eq "ccc_cbd_aaeeeeeedcdYdebcdacc``c`cccdc]ddddScceebV`bc^`Zd``d^__]TR_`a^BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB" || die "bad qual";

	print STDERR "PASSED\n";
}

sub _test() {
	_test_1();
	_test_2();
}

# Call _test() if the module is being run directly.
_test() unless caller();

1;
