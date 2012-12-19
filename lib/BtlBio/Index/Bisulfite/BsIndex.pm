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
# Bisulfite.pm
#
# Routines for in-silico bisulfite conversion of sequence data.
#
# Author: Ben Langmead
# Email: langmea@cs.jhu.edu
#

package BtlBio::Index::Bisulfite::BsIndex;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(bsc_convert_stream bsc_flush);

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../../..";
use BtlBio::Util::File;
use BtlBio::Alphabet::DNA;
use List::Util qw[min max];

##
# Given a reference to a sequence read from a fasta file, reverse-complement
# it, BSC-convert it, and write the result to the given filehandle.
#
sub bsc_flush {
	my ($rbuf, $ofh, $fw, $line_len) = @_;
	defined($rbuf) || die;
	defined($$rbuf) || die;
	return if length($$rbuf) == 0;
	my $len = length($$rbuf);
	$line_len = $line_len || 60;
	# Flush sequence so far
	revcomp_in_place($rbuf) unless $fw;
	$$rbuf =~ s/C/T/ig;
	$line_len > 0 || die;
	for(my $i = 0; $i < $len; $i += $line_len) {
		my $amt = min($line_len, $len-$i);
		print {$ofh} "".substr($$rbuf, $i, $amt)."\n";
	}
}

##
# Given an input filehandle that dispenses lines from one or more input FASTA
# files, write BSC-converted version to the given output filehandle.
#
# If 'watson' is true, write the bisulfite-converted Watson sequence (i.e. the
# sequence in the FASTA file).  Otherwise, write the bisulfite-converted Crick
# sequence (i.e. the reverse complement of each sequence in the FASTA file).
#
sub bsc_convert_stream($$$) {
	my ($ifh, $ofh, $watson) = @_;
	if($watson) {
		while(readline $ifh) {
			if(substr($_, 0, 1) eq ">") {
				# Name line
				print {$ofh} $_;
			} else {
				# Sequence line
				# Replace every C with T
				s/C/T/ig;
				print {$ofh} $_;
			}
		}
	} else {
		my $buf = "";
		my $line_len = 60;
		while(readline $ifh) {
			if(substr($_, 0, 1) eq ">") {
				# Name line
				# Flush previous sequence
				bsc_flush(\$buf, $ofh, 0, $line_len);
				$buf = "";
				print {$ofh} $_;
			} else {
				# Sequence line
				chomp;
				$line_len = max($line_len, length($_));
				$buf .= $_;
			}
		}
		bsc_flush(\$buf, $ofh, 0, $line_len);
	}
}

#
# TESTS
#

sub _test_bsc_flush_1() {
	print STDERR "Testing bsc_flush 1 ... ";
	my $fa = "AAACAGATCACCCGCTGAGCGGGTTATCTGTT";
	my $ostr = "";
	open(my $ofh, '>', \$ostr);
	bsc_flush(\$fa, $ofh, 1);
	close($ofh);
	$ostr =~ s/^\s+//;
	$ostr =~ s/\s+$//;
	my $ex = "AAATAGATTATTTGTTGAGTGGGTTATTTGTT";
	chomp($ostr);
	$ostr eq $ex || die "Expected on top, actual on bottom:\n$ex\n$ostr";
	print STDERR "PASSED\n";
}

sub _test_bsc_flush_2() {
	print STDERR "Testing bsc_flush 2 ... ";
	my $fa = "AAACAGATCACCCGCTGAGCGGGTTATCTGTT";
	my $ostr = "";
	open(my $ofh, '>', \$ostr);
	bsc_flush(\$fa, $ofh, 0);
	close($ofh);
	$ostr =~ s/^\s+//;
	$ostr =~ s/\s+$//;
	my $ex = "AATAGATAATTTGTTTAGTGGGTGATTTGTTT";
	chomp($ostr);
	$ostr eq $ex || die "Expected on top, actual on bottom:\n$ex\n$ostr";
	print STDERR "PASSED\n";
}

sub _test_bsc_convert_stream($$$) {
	my ($fa, $ex, $watson) = @_;
	my $ostr = "";
	open(my $ifh, '<', \$fa);
	open(my $ofh, '>', \$ostr);
	bsc_convert_stream($ifh, $ofh, $watson);
	close($ifh); close($ofh);
	chomp($ostr);
	$ostr =~ s/^\s+//;
	$ostr =~ s/\s+$//;
	$ex eq $ostr || die "\nExpected:\n$ex\n\nActual:\n$ostr";
}

sub _test_bsc_convert_stream_1() {
	print STDERR "Testing bsc_convert_stream 1 ... ";
	_test_bsc_convert_stream(
		">seq1\nAAACAGATCACCCGCTGAGCGGGTTATCTGTT\n>seq2\nCCCCCCCCCCCCCCCC",
		">seq1\nAAATAGATTATTTGTTGAGTGGGTTATTTGTT\n>seq2\nTTTTTTTTTTTTTTTT",
		1);
	print STDERR "PASSED\n";
}

sub _test_bsc_convert_stream_2() {
	print STDERR "Testing bsc_convert_stream 2 ... ";
	_test_bsc_convert_stream(
		">seq1\nAAACAGATCACCCGCTGAGCGGGTTATCTGTT\n>seq2\nCCCCCCCCCCCCCCCC",
		">seq1\nAATAGATAATTTGTTTAGTGGGTGATTTGTTT\n>seq2\nGGGGGGGGGGGGGGGG",
		0);
	print STDERR "PASSED\n";
}

sub _test_bsc_convert_stream_3() {
	print STDERR "Testing bsc_convert_stream 3 ... ";
	_test_bsc_convert_stream(
		">seq1\nCCCCC\n>seq2\nGGGGG",
		">seq1\nTTTTT\n>seq2\nGGGGG",
		1);
	_test_bsc_convert_stream(
		">seq1\nCCCCC\n>seq2\nGGGGG",
		">seq1\nGGGGG\n>seq2\nTTTTT",
		0);
	print STDERR "PASSED\n";
}

sub _test() {
	_test_bsc_flush_1();
	_test_bsc_flush_2();
	_test_bsc_convert_stream_1();
	_test_bsc_convert_stream_2();
	_test_bsc_convert_stream_3();
}

# Call _test() if the module is being run directly.
_test() unless caller();

1;
