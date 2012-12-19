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
# BtlBio::Read::CSFasta
#
# Routines converting between CSFASTA files and BtlBio::Read objects.
#
# Author: Ben Langmead
# Email: langmea@cs.jhu.edu
#

package BtlBio::Read::CSFasta;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(parse_csfasta_read);

use strict;
use warnings;
use Carp;
use FindBin qw($Bin);
use lib "$Bin/../..";
use BtlBio::Read::Read;

our $format_str = "CSFASTA";
our $format_ext = ".csfasta";

##
# Given a filehandle, parses a record from it and returns a BtlBio::Read
# object.  The mateid field is always set to 0.
#
sub parse_csfasta_read($$) {
	my ($fh_fa, $fh_qv) = @_;
	my ($orig_fa, $orig_qv) = ("", "");
	# Parse name lines
	my ($name_fa, $name_qv) = (undef, undef);
	while(readline $fh_fa) {
		my $line_fa = $_;
		my $line_qv = readline $fh_qv;
		defined($line_fa) == defined($line_qv) || die;
		next if /^\s*$/; # Skip empty lines
		$orig_fa .= $line_fa;
		$orig_qv .= $line_qv;
		chomp($line_fa); chomp($line_qv);
		substr($line_fa, 0, 1) eq ">" ||
			croak("Expected $format_str name line beginning with >, got:\n$orig_fa");
		substr($line_qv, 0, 1) eq ">" ||
			croak("Expected $format_str name line beginning with >, got:\n$orig_qv");
		$name_fa = substr($line_fa, 1);
		$name_qv = substr($line_qv, 1);
		$name_fa eq $name_qv ||
			croak("Names from .csfasta and _qv.QV files do not match: ".
			      "'$name_fa', '$name_qv'");
		last;
	}
	return undef unless defined($name_fa);
	# Parse sequence lines
	my $seq  = readline $fh_fa;
	my $qual = readline $fh_qv;
	defined($seq) ||
		croak("Expected $format_str sequence line following name line:\n$orig_fa");
	defined($qual) ||
		croak("Expected $format_str sequence line following name line:\n$orig_qv");
	$orig_fa .= $seq;
	$orig_qv .= $qual;
	chomp($seq); chomp($qual);
	my $rd = BtlBio::Read::Read->new(
		$name_fa,
		$seq,
		$qual,
		0,
		1, # color
		$orig_fa.$orig_qv);
	$rd->asciiize_quals();
	return $rd;
}

##
# Given a pair of filehandles, parses a pair of records from them and
# returns a pair of BtlBio::Read objects.  Reads parsed from the first
# filehandle argument have mateid set to 1, and reads parsed from the second
# have mateid set to 2.
#
sub parse_csfasta_pair($$) {
	my ($fh1_fa, $fh1_qv, $fh2_fa, $fh2_qv) = @_;
	my ($rd1, $rd2) = (
		parse_csfasta_read($fh1_fa, $fh1_qv),
		parse_csfasta_read($fh2_fa, $fh2_qv));
	defined($rd1) == defined($rd2) ||
		croak("Paired input $format_str files had different lengths");
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
	print STDERR "Testing parse_csfasta_read 1 ... ";
	my $ostr_fa = "";
	open(my $ofh_fa, '>', \$ostr_fa);
	my $ostr_qv = "";
	open(my $ofh_qv, '>', \$ostr_qv);
	print {$ofh_fa} qq!
>Read1
A1230
!;
	print {$ofh_qv} qq!
>Read1
10 20 30 10
!;
	open(my $ifh_fa, '<', \$ostr_fa);
	open(my $ifh_qv, '<', \$ostr_qv);
	my $rd1 = parse_csfasta_read($ifh_fa, $ifh_qv);
	my $rd2 = parse_csfasta_read($ifh_fa, $ifh_qv);
	defined($rd2) && die "Expected second call to parse_fasta_read to fail";
	
	$rd1->name   eq "Read1" || die "bad name: ".$rd1->name;
	$rd1->seq    eq "CGTA"  || die "bad seq";
	$rd1->primer eq "A"     || die "bad pimer";
	$rd1->qual   eq "+5?+"  || die "bad qual: '$rd1->{_qual}'";

	print STDERR "PASSED\n";
}

sub _test_2() {
	print STDERR "Testing parse_csfasta_read 2 ... ";
	my $ostr_fa = "";
	open(my $ofh_fa, '>', \$ostr_fa);
	my $ostr_qv = "";
	open(my $ofh_qv, '>', \$ostr_qv);
	my $fa = qq!

>Read1
C1020130201





>Read2
T10202222




!;
	my $qv = qq!

>Read1
10 7 2 10 1 3 33 2 10 1





>Read2
8 9 2 20 8 27 31 20




!;
	print {$ofh_fa} $fa;
	print {$ofh_qv} $qv;
	open(my $ifh_fa, '<', \$ostr_fa);
	open(my $ifh_qv, '<', \$ostr_qv);
	my $rd1 = parse_csfasta_read($ifh_fa, $ifh_qv);
	my $rd2 = parse_csfasta_read($ifh_fa, $ifh_qv);
	my $rd3 = parse_csfasta_read($ifh_fa, $ifh_qv);
	defined($rd3) && die "Expected third call to parse_csfasta_read to fail";
	
	$rd1->name eq "Read1" || die "bad name: ".$rd1->name;
	$rd1->seq  eq "CAGACTAGAC"    || die "bad seq";
	$rd1->primer eq "C"           || die "bad primer";
	#$rd1->qual eq "10 7 2 10 1 3 33 2 10 1" || die "bad qual";
	$rd1->qual eq "+(#+\"\$B#+\"" || die "bad qual";

	$rd2->name   eq "Read2" || die "bad name: ".$rd2->name;
	#$rd2->seq    eq "T10202222" || die "bad seq";
	$rd2->seq    eq "CAGAGGGG"  || die "bad seq";
	$rd2->primer eq "T"         || die "bad primer";
	$rd2->qual   eq ")*#5)<\@5" || die "bad qual";

	print STDERR "PASSED\n";
}

sub _test() {
	_test_1();
	_test_2();
}

# Call _test() if the module is being run directly.
_test() unless caller();

1;
