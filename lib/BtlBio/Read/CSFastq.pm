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
# BtlBio::Read::CSFastq
#
# Routines converting between CSFastq files and BtlBio::Read objects.
#
# Author: Ben Langmead
# Email: langmea@cs.jhu.edu
#

package BtlBio::Read::CSFastq;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(parse_csfastq_read);

use strict;
use warnings;
use Carp;
use FindBin qw($Bin);
use lib "$Bin/../..";
use BtlBio::Read::Read;

our $format_str = "CSFastq";
our $format_ext = ".fq";

##
# Given a filehandle, parses a record from it and returns a BtlBio::Read
# object.  The mateid field is always set to 0.
#
sub parse_csfastq_read($) {
	my ($fh) = @_;
	my $orig = "";
	# Parse name line
	my $name = undef;
	while(readline $fh) {
		# Skip empty lines
		next if /^\s*$/;
		$orig .= $_;
		chomp;
		substr($_, 0, 1) eq "@" ||
			croak("Expected $format_str name line beginning with @, got:\n$_");
		$name = substr($_, 1);
		last;
	}
	return undef unless defined($name);
	# Parse sequence line
	my $seq = readline $fh;
	defined($seq) ||
		croak("Expected $format_str sequence line following name line:\n$orig");
	$orig .= $seq;
	chomp($seq);
	# Parse second name line
	my $name2 = readline $fh;
	defined($name2) ||
		croak("Expected $format_str second name line following name & seq lines:\n$orig");
	$orig .= $name2;
	chomp($name2);
	# Parse quality line
	my $qual = readline $fh;
	defined($qual) ||
		croak("Expected $format_str quality line following name,seq,name2 lines:\n$orig");
	chomp($qual);
	return BtlBio::Read::Read->new(
		$name,
		$seq,
		$qual,
		0,
		1,  # color
		$orig);
}

##
# Given a pair of filehandles, parses a pair of records from them and
# returns a pair of BtlBio::Read objects.  Reads parsed from the first
# filehandle argument have mateid set to 1, and reads parsed from the second
# have mateid set to 2.
#
sub parse_csfastq_pair($$) {
	my ($fh1, $fh2) = @_;
	my ($rd1, $rd2) = (parse_csfastq_read($fh1), parse_csfastq_read($fh2));
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
	print STDERR "Testing parse_csfastq_read 1 ... ";
	my $ostr = "";
	open(my $ofh, '>', \$ostr);
	print {$ofh} qq!
\@Read1
ACGT
+
acbd
!;
	open(my $ifh, '<', \$ostr);
	my $rd1 = parse_csfastq_read($ifh);
	my $rd2 = parse_csfastq_read($ifh);
	defined($rd2) && die "Expected second call to parse_csfastq_read to fail";
	
	$rd1->name eq "Read1"   || die "bad name: ".$rd1->name;
	$rd1->seq  eq "ACGT"    || die "bad seq";
	$rd1->qual eq "acbd"    || die "bad qual";
	$rd1->colorspace        || die "expected colorspace";
	print STDERR "PASSED\n";
}

sub _test_2() {
	print STDERR "Testing parse_csfastq_read 2 ... ";
	my $ostr = "";
	open(my $ofh, '>', \$ostr);
	my $qseq = qq!

\@Read1
ACGT
+
acbd



\@Read2
TGACACA
+
IIIIIII


!;
	print {$ofh} $qseq;
	open(my $ifh, '<', \$ostr);
	my $rd1 = parse_csfastq_read($ifh);
	my $rd2 = parse_csfastq_read($ifh);
	my $rd3 = parse_csfastq_read($ifh);
	defined($rd3) && die "Expected third call to parse_csfastq_read to fail";
	
	$rd1->name eq "Read1"   || die "bad name: ".$rd1->name;
	$rd1->seq  eq "ACGT"    || die "bad seq";
	$rd1->qual eq "acbd"    || die "bad qual";
	$rd1->colorspace        || die "expected colorspace";

	$rd2->name eq "Read2"   || die "bad name: ".$rd2->name;
	$rd2->seq  eq "TGACACA" || die "bad seq";
	$rd2->qual eq "IIIIIII" || die "bad qual";
	$rd2->colorspace        || die "expected colorspace";
	print STDERR "PASSED\n";
}

sub _test_3() {
	print STDERR "Testing parse_csfastq_read 3 ... ";
	my $ostr = "";
	open(my $ofh, '>', \$ostr);
	print {$ofh} qq!
\@Read1
A1230
+
acbd
!;
	open(my $ifh, '<', \$ostr);
	my $rd1 = parse_csfastq_read($ifh);
	my $rd2 = parse_csfastq_read($ifh);
	defined($rd2) && die "Expected second call to parse_csfastq_read to fail";
	
	$rd1->name eq "Read1"   || die "bad name: ".$rd1->name;
	$rd1->seq  eq "CGTA"    || die "bad seq";
	$rd1->qual eq "acbd"    || die "bad qual";
	$rd1->primer eq "A"     || die "bad primer";
	$rd1->colorspace        || die "expected colorspace";
	print STDERR "PASSED\n";
}

sub _test() {
	_test_1();
	_test_2();
	_test_3();
}

# Call _test() if the module is being run directly.
_test() unless caller();

1;
