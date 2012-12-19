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
# BtlBio::Align::Bisulfite::BsFilter
#
# Object for filtering BsEvidence with a variety of criteria.
#
# Author: Ben Langmead
# Email: langmea@cs.jhu.edu
#

package BtlBio::Align::Bisulfite::BsFilter;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw();

use strict;
use warnings;
use Carp;
use FindBin qw($Bin);
use lib "$Bin/../../..";
use BtlBio::Align::Bisulfite::BsEvidence;

##
# Create a new bisulfite evidence object.  An evidence object is just an array
# ref with fields for:
#
sub new {
	my (
		$class,
		$mapq_min,
		$qual_min,
		$cyc_min,
		$cyc_max,
		$cyc_trim_begin,
		$cyc_trim_end,
		$rdlen_min) = @_;
	return bless {
		_mapq_min => $mapq_min, # evidence from alignments w/ MAPQ < this is ignored
		_qual_min => $qual_min, # minimum quality value to allow in
		_cyc_min  => $cyc_min,  # evidence from cycles < min are filtered
		_cyc_max  => $cyc_max,  # evidence from cycles > max are filtered
		_cyc_trim_begin => $cyc_trim_begin, # cycles to filter from beginning
		_cyc_trim_end   => $cyc_trim_end,   # cycles to filter from end
		_rdlen_min      => $rdlen_min       # minimum read length
	}, $class;
}

##
# Given a piece of evidence, return flags indicating which filters it passed.
#
sub filter($$) {
	my ($self, $ev) = @_;
	# Filter the evidence based on MAPQ?
	my $filt_mapq = ($ev->{_mapq} < $self->{_mapq_min}) || 0;
	# Filter the evidence based on quals?
	my $qual = $ev->{_qual_summ};
	my $filt_qual = ($qual < $self->{_qual_min}) || 0;
	# Filter the evidence based on cycle?
	my $filt_cyc = 0;
	$filt_cyc = 1 if $self->{_cyc_min} != -1 && $ev->{_cy} < $self->{_cyc_min};
	$filt_cyc = 1 if $self->{_cyc_max} != -1 && $ev->{_cy} > $self->{_cyc_max};
	$filt_cyc = 1 if $ev->{_cy} < $self->{_cyc_trim_begin};
	$filt_cyc = 1 if $ev->{_cy} >= $ev->{_allen} - $self->{_cyc_trim_end};
	# Filter the evidence based on nucleotide?
	my $al = $ev->{_watson_ev};
	my $filt_nuc = ($al ne "C" && $al ne "T") || 0;
	# Filter based on read length
	my $filt_rdl = ($ev->{_allen} < $self->{_rdlen_min}) || 0;
	defined($filt_mapq) || die;
	defined($filt_qual) || die;
	defined($filt_cyc) || die;
	defined($filt_rdl) || die;
	defined($filt_nuc) || die;
	return ($filt_mapq, $filt_qual, $filt_cyc, $filt_rdl, $filt_nuc);
}

#
# Simple tests
#

sub _test_1() {
	print STDERR "Testing BsFilter::filter 1 ... ";
	my $filt = BtlBio::Align::Bisulfite::BsFilter->new(
		10,       # MAPQ min
		10,       # qual min
		0,        # cyc min
		9999999,  # cyc max
		3,        # cyc trim begin
		3,        # cyc trim end
		100);     # min read length
	my $ev = BtlBio::Align::Bisulfite::BsEvidence->new(
		"ref1",   # ref name
		10,       # ref offset
		"read1",  # read ID
		"C",      # allele
		1,        # watson
		1,        # fw
		1,        # flags
		chr(30+33),# qual 1
		chr(40+33),# qual 2
		9,        # cyc
		101,      # alignment len
		-55,      # alignment score
		40);      # MAPQ
	my ($filt_mapq, $filt_qual, $filt_cyc, $filt_rdl, $filt_nuc) = $filt->filter($ev);
	defined($filt_mapq) || die;
	defined($filt_qual) || die;
	defined($filt_cyc) || die;
	defined($filt_rdl) || die;
	defined($filt_nuc) || die;
	!$filt_mapq || die;
	!$filt_qual || die;
	!$filt_cyc  || die;
	!$filt_rdl  || die;
	!$filt_nuc  || die;
	print STDERR "PASSED\n";
}

sub _test_2() {
	print STDERR "Testing BsFilter::filter 2 ... ";
	my $filt = BtlBio::Align::Bisulfite::BsFilter->new(
		10,       # MAPQ min
		10,       # qual min
		0,        # cyc min
		9999999,  # cyc max
		3,        # cyc trim begin
		3,        # cyc trim end
		100);     # minimum read length
	my $ev = BtlBio::Align::Bisulfite::BsEvidence->new(
		"ref1",   # ref name
		10,       # ref offset
		"read1",  # read ID
		"T",      # allele
		0,        # watson
		1,        # fw
		2,        # flags
		chr(7+33),# qual 1
		chr(4+33),# qual 2
		49,       # cyc
		50,       # alignment len
		-55,      # alignment score
		9);      # MAPQ
	my ($filt_mapq, $filt_qual, $filt_cyc, $filt_rdl, $filt_nuc) = $filt->filter($ev);
	defined($filt_mapq) || die;
	defined($filt_qual) || die;
	defined($filt_cyc) || die;
	defined($filt_rdl) || die;
	defined($filt_nuc) || die;
	$filt_mapq || die;
	$filt_qual || die;
	$filt_cyc  || die;
	$filt_rdl  || die;
	$filt_nuc  || die;
	print STDERR "PASSED\n";
}

sub _test_3() {
	print STDERR "Testing BsFilter::filter 3 ... ";
	my $filt = BtlBio::Align::Bisulfite::BsFilter->new(
		10,       # MAPQ min
		5,        # qual min
		0,        # cyc min
		9999999,  # cyc max
		1,        # cyc trim begin
		1,        # cyc trim end
		100);     # minimum read length
	my $ev = BtlBio::Align::Bisulfite::BsEvidence->new(
		"ref1",   # ref name
		10,       # ref offset
		"read1",  # read ID
		"G",      # allele
		0,        # watson
		1,        # fw
		0,        # flags
		chr(7+33),# qual 1
		chr(4+33),# qual 2
		100,      # cyc
		101,      # alignment len
		-55,      # alignment score
		10);      # MAPQ
	my ($filt_mapq, $filt_qual, $filt_cyc, $filt_rdl, $filt_nuc) = $filt->filter($ev);
	defined($filt_mapq) || die;
	defined($filt_qual) || die;
	defined($filt_cyc) || die;
	defined($filt_rdl) || die;
	defined($filt_nuc) || die;
	!$filt_mapq || die;
	!$filt_qual || die;
	 $filt_cyc  || die;
	!$filt_rdl  || die;
	!$filt_nuc  || die;
	print STDERR "PASSED\n";
}

sub _test_4() {
	print STDERR "Testing BsFilter::filter 4 ... ";
	my $filt = BtlBio::Align::Bisulfite::BsFilter->new(
		0,        # MAPQ min
		6,        # qual min
		1,        # cyc min
		40,       # cyc max
		1,        # cyc trim begin
		1,        # cyc trim end
		100);     # minimum read length
	my $ev = BtlBio::Align::Bisulfite::BsEvidence->new(
		"ref1",   # ref name
		10,       # ref offset
		"read1",  # read ID
		"G",      # allele
		0,        # watson
		1,        # fw
		1,        # flags
		chr(7+33),# qual 1
		chr(4+33),# qual 2
		41,       # cyc
		101,      # alignment len
		-55,      # alignment score
		0);       # MAPQ
	my ($filt_mapq, $filt_qual, $filt_cyc, $filt_rdl, $filt_nuc) = $filt->filter($ev);
	defined($filt_mapq) || die;
	defined($filt_qual) || die;
	defined($filt_cyc) || die;
	defined($filt_rdl) || die;
	defined($filt_nuc) || die;
	!$filt_mapq || die;
	$filt_qual || die;
	$filt_cyc  || die;
	!$filt_rdl  || die;
	!$filt_nuc  || die;
	print STDERR "PASSED\n";
}

sub _test() {
	_test_1();
	_test_2();
	_test_3();
	_test_4();
}

# Call _test() if the module is being run directly.
_test() unless caller();

1;
