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
# BtlBio::Alignment::Colorspace
#
# Routines for manipulating colorspace alignments.
#
# Author: Ben Langmead
# Email: langmea@cs.jhu.edu
#

package BtlBio::Alignment::Colorspace;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(get_cs_cq
             get_trimmed_cs_cq);

use strict;
use warnings;
use Carp;
use FindBin qw($Bin);
use lib "$Bin/../..";
use BtlBio::Alignment::Alignment;

sub get_cs_cq($) {
	my $al = shift;
	$al->color || croak("Must provide a colorspace alignment to get_cs_cq");
	defined($al->{_opts}{"CS"}) ||
		croak("No color sequence info (CS:Z) for aligned colorspace read");
	defined($al->{_opts}{"CQ"}) ||
		croak("No color quality info (CQ:Z) for aligned colorspace read");
	return ($al->opt("CS", "Z"), $al->opt("CQ", "Z"));
}

sub get_trimmed_cs_cq($) {
	my $al = shift;
	my ($cs, $cq) = get_cs_cq($al);
	# Is a primer base present?
	if(substr($cs, 0, 1) =~ /[A-Za-z]/ && substr($cs, 1, 1) =~ /[0-9]/) {
		$cs = substr($cs, 2);
		$cq = substr($cq, 1);
	}
	length($cs) == length($cq) ||
		croak("Error: expected CS:Z and CQ:Z strings to have same length ".
		      "after trimming, got:\n$cs\n$cq");
	return ($cs, $cq);
}

sub _test_get_cs_cq_1() {
	print STDERR "Testing get_cs_cq 1 ... ";
	my %h = ();
	$h{"CS"}{type}  = "Z";
	$h{"CS"}{value} = "A0123";
	$h{"CQ"}{type}  = "Z";
	$h{"CQ"}{value} = "IIII";
	my $al = BtlBio::Alignment::Alignment->new(
		"tmp",  # read name
		"AA",   # read sequence from alignment record
		"II",   # read qualities from alignment record
		1,      # colorspace?
		"ref",  # name of reference sequence ("text")
		0,      # offset into reference
		1,      # true -> aligned to Watson strand
		10,     # alignment score
		10,     # mapping quality
		"2",    # MD:Z string
		"2M",   # CIGAR string
		0, 0,   # mate 1/2
		16,     # SAM flags
		\%h,    # options hash ref
		"");    # original string
	my ($cs, $cq) = get_cs_cq($al);
	$cs eq "A0123" || die;
	$cq eq "IIII"  || die;
	print STDERR "PASSED\n";
}

sub _test() {
	_test_get_cs_cq_1();
}

# Call _test() if the module is being run directly.
_test() unless caller();

1;
