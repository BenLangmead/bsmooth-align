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
# BtlBio::Alignment::SAM
#
# Routines for moving SAM alignments into and out of
# BtlBio::Alignment::Alignment objects.
#
# Author: Ben Langmead
# Email: langmea@cs.jhu.edu
#

package BtlBio::Alignment::SAM;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(parse_sam_record);

use strict;
use warnings;
use Carp;
use FindBin qw($Bin);
use lib "$Bin/../..";
use BtlBio::Alignment::Alignment;

##
# Given an array ref of the fields in a SAM record, parse and return an
# Alignment object representing the alignment.
#
sub parse_sam_record {
	my ($line, $strict, $unaligned) = @_;
	chomp($line);
	my @fs = split(/\t/, $line, -1);
	scalar(@fs) >= 11 || confess("Bad SAM line; not 11 tokens: $line");
	my $flags = $fs[1];
	my $optflags = aln_string_to_opt_hash(join("\t", @fs[11..$#fs]));
	my $inpair = (($flags &   1) != 0);
	my $fw     = (($flags &  16) == 0);
	my $mate1  = (($flags &  64) != 0);
	my $mate2  = (($flags & 128) != 0);
	if(!$inpair) {
		!$mate1 || confess("Can't have both 0x1 and 0x40 set in SAM FLAGS");
		!$mate2 || confess("Can't have both 0x1 and 0x80 set in SAM FLAGS");
	} else {
		$mate1 != $mate2 || confess("If 0x1 is set in SAM FLAGS, either 0x40 or 0x80 must also be set");
	}
	my $md = $optflags->{"MD"}{value}; # could be undef
	my ($score, $color) = (0, 0);
	$score = $optflags->{"AS"}{value} if defined($optflags->{"AS"});
	$color = 1 if defined($optflags->{"CS"});
	my $qname = $fs[0];
	my $seq   = $fs[9];
	my $qual  = $fs[10];
	my $rname = $fs[2];  $rname = undef if $fs[2] eq "*";
	my $pos   = $fs[3];  $pos   = undef if $fs[2] eq "*";
	my $mapq  = $fs[4];  $mapq  = undef if $fs[2] eq "*";
	my $cigar = $fs[5];  $cigar = undef if $fs[2] eq "*";
	my $al = BtlBio::Alignment::Alignment->new(
		$qname,    # QNAME
		$seq,      # SEQ
		$qual,     # QUAL
		$color,    # colorspace?
		$rname,    # RNAME
		$pos,      # POS
		$fw,       # orientation
		$score,    # alignment score
		$mapq,     # MAPQ
		$md,       # MD:Z string
		$cigar,    # CIGAR
		$mate1,    # whether read is mate #1 (if this is known)
		$mate2,    # whether read is mate #2 (if this is known)
		$flags,    # SAM flags (if available)
		$optflags, # options
		$line,     # original record
		$unaligned);
	return $al;
}

#
# Simple tests
#

sub _test_1() {
	print STDERR "Testing parse_sam_record 1 ... ";
	my $line = "GTGTTATATGGCGA	4	*	0	0	*	*	0	0	GTGTTATATGGTGA	ABCDEFGHI12345	XM:i:0\n";
	my $al = parse_sam_record($line);
	!$al->aligned || die;
	$al->score == 0 || die;
	$al->flags == 4 || die "Expected 0, got ".$al->flags;
	print STDERR "PASSED\n";
}

sub _test_2() {
	print STDERR "Testing parse_sam_record 2 ... ";
	my $line = "GTGTTATATGGCGA	0	ref1	4	255	14M	*	0	0	GTGTTATATGGTGA	ABCDEFGHI12345	XA:i:0	MD:Z:14	NM:i:0";
	my $al = parse_sam_record($line);
	$al->aligned || die;
	$al->score == 0 || die;
	$al->flags == 0 || die;
	print STDERR "PASSED\n";
}

sub _test_3() {
	print STDERR "Testing parse_sam_record 3 ... ";
	my $line = "GTGTTATATGGCGA	0	ref1	4	255	14M	*	0	0	GTGTTATATGGTGA	ABCDEFGHI12345	AS:i:-40	XA:i:0	MD:Z:14	NM:i:0";
	my $al = parse_sam_record($line);
	$al->aligned || die;
	$al->score == -40 || die;
	$al->flags == 0 || die;
	$al->{_opts}{XA}{type} eq "i" || die;
	$al->{_opts}{XA}{value} == 0 || die;
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
