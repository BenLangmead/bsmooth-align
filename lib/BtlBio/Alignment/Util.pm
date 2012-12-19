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
# BtlBio::Alignment::Util
#
# 
#
# Author: Ben Langmead
# Email: langmea@cs.jhu.edu
#

package BtlBio::Alignment::Util;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(parse_cigar parse_md parse_cigar_ivals);

use strict;
use warnings;
use Carp;

##
# Parse a CIGAR string into a string of operators.  Operators are expanded into
# runs where appropriate.  = and X are collapsed into M.
#
sub parse_cigar($) {
	my ($cigar) = @_;
	my $ret = "";
	my $i = 0;
	my ($nms, $nds, $nis) = (0, 0, 0);
	while($i < length($cigar)) {
		substr($cigar, $i) =~ /^([0-9]+)/;
		defined($1) || croak("Could not parse number at pos $i: '$cigar'");
		my $runlen = $1;
		$i += length($1);
		$i < length($cigar) || confess("Bad cigar string: '$cigar'");
		my $op = substr($cigar, $i, 1);
		defined($op) || croak("Could not parse operation at pos $i: '$cigar'");
		if($op eq "X" || $op eq "=") {
			$op = "M";
		}
		$nms += $runlen if $op eq "M";
		$nds += $runlen if $op eq "D";
		$nis += $runlen if $op eq "I";
		$ret .= ($op x $runlen);
		$i++;
	}
	return ($ret, $nms, $nds, $nis);
}

##
# Parse a CIGAR string that may use the N operator.  Return a list of
# intervals covered by the alignment, where offset 0=the leftmost aligned
# position.
#
sub parse_cigar_ivals($) {
	my $cigar = shift;
	my @ivals = ();
	my ($left, $right) = (0, 0);
	my $i = 0;
	while($i < length($cigar)) {
		substr($cigar, $i) =~ /^([0-9]+)/;
		defined($1) || croak("Could not parse number at pos $i: '$cigar'");
		my $runlen = $1;
		$i += length($1); # skip over run length
		$i < length($cigar) || confess("Bad cigar string: '$cigar'");
		my $op = substr($cigar, $i, 1);
		defined($op) || croak("Could not parse operation at pos $i: '$cigar'");
		$op = "M" if $op eq "X" || $op eq "=";
		$right += $runlen if $op eq "M" || $op eq "D";
		if($op eq "N") {
			$right > $left || die "Expected $right > $left";
			push @ivals, [$left, $right];
			# Calculate new left
			$left = $right + $runlen;
			$right = $left;
		}
		$i++; # skip over operation
	}
	if($right > $left) {
		push @ivals, [$left, $right];
	}
	return \@ivals;
}

##
# Parse an MD:Z string into a string with length equal to query length.  Each
# position contains either a space, if the read matches the reference at that
# position, or a character, if the reference contains a character that doesn't
# match its opposite in the alignment.  In the latter case, the character in
# the string is the reference character.
#
sub parse_md($) {
	my ($md) = @_;
	my $i = 0;
	my $ret = "";
	while($i < length($md)) {
		# Starts with a number?
		my $ch = substr($md, $i, 1);
		if($ch =~ /[0-9]/) {
			# Parse the number off the beginning
			substr($md, $i) =~ /^([0-9]+)/;
			defined($1) || croak("Could not parse number at pos $i: '$md'");
			my $runlen = $1;
			$ret .= (" " x $runlen) if $runlen > 0;
			$i += length($runlen);
		} elsif($ch eq "^") {
			# Read gap
			$i++;
			substr($md, $i) =~ /^([A-Za-z]+)/;
			defined($1) || croak("Could not parse read gap at pos $i: '$md'");
			my $chrs = $1;
			$i += length($chrs);
			$ret .= $chrs;
		} else {
			# DNA character
			$ch =~ /[A-Z.]/i || croak("Bad char '$ch' at pos $i: '$md'");
			$ret .= $ch;
			$i++;
		}
	}
	return $ret;
}

sub _test_parse_cigar1() {
	print STDERR "Testing parse_cigar 1 ... ";
	my ($cs, $nms, $nds, $nis) = parse_cigar("9=9X9=");
	$cs eq ("M" x 27) || die;
	print STDERR "PASSED\n";
}

sub _test_parse_cigar2() {
	print STDERR "Testing parse_cigar 2 ... ";
	my ($cs, $nms, $nds, $nis) = parse_cigar("9=9X3D9=2I");
	$cs eq "MMMMMMMMMMMMMMMMMMDDDMMMMMMMMMII" || die;
	print STDERR "PASSED\n";
}

sub _test_parse_md1() {
	print STDERR "Testing parse_md 1 ... ";
	parse_md("15G4C4A4") eq "               G    C    A    " || die;
	print STDERR "PASSED\n";
}

sub _test_parse_md2() {
	print STDERR "Testing parse_md 2 ... ";
	parse_md("5AACC6") eq "     AACC      " || die;
	print STDERR "PASSED\n";
}

sub _test_parse_md3() {
	print STDERR "Testing parse_md 3 ... ";
	parse_md("25^T28") eq "                         T                            " || die;
	print STDERR "PASSED\n";
}

sub _parse_cigar_ivals_1() {
	print STDERR "Testing parse_cigar_ivals 1 ... ";
	my $ivals = parse_cigar_ivals("9=9X3D9=2I");
	scalar(@$ivals) == 1 || die "Expected 1 intervals, got ".scalar(@$ivals);
	$ivals->[0]->[0] == 0 || die;
	$ivals->[0]->[1] == 30 || die;
	print STDERR "PASSED\n";
}

sub _parse_cigar_ivals_2() {
	print STDERR "Testing parse_cigar_ivals 2 ... ";
	my $ivals = parse_cigar_ivals("9=9X3I9=2D");
	scalar(@$ivals) == 1 || die;
	$ivals->[0]->[0] == 0 || die;
	$ivals->[0]->[1] == 29 || die;
	print STDERR "PASSED\n";
}

sub _parse_cigar_ivals_3() {
	print STDERR "Testing parse_cigar_ivals 3 ... ";
	my $ivals = parse_cigar_ivals("9=9X3N9=2I");
	scalar(@$ivals) == 2 || die;
	$ivals->[0]->[0] == 0 || die;
	$ivals->[0]->[1] == 18 || die;
	$ivals->[1]->[0] == 21 || die;
	$ivals->[1]->[1] == 30 || die;
	print STDERR "PASSED\n";
}

sub _test() {
	_test_parse_cigar1();
	_test_parse_cigar2();
	_test_parse_md1();
	_test_parse_md2();
	_test_parse_md3();
	_parse_cigar_ivals_1();
	_parse_cigar_ivals_2();
	_parse_cigar_ivals_3();
}

# Call _test() if the module is being run directly.
_test() unless caller();

1;
