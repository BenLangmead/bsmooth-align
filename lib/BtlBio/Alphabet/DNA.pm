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
# BtlBio::Alphabet::DNA
#
# Routines for complementing, reverse complementing, dealing with ambiguous
# IUPAC nucleotides and converting to and from colors.
#
# Author: Ben Langmead
# Email: langmea@cs.jhu.edu
#

package BtlBio::Alphabet::DNA;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(revcompMap
             comp
             plus
             revcomp
             revcomp_in_place
             unambig
             normalize_nucs
             normalize_nucs_inplace);

use strict;
use warnings;
use Carp;

my %revcompMap = (
	"A" => "T",
	"T" => "A",
	"C" => "G",
	"G" => "C",
	"R" => "Y",
	"Y" => "R",
	"M" => "K",
	"K" => "M",
	"S" => "S",
	"W" => "W",
	"B" => "V",
	"V" => "B",
	"H" => "D",
	"D" => "H",
	"N" => "N",
	"." => ".",
	"a" => "t",
	"t" => "a",
	"c" => "g",
	"g" => "c",
	"r" => "y",
	"y" => "r",
	"m" => "k",
	"k" => "m",
	"s" => "s",
	"w" => "w",
	"b" => "v",
	"v" => "b",
	"h" => "d",
	"d" => "h",
	"n" => "n",
);

##
# Map from dinucleotides to the colors used to encode them.  Note that both
# nucleotides must be uppercase.
#
my %colorMap = (
	"AA" => 0,
	"AC" => 1,
	"AG" => 2,
	"AT" => 3,
	"CA" => 1,
	"CC" => 0,
	"CG" => 3,
	"CT" => 2,
	"GA" => 2,
	"GC" => 3,
	"GG" => 0,
	"GT" => 1,
	"TA" => 3,
	"TC" => 2,
	"TG" => 1,
	"TT" => 0,
	"NA" => ".",
	"NC" => ".",
	"NG" => ".",
	"NT" => ".",
	"AN" => ".",
	"CN" => ".",
	"GN" => ".",
	"TN" => ".",
	"NN" => ".",
);

my %dnaEncode = (
	"0" => "A",
	"1" => "C",
	"2" => "G",
	"3" => "T",
	"." => "N",
);

sub colorize {
	my ($s, $dna) = @_;
	my $cstr = "";
	for (0..length($s)-2) {
		my $dinuc = uc substr($s, $_, 2);
		defined($colorMap{$dinuc}) || die "Bad dinuc: $dinuc";
		my $c .= $colorMap{$dinuc};
		$c = $dnaEncode{$c} if $dna;
		$cstr .= $c;
	}
	return $cstr;
}

my %compat = (
	"A" => "A",
	"T" => "T",
	"C" => "C",
	"G" => "G",
	"R" => "AG",
	"Y" => "CT",
	"M" => "AC",
	"K" => "GT",
	"S" => "CG",
	"W" => "AT",
	"B" => "CGT",
	"V" => "ACG",
	"H" => "ACT",
	"D" => "AGT",
	"N" => "N"
);

my %incompat = (
	"A" => "CGT",
	"T" => "ACG",
	"C" => "AGT",
	"G" => "ACT",
	"R" => "CT",
	"Y" => "AG",
	"M" => "GT",
	"K" => "AC",
	"S" => "AT",
	"W" => "CG",
	"B" => "A",
	"V" => "T",
	"H" => "G",
	"D" => "C",
	"N" => "N"
);

my %unambigSet = (
	"A" => 1, "a" => 1,
	"C" => 1, "c" => 1,
	"G" => 1, "g" => 1,
	"T" => 1, "t" => 1
);

##
# Return the complement, incl. if it's IUPAC.
#
sub comp($) {
	my $ret = $revcompMap{$_[0]} || die "Can't reverse-complement '$_[0]'";
	return $ret;
}

##
# Return the complement, incl. if it's IUPAC.
#
sub revcomp {
	my ($ret, $color) = @_;
	$ret = reverse $ret;
	unless($color) {
		for(my $i = 0; $i < length($ret); $i++) {
			substr($ret, $i, 1) = comp(substr($ret, $i, 1));
		}
	}
	return $ret;
}

##
# Reverse-complement a string without making a copy.  Caller passes in a ref
# to the string.
#
sub revcomp_in_place {
	my ($sref, $color) = @_;
	my $half = int(length($$sref)/2);
	for(my $i = 0; $i < $half; $i++) {
		my $tmp = substr($$sref, $i, 1);
		my $swapi = -$i-1;
		if($color) {
			substr($$sref, $i, 1) = substr($$sref, $swapi, 1);
			substr($$sref, $swapi, 1) = $tmp;
		} else {
			substr($$sref, $i, 1) = comp(substr($$sref, $swapi, 1));
			substr($$sref, $swapi, 1) = comp($tmp);
		}
	}
	if((length($$sref) & 1) != 0) {
		if($color) {
			substr($$sref, $half, 1) = substr($$sref, $half, 1);
		} else {
			substr($$sref, $half, 1) = comp(substr($$sref, $half, 1));
		}
	}
}

##
# Return true iff it's unambiguous.
#
sub unambig($) {
	return $unambigSet{$_[0]};
}

##
# Manipulate DNA in an integer-indexed fashion.
#
sub plus($$) {
	my ($c, $amt) = @_;
	my %ctoi = ("A" => 0, "C" => 1, "G" => 2, "T" => 3);
	my %itoc = (0 => "A", 1 => "C", 2 => "G", 3 => "T");
	$c = uc $c;
	defined($ctoi{$c}) || die;
	return $itoc{($ctoi{$c}+$amt) % 4};
}

##
# Turn all non-standard characters to Ns.
#
sub normalize_nucs_inplace($) {
	$$_[0] =~ s/[^ACGTacgt]/N/g;
}

##
# Turn all non-standard characters to Ns.
#
sub normalize_nucs($) {
	my $s = shift;
	$s =~ s/[^ACGTacgt]/N/g;
	return $s;
}

1;
