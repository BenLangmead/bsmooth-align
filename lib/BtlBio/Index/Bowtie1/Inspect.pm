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
# BtlBio::Index::Bowtie1::Inspect
#
# For inspecting Bowtie 1 indexes.
#
# Author: Ben Langmead
# Email: langmea@cs.jhu.edu
#

package BtlBio::Index::Bowtie1::Inspect;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(bowtie1_index_restore_string
             bowtie1_index_restore_hash
             bowtie1_index_restore_length_hash
             bowtie1_index_quick_length_hash);

use strict;
use warnings;
use Carp;
use FindBin qw($Bin);
use lib "$Bin/../../..";
use BtlBio::Util::Tool;

##
# Run bowtie-inspect on the given index and return the FASTA output.
#
sub bowtie1_index_restore_string {
	my ($base, $btinspect_exe, $slow) = @_;
	$btinspect_exe = find_tool("bowtie-inspect", $btinspect_exe);
	$slow = $slow || 0;
	$slow = $slow ? "-e" : "";
	my $cmd = "$btinspect_exe $slow $base";
	return `$cmd`;
}

##
# Run bowtie-inspect on the given index and return a hashed version of the Fasta
#
sub bowtie1_index_restore_hash {
	my ($base, $btinspect_exe, $slow) = @_;
	$btinspect_exe = find_tool("bowtie-inspect", $btinspect_exe);
	$slow = $slow || 0;
	$slow = $slow ? "-e" : "";
	my $cmd = "$btinspect_exe $slow $base";
	open(CMD, "$cmd |") || croak("Could not open command '$cmd'");
	my %hash = ();
	my $name = "";
	while(<CMD>) {
		chomp;
		if(/^>/) {
			$name = substr($_, 1);
			$hash{$name} = "";
		} else {
			$hash{$name} .= $_;
		}
	}
	return \%hash;
}

##
# Run bowtie-inspect on the given index and return a hashed version of the Fasta
#
sub bowtie1_index_restore_length_hash {
	my ($base, $btinspect_exe, $slow) = @_;
	$btinspect_exe = find_tool("bowtie-inspect", $btinspect_exe);
	$slow = $slow || 0;
	$slow = $slow ? "-e" : "";
	my $cmd = "$btinspect_exe $slow $base";
	open(CMD, "$cmd |") || croak("Could not open command '$cmd'");
	my %hash = ();
	my $name = "";
	while(<CMD>) {
		chomp;
		if(/^>/) {
			$name = substr($_, 1);
			$hash{$name} = 0;
		} else {
			$hash{$name} += length($_);
		}
	}
	return \%hash;
}

##
# Run bowtie-inspect -s (summary mode) on the given index and parse its output
# to get the lengths of each of the reference strings.
#
sub bowtie1_index_quick_length_hash {
	my ($base, $btinspect_exe, $slow) = @_;
	$btinspect_exe = find_tool("bowtie-inspect", $btinspect_exe);
	my $cmd = "$btinspect_exe -s $base";
	open(CMD, "$cmd |") || croak("Could not open command '$cmd'");
	my %hash = ();
	while(<CMD>) {
		chomp;
		next unless /^Sequence/;
		my @s = split(/\t/, $_, -1);
		my ($tmp, $fullname, $len) = @s;
		defined($len) || die;
		my $name = $fullname;
		$name =~ s/\s.*//;
		$hash{$name} = $len;
	}
	return \%hash;
}

1;
