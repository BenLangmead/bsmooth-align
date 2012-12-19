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
# BtlBio::Format::FastaIndexed
#
# Routines for parsing indexed FASTA files, as well as for saving and loading
# record indexes.  A record index contins a map from reference sequence names
# to their (a) lengths in characters, and (b) offset into the file, such that
# one could scan there with seek and arrive at the beginning of the name line.
#
# Author: Ben Langmead
# Email: langmea@cs.jhu.edu
#

package BtlBio::Format::FastaIndexed;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(fasta_idx_save_recidx
             fasta_idx_load_recidx
             fasta_idx_index_handle
             fasta_idx_index
             fasta_idx_get_index
             read_fasta_idx_handle_to_hash
             find_fasta_idx_handle_to_hash
             find_fasta_idx_to_hash
             read_fasta_idx_to_hash);

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../..";
use BtlBio::Util::File;

##
# Save a record index to given filehandle.
#
sub fasta_idx_save_recidx {
	my ($rec_idx, $fh) = @_;
	for my $k (keys %$rec_idx) {
		print {$fh} "$k\t".$rec_idx->{$k}{fn}.
		              "\t".$rec_idx->{$k}{off}.
		              "\t".$rec_idx->{$k}{len}.
		              "\n";
	}
	return $rec_idx;
}

##
# Load a record index from given filehandle.
#
sub fasta_idx_load_recidx {
	my ($rec_idx, $fh) = @_;
	while(readline $fh) {
		chomp;
		my @s = split(/\t/, $_, -1);
		scalar(@s) == 4 || die "Expected 4 tokens, got:\n$_";
		$rec_idx->{$s[0]}{fn}  = $s[1];
		$rec_idx->{$s[0]}{off} = $s[2];
		$rec_idx->{$s[0]}{len} = $s[3];
	}
	return $rec_idx;
}

##
# Given a hash and a filehandle, read the filehandle as a fasta file and
# install the record index into the hash.
#
sub fasta_idx_index_handle {
	my ($rec_idx, $fh, $fn) = @_;
	my $len = 0;
	my $name = "";
	my $off = undef;
	while(1) {
		my $cur_off = tell $fh;
		my $line = readline $fh;
		last unless defined($line);
		chomp($line);
		if(substr($line, 0, 1) eq ">") {
			# Name line
			if($name ne "") {
				defined($off) || die;
				$name =~ s/\s.*//;
				$rec_idx->{$name}{fn} = $fn;
				$rec_idx->{$name}{off} = $off;
				$rec_idx->{$name}{len} = $len;
			}
			$name = substr($line, 1);
			$off = $cur_off;
			$len = 0;
		} else {
			# Sequence line
			$len += length($line);
		}
	}
	if($name ne "") {
		defined($off) || die;
		$name =~ s/\s.*//;
		$rec_idx->{$name}{fn} = $fn;
		$rec_idx->{$name}{off} = $off;
		$rec_idx->{$name}{len} = $len;
	}
	return $rec_idx;
}

##
# Given a hash and a filename, make a record index, install it in the hash,
# and possibly save it to a file as well.
#
sub fasta_idx_index {
	my ($rec_idx, $fn, $save) = @_;
	my $fh = openex($fn);
	fasta_idx_index_handle($rec_idx, $fh, $fn);
	close($fh);
	if($save) {
		open(my $ofh, ">$fn.recidx") ||
			die "Could not open '$fn.recidx' for writing record index";
		fasta_idx_save_recidx($rec_idx, $ofh);
		close($ofh);
	}
	return $rec_idx;
}

##
# Get a record index for the fasta file at the given path.  First, we look to
# see if there's an ABC.recidx file, where ABC = the fasta file's name.  If
# that doesn't work, we 
#
sub fasta_idx_get_index {
	my ($fn, $save) = @_;
	my %rec_idx = ();
	if(-f "$fn.recidx") {
		my $fh = openex("$fn.recidx");
		fasta_idx_load_recidx(\%rec_idx, $fh);
		close($fh);
	} else {
		fasta_idx_index(\%rec_idx, $fn, $save);
	}
	return \%rec_idx;
}

##
# Read a Fasta file and store name/sequence key/value pairs in the given hash.
# Indexed fasta treated exactly the same as normal fasta here.
#
# Args:
#
# 1. Fasta filehandle open for reading
# 2. Hash ref to destination hash
#
# Returns: nothing
#
sub read_fasta_idx_handle_to_hash {
	my ($fh, $h, $conf) = @_;
	my $name = "";
	while(readline $fh) {
		chomp;
		if(substr($_, 0, 1) eq ">") {
			if($conf->{remove_empties} && $name ne "" && $h->{$name} eq "") {
				delete $h->{$name};
			}
			$name = substr($_, 1);
			if($conf->{truncate_names}) {
				$name =~ s/\s.*//;
			}
			if(defined($h->{$name})) {
				die "Already saw FASTA record with name '$name'";
			}
			$h->{$name} = "";
		} else {
			$name ne "" || die "No name prior to sequence line:\n$_";
			$h->{$name} .= $_;
		}
	}
	if($conf->{remove_empties} && $name ne "" && $h->{$name} eq "") {
		delete $h->{$name};
	}
}

##
# Search through a Fasta file for a particular reference sequence and store
# a name/sequence key/value pair for it in the given hash.
#
# Args:
#
# 1. Seekable filehandle for Fasta file to read
# 2. Name of sequence
# 3. Hash ref to destination hash
#
# Returns: nothing
#
sub find_fasta_idx_handle_to_hash {
	my ($fh, $rec_idx, $nm, $h, $conf) = @_;
	defined($rec_idx->{$nm}) || return 0;
	my $off = $rec_idx->{$nm}{off};
	seek $fh, $off, 0;
	my $first = 1;
	my $name = "";
	while(readline $fh) {
		chomp;
		if(substr($_, 0, 1) eq ">") {
			return 1 unless $first;
			$name = substr($_, 1);
			$name =~ s/\s.*// if $conf->{truncate_names};
			$h->{$name} = "";
			$first = 0;
		} else {
			$first && die;
			$h->{$name} .= $_;
		}
		
	}
	return 1;
}

##
# Search through a Fasta file for a particular reference sequence and store
# a name/sequence key/value pair for it in the given hash.
#
# Args:
#
# 1. File name for Fasta file to read
# 2. Name of sequence
# 3. Hash ref to destination hash
#
# Returns: nothing
#
sub find_fasta_idx_to_hash {
	my ($fn, $rec_idx, $nm, $h, $conf) = @_;
	my $fh = openex($fn);
	my $ret = find_fasta_idx_handle_to_hash($fh, $rec_idx, $nm, $h, $conf);
	close($fh);
	return $ret;
}

##
# Read a Fasta file and store name/sequence key/value pairs in the given hash.
#
# Args:
#
# 1. File name for Fasta file to read
# 2. Hash ref to destination hash
#
# Returns: nothing
#
sub read_fasta_idx_to_hash {
	my ($fn, $h, $conf) = @_;
	my $fh = openex($fn);
	read_fasta_idx_handle_to_hash($fh, $h, $conf);
	close($fh);
}

sub _test_1() {
	require BtlBio::Util::Temp; BtlBio::Util::Temp->import();
	print STDERR "Testing read_fasta_idx_to_hash 1 ... ";
	my ($ofh, $ofn) = temp_filehandle(".fa");
	print $ofh ">ref1 blah\n";
	print $ofh "ACGATCGATCGTACGTAG\n\n";
	print $ofh "A\n";
	print $ofh ">ref2 blah\n";
	print $ofh "A\n\n";
	print $ofh "A\n";
	print $ofh "C\n";
	close($ofh);
	my %h = ();
	read_fasta_idx_to_hash($ofn, \%h, {});
	scalar(keys %h) == 2 || die;
	defined($h{"ref1 blah"}) || die;
	defined($h{"ref2 blah"}) || die;
	$h{"ref1 blah"} eq "ACGATCGATCGTACGTAGA" || die;
	$h{"ref2 blah"} eq "AAC" || die;
	print STDERR "PASSED\n";
	unlink("$ofn");
	unlink("$ofn.recidx");
}

sub _test_2() {
	require BtlBio::Util::Temp; BtlBio::Util::Temp->import();
	print STDERR "Testing find_fasta_to_hash 1 ... ";
	my ($ofh, $ofn) = temp_filehandle(".fa");
	print $ofh ">ref1 blah\n";
	print $ofh "ACGATCGATCGTACGTAG\n\n";
	print $ofh "A\n";
	print $ofh ">ref2 blah\n";
	print $ofh "A\n\n";
	print $ofh "A\n";
	print $ofh "C\n";
	close($ofh);
	my $rec_idx = fasta_idx_get_index($ofn, 1);
	scalar(keys %$rec_idx) == 2 || die "Expected 2 keys, got ".scalar(keys %$rec_idx);
	my %h = ();
	find_fasta_idx_to_hash($ofn, $rec_idx, "ref1", \%h, { truncate_names => 1 });
	scalar(keys %h) == 1 || die "Expected 1 key got ".scalar(keys %h);
	defined($h{"ref1"}) || die;
	$h{"ref1"} eq "ACGATCGATCGTACGTAGA" || die;
	%h = ();
	$rec_idx = fasta_idx_get_index($ofn, 1);
	find_fasta_idx_to_hash($ofn, $rec_idx, "ref2", \%h, { truncate_names => 1 });
	defined($h{"ref2"}) || die;
	scalar(keys %h) == 1 ||  die "Expected 1 key got ".scalar(keys %h);
	$h{"ref2"} eq "AAC" || die;
	print STDERR "PASSED\n";
	unlink("$ofn");
	unlink("$ofn.recidx");
}

sub _test() {
	_test_1();
	_test_2();
}

# Call _test() if the module is being run directly.
_test() unless caller();

1;
