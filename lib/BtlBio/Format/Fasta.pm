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
# BtlBio::Format::Fasta
#
# Parsing FASTA files.
#
# Author: Ben Langmead
# Email: langmea@cs.jhu.edu
#

package BtlBio::Format::Fasta;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(find_fasta_handle_to_hash
             find_fasta_to_hash
             read_fasta_handle_to_hash
             read_fasta_lens_handle_to_hash
             read_fasta_lens_to_hash
             read_fasta_to_hash);

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../..";
use BtlBio::Util::File;

##
# Read a Fasta file and store name/sequence key/value pairs in the given hash.
#
# Args:
#
# 1. Fasta filehandle open for reading
# 2. Hash ref to destination hash
#
# Returns: nothing
#
sub read_fasta_handle_to_hash {
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
# 1. File name for Fasta file to read
# 2. Name of sequence
# 3. Hash ref to destination hash
#
# Returns: nothing
#
sub find_fasta_handle_to_hash {
	my ($fh, $nm, $h, $conf) = @_;
	my $name = "";
	my $collect = 0;
	while(readline $fh) {
		chomp;
		if(substr($_, 0, 1) eq ">") {
			return 1 if $collect;
			$name = substr($_, 1);
			if($conf->{truncate_names}) {
				$name =~ s/\s.*//;
			}
			$collect = 1 if $name eq $nm;
			$h->{$name} = "" if $collect;
		} else {
			$name ne "" || die "No name prior to sequence line:\n$_";
			$h->{$name} .= $_ if $collect;
		}
	}
	return $collect; # didn't find it
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
sub find_fasta_to_hash {
	my ($fn, $nm, $h, $conf) = @_;
	my $fh = openex($fn);
	my $ret = find_fasta_handle_to_hash($fh, $nm, $h, $conf);
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
sub read_fasta_to_hash {
	my ($fn, $h, $conf) = @_;
	my $fh = openex($fn);
	read_fasta_handle_to_hash($fh, $h, $conf);
	close($fh);
}

##
# Read a Fasta file and store name/sequence-length key/value pairs in the
# given hash.
#
# Args:
#
# 1. Fasta filehandle open for reading
# 2. Hash ref to destination hash
#
# Returns: nothing
#
sub read_fasta_lens_handle_to_hash {
	my ($fh, $h, $conf) = @_;
	my $name = "";
	while(readline $fh) {
		chomp;
		if(substr($_, 0, 1) eq ">") {
			if($conf->{remove_empties} && $name ne "" && $h->{$name} == 0) {
				delete $h->{$name};
			}
			$name = substr($_, 1);
			if($conf->{truncate_names}) {
				$name =~ s/\s.*//;
			}
			if(defined($h->{$name})) {
				die "Already saw FASTA record with name '$name'";
			}
			$h->{$name} = 0;
		} else {
			$name ne "" || die "No name prior to sequence line:\n$_";
			$h->{$name} += length($_);
		}
	}
	if($conf->{remove_empties} && $name ne "" && $h->{$name} == 0) {
		delete $h->{$name};
	}
}

##
# Read a Fasta file and store name/sequence-length key/value pairs in the
# given hash.
#
# Args:
#
# 1. File name for Fasta file to read
# 2. Hash ref to destination hash
#
# Returns: nothing
#
sub read_fasta_lens_to_hash {
	my ($fn, $h, $conf) = @_;
	my $fh = openex($fn);
	read_fasta_lens_handle_to_hash($fh, $h, $conf);
	close($fh);
}

sub _test_1() {
	require BtlBio::Util::Temp; BtlBio::Util::Temp->import();
	print STDERR "Testing read_fasta_lens_to_hash 1 ... ";
	my ($ofh, $ofn) = temp_filehandle(".fa");
	print $ofh ">ref1 blah\n";
	print $ofh "ACGATCGATCGTACGTAG\n\n";
	print $ofh "A\n";
	close($ofh);
	my %h = ();
	read_fasta_lens_to_hash($ofn, \%h, {});
	defined($h{"ref1 blah"}) || die;
	$h{"ref1 blah"} == 19 || die;
	print STDERR "PASSED\n";
	unlink("$ofn");
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
	my %h = ();
	find_fasta_to_hash($ofn, "ref1", \%h, { truncate_names => 1 });
	scalar(keys %h) == 1 || die;
	defined($h{"ref1"}) || die;
	$h{"ref1"} eq "ACGATCGATCGTACGTAGA" || die;
	%h = ();
	find_fasta_to_hash($ofn, "ref2", \%h, { truncate_names => 1 });
	defined($h{"ref2"}) || die;
	scalar(keys %h) == 1 || die;
	$h{"ref2"} eq "AAC" || die;
	print STDERR "PASSED\n";
	unlink("$ofn");
}

sub _test() {
	_test_1();
	_test_2();
}

# Call _test() if the module is being run directly.
_test() unless caller();

1;
