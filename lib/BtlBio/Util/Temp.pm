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
# BtlBio::Util::Temp
#
# Some functions helpful for making temporary files, including temporary files
# pre-filled with DNA sequences and sequence reads.
#
# Author: Ben Langmead
# Email: langmea@cs.jhu.edu
#

package BtlBio::Util::Temp;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(temp_fasta_filename temp_fastq_filename temp_name temp_filehandle);

use strict;
use warnings;

# This function generates random strings of a given length
our $cnt = 1;
sub _rand_str($) {
	my $len = shift;
	if($cnt == 1) {
		srand(time ^ ($$ + ($$ << 15)));
	}
	$cnt++;
	my @chars = ('a'..'z', 'A'..'Z', '0'..'9', '_');
	my $str = "";
	foreach (1..$len) {
		$str .= $chars[int(rand(scalar(@chars)))];
	}
	return $str;
}

##
# Takes arrays of sequences and perhaps a parallel array of names, and writes
# everything to a temporary FASTA file, and returns its name.
#
sub temp_fasta_filename {
	my ($seqs, $names) = @_;
	defined($seqs) || die "sequences argument not defined";
	my ($fh, $tmpfn) = temp_filehandle(".fa");
	for(0..scalar(@$seqs)-1) {
		my $name = "noname$_";
		if(defined($names) && $_ < scalar(@$names)) {
			$name = $names->[$_];
		}
		print {$fh} ">$name\n".$seqs->[$_]."\n";
	}
	close($fh);
	return $tmpfn;
}

##
# Takes arrays of sequences and perhaps a parallel array of names, and writes
# everything to a temporary FASTA file, and returns its name.
#
sub temp_fastq_filename {
	my ($seqs, $quals, $names) = @_;
	defined($seqs) || die "sequences argument not defined";
	my ($fh, $tmpfn) = temp_filehandle(".fq");
	for(0..scalar(@$seqs)-1) {
		my $name = "noname$_";
		my $qual = undef;
		if(defined($names) && $_ < scalar(@$names)) {
			$name = $names->[$_];
		}
		if(defined($quals) && $_ < scalar(@$quals)) {
			$qual = $quals->[$_];
		}
		$qual = $qual || "I" x length($seqs->[$_]);
		print {$fh} "\@$name\n".$seqs->[$_]."\n+\n$qual\n";
	}
	close($fh);
	return $tmpfn;
}

##
# Return the name of a temporary file.
#
sub temp_name {
	my ($suffix, $dir) = @_;
	$suffix = $suffix || ".tmp";
	my $temp_dir = $dir || $ENV{TMPDIR} || "/tmp/";
	$temp_dir .= "/" unless substr($temp_dir, -1) eq "/";
	return "$temp_dir"._rand_str(10)."$suffix";
}

##
# Return a filehandle open for writing, along with the temporary file's name.
#
sub temp_filehandle {
	my ($suffix) = @_;
	my $fn = temp_name($suffix);
	my $fh;
	open($fh, ">$fn") || croak("Could not open temporary file '$fn' for writing");
	return ($fh, $fn);
}

1;
