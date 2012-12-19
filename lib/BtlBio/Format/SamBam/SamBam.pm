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
# BtlBio::Format::SamBam::SamBam
#
# Input and output streams to/from .sam and .bam files.
#
# Author: Ben Langmead
#   Date: 5/16/2012
#

package BtlBio::Format::SamBam::SamBam;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(sam_to_bam_filename
             sam_to_sorted_bam_filename
             sam_input_stream
             command_input_stream
             command_multi_input_stream
             pileup_input_stream
             mpileup_input_stream
             mpileup_multi_input_stream
             output_stream);

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../../..";
use BtlBio::Util::Tool;
use BtlBio::Util::File;
use BtlBio::Util::Temp;

##
# Return path to BWA executable.
#
sub find_exe() {
	return find_tool("samtools");
}

##
# Given a filename for a .sam file, make a corresponding .bam file and return
# it.
#
sub sam_to_bam_filename {
	my ($fn, $args_view, $ofn, $samtools_exe, $verbose) = @_;
	$args_view = $args_view || "";
	$ofn = $ofn || temp_name(".bam");
	$samtools_exe = $samtools_exe || find_exe();
	my $cmd = "$samtools_exe view $args_view -bS $fn > $ofn";
	print STDERR "samtools view command: $cmd\n" if $verbose;
	system($cmd);
	$? == 0 || die "Bad exitlevel $? from '$cmd'";
	return $ofn;
}

##
# Given a filename for a .sam file, make a corresponding sorted .bam file and
# return it.
#
sub sam_to_sorted_bam_filename {
	my ($fn, $args_view, $ofn, $samtools_exe, $verbose) = @_;
	$args_view = $args_view || "";
	$ofn = $ofn || temp_name(".sorted.bam");
	$ofn =~ s/\.bam$//;
	$samtools_exe = $samtools_exe || find_exe();
	my $tmp = temp_name(".bam");
	my $cmd = "$samtools_exe view $args_view -bS $fn > $tmp";
	print STDERR "samtools view command: $cmd\n" if $verbose;
	system($cmd);
	$? == 0 || die "Bad exitlevel $? from '$cmd'";
	$cmd = "$samtools_exe sort $tmp $ofn";
	print STDERR "samtools sort command: $cmd\n" if $verbose;
	system($cmd);
	$? == 0 || die "Bad exitlevel $? from '$cmd'";
	return "$ofn.bam";
}

##
# Given a (possibly compressed) .sam or .bam file, open a stream that reads
# the file as SAM and return it.
#
sub sam_input_stream {
	my ($fn, $samtools_exe, $format) = @_;
	$samtools_exe = $samtools_exe || find_exe();
	-f $fn || die "File '$fn' given to sam_input_stream but doesn't exist";
	my $is_sam = (defined($format) && $format eq "sam");
	$is_sam = ($fn =~ /\.sam(\.gz|\.bz2)?$/) unless defined($format);
	my $is_bam = (defined($format) && $format eq "bam");
	$is_bam = ($fn =~ /\.bam$/) unless defined($format);
	if($is_sam) {
		return (openex($fn), "sam");
	} elsif($is_bam) {
		my $fh;
		open($fh, "$samtools_exe view -h $fn |") ||
			die "Could not open samtools pipe for '$fn'";
		return ($fh, "bam");
	} else {
		die "Not sure how to open SAM stream from file '$fn'";
	}
	return (undef, undef);
}

##
# Given a (possibly compressed) .sam or .bam file, open a stream that reads
# the file using the given samtools command and return the open filehandle.
#
sub command_input_stream($$$$) {
	my ($fn, $stcmd, $args, $samtools_exe) = @_;
	$samtools_exe = $samtools_exe || find_exe();
	-f $fn || die "File '$fn' given to pileup_input_stream but doesn't exist";
	my $inpipe = "";
	my $infile = $fn;
	my $inflags = "";
	if($fn =~ /\.sam\.gz$/) {
		$inpipe = "gzip -dc $fn |";
		$infile = "-";
		$inflags = "-S";
	}
	elsif($fn =~ /\.sam\.bz2$/) {
		$inpipe = "bzip2 -dc $fn |";
		$infile = "-";
		$inflags = "-S";
	}
	elsif($fn =~ /\.sam$/) {
		$inflags = "-S";
	}
	else {
		$fn =~ /\.bam$/ || die;
	}
	my $fh;
	my $cmd = "$inpipe $samtools_exe $stcmd $inflags $args $infile";
	open($fh, "$cmd |") || die "Could not open pipe '$cmd' for reading";
	return $fh;
}

##
# Given a set of .sam or .bam files, open a stream that reads the files using
# the given samtools command and return the open filehandle.
#
sub command_multi_input_stream($$$$) {
	my ($fns, $stcmd, $args, $samtools_exe) = @_;
	$samtools_exe = $samtools_exe || find_exe();
	my ($nsam, $nbam) = (0, 0);
	for my $fn (@$fns) {
		if($fn =~ /\.sam\.gz$/ || $fn =~ /\.sam\.bz2$/) {
			die "Can't have compressed SAM inputs for a multi-input samtools command";
		}
		if($fn =~ /\.sam$/) {
			$nsam++;
		}
		if($fn =~ /\.bam$/) {
			$nbam++;
		}
		-f $fn || die "File '$fn' given to pileup_input_stream but doesn't exist";
	}
	$nsam == 0 || $nbam == 0 || die "Can't have both SAM and BAM inputs to command_multi_input_stream";
	my $fnstr = join(" ", @$fns);
	my $inpipe = ""; # unused
	my $infile = $fnstr;
	my $inflags = "";
	$inflags = "-S" if $nsam > 0;
	my $fh;
	my $cmd = "$inpipe $samtools_exe $stcmd $inflags $args $infile";
	open($fh, "$cmd |") || die "Could not open pipe '$cmd' for reading";
	return $fh;
}

##
# Given a (possibly compressed) .sam or .bam file, open a stream that reads
# the file as a pileup and return the open filehandle.
#
sub pileup_input_stream($$$) {
	return command_input_stream($_[0], "pileup", $_[1], $_[2]);
}

##
# Given a (possibly compressed) .sam or .bam file, open a stream that reads
# the file as a mpileup and return the open filehandle.
#
sub mpileup_input_stream {
	return command_input_stream($_[0], "mpileup", $_[1], $_[2]);
}

##
# Given a set of .sam or .bam files, open a stream that reads all of them in
# tandem as an mpileup and return the open filehandle.
#
sub mpileup_multi_input_stream($$$) {
	return command_multi_input_stream($_[0], "mpileup", $_[1], $_[2]);
}

##
# Create a new output stream for either SAM or BAM output.  For SAM output,
# the output stream is a file.  The BAM output, the output stream is piped
# through samtools view -b then redirected to a file.
#
sub output_stream($$$) {
	my ($fn, $type, $samtools_exe) = @_;
	$samtools_exe = $samtools_exe || find_exe();
	my $fh;
	if($type eq "sam") {
		open($fh, ">$fn") || die "Could not open '$fn' for writing";
	} else {
		$type eq "bam" || die;
		my $pipe = "| $samtools_exe view -bS - > $fn";
		open($fh, $pipe) || die "Could not open pipe '$pipe' for writing";
	}
	return $fh;
}

1;
