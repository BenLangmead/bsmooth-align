#!/usr/bin/env perl

#
# Copyright 2012, Ben Langmead <blangmea@jhsph.edu>
#
# This file is part of The BSmooth Alignment Pipeline.
#
# The BSmooth Alignment Pipeline is free software: you can redistribute it
# and/or modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
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
# bs_sam_extract.pl
#
# Extract evidence from BS-SAM files.  A BS-SAM file is output by a tool like
# Merman that aligns in a fully bisulfite-aware fashion.  Tools that aren't
# bisulfite aware can still be harnessed using the bswc_* scripts included
# here, and the SAM files output by those tools are referred to as BSWC_SAM
# files.
#

use strict;
use warnings;
use Getopt::Long;
use File::Path;
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use BtlBio::Alignment::Alignment;
use BtlBio::Alignment::SAM;
use BtlBio::Align::Bisulfite::BsAlign;
use BtlBio::Align::Bisulfite::BsEvidence;
use BtlBio::Util::File;
use BtlBio::Util::Tool;
use BtlBio::Util::Temp;
use BtlBio::Format::SamBam::SamBam;
use POSIX;
use Carp;

my $samtools_exe = "samtools";  # path to samtools bin; only for --keep-bam
my $temp_dir = $ENV{TMPDIR};    # temp dir for intermediates; reads, alns
my $bscpg = 1;                  # 1 -> do CpG->YpG conversion
my $bsc = 0;                    # 1 -> do C->Y conversion
my $ev_out_dir = ".";           # directory to write evidence
my $echo_sam = 0;               # print out SAM as it comes in

my $usage = qq!

Given .sam files output by a bs_*.pl wrapper script, extract and output
methylation evidence.

Usage:
  bs_sam_extract.pl [options*] <sam/bam>

  <sam/bam>
    A list of .sam and/or .bam files containing alignments generated by a
    bs_*.pl script (NOT bswc_*.pl script - use bscw_sam_extract.pl for that).

Options (defaults in parentheses):
  --echo-sam             Print SAM alignments to STDOUT as they're read in
  --bscpg                CpG Cs become Ys, other Cs become Ts; extract evidence
                         from CpG Cs (on)
   OR:
  --bsc                  All Cs become Ys; extract evidence from all Cs (off)
   *** --bsc/--bscpg must be set to match the exact ***
   *** settings used when creating the alignments   ***
!;

##
# Print and run given command.  Die if it returns non-zero.
#
sub run($) {
	my $cmd = shift;
	print STDERR "Running '$cmd'...\n";
	system($cmd);
	$? == 0 || croak("Command '$cmd' failed with exitlevel $?\n");
}

##
# Die and print the usage message.
#
sub dieusage($$) {
	my ($msg, $level) = @_;
	print STDERR "$usage\n--\n";
	print STDERR "Error $level:\n";
	print STDERR "$msg\n";
	exit $level;
}

GetOptions (
	"samtools=s"      => \$samtools_exe,
	"temp=s"          => \$temp_dir,
	"bscpg"           => sub { $bscpg = 1; $bsc = 0; },
	"bsc"             => sub { $bscpg = 0; $bsc = 1; },
	"echo-sam"        => \$echo_sam,
	"output-dir=s"    => \$ev_out_dir
) || dieusage("Bad option", 0);

mkpath($ev_out_dir);

# Parse out arguments from among the double-dashes
my @sam_fns = ();
for my $a (@ARGV) {
	next if $a eq "--";
	push @sam_fns, $a;
}

# Now we have filehandles from which to read the two streams of alignments,
# one with alignments to BSW/BSWR, one with alignments to BSC/BSCR
my $nheads = 0;
my $nrecs = 0;
my $nrec_ival = 1000;
my %ev_out_fhs = ();
##
# Given a piece of methylation evidence, simply print it to STDOUT.
#
my $nev = 0;
open(MARK, ">$ev_out_dir/.".($bscpg ? "cpg" : "c").".ev") || die;
print MARK "\n"; close(MARK);
sub evidenceSink($) {
	my $ev = shift;
	$nev++;
	my $tname = $ev->{_tname};
	$tname =~ s/\s.*//;
	$tname = sanitize_filename($tname);
	my $outfn = "$ev_out_dir/$tname.ev.tsv";
	if(!defined($ev_out_fhs{$outfn})) {
		open($ev_out_fhs{$outfn}, ">$outfn") ||
			die "Could not open '$outfn' for writing";
		my $outbinfn = "$ev_out_dir/.$tname.ev.tsv.bin";
		open(TMP, ">$outbinfn") ||
			die "Could not open '$outbinfn' for writing";
		print TMP "$ev->{_tname}\n"; # original, not sanitized
		close(TMP);
	}
	print {$ev_out_fhs{$outfn}} $ev->to_record."\n";
}
print STDERR "Processing SAM output from aligners...\n";
for my $fn (@sam_fns) {
	my ($al_fh, $format) = sam_input_stream($fn, $samtools_exe);
	while(1) {
		# Output is in SAM format
		my $line = readline $al_fh;
		last unless defined($line);
		print $line if $echo_sam;
		# Pass SAM line through to .sam/.bam files if --keep-sam/--keep-bam are
		# specified
		my $firstc = substr($line, 0, 1);
		if($firstc eq "\@") {
			# Header line
			$nheads++;
			next; # skip
		}
		# Create Alignment objects for both
		my $al = parse_sam_record($line);
		# Are there alignments on both strands?
		bs_analyze_alignment($al, \&evidenceSink) if $al->aligned;
		if((++$nrecs % $nrec_ival) == 0) {
			print STDERR "  Processed $nrecs BS-SAM records...\n";
		}
	}
	close($al_fh);
}
print STDERR "Summary:\n";
print STDERR "  SAM header lines parsed: $nheads\n";
print STDERR "  SAM record pairs parsed: $nrecs\n";
print STDERR "  Pieces of evidence extracted: $nev\n";

for my $ofh (%ev_out_fhs) { close($ofh); }
