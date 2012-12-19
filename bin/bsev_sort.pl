#!/usr/bin/env perl

#
# Copyright 2012, Ben Langmead <langmea@cs.jhu.edu>
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
# Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
#
# This file is part of BSmooth.
#
# Takes a directory of BSC evidence files and outputs a sorted & compressed
# version.
# 

use strict;
use warnings;
use File::Path;
use File::Basename;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use BtlBio::Util::File;
use BtlBio::Util::Temp;
use BtlBio::Align::Bisulfite::BsEvidence;
use BtlBio::Align::Bisulfite::BsFilter;
use BtlBio::Align::Bisulfite::BsSummary;

my $output_dir = "";    # place sorted output here
my $ev_dir = "";        # evidence directories
my $temp_dir = "";      # do intermediate things here
my $threads = 0;        # # threads
my $compress = "";      # type of compression for final output
my @name_map_fn = ();   # files with name mappings

my $usage = qq!
Usage:
  bsev_sort.pl --ev=<ev-path> --out=<out-path> [options*]

Take a directory of read-level measurements output by one of the BSmooth
alignment scripts and sort it with respect to the reference genome.

Required arguments:
  
  <ev-path>
    Path to directory containing read-level measurements
  <out-path>
    Path to directory where sorted output should be placed

Options (defaults in parentheses):
 
  --num-threads=<int>      Use <int> threads/CPUs simultaneously (1)
  --name-map=<path>        Many-to-one map from raw reference names to bin
                           names.  Output will be binned by bin name.
  --gzip                   gzip compress read-level measurement output (off)
  --bzip2                  bzip2 compress read-level measurement output (off)
  --temp=<path>            Set path used for temporary FASTA files (TMPDIR)
  --help/--usage           Print this message and quit

!;

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

$SIG{INT} = sub { confess(); };

GetOptions (
	"evidence-directory=s" => \$ev_dir,
	"output-directory=s" => \$output_dir,
	"temp-directory|temporary-directory=s" => \$temp_dir,
	"num-threads|nthreads=i" => \$threads,
	"compress=s" => \$compress,
	"name-map=s" => \@name_map_fn,
	"gzip"  => sub { $compress = "gzip"; },
	"bzip2" => sub { $compress = "bzip2"; },
	"help|usage" => sub { print $usage; exit 0; }
) || dieusage("Bad option", 0);

$ev_dir     ne "" || dieusage("Must specify --evidence=<path>", 1);
$output_dir ne "" || dieusage("Must specify --output=<path>", 2);

print STDERR "---- Settings ----\n";
print STDERR "Evidence directory: $ev_dir\n";
print STDERR "Compression type: $compress\n" if $compress ne "";
print STDERR "# threads: $threads\n";
print STDERR "Name maps:\n";
for my $f (@name_map_fn) { print STDERR "  $f\n"; }
print STDERR "Output directory: $output_dir\n";
print STDERR "------------------\n";

my $name_map_str = "";
for my $f (@name_map_fn) { $name_map_str .= " --name-map $f "; }

# Establish whether the evidence is C or CpG
my $cpg = -f "$ev_dir/.cpg.ev";

print STDERR "Constructing sort command\n";
# --bin-key=1:  bin by reference sequence
# --sort-key=2: sort by reference offset
# --numeric=2:  reference offset is numeric
$temp_dir = " --intermediate=".temp_name(undef, $temp_dir) if $temp_dir ne "";
$threads = $threads > 0 ? " --num-threads=$threads" : "";
my $compress_str = " ";
$compress_str .= "--compress=$compress" if $compress ne "";
my $input_str = " ";
$input_str .= "--input=$ev_dir ";
my $cmd = "perl $Bin/disk_sort.pl".
          " --numeric=2".
          " --bin-key=1".
          " --sort-key=1".
          " --sort-key=2".
          " --output=$output_dir".
          " --suffix=.tsv".
          $input_str.
          $temp_dir.
          $compress_str.
          $name_map_str.
          $threads;

print STDERR "Command:\n$cmd\n";
system($cmd);
$? == 0 || die "Bad exitlevel from sort command: $?";

my $mark_fn = $cpg ? ".cpg.sorted.ev" : ".c.sorted.ev";
open(MARK, ">$output_dir/$mark_fn");
print MARK "\n";
close(MARK);
