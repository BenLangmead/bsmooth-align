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
# bowtie2_bs_build.pl
#
# Build a pair of Bowtie 2 (Bowtie 2, not Bowtie 1) indexes corresponding to the
# in-silico bisulfite-treated strands of the reference.
#

use strict;
use warnings;
use Getopt::Long;
use File::Path;
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use BtlBio::Index::Bisulfite::BsIndex;
use BtlBio::Util::File;
use BtlBio::Util::Tool;

my $keep_fasta = 0;
my $skip_index = 0;
my $temp_dir = $ENV{TMPDIR} || "/tmp";
my $output_dir = ".";
my $bowtie_build_exe = "bowtie2-build";
my $idx_name = "";

my $usage = qq!
Usage:
  bowtie2_bs_build.pl --name=<string> [options*] <fasta-files>

Build a pair of Bowtie 2 (Bowtie 2, not Bowtie 1) indexes corresponding to the
in-silico bisulfite-treated strands of the reference.

Required arguments:
  --name=<string>         Name prefix for index files

Options (defaults in parentheses):
  --bowtie2-build=<path>  Path to bowtie2-build executable (from PATH)
  --keep-fasta            Keep bisulfite-converted FASTA files (delete)
  --skip-index            Generate \& keep FASTA; don't build index (off)
  --temp=<path>           Set path used for temporary FASTA files (TMPDIR)
  --output=<path>         Directory to store index, and/or fasta if
                          --skip-index or --keep-fasta are specified (.)
  --help/--usage          Print this message and quit

Note: You will need to have enough free space on the partitions 
!;

sub dieusage($$) {
	my ($msg, $level) = @_;
	print STDERR "Error $level:\n";
	print STDERR "$msg\n";
	print STDERR "---\n";
	print STDERR "$usage\n";
	exit $level;
}

GetOptions (
	"name=s"          => \$idx_name,
	"bowtie2-build=s" => \$bowtie_build_exe,
	"keep-fasta"      => \$keep_fasta,
	"skip-index"      => sub { $keep_fasta = 1; $skip_index = 1; },
	"temp=s"          => \$temp_dir,
	"output=s"        => \$output_dir,
	"help|usage"      => sub { print $usage; exit 0; }
) || dieusage("Bad option", 0);

if($idx_name eq "") {
	print "$usage\n---\nError: --name not specified\n"; exit 1;
}
if(scalar(@ARGV) == 0) {
	print "$usage\n---\nError: no <fasta-files> specified\n"; exit 1;
}

$bowtie_build_exe = find_tool("bowtie2-build", $bowtie_build_exe);
$idx_name ne "" || dieusage("Must specify --name", 20);

print STDERR "bowtie2-build executable: $bowtie_build_exe\n";
print STDERR "Index name: $idx_name\n";
print STDERR "Keep fasta files: $keep_fasta\n";
print STDERR "Skip building index: $skip_index\n";
print STDERR "Temporary directory: $temp_dir\n";
print STDERR "Output directory: $output_dir\n";
print STDERR "FASTA input files:\n";
my @fa = @ARGV;
for (@fa) { print STDERR "  $_\n"; }
scalar(@fa) > 0 ||
	die "Must specify one or more FASTA input files as arguments!";

# Convert fasta files
mkpath($temp_dir);
my $fasta_dir = $temp_dir;
$fasta_dir = $output_dir if $keep_fasta || $skip_index;
mkpath($fasta_dir);
mkpath($output_dir);

for(my $cr = 0; $cr < 2; $cr++) {
	my $strand_name = ($cr ? "Crick" : "Watson");
	my $wa = ($cr ? 0 : 1);
	print STDERR "-- Converting FASTA files: $strand_name --\n";
	print STDERR "Output dir: $fasta_dir\n";
	my @conv_fa = ();
	my $idx_fn = "$output_dir/$idx_name";
	my $desig = (lc $strand_name);
	$idx_fn .= ".$desig";
	for my $f (@fa) {
		my $base = fileparse($f);
		print STDERR "Converting $base...\n";
		my $ofn = "$fasta_dir/$base.$desig.fa";
		my $ifh = openex($f);
		open(my $ofh, ">$ofn") || die "Could not open '$ofn' for writing";
		bsc_convert_stream($ifh, $ofh, $wa);
		close($ifh); close($ofh);
		push @conv_fa, $ofn;
	}
	
	# Build index
	print STDERR "-- Building index: $strand_name --\n";
	my $fnstr = join(",", @conv_fa);
	my $cmd = "$bowtie_build_exe $fnstr $idx_fn";
	print STDERR "$cmd\n";
	system($cmd);
	$? == 0 || die "bowtie2-build command failed (see above)";
	
	# Purge temporary FASTA files
	unless($keep_fasta) {
		print STDERR "-- Deleting tempoarary FASTA: $strand_name --\n";
		for my $f (@conv_fa) {
			print STDERR "Deleting $f ...\n";
			unlink($f);
		}
	}
}
