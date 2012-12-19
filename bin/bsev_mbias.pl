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
# Takes a directory of BSmooth read-level measurements and generates an M-bias
# table and plot.
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

my $output_dir = ".";    # place sorted output here
my $ev_dir = ();         # evidence directories
my $mapq_min = undef;    # cutoff for repetitive mapq
my $nthreads = 1;        # # threads to use
my $m1 = 1;              # include measurements from mate #1?
my $m2 = 1;              # include measurements from mate #2?
my $unp = 1;             # include measurements from unpaired reads?
my @incl_bins = ();      # include only these bins
my @excl_bins = ();      # exclude these bins
my %inclSet = ();        # set version of incl_bins
my %exclSet = ();        # set version of excl_bins
my $rev_cy = 0;          # reverse cycle?
my $strat_len = 0;       # >0, then stratify by length
my $strat_score = 0;     # >0, then stratify by score
my $readlen_min = undef; # minimum read len to accept
my $readlen_max = undef; # maximum read len to accept

my $usage = qq!
Usage:
  bsev_mbias.pl --evidence=<ev-path> --output=<out-path> [options*]

Take a directory of read-level measurements output by one of the BSmooth
alignment scripts and generate M-bias tables and plots.

Required arguments:
  
  <ev-path>
    Path to directory containing read-level measurements
  <out-path>
    Path to directory where M-bias tables and plots should be placed.

Options (defaults in parentheses):
 
  --num-threads=<int>      Use <int> threads/CPUs simultaneously (1)
  --help/--usage           Print this message and quit
  --reverse-cycle          Invert cycles to be w/r/t end of read

 Stratification:
 
  --strat-length=<int>     Make additional set of M-bias tables stratified by
                           length.  Each stratum contains <int> distinct read
                           lengths.
  --strat-score=<int>      Make additional set of M-bias tables stratified by
                           alignment score.  Each stratum contains <int>
                           distinct scores.

 Filtering:

  --include=<name>         Include chromosome <name>, exclude others not
                           specified with --include
  --exclude=<name>         Exclude chromosome <name>
  
  --mapq-min=<int>         Measurements with MAPQ < <int> not factored in (8)
  --readlen-min=<int>      Filter measurements from reads shorter than <int>
  --readlen-max=<int>      Filter measurements from reads longer than <int>
  --exclude-mate1          Measurements from mate 1 not factored in
  --exclude-mate2          Measurements from mate 2 not factored in
  --exclude-unpaired       Exclude measurements from mates that align in
                           unpaired but not in paired-end fashion

Output:

 Output is a collection of M-bias tables.  Each M-bias table summarizes the
 number of the read-level measurements that are methylated and unmethylated,
 stratified by read position.  This can be a helpful tool for identifying
 artifacts, especially at the extreme ends of the reads.
 
 In an M-bias table, each row is a read position.  The tab-separated columns
 are:
 
 1. Read position (0-based)
 2. # unmethylated read-level measurements
 3. # methylated read-level measurements
 4. # read-level measurements with nucleotide other than C or T

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
	"evidence-directory=s"                => \$ev_dir,
	"output-directory=s"                  => \$output_dir,
	"strat-length=i"                      => \$strat_len,
	"strat-cost=i"                        => \$strat_score,
	"reverse-cycle"                       => \$rev_cy,
	"mapq-min=i"                          => \$mapq_min,
	"readlen-min=i"                       => \$readlen_min,
	"readlen-max=i"                       => \$readlen_max,
	"include-chromosome|inc-chromosome=s" => \@incl_bins,
	"exclude-chromosome|ex-chromosome=s"  => \@excl_bins,
	"nthreads|num-threads=i"              => \$nthreads,
	"exclude-mate1"                       => sub { $m1 = 0; },
	"exclude-mate2"                       => sub { $m2 = 0; },
	"exclude-unpaired"                    => sub { $unp = 0; },
	"help|usage"                          => sub { print $usage; exit 0; }
) || die "Bad option";

print STDERR "---- Settings ----\n";
print STDERR "Evidence directory: $ev_dir\n";
print STDERR "Minimum MAPQ: ".(defined($mapq_min) ? $mapq_min : "(none)")."\n";
print STDERR "Minimum read length: ".(defined($readlen_min) ? $readlen_min : "(none)")."\n";
print STDERR "Maximum read length: ".(defined($readlen_max) ? $readlen_max : "(none)")."\n";
print STDERR "Reverse cycle: ".($rev_cy ? "yes" : "no")."\n";
if(scalar(@incl_bins) > 0) {
	print STDERR "Chromosomes to include:\n";
	for my $chr (@incl_bins) { print STDERR "  $chr\n"; }
}
if(scalar(@excl_bins) > 0) {
	print STDERR "Chromosomes to exclude:\n";
	for my $chr (@excl_bins) { print STDERR "  $chr\n"; }
}
print STDERR "Output directory: $output_dir\n";
print STDERR "------------------\n";

mkpath($output_dir);

# M/U evidence table (not stratified)
my %ev = ();

# M/U evidence table stratified by alignment length and alignment cost
my %evStrat = ();

# M/U evidence table stratified by alignment length
my %evLenStrat = ();

# M/U evidence table stratified by alignment cost
my %evCostStrat = ();

# M/U evidence table stratified by alignment cost, including only full-length
# alignments
my %evCostFullStrat = ();

# By default, stratify length by 10 and cost by 50
my $len_strat = $strat_len;
my $len_color_strat = 10;
my $cost_strat = $strat_score;

# Lowest and highest cycles observed
my $lo_cyc = 999999;
my $hi_cyc = -1;

# Lowest and highest length strata observed
my $lo_len = 999999;
my $hi_len = -1;

# Lowest and highest cost strata observed
my $lo_cost = 999999;
my $hi_cost = -1;

# Non-zero -> plot all strata, instead of just the summaries
my $allPlots = 0;
# Non-zero -> make stratified plots
my $strataPlots = $strat_len > 0 && $strat_score > 0;

my $ylo = 0.0;
my $yhi = 1.0;
my $sampName = "";

my $pdfWidth = 10;
my $pdfHeight= 7;

# How long is a "full-length" alignment?
my $fullLength = -1;

##
# Given a len stratum id, return the low and high bounds for what
# goes into the stratum
#
sub lenRange($) {
	my $len_st = shift;
	return ($len_st*$len_strat, $len_st*$len_strat+$len_strat);
}

##
# Given a stratum id, return the low and high bounds for the associated
# len stratum
#
sub stratToLenRange($) {
	my $strat = shift;
	my ($len_st, $cost_st) = split /_/, $strat;
	return lenRange($len_st);
}

##
# Given a length stratum id, return a string describing the stratum.
#
sub lenStratToStr($) {
	my $st = shift;
	return "".($st*$len_strat).":".($st*$len_strat+$len_strat-1);
}

##
# Given a cost stratum id, return the low and high bounds for what
# goes into the stratum
#
sub costRange($) {
	my $cost_st = shift;
	return ($cost_st*$cost_strat, $cost_st*$cost_strat+$cost_strat);
}

##
# Given a stratum id, return the low and high bounds for the associated
# cost stratum
#
sub stratToCostRange($) {
	my $strat = shift;
	my ($len_st, $cost_st) = split /_/, $strat;
	return costRange($cost_st);
}

##
# Given a cost stratum id, return a string describing the stratum.
#
sub costStratToStr($) {
	my $st = shift;
	return "".($st*$cost_strat).":".($st*$cost_strat+$cost_strat-1);
}

# Set up the set of evidence files to exclude
for (@excl_bins) { $exclSet{$_} = 1; }

# Set up the set of evidence files to include
for (@incl_bins) { $inclSet{$_} = 1; }
my $has_incl = scalar(@incl_bins) > 0;

my $nskip = 0;
my $nskip_rdlen = 0;
my $nskip_mate = 0;
my $nskip_chr = 0;
my $nskip_mapq = 0;

# For each evidence input file
my $nev = 0;
my $nev_ival = 50000;
for my $f (<$ev_dir/*>) {
	my $ev_fh = openex($f);
	print STDERR "Reading evidence from $f...\n";
	# Iterate through lines
	while(readline $ev_fh) {
		chomp;
		my ($ev_tname, $ev_toff, $ev_readid, $ev_al, $ev_watson, $ev_fw,
		    $ev_flags, $ev_qu1, $ev_qu2, $ev_cy, $ev_allen,
		    $ev_score, $ev_mapq) = split(/\t/, $_, -1);
		if((++$nev % $nev_ival) == 0) {
			my $pct       = 100.0 * $nskip       / ($nev || 1);
			my $pct_chr   = 100.0 * $nskip_chr   / ($nev || 1);
			my $pct_mate  = 100.0 * $nskip_mate  / ($nev || 1);
			my $pct_rdlen = 100.0 * $nskip_rdlen / ($nev || 1);
			my $pct_mapq  = 100.0 * $nskip_mapq  / ($nev || 1);
			printf STDERR "Parsed $nev read-level measurements; skipped: ".
			              "($nskip_chr, %0.2f%% for chromosome, $nskip_mate, ".
			              "%0.2f%% for mate, $nskip_rdlen, %0.2f%% for read ".
			              "length, $nskip_mapq, %0.2f%% for MAPQ)...\n",
			              $pct_chr, $pct_mate, $pct_rdlen, $pct_mapq;
		}
		if($has_incl) {
			if(!defined($inclSet{$ev_tname})) {
				$nskip_chr++; $nskip++; next;
			}
		}
		if(defined($exclSet{$ev_tname})) {
			$nskip_chr++; $nskip++; next;
		}
		#   0x1 template having multiple segments in sequencing
		#   0x2 each segment properly aligned according to the aligner
		#   0x4 segment unmapped
		#   0x8 next segment in the template unmapped
		#  0x10 SEQ being reverse complemented
		#  0x20 SEQ of the next segment in the template being reversed 0x40 the first segment in the template
		#  0x80 the last segment in the template
		# 0x100 secondary alignment
		# 0x200 not passing quality controls 0x400 PCR or optical duplicate
		
		# Parse some SAM flags
		my $ev_paired      = ($ev_flags & 0x1 ) != 0;
		my $ev_mate_unmap  = ($ev_flags & 0x8 ) != 0;
		my $ev_mate1       = ($ev_flags & 0x40) != 0;
		my $ev_mate2       = ($ev_flags & 0x80) != 0;
		!$ev_mate1 || !$ev_mate2 || die;
		if($ev_paired && $ev_mate1 && !$m1) { $nskip_mate++; $nskip++; next; }
		if($ev_paired && $ev_mate2 && !$m2) { $nskip_mate++; $nskip++; next; }
		if($ev_paired && $ev_mate_unmap && !$unp) { $nskip_mate++; $nskip++; next; }
		defined($ev_mapq) || die "Not enough tokens in methylation evidence record:\n$_";
		$ev_toff  == int($ev_toff ) || die;
		$ev_mapq  == int($ev_mapq ) || die;
		$ev_cy    == int($ev_cy   ) || die;
		$ev_allen == int($ev_allen) || die;
		$ev_score == int($ev_score) || die;
		$ev_score = abs($ev_score);
		$ev_cy < $ev_allen || 
			die "Sequencing cycle was greater than alignment length";
		if(defined($readlen_min) && $ev_allen < $readlen_min) {
			$nskip_rdlen++; $nskip++; next;
		}
		if(defined($readlen_max) && $ev_allen > $readlen_max) {
			$nskip_rdlen++; $nskip++; next;
		}
		if(defined($mapq_min) && $ev_mapq < $mapq_min) {
			$nskip_mapq++; $nskip++; next;
		}
		if($rev_cy) {
			$ev_cy = $ev_allen - $ev_cy - 1;
		}
		if(!$ev_watson) {
			$ev_al = "T" if $ev_al eq "A";
			$ev_al = "C" if $ev_al eq "G";
		}
		$ev_al = uc $ev_al;
		$ev_al = "O" if $ev_al ne "C" && $ev_al ne "T";
		$lo_cyc = $ev_cy if $ev_cy < $lo_cyc;
		$hi_cyc = $ev_cy if $ev_cy > $hi_cyc;
		$ev{$ev_cy}{$ev_al}++;
		my $cost_st = 0;
		if(defined($ev_score) && $cost_strat > 0) {
			$cost_st = int($ev_score/$cost_strat);
			$hi_cost = $cost_st if $cost_st > $hi_cost;
			$lo_cost = $cost_st if $cost_st < $lo_cost;
		}
		my $len_st = 0;
		if(defined($ev_allen) && $len_strat > 0) {
			$len_st = int($ev_allen/$len_strat);
			$hi_len  = $len_st  if $len_st  > $hi_len;
			$lo_len  = $len_st  if $len_st  < $lo_len;
		}
		my $stratum = "${len_st}_${cost_st}";
		if($strat_len > 0) {
			$evLenStrat{$len_st}{$ev_cy}{$ev_al}++;
		}
		if($strat_score > 0) {
			$evCostStrat{$cost_st}{$ev_cy}{$ev_al}++;
			$evCostFullStrat{$cost_st}{$ev_cy}{$ev_al}++ if $ev_allen >= $fullLength;
		}
		if($strat_len > 0 || $strat_score > 0) {
			$evStrat{$stratum}{$ev_cy}{$ev_al}++;
		}
	}
	close($ev_fh);
}
print STDERR "$nev read-level measurements tallied\n";
print STDERR "  $nskip skipped due to filters\n";
print STDERR "    $nskip_chr skipped due to chromosome filter\n";
print STDERR "    $nskip_mate skipped due to mate filter\n";
print STDERR "    $nskip_rdlen skipped due to read length filter\n";
print STDERR "    $nskip_mapq skipped due to MAPQ filter\n";
if($strat_len > 0) {
	print STDERR "Found the following length strata:\n";
	for (keys %evLenStrat) { print STDERR "  $_\n"; }
}
if($strat_score > 0) {
	print STDERR "Found the following cost strata:\n";
	for (keys %evCostStrat) { print STDERR "  $_\n"; }
	print STDERR "Found the following cost strata for full-length alignments:\n";
	for (keys %evCostFullStrat) { print STDERR "  $_\n"; }
}
if($strat_len > 0 || $strat_score > 0) {
	print STDERR "Found the following strata:\n";
	for (keys %evStrat) { print STDERR "  $_\n"; }
}

# Open one output file per srtatum
my %stfh = ();
my %lenfh = ();
my %costfh = ();
my %costFullfh = ();
my @lenStratFn = ();
my @lenStratLabs = ();
my @costStratFn = ();
my @costStratLabs = ();
my @costFullStratFn = ();
my @costFullStratLabs = ();

my $ext = ".tsv";
$ext = ".rev$ext" if $rev_cy;

print STDERR "Writing tables...\n";
if($strat_len > 0 || $strat_score > 0) {
	for (keys %evStrat) {
		defined($stfh{$_}) && die;
		my $fn = "mbias_$_$ext";
		open($stfh{$_}, ">$output_dir/$fn") || die "Could not open '$fn' for writing";
		print {$stfh{$_}} "Offset\tT\tC\tOther\n";
	}
}
if($strataPlots) {
	if($strat_len > 0) {
		for (sort {$a <=> $b} keys %evLenStrat) {
			defined($lenfh{$_}) && die;
			my $fn = "mbias_len_$_$ext";
			open($lenfh{$_}, ">$output_dir/$fn") || die "Could not open '$fn' for writing";
			print {$lenfh{$_}} "Offset\tT\tC\tOther\n";
			push @lenStratFn, $fn;
			push @lenStratLabs, lenStratToStr($_);
		}
	}
	if($strat_score > 0) {
		for (sort {$a <=> $b} keys %evCostStrat) {
			defined($costfh{$_}) && die;
			my $fn = "mbias_cost_$_$ext";
			open($costfh{$_}, ">$output_dir/$fn") || die "Could not open '$fn' for writing";
			print {$costfh{$_}} "Offset\tT\tC\tOther\n";
			push @costStratFn, $fn;
			push @costStratLabs, costStratToStr($_);
		}
		for (sort {$a <=> $b} keys %evCostFullStrat) {
			defined($costFullfh{$_}) && die;
			my $fn = "mbias_cost_full_length_$_$ext";
			open($costFullfh{$_}, ">$output_dir/$fn") || die "Could not open '$fn' for writing";
			print {$costFullfh{$_}} "Offset\tT\tC\tOther\n";
			push @costFullStratFn, $fn;
			push @costFullStratLabs, costStratToStr($_);
		}
	}
}
for (keys %evStrat)     { defined($stfh{$_})   || die; }
if($strataPlots) {
	for (keys %evLenStrat)  { defined($lenfh{$_})  || die; }
	for (keys %evCostStrat) { defined($costfh{$_}) || die; }
	for (keys %evCostFullStrat) { defined($costFullfh{$_}) || die; }
}

my $oafn = "mbias$ext";
open(OALL, ">$output_dir/$oafn") || die;
print OALL "Offset\tT\tC\tOther\n";

# Output stratified tables and whole table
for my $c ($lo_cyc..$hi_cyc) {
	if($strataPlots) {
		if($strat_len > 0 || $strat_score > 0) {
			for my $st (keys %evStrat) {
				# Compose line for this cycle and this stratum
				defined($stfh{$st}) || die;
				my $u = $evStrat{$st}{$c}{T} || 0;
				my $m = $evStrat{$st}{$c}{C} || 0;
				my $o = $evStrat{$st}{$c}{O} || 0;
				print {$stfh{$st}} "".($c+1)."\t$u\t$m\t$o\n";
			}
		}
		if($strat_len > 0) {
			for my $ln (keys %evLenStrat) {
				# Compose line for this cycle and this stratum
				defined($lenfh{$ln}) || die;
				my $u = $evLenStrat{$ln}{$c}{T} || 0;
				my $m = $evLenStrat{$ln}{$c}{C} || 0;
				my $o = $evLenStrat{$ln}{$c}{O} || 0;
				print {$lenfh{$ln}} "".($c+1)."\t$u\t$m\t$o\n";
			}
		}
		if($strat_score > 0) {
			for my $co (keys %evCostStrat) {
				# Compose line for this cycle and this stratum
				defined($costfh{$co}) || die;
				my $u = $evCostStrat{$co}{$c}{T} || 0;
				my $m = $evCostStrat{$co}{$c}{C} || 0;
				my $o = $evCostStrat{$co}{$c}{O} || 0;
				print {$costfh{$co}} "".($c+1)."\t$u\t$m\t$o\n";
			}
			for my $co (keys %evCostFullStrat) {
				# Compose line for this cycle and this stratum
				defined($costFullfh{$co}) || die;
				my $u = $evCostFullStrat{$co}{$c}{T} || 0;
				my $m = $evCostFullStrat{$co}{$c}{C} || 0;
				my $o = $evCostFullStrat{$co}{$c}{O} || 0;
				print {$costFullfh{$co}} "".($c+1)."\t$u\t$m\t$o\n";
			}
		}
	}
	# Compose line for this cycle
	my $u = $ev{$c}{T} || 0;
	my $m = $ev{$c}{C} || 0;
	my $o = $ev{$c}{O} || 0;
	print OALL "".($c+1)."\t$u\t$m\t$o\n";
}
close(OALL);
if($strataPlots) {
	for (keys %evStrat)     { defined($stfh{$_}) || die; close($stfh{$_}); }
	for (keys %evLenStrat)  { close($lenfh{$_}); }
	for (keys %evCostStrat) { close($costfh{$_}); }
	for (keys %evCostFullStrat) { close($costFullfh{$_}); }
}

if($strataPlots) {
	# Output a table showing stratum membership
	my $oamemfn = "$output_dir/ev_stratum$ext";
	open(OAMEM, ">$oamemfn") || die "Could not open '$oamemfn' for writing";
	# Print column headers
	for my $co ($lo_cost..$hi_cost) {
		print OAMEM costStratToStr($co);
		print OAMEM "\t" if $co < $hi_cost;
	}
	print OAMEM "\n";
	for my $ln ($lo_len..$hi_len) {
		print OAMEM lenStratToStr($ln)."\t";
		for my $co ($lo_cost..$hi_cost) {
			my $stratum = "${ln}_${co}";
			my $ev = 0;
			for ($lo_cyc..$hi_cyc) {
				my $u = $evStrat{$stratum}{$_}{T} || 0;
				my $m = $evStrat{$stratum}{$_}{C} || 0;
				$ev += ($u + $m);
			}
			my $u0 = $evStrat{$stratum}{$lo_cyc}{T} || 0;
			my $m0 = $evStrat{$stratum}{$lo_cyc}{C} || 0;
			my $rds = $u0 + $m0;
			print OAMEM "$rds";
			print OAMEM "\t" if $co < $hi_cost;
		}
		print OAMEM "\n";
	}
	close(OAMEM);
}

print STDERR "Done\n";
