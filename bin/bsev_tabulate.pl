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
# bsev_tabulate.pl
#
# Takes a directory of BSC evidence files and outputs tables describing
# methylation evidence at each cytosine in the genome.  We accomplish this by
# walking through the reference looking for Cs/CpGs and pulling evidence from
# the corresponding unsorted evidence file.
#
# Many reference sequences
#
# This is trickier when there are many reference sequences, since in the
# previous step, we may not have wanted to 
#
# The tables are:
# 1. $NAME.c_strand_tab.tsv
# 2. $NAME.cpg_strand_tab.tsv
# 3. $NAME.cpg_tab.tsv
#
# Where $NAME is set by the user with the --name option.
#
# $NAME.c_strand_tab.tsv contains a table where each row is a position in the
# genome that corresponds to a cytosine on either the watson strand (i.e. C) or
# the crick strand (i.e. G).
#
# TODO: Need to be able to aggregate summaries for # loci, # pieces of ev, etc.
# across threads.
#
# TODO: Need mate-specific cycle filters.
#

use strict;
use warnings;
use File::Path;
use File::Basename;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use File::Spec;
use BtlBio::Util::File;
use BtlBio::Util::Temp;
#use BtlBio::Align::Bisulfite::BsEvidence;
#use BtlBio::Align::Bisulfite::BsFilter;
use BtlBio::Format::FastaIndexed;
use BtlBio::Align::Bisulfite::BsSummary;
use BtlBio::Contrib::ForkManager;
use Carp;

my %ref = ();           # all the reference sequences
my $output_dir_cpg     = undef; # output directory for merged CpG ev
my $output_dir_cpg_str = undef; # output directory for stranded CpG ev
my $output_dir_c       = undef; # output directory for C ev
my $mapq_min = 20;      # evidence from alignments w/ MAPQ < this is ignored
my $qual_strat = 10;    # stratify base quals into intervals of this size
my $max_strat = 40;     # qualities greater than this are rounded to this
my $qual_min = 0;       # minimum quality value to allow in
my $readlen_min = -1;    # evidence from reads shorter are filtered
my $readlen_max = -1;    # evidence from reads longer are filtered
my $cyc_min = -1;        # evidence from cycles < min are filtered
my $cyc_max = -1;        # evidence from cycles > max are filtered
my $cyc_trim_begin = 0;  # ev from poss more than this far from start are filtered
my $cyc_trim_end = 0;    # ev from poss more than this far from end are filtered
my $cyc_trim_begin1 = 0; # ev from poss more than this far from start on mate 1 are filtered
my $cyc_trim_end1 = 0;   # ev from poss more than this far from end on mate 1 are filtered
my $cyc_trim_begin2 = 0; # ev from poss more than this far from start on mate 2 are filtered
my $cyc_trim_end2 = 0;   # ev from poss more than this far from end on mate 2 are filtered
my $stop_after_poss = 0;# stop after this many fasta positions
my $nthreads = 1;       # number of threads to use
my $conv_64_33 = 0;     # convert quals from 64-based to 33-based
my $rec_idx_load = "";  # save record index here
my $rec_idx_save = "";  # load record index from here
my $compress = "none";  # compression type
my @ev_dirs = ();       # evidence directories
my @ref_fns = ();       # fasta files with references
my @name_map_fn = ();   # files with reference name mappings
my %name_map = ();      # map from long to short names

my $usage = qq!
Usage:
  bsev_tabulate.pl [options*] -- <ev-dirs> -- <refs>

Take a directory of *sorted* read-level measurements, and a set of FASTA files
with all the reference sequences, and output a set of summary tables with one
row per methylation locus.

Required arguments:
  
  <ev-dirs>
    Paths to directories containing sorted read-level measurements.  If many
    are specified, they are effectively merged.
  <refs>
    Fasta files containing all original reference sequences.

Options (defaults in parentheses):

  --cpg=<path>              Output CpG summaries with merged strands to <path>
  --cpg-str=<path>          Output stranded CpG summaries to <path>
  --c=<path>                Output all-C summaries to <path> (only relevant if
                            all-C measurements were output by alignment script)
  --num-threads=<int>       Use <int> threads/CPUs simultaneously (1)
  --name-map=<path>         Many-to-one map from raw reference names to bin
                            names.  Output will be binned by bin name.
  --gzip                    gzip compress table output (off)
  --bzip2                   bzip2 compress tabl output (off)

 Filtering:
  --mapq-min=<int>          Filter measurements with MAPQ < <int>
  --qual-min=<int>          Filter measurements with base qual < <int>
  --readlen-min=<int>       Filter measurements from reads shorter than <int>
  --readlen-max=<int>       Filter measurements from reads longer than <int>
  --cyc-min=<int>           Filter measurements w/ sequencing cycle < <int>
  --cyc-max=<int>           Filter measurements w/ sequencing cycle > <int>
  --cyc-trim-begin=<int>    Filter measurements <= <int> positions from the 
                            beginning of the read
  --cyc-trim-end=<int>      Filter measurements <= <int> positions from the end
                            of the read
  --cyc-trim-begin-mate1=<int> Like --cyc-trim-begin but just for mate 1
  --cyc-trim-end-mate1=<int>   Like --cyc-trim-end but just for mate 1
  --cyc-trim-begin-mate2=<int> Like --cyc-trim-begin but just for mate 2
  --cyc-trim-end-mate2=<int>   Like --cyc-trim-end but just for mate 2
 
 Misc:
  --temp=<path>             Set path used for temporary FASTA files (TMPDIR)
  --help/--usage            Print this message and quit

Output format:

 One genome locus is described per line.  Each line has 12 tab-delimeted
 fields.  The fields are:
 
 1.  Reference name
 2.  Reference offset (1-based)
 3.  Locus type.  W=cytosine on the Watson strand, C=cytosine on the Crick
     strand, M="merged" - summarizes measurements from both Watson and Crick
     overlapping a CpG
 4.  Base qualities for methylated read-level measurements.  Length of this
     string equals the number of methylated read-level measurements.  Qualities
     are Phred+33 encoded.
 5:  Number of distinct read positions where methylated read-level measurements
     were observed
 6.  Base qualities for unmethylated read-level measurements.  Length of this
     string equals the number of unmethylated read-level measurements.
     Qualities are Phred+33 encoded.
 7.  Number of distinct read positions where unmethylated read-level
     measurements were observed
 8.  # read-level measurements filtered out because they came from a disallowed
     sequencing cycle
 9.  # read-level measurements filtered out because they came from a disallowed
     read length
 10. # read-level measurements filtered out because they did not indicate
     presense or absense of methylation (e.g. a 'G' on the Watson strand)
 11. # read-level measurements filtered out because they came from a read with
     a disallowed mapping quality (MAPQ)
 12. # read-level measurements filtered out because they came from a read
     position with a disallowed base quality
 
 Note: filters are applied in the order listed above (cycle, read length,
 allele, MAPQ, base quality), and a single filtered-out read-level measurement
 will contribute to only one of the filter counters.

!;

##
# Die and print the usage message.
#
sub dieusage($$) {
	my ($msg, $level) = @_;
	print STDERR "$usage\n--\nError $level:\n$msg\n";
	exit $level;
}

$SIG{INT} = sub { confess(); };

GetOptions (
	"evidence-directory=s"            => \@ev_dirs,
	"reference=s"                     => \@ref_fns,
	"cpg-output-directory|cpg=s"          => \$output_dir_cpg,
	"cpg-stranded-output-directory|cpg-str=s" => \$output_dir_cpg_str,
	"c-output-directory|c=s"          => \$output_dir_c,
	"name-map=s"                      => \@name_map_fn,
	"mapq-min=i"                      => \$mapq_min,
	"qual-min=i"                      => \$qual_min,
	"readlen-min=i"                   => \$readlen_min,
	"readlen-max=i"                   => \$readlen_max,
	"cyc-min=i"                       => \$cyc_min,
	"cyc-max=i"                       => \$cyc_max,
	"cyc-trim-begin=i"                => \$cyc_trim_begin,
	"cyc-trim-end=i"                  => \$cyc_trim_end,
	"cyc-trim-begin-mate1=i"          => \$cyc_trim_begin1,
	"cyc-trim-end-mate1=i"            => \$cyc_trim_end1,
	"cyc-trim-begin-mate2=i"          => \$cyc_trim_begin2,
	"cyc-trim-end-mate2=i"            => \$cyc_trim_end2,
	"stop-after=i"                    => \$stop_after_poss,
	"conv-64-to-33"                   => \$conv_64_33,
	"save-record-index=s"             => \$rec_idx_save,
	"load-record-index=s"             => \$rec_idx_load,
	"num-threads|threads|num-cpus|cpus=i" => \$nthreads,
	"compress=s"                      => \$compress,
	"gzip"                            => sub { $compress = "gzip"; },
	"bzip2"                           => sub { $compress = "bzip2"; },
	"help|usage"                      => sub { print $usage; exit 0; }
) || dieusage("Bad option", 0);

# Parse out arguments from among the double-dashes
my $dds = 0;
for my $a (@ARGV) {
	if($a eq "--") { $dds++; next; }
	if($dds == 0) {
		push @ev_dirs, $a;
	} elsif($dds == 1) {
		push @ref_fns, $a;
	} else {
		die "Error: Too many --s in argument list";
	}
}

scalar(@ev_dirs) > 0 || dieusage("Must specify at least one <ev-dir> ".
                                 "(directory containing sorted read-level ".
                                 "measurements)", 0);

print STDERR "---- Settings ----\n";
print STDERR "Minimum MAPQ: $mapq_min\n";
print STDERR "Name-map files:\n";
for(0..$#name_map_fn) {print STDERR "  $name_map_fn[$_]\n"}
print STDERR "Evidence directories:\n";
for(0..$#ev_dirs) {print STDERR "  $ev_dirs[$_]\n"}
print STDERR "Reference files:\n";
for(0..$#ref_fns) {print STDERR "  $ref_fns[$_]\n"}
print STDERR "Convert base-64 quals to base-33?: $conv_64_33\n";
print STDERR "Comrpession type: $compress\n";
print STDERR "Filter measurements from reads shorter than: $readlen_min\n" if $readlen_min != -1;
print STDERR "Filter measurements from reads longer than: $readlen_max\n" if $readlen_max != -1;
print STDERR "Trim positions from beginning: $cyc_trim_begin\n";
print STDERR "Trim positions from end: $cyc_trim_end\n";
print STDERR "Trim positions from beginning of mate 1: $cyc_trim_begin1\n";
print STDERR "Trim positions from end of mate 1: $cyc_trim_end1\n";
print STDERR "Trim positions from beginning of mate 2: $cyc_trim_begin2\n";
print STDERR "Trim positions from end of mate 2: $cyc_trim_end2\n";
print STDERR "Comrpession type: $compress\n";
print STDERR "Output directory for merged CpG ev: ".($output_dir_cpg || "(none)")."\n";
print STDERR "Output directory for stranded CpG ev: ".($output_dir_cpg_str || "(none)")."\n";
print STDERR "Output directory for C ev: ".($output_dir_c || "(none)")."\n";
print STDERR "------------------\n";

my $do_cpg  = defined($output_dir_cpg);
my $do_cpgs = defined($output_dir_cpg_str);
my $do_c    = defined($output_dir_c);

$do_cpg || $do_cpgs || $do_c || die "Must specify at least one of: --cpg, --cpg-str, --c";

my @m_heads = ();
my @u_heads = ();
if(0) {
	for(my $i = 0; $i <= $max_strat; $i += $qual_strat) {
		push @m_heads, "M$i";
		push @u_heads, "U$i";
	}
} else {
	push @m_heads, "Mstr";
	push @u_heads, "Ustr";
}

my @heads = (
	"ref",          # 1.
	"off",          # 2.
	"strand",       # 3.
	@m_heads,
	"Mcy",          # 9.
	@u_heads,
	"Ucy",          # 15.
	"filt_cycle",   # 16.
	"filt_readlen", # 17.
	"filt_allele",  # 18.
	"filt_mapq",    # 19.
	"filt_baseq");  # 20.

print STDERR "Parsing name-map files...\n";
for(@name_map_fn) {
	open(NM, $_) || die "Could not open name-map file '$_' for reading";
	while(<NM>) {
		chomp;
		next if /^\s*$/;
		next if /^\s*#/;
		my @s = split(/\t/);
		$name_map{$s[0]} = $s[1];
	}
	close(NM);
}

# Map from bins to all the evidence files for the bin
my %bin_to_fn = ();

# Find all the evidence files and their respective bin
print STDERR "Walking reference...\n";
#my $cpg = undef;
my $cpg = 1;
for my $evd (@ev_dirs) {
	if(! -f "$evd/.c.sorted.ev" && ! -f "$evd/.cpg.sorted.ev") {
		die "Error: '$evd' is not a sorted Bsmooth evidence directory\n";
	}
	# Establish whether the evidence is C or CpG
	#$cpg = -f "$evd/.cpg.sorted.ev";
	# Compile a map from bins to all the evidence files for the bin
	for my $evfn (<$evd/*>) {
		my ($base, $dir) = fileparse($evfn);
		my $evbinfn = "$dir/.$base.full_name";
		if(! -f $evbinfn) {
			die "Error: evidence file '$evfn' missing bin file '$evbinfn'\n";
		}
		# Determine bin
		open(BIN, $evbinfn) || die "Could not open '$evbinfn' for reading";
		my $bin = <BIN>; chomp($bin);
		close(BIN);
		# Add evidence file to bin
		defined($bin_to_fn{$bin}) && die;
		$bin_to_fn{$bin} = $evfn;
	}
}
mkpath($output_dir_cpg)     if $do_cpg;
mkpath($output_dir_cpg_str) if $do_cpgs;
mkpath($output_dir_c)       if $do_c;

# Now we know, for each reference sequence, what are all the files that have
# evidence for that sequence.  Now we scan through the reference sequences and
# emit tables for each.

my %revcomp_map = (
	"A" => "T", "a" => "t",
	"C" => "G", "c" => "g",
	"G" => "C", "g" => "c",
	"T" => "A", "t" => "a"
);

my $cthreads = 1;
my $tid = 1;
my $ctid = 0;
my @child_pids = ();
my %rec_idx = ();
print STDERR "Making record index of FASTA inputs...\n";
if($rec_idx_load ne "") {
	my $fh = openex($rec_idx_load);
	fasta_idx_load_recidx(\%rec_idx, $fh);
	print STDERR "  Loaded record index with ".(scalar(keys %rec_idx))." entries\n";
	close($fh);
} else {
	for my $fn (@ref_fns) {
		$fn = File::Spec->rel2abs($fn);
		fasta_idx_index(\%rec_idx, $fn, 0);
	}
	print STDERR "  Made record index with ".(scalar(keys %rec_idx))." entries\n";
}
if($rec_idx_save ne "") {
	open(my $ofh, ">$rec_idx_save") || die "Could not open '$rec_idx_save' for writing";
	fasta_idx_save_recidx(\%rec_idx, $ofh);
	close($ofh);
}

# Now read through all the references once and emit table rows as we go
my ($ofh_c_w, $ofh_c_c, $ofh_cpg_w, $ofh_cpg_c, $ofh_cpg_m) =
	(undef, undef, undef, undef, undef);
# Holds all the various aspects of the piece of evidence currently under consideration
my $last_name = "";
my $sorted_fh = undef;
my $sorted_line = undef;
my $pos = 0;
my $pos_ival = 500000;
my $nloci = 0;
my ($nloci_cov, $nloci_uncov) = (0, 0);
my $nev = 0;
my $nev_summ = 0;
my ($nev_w, $nev_c) = (0, 0);
my $nfev = 0;
my ($nfev_w, $nfev_c) = (0, 0);
my ($nrepcnt_c, $nrepcnt_cpg, $nrepcnt_cpg_str) = (0, 0, 0);
my $first = 1;
my $cseq = 0;
my $mine = 0;
my $qb = $conv_64_33 ? 64 : 33;
my ($nfilt_cyc, $nfilt_rdl, $nfilt_nuc, $nfilt_mapq, $nfilt_qual) = (0, 0, 0, 0, 0);
my $pm = undef;
my $childFailed = 0;
my $childFailedPid = 0;
my $summary_dir = temp_name();
(-d $summary_dir || -f $summary_dir) && die;
mkpath($summary_dir);
if($nthreads > 1) {
	$pm = new Parallel::ForkManager($nthreads);
	# Set callback for when a child finishes up so we can get its exit code
	$pm->run_on_finish(sub {
		my ($pid, $exit_code, $ident) = @_;
		if($nthreads > 0 && $exit_code != 0) {
			$childFailed = $exit_code;
			$childFailedPid = $pid;
		}
	});
}
# Get reference record index
my @fns = values %bin_to_fn;
@fns = sort { -s $b <=> -s $a } @fns;
my $expected_nfields = undef;
for my $fn (@fns) {
	$ctid++;
	next if ($nthreads > 1 && $pm->start);
	$tid = $ctid;
	open(SUMM, ">$summary_dir/$tid") || die "Could not open $summary_dir/$tid for writing";
	my $fn_base = fileparse($fn);
	$fn_base =~ s/\.gz$//;
	$fn_base =~ s/\.bz2$//;
	$fn_base =~ s/\.tsv$//;
	my $fh = openex($fn);       # evidence filehandle
	my $fh_next = readline $fh; # peek ahead
	defined($fh_next) || die;
	chomp($fh_next);
	$ofh_c_w = undef;
	$ofh_c_c = undef;
	my $open_pre = ">";
	$open_pre = "| gzip  -c > " if $compress eq "gzip";
	$open_pre = "| bzip2 -c > " if $compress eq "bzip2";
	my $open_suff = "";
	$open_suff = ".gz"  if $compress eq "gzip";;
	$open_suff = ".bz2" if $compress eq "bzip2";
	if($do_c) {
		open($ofh_c_w, "$open_pre$output_dir_c/$fn_base.c.tsv$open_suff") ||
			die "Could not open '$output_dir_c/$fn_base.c.tsv$open_suff' for writing";
		$ofh_c_c = $ofh_c_w;
		print {$ofh_c_w} join("\t", @heads)."\n";
	}
	if($do_cpgs) {
		open($ofh_cpg_w, "$open_pre$output_dir_cpg_str/$fn_base.stranded.cpg.tsv$open_suff") ||
			die "Could not open '$output_dir_cpg_str/$fn_base.stranded.cpg.tsv$open_suff' for writing";
		$ofh_cpg_c = $ofh_cpg_w;
		print {$ofh_cpg_w} join("\t", @heads)."\n";
	}
	if($do_cpg) {
		open($ofh_cpg_m, "$open_pre$output_dir_cpg/$fn_base.cpg.tsv$open_suff") ||
			die "Could not open '$output_dir_cpg/$fn_base.cpg.tsv$open_suff' for writing";
		print {$ofh_cpg_m} join("\t", @heads)."\n";
	}
	my $fa_fh = undef;
	my $cur_ref_off = undef;
	my $cur_ref_len = -1;
	my ($c_prev, $c, $c_next) = (-1, -1, -1);
	my $ref_carry = ""; # carryover
	my ($ev_tname_last, $ev_toff_last) = ("", -1);
	my $summ      = BtlBio::Align::Bisulfite::BsSummary->new();
	my $summ_prev = BtlBio::Align::Bisulfite::BsSummary->new();
	$summ->{_empty} || die;
	$summ_prev->{_empty} || die;
	# Next batch of read-level measurements
	my $iter = 0;
	my $last_tname = "";
	my $last_toff = -1;
	##
	# This loop is extremely complicated because it has to deal with (a) chunks
	# of measurements, (b) the boundaries between references, (c) the end of
	# the measurement file.  And it has to do all this while walking through
	# the reference sequence in tandem.
	#
	while(defined($fh_next)) {
		$summ_prev = $summ;
		$summ = BtlBio::Align::Bisulfite::BsSummary->new();
		my ($tname, $toff) = (undef, undef);
		# Set to 1 if current summary is last one in reference sequence (incl.
		# if it's the last in the file).
		my $tname_end = 0;
		while(defined($fh_next)) {
			# This loop (a) grabs the next chunk of measurements, and (b) sets
			# the tname_end variable if it's the last chunk for this reference
			my $line = $fh_next;
			$fh_next = readline $fh;
			my $last_in_batch = 0;
			my ($ev_tname, $ev_toff, $ev_readid, $ev_al, $ev_watson, $ev_fw,
			    $ev_flags, $ev_qu1, $ev_qu2, $ev_cy, $ev_allen,
			    $ev_score, $ev_mapq) = split(/\t/, $line, -1);
			$tname_end && die;
			if(defined($fh_next)) {
				chomp($fh_next);
				my ($ev_tname_next, $ev_toff_next) = split(/\t/, $fh_next, -1);
				$tname_end = $ev_tname ne $ev_tname_next;
				$last_in_batch = $tname_end || $ev_toff_next != $ev_toff;
			} else {
				$tname_end = 1;
				$last_in_batch = 1;
			}
			$ev_tname ne "" || die;
			$ev_toff > 0 || die "Expected 1-based ev_toff, got $ev_toff";
			$ev_toff--;
			($tname, $toff) = ($ev_tname, $ev_toff);
			$ev_watson == 0 || $ev_watson == 1 || die;
			$ev_fw     == 0 || $ev_fw     == 1 || die;
			int($ev_toff  ) == $ev_toff   || die;
			int($ev_cy    ) == $ev_cy     || die;
			int($ev_flags ) == $ev_flags  || die;
			int($ev_score ) == $ev_score  || die;
			int($ev_allen ) == $ev_allen  || die;
			int($ev_mapq  ) == $ev_mapq   || die;
			defined($ev_mapq) || croak("Not enough tokens in read-level measurement:\n$line");
			my $ev_watson_al = $ev_al;
			if(!$ev_watson) {
				$ev_watson_al = $revcomp_map{$ev_al} || $ev_al;
			}
			my ($xqu1, $xqu2) = (ord($ev_qu1) - $qb, ord($ev_qu2) - $qb);
			$xqu2 = -1 if $ev_qu2 eq "-1";
			$xqu1 = 0 if $xqu1 < 0;
			$xqu2 = $xqu1 if $xqu2 < 0;
			my $ev_qual_summ = ($xqu1 + $xqu2)/2;
			$nev++;
			if($ev_watson) {
				$nev_w++;
			} else {
				$nev_c++;
			}
			my $filt_mapq = $ev_mapq < $mapq_min; # Filter the evidence based on MAPQ?
			my $filt_qual = $ev_qual_summ < $qual_min; # Filter the evidence based on quals?
			my $filt_cyc = 0;     # Filter the evidence based on cycle?
			my $filt_readlen = 0; # Filter the evidence based on read len?
			$filt_readlen = 1 if $readlen_min != -1 && $ev_allen < $readlen_min;
			$filt_readlen = 1 if $readlen_max != -1 && $ev_allen > $readlen_max;
			$filt_cyc = 1 if $cyc_min != -1 && $ev_cy < $cyc_min;
			$filt_cyc = 1 if $cyc_max != -1 && $ev_cy > $cyc_max;
			$filt_cyc = 1 if $ev_cy < $cyc_trim_begin;
			$filt_cyc = 1 if $ev_cy >= $ev_allen - $cyc_trim_end;
			my $ev_paired      = ($ev_flags & 0x1 ) != 0;
			my $ev_mate_unmap  = ($ev_flags & 0x8 ) != 0;
			my $ev_mate1       = ($ev_flags & 0x40) != 0;
			my $ev_mate2       = ($ev_flags & 0x80) != 0;
			if($ev_paired && $ev_mate1) {
				$filt_cyc = 1 if $ev_cy < $cyc_trim_begin1;
				$filt_cyc = 1 if $ev_cy >= $ev_allen - $cyc_trim_end1;
			}
			if($ev_paired && $ev_mate2) {
				$filt_cyc = 1 if $ev_cy < $cyc_trim_begin2;
				$filt_cyc = 1 if $ev_cy >= $ev_allen - $cyc_trim_end2;
			}
			# Filter the evidence based on nucleotide?
			my $al = $ev_watson_al;
			my $filt_nuc = ($al ne "C" && $al ne "T");
			my $filt = $filt_mapq || $filt_qual || $filt_cyc || $filt_readlen || $filt_nuc;
			if($filt_cyc) {
				$summ->{_filt_cyc}++; $nfilt_cyc++;
			} elsif($filt_readlen) {
				$summ->{_filt_rdl}++; $nfilt_rdl++;
			} elsif($filt_nuc) {
				$summ->{_filt_nuc}++; $nfilt_nuc++;
			} elsif($filt_mapq) {
				$summ->{_filt_mapq}++; $nfilt_mapq++;
			} elsif($filt_qual) {
				$summ->{_filt_qual}++; $nfilt_qual++;
			}
			if(!$filt) {
				# Evidence not filtered
				$ev_qual_summ >= 0 || die;
				my $strat = ($ev_qual_summ/$qual_strat);
				$strat = $max_strat if $strat > $max_strat;
				($strat >= 0 && $strat <= $max_strat) || die;
				$al eq "C" || $al eq "T" || die;
				$summ->{_qualstr}{$al} = "" unless defined($summ->{_qualstr}{$al});
				$summ->{_qualstr}{$al} .= chr($ev_qual_summ+33);
				$summ->{_mapqstr}{$al} = "" unless defined($summ->{_mapqstr}{$al});
				$summ->{_mapqstr}{$al} .= chr($ev_mapq+33);
				$summ->{_count}{$al}++;
				$summ->{_qual}{$al}[$strat]++; # Tally the base quality stratum
				$summ->{_cycs}{$al}{$ev_cy}++; # Tally the cycle
				$summ->{_ncycs}{$al}++ if $summ->{_cycs}{$al}{$ev_cy} == 1;
				$summ->{_ncycs}{$al} <= length($summ->{_qualstr}{$al}) || die;
			} else {
				$nfev++;
				if($ev_watson) {
					$nfev_w++;
				} else {
					$nfev_c++;
				}
			}
			$summ->{_empty} = 0;
			$summ->rep_ok($qual_strat, $max_strat) || die;
			last if $last_in_batch;
		} # while(readline $fh)
		$summ->{_empty} && die;
		$nev_summ += $summ->num_measurements();
		# Did we move on to a new reference sequence?
		if($tname ne $last_tname) {
			# Look up where this reference appears in the input
			defined($rec_idx{$tname}) ||
				die "No entry for reference '$tname' in record index";
			my $fa_fn    = $rec_idx{$tname}{fn};
			my $seek_off = $rec_idx{$tname}{off};
			close($fa_fh) if defined($fa_fh);
			$fa_fh = openex($fa_fn);
			seek $fa_fh, $seek_off, 0;
			my $name_line = readline $fa_fh; chomp($name_line);
			substr($name_line, 0, 1) eq ">" || die;
			my $name = substr($name_line, 1);
			$name =~ s/\s.*//;
			$name eq $tname ||
				die "Expected name '$tname', got '$name' after seek to $fa_fn:$seek_off";
			$cur_ref_len = $rec_idx{$tname}{len};
			$cur_ref_len > 0 || die;
			$cur_ref_off = -1;
			# Read the first char, stick in c_next
			($c, $c_prev) = ("", "");
			while(1) {
				my $n = read $fa_fh, $c_next, 1;
				$n == 1 || die "Return from read was $n for ref:$tname, off:$cur_ref_off, len:$cur_ref_len, chars are $c_prev$c$c_next";
				last unless $c_next =~ /\s/;
				$c_next eq ">" && die;
			}
			$toff < $cur_ref_len || die;
			int($cur_ref_len) == $cur_ref_len || die;
			$c_prev = -1;
			$ref_carry = "";
		}
		$toff < $cur_ref_len || die "tname:$tname, toff:$toff, cur_ref_len:$cur_ref_len";
		#
		# toff_cur contains the 0-based offset of a piece of evidence
		#
		# As long as we're not yet to the next piece of evidence
		for(my $end = 0; $end <= 1; $end++) {
			last if $end && !$tname_end;
			$summ_prev = $summ if $end;
			while($cur_ref_off+1 <= $toff ||
			     ($end && $cur_ref_off+1 < $cur_ref_len))
			{
				# Move ahead in the reference
				$cur_ref_off++;
				$cur_ref_off < $cur_ref_len || die "ref:$tname, cur_ref_off:$cur_ref_off, cur_ref_len:$cur_ref_len, toff:$toff";
				$c_prev = $c;
				$c = $c_next;
				my $off0 = $cur_ref_off;
				my $off1 = $off0+1;
				while(1) {
					# Off the end?
					if($cur_ref_off == $cur_ref_len-1) {
						$c_next = "EOF";
						last;
					} else {
						my $n = read $fa_fh, $c_next, 1;
						$n == 1 || die "Return from read was $n for ref:$tname, off:$cur_ref_off, len:$cur_ref_len, chars are $c_prev$c$c_next";
						last unless $c_next =~ /\s/;
						$c_next eq ">" && die;
					}
				}
				if((++$pos % $pos_ival) == 0) {
					print STDERR "  TID $tid: $pos poss ($nloci mloci (cov:$nloci_cov, uncov:$nloci_uncov), $nev measures ($nev_w W/$nev_c C), $nfev filtered ($nfev_w W/$nfev_c C), now in $tname)...\n";
				}
				# cur_ref_off is the 0-based offset of the current reference
				# character.  toff_cur is the 0-based offset of the current set
				# of evidence
				if($toff == $cur_ref_off) {
					$end && die;
					# Last character we got is on top of 
					if($c eq "C" && $do_c) {
						$nloci++; $nloci_cov++;
						my ($fs, $repcnt) = $summ->summary_tsv_nostrata();
						$expected_nfields = scalar(@$fs) unless defined($expected_nfields);
						scalar(@$fs) == $expected_nfields ||
							die "Bad number of fields (expected $expected_nfields):\n".
							    join("\t", @$fs)."\n";
						$nrepcnt_c += $repcnt;
						print {$ofh_c_w} "$tname\t$off1\tW\t".join("\t", @$fs)."\n";
					}
					if($c eq "G" && $do_c) {
						$nloci++; $nloci_cov++;
						my ($fs, $repcnt) = $summ->summary_tsv_nostrata();
						$expected_nfields = scalar(@$fs) unless defined($expected_nfields);
						scalar(@$fs) == $expected_nfields ||
							die "Bad number of fields (expected $expected_nfields):\n".
							    join("\t", @$fs)."\n";
						$nrepcnt_c += $repcnt;
						print {$ofh_c_c} "$tname\t$off1\tC\t".join("\t", @$fs)."\n";
					}
					if(($do_cpg || $do_cpgs) && $c eq "G" && $c_prev eq "C") {
						# Evidence covers the G in a CpG
						$nloci += 2; $nloci_cov++;
						my ($fs, $repcnt) = $summ->summary_tsv_nostrata();
						$expected_nfields = scalar(@$fs) unless defined($expected_nfields);
						scalar(@$fs) == $expected_nfields ||
							die "Bad number of fields (expected $expected_nfields):\n".
							    join("\t", @$fs)."\n";
						$nrepcnt_cpg_str += $repcnt if $do_cpgs;
						my $cline = "$tname\t$off1\tC\t".join("\t", @$fs)."\n";
						if($last_tname eq $tname && $last_toff == $cur_ref_off-1) {
							# Previous evidence covers the C in the CpG
							my ($fs, $repcnt) = $summ_prev->summary_tsv_nostrata();
							$expected_nfields = scalar(@$fs) unless defined($expected_nfields);
							scalar(@$fs) == $expected_nfields ||
								die "Bad number of fields (expected $expected_nfields):\n".
								    join("\t", @$fs)."\n";
							print {$ofh_cpg_w} "$tname\t$off0\tW\t".join("\t", @$fs)."\n" if $do_cpgs;
							$summ_prev->merge_inplace($summ, 1);
							($fs, $repcnt) = $summ_prev->summary_tsv_nostrata();
							$expected_nfields = scalar(@$fs) unless defined($expected_nfields);
							scalar(@$fs) == $expected_nfields ||
								die "Bad number of fields (expected $expected_nfields):\n".
								    join("\t", @$fs)."\n";
							if($do_cpg) {
								$nrepcnt_cpg_str += $repcnt;
								$nrepcnt_cpg += $repcnt;
								print {$ofh_cpg_m} "$tname\t$off0\tM\t".join("\t", @$fs)."\n";
								$nloci_cov++;
							}
						} else {
							# No evidence covering the C in the CpG
							if($do_cpg || $do_cpgs) {
								print {$ofh_cpg_w} "$tname\t$off0\tW\t".bs_summary_blank_tsv_nostrata()."\n"  if $do_cpgs;
								$nrepcnt_cpg += $repcnt;
								print {$ofh_cpg_m} "$tname\t$off0\tM\t".join("\t", @$fs)."\n"  if $do_cpg;
								$nloci_uncov++;
							}
						}
						print {$ofh_cpg_c} $cline if $do_cpgs;
					}
				} else {
					# There's no evidence on this position
					if($c eq "C" && $do_c) {
						$nloci++; $nloci_uncov++;
						print {$ofh_c_w} "$tname\t$off1\tW\t".bs_summary_blank_tsv_nostrata()."\n";
					}
					if($c eq "G" && $do_c) {
						$nloci++; $nloci_uncov++;
						print {$ofh_c_c} "$tname\t$off1\tC\t".bs_summary_blank_tsv_nostrata()."\n";
					}
					if(($do_cpgs || $do_cpg) && $c eq "G" && $c_prev eq "C") {
						# Was Watson evidence empty?
						$nloci++;
						if($last_tname eq $tname && $last_toff == $cur_ref_off-1) {
							my ($fs, $repcnt) = $summ_prev->summary_tsv_nostrata();
							$expected_nfields = scalar(@$fs) unless defined($expected_nfields);
							scalar(@$fs) == $expected_nfields ||
								die "Bad number of fields (expected ".
								    "$expected_nfields):\n".join("\t", @$fs).
								    "\n";
							$nrepcnt_cpg += $repcnt if $do_cpg;
							$nrepcnt_cpg_str += $repcnt if $do_cpgs;
							print {$ofh_cpg_w} "$tname\t$off0\tW\t".join("\t", @$fs)."\n" if $do_cpgs;
							print {$ofh_cpg_m} "$tname\t$off0\tM\t".join("\t", @$fs)."\n" if $do_cpg;
							$nloci_cov++;
						} else {
							print {$ofh_cpg_w} "$tname\t$off0\tW\t".bs_summary_blank_tsv_nostrata()."\n" if $do_cpgs;
							print {$ofh_cpg_m} "$tname\t$off0\tM\t".bs_summary_blank_tsv_nostrata()."\n" if $do_cpg;
							$nloci_uncov++;
						}
						print {$ofh_cpg_c} "$tname\t$off1\tC\t".bs_summary_blank_tsv_nostrata()."\n" if $do_cpgs;
					}
				}
			}
			($last_tname, $last_toff) = ($tname, $toff);
		}
	} # while(!$done)
	$nev == $nev_summ || die "Number of measurements=$nev; only $nev_summ made it into summaries";
	close($ofh_c_w) if defined($ofh_c_w);
	close($ofh_cpg_w) if defined($ofh_cpg_w);
	close($ofh_cpg_m) if defined($ofh_cpg_m);
	close($fh);
	close($fa_fh) if defined($fa_fh);
	print SUMM join("\t",
		($pos, $nloci, $nloci_cov, $nloci_uncov, $nev, $nev_w, $nev_c, $nfev,
		 $nfev_w, $nfev_c, $nfilt_cyc, $nfilt_rdl, $nfilt_nuc, $nfilt_mapq,
		 $nfilt_qual, $nrepcnt_c, $nrepcnt_cpg, $nrepcnt_cpg_str))."\n";
	close(SUMM);
	$pm->finish if $nthreads > 1; # end of fork
}
if($nthreads > 1) {
	$pm->wait_all_children;
	if($childFailed) {
		croak("Child with PID $childFailedPid exited abnormally\n");
	} else { print STDERR "All children succeeded\n"; }
}
sleep(2);
print STDERR "Combining summaries into master summary\n";
my $tpos = 0;
my $tnloci = 0;
my ($tnloci_cov, $tnloci_uncov) = (0, 0);
my $tnev = 0;
my ($tnev_w, $tnev_c) = (0, 0);
my $tnfev = 0;
my ($tnfev_w, $tnfev_c) = (0, 0);
my ($tnfilt_cyc, $tnfilt_rdl, $tnfilt_nuc, $tnfilt_mapq, $tnfilt_qual) = (0, 0, 0, 0, 0);
my ($tnrepcnt_c, $tnrepcnt_cpg, $tnrepcnt_cpg_str) = (0, 0, 0);
while(<$summary_dir/*>) {
	open(SUMM, $_) || die "Could not open '$_' for reading";
	my $summ = <SUMM>;
	chomp($summ);
	my ($pos, $nloci, $nloci_cov, $nloci_uncov, $nev, $nev_w, $nev_c, $nfev,
	    $nfev_w, $nfev_c, $nfilt_cyc, $nfilt_rdl, $nfilt_nuc, $nfilt_mapq,
	    $nfilt_qual, $nrepcnt_c, $nrepcnt_cpg, $nrepcnt_cpg_str) =
	    split(/\t/, $summ, -1);
	$tpos             += $pos;
	$tnloci           += $nloci;
	$tnloci_cov       += $nloci_cov;
	$tnloci_uncov     += $nloci_uncov;
	$tnev             += $nev;
	$tnev_w           += $nev_w;
	$tnev_c           += $nev_c;
	$tnfev            += $nfev;
	$tnfev_w          += $nfev_w;
	$tnfev_c          += $nfev_c;
	$tnfilt_cyc       += $nfilt_cyc;
	$tnfilt_rdl       += $nfilt_rdl;
	$tnfilt_nuc       += $nfilt_nuc;
	$tnfilt_mapq      += $nfilt_mapq;
	$tnrepcnt_c       += $nrepcnt_c;
	$tnrepcnt_cpg     += $nrepcnt_cpg;
	$tnrepcnt_cpg_str += $nrepcnt_cpg_str;
	close(SUMM);
}
print STDERR "Removing summary directory '$summary_dir'\n";
rmtree($summary_dir);

my $tnfilt = $tnfilt_cyc + $tnfilt_rdl + $tnfilt_nuc + $tnfilt_mapq + $tnfilt_qual;
my $tnloci_cov_pct   = (100.00 * $tnloci_cov   / (1.0 * $tnloci || 1));
my $tnloci_uncov_pct = (100.00 * $tnloci_uncov / (1.0 * $tnloci || 1));
my $tnev_w_pct  = (100.00 * $tnev_w  / (1.0 * $tnev  || 1));
my $tnev_c_pct  = (100.00 * $tnev_c  / (1.0 * $tnev  || 1));
my $tnfev_w_pct = (100.00 * $tnfev_w / (1.0 * $tnfev || 1));
my $tnfev_c_pct = (100.00 * $tnfev_c / (1.0 * $tnfev || 1));
my $tnfev_pct   = (100.00 * $tnfev   / (1.0 * $tnev || 1));
my $tnnfev_pct  = (100.00 * ($tnev-$tnfev) / (1.0 * $tnev || 1));
my $tnfilt_pct  = (100.00 * $tnfilt  / (1.0 * $tnev || 1));
my $tnfilt_cyc_pct  = (100.00 * $tnfilt_cyc   / (1.0 * $tnfilt || 1));
my $tnfilt_rdl_pct  = (100.00 * $tnfilt_rdl   / (1.0 * $tnfilt || 1));
my $tnfilt_nuc_pct  = (100.00 * $tnfilt_nuc   / (1.0 * $tnfilt || 1));
my $tnfilt_mapq_pct = (100.00 * $tnfilt_mapq  / (1.0 * $tnfilt || 1));
my $tnfilt_qual_pct = (100.00 * $tnfilt_qual  / (1.0 * $tnfilt || 1));
my $tnrepcnt_c_pct       = (100.00 * $tnrepcnt_c       / (1.0 * ($tnev) || 1));
my $tnrepcnt_cpg_pct     = (100.00 * $tnrepcnt_cpg     / (1.0 * ($tnev) || 1));
my $tnrepcnt_cpg_str_pct = (100.00 * $tnrepcnt_cpg_str / (1.0 * ($tnev) || 1));

printf STDERR "Totals:\n";
printf STDERR "  $tpos reference positions\n";
printf STDERR "  $tnloci methylation loci\n";
printf STDERR "    $tnloci_cov (%.2f%%) covered\n", $tnloci_cov_pct;
printf STDERR "    $tnloci_uncov (%.2f%%) uncovered\n", $tnloci_uncov_pct;
printf STDERR "  $tnev read-level measurements\n";
printf STDERR "    $tnev_w (%.2f%%) watson\n", $tnev_w_pct;
printf STDERR "    $tnev_c (%.2f%%) crick\n", $tnev_c_pct;
printf STDERR "    --\n";
printf STDERR "    %d (%.2f%%) read-level measurements surviving filters\n", $tnev-$tnfev, $tnnfev_pct;
printf STDERR "    $tnfev (%.2f%%) read-level measurements filtered out\n", $tnfev_pct;
printf STDERR "      $tnfev_w (%.2f%%) watson\n", $tnfev_w_pct;
printf STDERR "      $tnfev_c (%.2f%%) crick\n", $tnfev_c_pct;
printf STDERR "      --\n";
printf STDERR "      $tnfilt_cyc (%.2f%%) due to cycle filter\n", $tnfilt_cyc_pct;
printf STDERR "      $tnfilt_rdl (%.2f%%) due to read-length filter\n", $tnfilt_rdl_pct;
printf STDERR "      $tnfilt_nuc (%.2f%%) due to allele filter\n", $tnfilt_nuc_pct;
printf STDERR "      $tnfilt_mapq (%.2f%%) due to MAPQ filter\n", $tnfilt_mapq_pct;
printf STDERR "      $tnfilt_qual (%.2f%%) due to alignment score filter\n", $tnfilt_qual_pct;
printf STDERR "    --\n";
printf STDERR "    $tnrepcnt_c (%.2f%%) reported in C/G tables\n", $tnrepcnt_c_pct;
printf STDERR "    --\n";
printf STDERR "    $tnrepcnt_cpg (%.2f%%) reported in CpG tables\n", $tnrepcnt_cpg_pct;
printf STDERR "    --\n";
printf STDERR "    $tnrepcnt_cpg_str (%.2f%%) reported in stranded CpG tables\n", $tnrepcnt_cpg_str_pct;
