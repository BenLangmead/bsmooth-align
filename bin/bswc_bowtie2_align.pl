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
# bowtie2_bs_align.pl
#
# Given a pair of Bowtie 2 (Bowtie 2, not Bowtie 1) indexes corresponding to the
# in-silico bisulfite-treated strands of the reference, and a set of reads,
# align the reads to the reference in a bisulfite-aware, unbiased manner.
#
# Input        Convert     -----------  Align   ------------    Merge SAM,
# (any format) (to fastq) | Converted |        |   Watson   |   extract
#   ----------  --------> |   Reads   | -----> | Alignments |-- evidence
#  | Unpaired |           -----------          ------------    \___________>
#  |  Reads   |            -----------          ------------    /
#   ----------  --------> | Converted | -----> |   Crick    |--
#                         |   Reads   |        | Alignments |
#                          -----------          ------------
#
# Intermediate reads files versus named pipes:
#
# We'd like to avoid making intermediate files, but we might need to make
# intermediate files for: (a) the input reads, so that they can be streamed
# separately through the Watson and Crick processes, (b) the in-silico
# bisulfite converted reads, and (c) the alignments, either pre- or
# post-merging, that come from the Watson and Crick aligners.  One option for
# cases (a) and (b) is to store all the reads in intermediate files, then
# proceed to the next step.
#
# Intermediate alignment files:
#
# We'd like to avoid creating large intermediate files, but we can't avoid
# making at least one intermediate .sam file if we run the Watson/Crick bowtie2
# processes in series.  If we run them in parallel, we can avoid the
# intermediate .sam file by aggregating output from the two 'bowtie2's and
# immediately emitting evidence records, without storing alignments.
#
# Handling 2/4 strands & pairs:
#
# 2 strands (i.e. the read can be from BSW or BSC) is straightforward:
# in-silico BS-converted BSW and BSC strands have their own indexes, and each
# read is bisulfite converted (BS) once and aligned to BSREFW and BSREFC with
# --norc.  Handling 4 strands (i.e. the read can be from BCW, BSC, BSWR or
# BSCR) is tricker.  It's NOT correct to simply take away --norc.  Rather, for
# reads that originate from one of the R strands, we have to take the reverse
# complement of the read, then in-silico bisulfite convert, then align with
# --norc.  Alternatively, we could just do a reverse-complement version of
# in-silico BS conversion where we turn Gs into As instead of turning Cs into
# Ts, then align with --nofw.  They're equivalent.
#
# Paired-end makes things a somewat trickier still.  Note the *fragment* can
# come from BSW, BSWR, BSC or BSCR, and the relative orientation of the two
# ends of the fragment can vary.  Consider the following case breakdown:
#
# Case 1. Fragment is sequenced --ff:
#  Case 1a. -ff and sequenced from BSW/BSC
#   Align ISC(M1), ISC(M2) to BSREFW/BSREFC with --ff
#  Case 1b. -ff and sequenced from BSWR/BSCR
#   Align ISC(RC(M2)), ISC(RC(M1)) to BSREFW/BSREFC with --ff
# Case 2. Fragment is sequenced --fr:
#  Case 2a. -fr and sequenced from BSW/BSC
#   Align ISC(M1), ISC(RC(M2)) to BSREFW/BSREFC with --ff
#  Case 2b. -fr and sequenced from BSWR/BSCR
#   Align ISC(M2), ISC(RC(M1)) to BSREFW/BSREFC with --ff
# Case 3. Fragment is sequenced --rf:
#  Case 3a. -rf and sequenced from BSW/BSC
#   Align ISC(RC(M1)), ISC(M2) to BSREFW/BSREFC with --ff
#  Case 3b. -rf and sequenced from BSWR/BSCR
#   Align ISC(RC(M2)), ISC(M1) to BSREFW/BSREFC with --ff
# Case 4. Fragment is sequenced --rr:
#  Case 4a. -rr and sequenced from BSW/BSC
#   Align ISC(RC(M1)), ISC(RC(M2)) to BSREFW/BSREFC with --ff
#  Case 4b. -rr and sequenced from BSWR/BSCR
#   Align ISC(M2), ISC(M1) to BSREFW/BSREFC with --ff
#

use strict;
use warnings;
use Getopt::Long;
use File::Path qw(rmtree mkpath);
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use BtlBio::Index::Bowtie2::Index;
use BtlBio::Alignment::Alignment;
use BtlBio::Alignment::SAM;
use BtlBio::Align::Bisulfite::BsAlign;
use BtlBio::Align::Bisulfite::BsEvidence;
use BtlBio::Read::Read;
use BtlBio::Read::Qseq;
use BtlBio::Read::Fastq;
use BtlBio::Read::Fasta;
use BtlBio::Util::File;
use BtlBio::Util::Tool;
use BtlBio::Util::Temp;
use BtlBio::Format::Fasta;
use BtlBio::Alphabet::DNA;
use POSIX;
use Storable qw(dclone);
use Carp;
use List::Util qw(max min);

my $ALIGNER = "bowtie2";        # aligner binary name
my $aligner_exe = $ALIGNER;     # path to aligner bin
my $samtools_exe = "samtools";  # path to samtools bin; only for --keep-bam
my @aligner_args = ();          # arguments to pass to aligner
my %ref = ();                   # all the reference sequences
my $bscpg = 1;                  # 1 -> do CpG->YpG conversion
my $bsc = 0;                    # 1 -> do C->Y conversion
my $temp_dir = $ENV{TMPDIR};    # temp dir for intermediates; reads, alns
$temp_dir = $temp_dir || "/tmp";# temp dir for intermediates; reads, alns
my $idx_name = "";              # aligner index
my $sam_fn = undef;             # put Watson/Crick SAM results here
my $bam_fn = undef;             # put Watson/Crick BAM results here
my $format = "fastq";           # format of input files
my $simul_align = 1;            # run 2 aligners in tandam for watson/crick
my $use_named_pipe = 1;         # use named pipes, not intermediate aln files
my $echo_sam = 0;               # print out SAM as it comes in
my $ev_out_dir = undef;         # directory to write evidence
my $strand4 = 0;                # using 4-strand protocol?
my $keep_dir = undef;           # where to keep intermediate reads/alns
my $stop_after_aln = 0;         # stop after alignment step, leaving .sam/.bam
my $m1fw = 1;                   # first mate aligns forward
my $m2fw = 0;                   # second mate aligns backward
my $stop_after = 0;             # stop after this many reads
my @name_map_fn = ();           # list of files with long/short name maps
my %name_map = ();              # map from long to short names
my $metrics_fn = "";            # metrics basename; .crick.tsv/.watson.tsv added
my $compress = "";              # manner of compression ev output
my $no_unpaired = 0;            # disallow unpaired evidence when input is paired?
my $stderr_pref = undef;        # prefix for stderr files
my $skip_reads = 0;             # skip this many reads off the bat
my @mate1_fns = ();             # read files with unpaired reads/mate #1s
my @mate2_fns = ();             # read files with mate #2s
my @ref_fns = ();               # fasta files with references

my $usage = qq!
Usage:
  bowtie2_bs_align.pl [options*] \
    -- <index> -- <refs> -- <bowtie2-args> -- <mate1s/reads> -- [<mate2s>]

Align reads to reference with Bowtie 2 and using in-silico bisulfite conversion
to impart bisulfite awareness and remove methylation bias.  Cs in the genome
are treated as Ys, matching Cs or Ts with no penalty.  Bowtie 2 is run twice:
once aligning to the bisulfite-treated Watson strand and once aligning to the
bisulfite-treated Crick strand.  The two Bowtie 2 processes run simultaneously.
Output is a directory of read-level measurement files, which can be sorted with
with the bsev_sort.pl script; the sorted directory can then be tabulated with
the bsev_tabulate.pl script.  Can also output alignments.

Required arguments:
  
  Required arguments are separated from options list and from each other by
  double dashes (--).  This allows more flexible use of wildcards (e.g. *.fa).
  
  <index>
    Name prefix for index files.  E.g. if index is in /tmp/idx.*.watson.bt2
    and /tmp/idx.*.crick.bt2, specify "/tmp/idx"
  <refs>
    Fasta files containing all original reference sequences
  <bowtie2-args>
    Arguments to pass to both the Watson and Crick Bowtie 2 jobs.  Do not
    specify options that affect Bowtie 2's (a) output format, (b) permitted
    alignment orientations (i.e. --ff, --fr, --rf), (c) metrics output (--met,
    --met-file).  These are set by Bsmooth and overriding them could cause
    problems.
  <mate1s/reads>
    Files containing unpaired reads or mate #1s if reads are paired.
  [<mate2s>]
    Files containing mate #2s if reads are paired, ignored otherwise.  If reads
    are paired, #1 files and #2 files must be specified in same order.

Options (defaults in parentheses):
 
 Alignment:
  --bowtie2=<path>         Path to bowtie2 executable (def: from PATH)
  --four-strand            Align as though read/fragment could have originated
                           from BSW, BSWR, BSC, BSCR (def: just BSW, BSC)
  --metrics=<path>         Send Bowtie 2 perf metrics to <path>.<strand>.tsv
  
  --bscpg                  Extract read-level measurements from CpG Cs (on)
   OR:
  --bsc                    Extract read-level measurements from all Cs (off)
 
 Paired-end:
  --fr/--ff/--rf/--rr      For paired-end: mates are in FR/FF/RF/RR orientation
  --no-unpaired            Ignore read-level measurements from unpaired mate
                           alignments
 
 Input:
  --fastq                  Reads files are FASTQ format (default)
  --qseq                   Reads files are in Illumina _qseq.txt format
  --fasta                  Reads files are in FASTA format
  -s/--skip-reads <int>    Skip first <int> input reads
  -u/--stop-after <int>    Only align first <int> input reads
 
 Output:
  --out=<path>             Directory to store read-level measurement records
                           extracted from the alignments.
  --sam=<prefix>           Store Watson alignment .sam in <prefix>.watson.sam,
                           Crick in <prefix>.crick.sam, unaligned in
                           <prefix>.unaligned.sam
  --bam=<prefix>           Like --keep-sam but stores .bam; requires samtools
  
  --stderr=<prefix>        Redirect aligner stderr to <prefix>.watson.stderr/
                           <prefix>.crick.stderr
  --gzip                   gzip compress read-level measurement output
  --bzip2                  bzip2 compress read-level measurement output
 
 Misc:
  --name-map=<path>        Many-to-one map from raw reference names to bin
                           names.  Output will be binned by bin name.
  --samtools=<path>        samtools binary for converting SAM to BAM when
                           --keep-bam is used
  --stop-after-alignment   Don't output read-level measurements
  --temp=<path>            Set path used for temporary FASTA files (TMPDIR)
  --keep=<path>            Keep intermediate/temporary files at <path>
  --help/--usage           Print this message and quit

Output:

 Read-level measurements:

  When --out is specified, read-level measurements (RLM) are extracted and
  output.  In these output files, each RLM is on a separate line, with the
  following tab-separated fields.

   1.  Reference name
   2.  Reference offset (1-based)
   3.  Read name (from which RLM was extracted)
   4.  Allele (A/C/G/T)
   5.  Strand: 1=Watson, 0=Crick
   6.  Strand of strand: 1=Forward, 2=Reverse complement
   7.  Base quality 1
   8.  Base quality 2
   9.  Read position, 0-based offset from 5' end
   10. Read length
   11. Alignment score (AS:i field from alignment)
   12. Mapping quality (MAPQ) of alignment

 Alignments:
 
  When --keep-sam and/or --keep-bam are used to store alignment output from
  BSmooth.  Alignments are divided into three output files, one for alignments
  to the Watson strand, one for alignments to the Crick strand, and one for
  reads that did not align.  These output files are in SAM (or BAM) format.
  BSmooth adds one (for unaligned reads) or two (for aligned reads) extra
  fields to convey additional information:
  
   XB:Z:   This flag is present for aligned reads only.  It describes which of
           the four possible bisulfite strands the read aligned to.
           W = Watson, C = Crick, WR = reverse complement of Watson, CR =
           reverse complement of Crick.  Note: WR and C do not have the same
           sequence when bisulfite treatmet is used\!
   YO:Z:   The original read sequence.  The SAM SEQ field contains the
           C-to-T-converted sequence.

!;

##
# Print and run given command.  Die if it returns non-zero.
#
my %signo = ();
my @signame = ();
{
	# Get signal info
	use Config;
	my $i = 0;
	for my $name (split(' ', $Config{sig_name})) {
		$signo{$name} = $i;
		$signame[$i] = $name;
		$i++;
	}
}

sub run($) {
	my $cmd = shift;
	print STDERR "Running '$cmd'...\n";
	system($cmd);
	if($? == -1) {
	    confess("Failed to execute '$cmd': $!");
	} elsif($? & 127) {
		my $signm = "(unknown)";
		$signm = $signame[$? & 127] if defined($signame[$? & 127]);
		my $ad = "";
		$ad = "(core dumped)" if (($? & 128) != 0);
	    confess(sprintf "'$cmd' died with signal %d (%s) $ad", ($? & 127), $signm);
	} elsif($? != 0) {
	    confess(sprintf "'$cmd' exited with value %d", $? >> 8);
	}
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

$SIG{INT} = sub { confess(); };

my $index_ext = "bt2";

##
# Remove any trailing watson/crick.*.<ext>
#
sub strip_idx_name($) {
	my $idx = shift;
	$idx =~ s/\.rev\.[1-4]\.$index_ext$//;
	$idx =~ s/\.[1-4]\.$index_ext$//;
	$idx =~ s/\.watson$//;
	$idx =~ s/\.crick$//;
	return $idx;
}

GetOptions (
	"index=s"         => \$idx_name,
	"bowtie2=s"       => \$aligner_exe,
	"samtools=s"      => \$samtools_exe,
	"echo-sam"        => \$echo_sam,
	"bowtie2-args=s"  => \@aligner_args,
	"aligner-args=s"  => \@aligner_args,
	"no-named-pipes"  => sub { $use_named_pipe = 0; },
	"serial"          => sub { $simul_align = 0; $use_named_pipe = 0; },
	"keep-sam|sam=s"  => \$sam_fn,
	"keep-bam|bam=s"  => \$bam_fn,
	"bscpg"           => sub { $bscpg = 1; $bsc = 0; },
	"bsc"             => sub { $bscpg = 0; $bsc = 1; },
	"read=s"          => \@mate1_fns,
	"mate1=s"         => \@mate1_fns,
	"mate2=s"         => \@mate2_fns,
	"four-strand"     => \$strand4,
	"fastq"           => sub { $format = "fastq" },
	"qseq"            => sub { $format = "qseq"  },
	"fasta"           => sub { $format = "fasta" },
	"temp-directory|temporary-directory=s" => \$temp_dir,
	"output-directory|evidence-directory=s" => \$ev_out_dir,
	"keep=s"          => \$keep_dir,
	"stop-after-alignment" => \$stop_after_aln,
	"s|skip-reads=i"  => \$skip_reads,
	"u|stop-after=i"  => \$stop_after,
	"name-map=s"      => \@name_map_fn,
	"metrics=s"       => \$metrics_fn,
	"compress=s"      => \$compress,
	"stderr-prefix=s" => \$stderr_pref,
	"gzip"            => sub { $compress = "gzip";  },
	"bzip2"           => sub { $compress = "bzip2"; },
	"ff"              => sub { $m1fw = 1; $m2fw = 1; },
	"rr"              => sub { $m1fw = 0; $m2fw = 0; },
	"fr"              => sub { $m1fw = 1; $m2fw = 0; },
	"rf"              => sub { $m1fw = 0; $m2fw = 1; },
	"no-unpaired"     => \$no_unpaired,
	"help|usage"      => sub { print $usage; exit 0; }
) || dieusage("Bad option", 0);

# Parse out arguments from among the double-dashes
my $dds = 0;
for my $a (@ARGV) {
	if($a eq "--") {
		$dds++; next;
	}
	if($dds == 0) {
		$idx_name = $a;
	} elsif($dds == 1) {
		push @ref_fns, $a;
	} elsif($dds == 2) {
		push @aligner_args, $a;
	} elsif($dds == 3) {
		push @mate1_fns, $a;
	} else {
		$dds > 0 || die "";
		push @mate2_fns, $a;
	}
}

if($idx_name eq "") {
	print "$usage\n---\nError: <index> not specified\n"; exit 1;
}
if(scalar(@ref_fns) == 0) {
	print "$usage\n---\nError: <refs> not specified\n"; exit 1;
}
if(scalar(@mate1_fns) == 0) {
	print "$usage\n---\nError: <mate1s/reads> not specified\n"; exit 1;
}
mkpath($ev_out_dir) if $ev_out_dir ne ".";
defined($sam_fn) || defined($bam_fn) || defined($ev_out_dir) ||
	die "No output types specified.  Please specify one or more of: ".
	    "--out, --sam, --bam.";

scalar(@mate1_fns) > 0 || dieusage("No reads specified", 1);
!defined($keep_dir) || ! -d $keep_dir ||
	die "--keep directory '$keep_dir' already exists; remove it or specify different --keep directory";
$temp_dir = $keep_dir || temp_name(undef, $temp_dir);
mkpath($temp_dir);

my $paired = scalar(@mate2_fns) > 0;
my ($inter_wat, $inter_cri) =
	("$temp_dir/bsreads_watson.tab6", "$temp_dir/bsreads_crick.tab6");

print STDERR "Original index name: $idx_name\n";
$idx_name = strip_idx_name($idx_name);
print STDERR "Stripped index name: $idx_name\n";
print STDERR "Reference files:\n";
for(0..$#ref_fns) {print STDERR "  $ref_fns[$_]\n"}
print STDERR "Name-map files:\n";
for(0..$#name_map_fn) {print STDERR "  $name_map_fn[$_]\n"}
print STDERR "$ALIGNER executable: $aligner_exe\n";
print STDERR "$ALIGNER arguments: ".join(" ", @aligner_args)."\n";
print STDERR "Methylation output directory: $ev_out_dir\n";
print STDERR "Keep Watson/Crick .sam: ".($sam_fn || "no")."\n";
print STDERR "Keep Watson/Crick .bam: ".($bam_fn || "no")."\n";
print STDERR "samtools executable override: $samtools_exe\n" if $bam_fn;
print STDERR "Read-level measurements: ".($bscpg ? "just CpGs" : "all Cs")."\n";
print STDERR "Four-strand?: ".($strand4 ? "yes":"no")."\n";
print STDERR "Paired-end?: ".($paired ? "yes" : "no")."\n";
if($paired) {
	print STDERR "Mate #1 inputs:\n";
	for(0..$#mate1_fns) {print STDERR "  $mate1_fns[$_]\n"}
	print STDERR "Mate #2 inputs:\n";
	for(0..$#mate2_fns) {print STDERR "  $mate2_fns[$_]\n"}
} else {
	print STDERR "Unpaired inputs:\n";
	for(0..$#mate1_fns) {print STDERR "  $mate1_fns[$_]\n"}
}
print STDERR "Ignore evidence from unparied alignment of paired mate: $no_unpaired\n";
print STDERR "Input format: $format\n";
print STDERR "Skip first # reads: ".($skip_reads > 0 ? $skip_reads : "(no skipping)")."\n";
print STDERR "Stop after # reads: ".($stop_after > 0 ? $stop_after : "(no limit)")."\n";
print STDERR "Temporary directory: $temp_dir\n";
print STDERR "Intermediate Watson file name: $inter_wat\n";
print STDERR "Intermediate Crick file name: $inter_cri\n";

print STDERR "Locating tools...\n";
$aligner_exe = find_tool($ALIGNER, $aligner_exe);
$samtools_exe = find_tool("samtools", $samtools_exe) if $bam_fn;
print STDERR "  Found $ALIGNER: $aligner_exe\n";
print STDERR "  Found samtools: $samtools_exe\n" if $bam_fn;

print STDERR "Checking index...\n";
$idx_name ne "" || dieusage("Must specify --index", 20);
# Is index present?
bowtie2_index_exists("$idx_name.watson") ||
	die "Index with basename does not exist: '$idx_name.watson'";
bowtie2_index_exists("$idx_name.crick") ||
	die "Index with basename does not exist: '$idx_name.crick'";

print STDERR "Checking reads...\n";
# Do --mate1 and --mate2 match up?
scalar(@mate1_fns) > 0 ||
	die "Must specify reads using --reads/--mate1/--mate2";
!$paired || scalar(@mate1_fns) == scalar(@mate2_fns) ||
	die "Must specify same number of #1/#2 mate files with --mate1/--mate2";
# Do all the files exist?
for(@mate1_fns) { $_ eq "-" || -f $_ || die "Could not find reads file '$_'"; }
for(@mate2_fns) { $_ eq "-" || -f $_ || die "Could not find reads file '$_'"; }

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

my $nreads = 0;

##
# Write in-silico-converted read to filehandle in "tabbed" format.
#
sub handle_read($$$$) {
	my ($read, $ofh, $rc, $mate) = @_;
	my ($name, $seq) = ($read->{_name}, $read->{_seq});
	$read = bsc_ize_read(dclone($read), $rc, $mate);
	print {$ofh} $read->to_tab6;
	delete $read->{_bsc_orig_seq};
	($read->{_name}, $read->{_seq}) = ($name, $seq);
}

##
# Write in-silico-converted pair to filehandle in "tabbed" format.
#
sub handle_pair($$$$$$$) {
	# confusing that we've labeled the variables "1" and "2" when the
	# user might have passed 2, 1 instead.
	my ($m1, $m2, $ofh, $rc1, $rc2, $mate1, $mate2) = @_;
	my ($name1, $seq1) = ($m1->{_name}, $m1->{_seq});
	my ($name2, $seq2) = ($m2->{_name}, $m2->{_seq});
	$m1 = bsc_ize_read(dclone($m1), $rc1, $mate1);
	$m2 = bsc_ize_read(dclone($m2), $rc2, $mate2);
	print {$ofh} reads_to_tab6($m1, $m2);
	delete $m1->{_bsc_orig_seq};
	delete $m2->{_bsc_orig_seq};
	($m1->{_name}, $m1->{_seq}) = ($name1, $seq1);
	($m2->{_name}, $m2->{_seq}) = ($name2, $seq2);
}

my $parse_func = \&parse_fastq_read;
$parse_func    = \&parse_qseq_read  if $format =~ /qseq/i;
$parse_func    = \&parse_fasta_read if $format =~ /fasta/i;

##
# Given a filename and a pair of output streams, one for sending converted
# reads to the Watson aligner and one for the Crick aligner, convert all reads
# in the file and send them off.  This function assumes the two aligners are
# running simultaneosuly.
#
sub cat_file_unpaired($$$) {
	my ($ifn, $ofh_wat, $ofh_cri) = @_;
	return if $stop_after > 0 && $nreads >= $stop_after;
	print STDERR "  Writing reads from file '$ifn' ...\n";
	my $ifh;
	if($ifn eq "-") { $ifh = *STDIN; } else { $ifh = openex($ifn); }
	while(my $read = $parse_func->($ifh)) {
		if($skip_reads > 0) { $skip_reads--; next; }
		handle_read($read, $ofh_wat, 0, 0);
		handle_read($read, $ofh_cri, 0, 0);
		if($strand4) {
			handle_read($read, $ofh_wat, 1, 0);
			handle_read($read, $ofh_cri, 1, 0);
		}
		$nreads++;
		last if $stop_after > 0 && $nreads >= $stop_after;
	}
	close($ifh) unless $ifn eq "-";
}

##
# Given a pair of filenames and a set of four output streams, two for sending
# pairs to the Watson aligner and two for sending pairs to the Crick aligner,
# convert all pairs in the files and send them off.  This function assumes the
# two aligners are running simultaneosuly.
#
sub cat_file_paired($$$$$$) {
	my (
		$ifn1,      # file with mate #1s
		$ifn2,      # file with mate #2s
		$m1fw,      # true iff fragments are --ff or --fr
		$m2fw,      # true iff fragments are --rf or --rr
		$ofh_wat,   # output fh for mate #1s sent to Watson aligner
		$ofh_cri)   # output fh for mate #1s sent to Crick aligner
		= @_;
	
	return if $stop_after > 0 && $nreads >= $stop_after;
	print STDERR "  Writing reads from files '$ifn1'/'$ifn2' ...\n";
	my ($ifh1, $ifh2) = (undef, undef);
	if($ifn1 eq "-") { $ifh1 = *STDIN; } else { $ifh1 = openex($ifn1); }
	if($ifn2 eq "-") { $ifh2 = *STDIN; } else { $ifh2 = openex($ifn2); }
	while(1) {
		my $read1 = $parse_func->($ifh1);
		my $read2 = $parse_func->($ifh2);
		defined($read1) == defined($read2) || confess();
		last unless defined($read1);
		if($skip_reads > 0) { $skip_reads--; next; }
		handle_pair($read1, $read2, $ofh_wat, $m1fw ? 0:1, $m2fw ? 0:1, 1, 2);
		handle_pair($read1, $read2, $ofh_cri, $m1fw ? 0:1, $m2fw ? 0:1, 1, 2);
		if($strand4) {
			handle_pair($read2, $read1, $ofh_wat, $m2fw, $m1fw, 2, 1);
			handle_pair($read2, $read1, $ofh_cri, $m2fw, $m1fw, 2, 1);
		}
		$nreads++;
		last if $stop_after > 0 && $nreads >= $stop_after;
	}
	close($ifh1) unless $ifn1 eq "-";
	close($ifh2) unless $ifn2 eq "-";
}

# We can run two instances of Bowtie simultaneously, or we can run them one
# after the other.
if($use_named_pipe) {
	print STDERR "Creating named pipes...\n";
	mkfifo($inter_wat, 0700);
	mkfifo($inter_cri, 0700);
}

# fork off a process whose job it is to feed Bowtie
my $pid_feed = 0;
$pid_feed = fork() if $use_named_pipe;
if($pid_feed == 0) {
	# Open named pipes for writing to Watson and Crick processes
	open(my $ofh_wat, ">$inter_wat") || die "Can't open '$inter_wat' for writing";
	open(my $ofh_cri, ">$inter_cri") || die "Can't open '$inter_cri' for writing";
	for (0..$#mate1_fns) {
		my ($m1fn, $m2fn) = ($mate1_fns[$_], $mate2_fns[$_]);
		cat_file_paired($m1fn, $m2fn, $m1fw, $m2fw, $ofh_wat, $ofh_cri) if $paired;
		cat_file_unpaired($m1fn, $ofh_wat, $ofh_cri) if !$paired;
	}
	for($ofh_wat, $ofh_cri) { close($_); }
	exit 0 if $use_named_pipe;
}
print STDERR "Forked off paired feeder process (pid=$pid_feed)...\n" if $pid_feed != 0;

my ($al_wat_fn, $al_cri_fn) = (undef, undef);
my ($al_wat_fh, $al_cri_fh) = (undef, undef);

# Open a bowtie process that reads from the intermediate files/pipes and
# outputs alignments in .sam format.

my $aligner_req_args = "--quiet --ff --norc --sam-no-qname-trunc --reorder";

my ($stderr_redir_wat, $stderr_redir_cri) = ("", "");
my $cmd_wat = "$aligner_exe -x $idx_name.watson ".join(" ", @aligner_args)." $aligner_req_args --tab6 $inter_wat";
my $cmd_cri = "$aligner_exe -x $idx_name.crick " .join(" ", @aligner_args)." $aligner_req_args --tab6 $inter_cri";
if($metrics_fn ne "") {
	$cmd_wat .= " --metrics=1 --metrics-file=$metrics_fn.watson.tsv";
	$cmd_cri .= " --metrics=1 --metrics-file=$metrics_fn.crick.tsv";
}
if(defined($stderr_pref)) {
	$stderr_redir_wat = "2>$stderr_pref.watson.stderr";
	$stderr_redir_cri = "2>$stderr_pref.crick.stderr";
	$cmd_wat .= " $stderr_redir_wat";
	$cmd_cri .= " $stderr_redir_cri";
}

if($simul_align) {
	# We're going to run two aligners simultaneously
	$cmd_wat .= " |"; # Will read your output directly
	$cmd_cri .= " |"; # Yours too
	print STDERR "Opening simultaneous alignment pipes:\n";
	print STDERR "  $cmd_wat\n$cmd_cri\n";
	open($al_wat_fh, $cmd_wat) || die "Can't open pipe '$cmd_wat' for reading";
	open($al_cri_fh, $cmd_cri) || die "Can't open pipe '$cmd_cri' for reading";
} else {
	# Set $al_wat_fn and $al_cri_fn to temporary filenames
	($al_wat_fn, $al_cri_fn) =
		("$temp_dir/als_watson.sam", "$temp_dir/als_crick.sam");
	# We're going to run the aligners one after the other
	$cmd_wat .= " -S $al_wat_fn"; # Will read output later
	$cmd_cri .= " -S $al_cri_fn"; # Will read output later
	print STDERR "Running Watson aligner...\n";
	run($cmd_wat); # Run aligner on BSW/BSWR strands
	print STDERR "Running Crick aligner...\n";
	run($cmd_cri); # Run aligner on BSC/BSCR strands
	open($al_wat_fh, $al_wat_fn) || die "Can't open '$al_wat_fn' for reading";
	open($al_cri_fh, $al_cri_fn) || die "Can't open '$al_cri_fn' for reading";
}

my ($sam_fh_wat, $sam_fh_cri, $sam_fh_un) = (undef, undef, undef);
my ($bam_fh_wat, $bam_fh_cri, $bam_fh_un) = (undef, undef, undef);
if(defined($sam_fn)) {
	# Open .sam output pipes
	my ($sam_fn_wat, $sam_fn_cri, $sam_fn_un) =
		("$sam_fn.watson.sam", "$sam_fn.crick.sam", "$sam_fn.unaligned.sam");
	open($sam_fh_wat, ">$sam_fn_wat") || die "Can't open '$sam_fn_wat' for writing";
	open($sam_fh_cri, ">$sam_fn_cri") || die "Can't open '$sam_fn_cri' for writing";
	open($sam_fh_un,  ">$sam_fn_un")  || die "Can't open '$sam_fn_un' for writing";
}
if(defined($bam_fn)) {
	# Open .bam output pipes
	my ($bam_fn_wat, $bam_fn_cri, $bam_fn_un) =
		("$bam_fn.watson.bam", "$bam_fn.crick.bam", "$bam_fn.unaligned.bam");
	open($bam_fh_wat, "| $samtools_exe view -Sb - > $bam_fn_wat") || die "Can't open '$bam_fn_wat' for writing";
	open($bam_fh_cri, "| $samtools_exe view -Sb - > $bam_fn_cri") || die "Can't open '$bam_fn_cri' for writing";
	open($bam_fh_un,  "| $samtools_exe view -Sb - > $bam_fn_un")  || die "Can't open '$bam_fn_un' for writing";
}

# Now we have filehandles from which to read the two streams of alignments,
# one with alignments to BSW/BSWR, one with alignments to BSC/BSCR
my $nskip_unal = 0;    # # times there were 0 alignments
my $nskip_mult = 0;    # # times there was >1 alignment
my $nheads = 0;        # # header lines
my $nrecs = 0;         # # reads (pairs or ends, as appropriate)
my $nrecs_p = 0;       # # ends
my $nrec_ival = 5000;  # print msg every this many records
my %ev_out_fhs = ();   # filehandles for evidence files

##
# Given a piece of methylation evidence, simply print it to STDOUT.
#
my $nev = 0;
# Mark the evidence directory with a file indicating what type of evidence is
# there (either CpG-level or C-level)
if(defined($ev_out_dir)) {
	open(MARK, ">$ev_out_dir/.".($bscpg ? "cpg" : "c").".ev") || die;
	print MARK "\n"; close(MARK);
}

sub evidenceSink($) {
	my $ev = shift;
	$nev++;
	# Bin by reference name
	my $tname = $ev->{_tname};
	my $full_name = $tname;
	$tname =~ s/\s.*//;
	if(defined($name_map{$tname})) {
		$tname = $name_map{$tname};
		$full_name = $tname;
	}
	$tname = sanitize_filename($tname);
	my $outfn = "$ev_out_dir/$tname.ev.tsv";
	# Open output file for this reference, if there isn't one open already
	if(!defined($ev_out_fhs{$outfn})) {
		my $outpipe = ">$outfn";
		$outpipe = "| gzip  -c > $outfn.gz"  if $compress eq "gzip";
		$outpipe = "| bzip2 -c > $outfn.bz2" if $compress eq "bzip2";
		open($ev_out_fhs{$outfn}, $outpipe) ||
			confess("Could not open output pipe '$outpipe' for writing");
		my $outbinfn = "$ev_out_dir/.$tname.ev.tsv.bin";
		open(TMP, ">$outbinfn") ||
			confess("Could not open '$outbinfn' for writing");
		print TMP "$full_name\n"; # original, not sanitized
		close(TMP);
	}
	# Print the record
	print {$ev_out_fhs{$outfn}} $ev->to_record."\n";
}
my $ev_sink = undef;
$ev_sink = \&evidenceSink if defined($ev_out_dir);

sub same_name($$) {
	my ($al1, $al2) = @_;
	my ($nm1, $nm2) = ($al1->{_rdname}, $al2->{_rdname});
	$nm1 =~ s/:.*$//;
	$nm2 =~ s/:.*$//;
	$nm1 eq $nm2 || die;
}

sub order_mates {
	my ($a, $b, $m1fw, $m2fw, $rc) = @_;
	my ($m1, $m2);
	if(substr($a->{_rdname}, -2, 1) eq "1") {
		($m1, $m2) = ($a, $b);
	} else {
		($m1, $m2) = ($b, $a);
	}
	# Now m1 has mate 1 and m2 has mate 2
	if(defined($rc)) {
		# Sanity check the orientations of the mates
		$m1fw = ($rc ? !$m1fw : $m1fw);
		$m2fw = ($rc ? !$m2fw : $m2fw);
		$m1fw = ($m1fw ? "+" : "-");
		$m2fw = ($m2fw ? "+" : "-");
		substr($m1->{_rdname}, -1) eq $m1fw ||
			confess("Expected '$m1fw' (m1fw=$m1fw, m2fw=$m2fw, rc=$rc), got $m1->{_rdname}");
		substr($m2->{_rdname}, -1) eq $m2fw ||
			confess("Expected '$m2fw' (m1fw=$m1fw, m2fw=$m2fw, rc=$rc), got $m2->{_rdname}");
	}
	return ($m1, $m2);
}

sub next_lines($$) {
	my ($al_wat_fh, $al_cri_fh) = @_;
	my $line_wat = readline $al_wat_fh;
	my $line_cri = readline $al_cri_fh;
	defined($line_wat) == defined($line_cri) ||
		confess("SAM output from Watson/Crick ended on different lines\n".
		        "Watson: ".(defined($line_wat) ? $line_wat : "(undefined)")."\n".
		        "Crick:  ".(defined($line_cri) ? $line_cri : "(undefined)")."\n");
	if(defined($line_wat)) {
		print $line_wat if $echo_sam;
		print $line_cri if $echo_sam;
	}
	return ($line_wat, $line_cri);
}

##
# Given an array ref pointing to a list of alignments, or alignment pairs,
# return an alignment with maximal MAPQ.  Break ties randomly.
#
sub select_al($$$) {
	my ($als, $wats, $strand4) = @_;
	scalar(@$als) > 0 || die;
	my $paired = scalar(@{$als->[0]}) > 1;
	if(scalar(@$als) == 1) {
		my $new_mapq1 = $als->[0]->[0]->{_mapq};
		my $new_mapq2 = undef;
		$new_mapq2 = $als->[0]->[1]->{_mapq} if $paired;
		return (0, $new_mapq1, $new_mapq2);
	}
	my @max_mapq_l = ();
	my @strands = ();
	my $max_mapq = -1;
	my $n = 0;
	for my $i (0..scalar(@$als)-1) {
		$n++;
		my $mapq = $als->[$i]->[0]->{_mapq};
		$mapq   += $als->[$i]->[1]->{_mapq} if $paired;
		if($mapq > $max_mapq) {
			@max_mapq_l = ($i);
			$max_mapq = $mapq;
		} elsif($mapq == $max_mapq) {
			# Tied - we'll select randomly later
			push @max_mapq_l, $i;
		}
	}
	my $seli = int(rand(scalar(@max_mapq_l)));
	my $new_mapq1 = $als->[$seli]->[0]->{_mapq};
	my $new_mapq2 = undef;
	$new_mapq2 = $als->[$seli]->[1]->{_mapq} if $paired;
	if($strand4 && $n == 2) {
		# Check for the special case where the read aligns to W and CR or to C
		# and WR, in which case the ambiguity may be due to the read being
		# highly methylated.  We could further confirm this by seeing if the
		# two alignments are in corresponding spots.  We might report "strand
		# MAPQ" as a separate number.
		my $fw1 = (substr($als->[0]->[0]->{_rdname}, -1) eq "+") ? 1 :0;
		my $fw2 = (substr($als->[1]->[0]->{_rdname}, -1) eq "+") ? 1 : 0;
		if($fw1 != $fw2 && $wats->[0] != $wats->[1]) {
			$new_mapq1 = int($als->[$seli]->[0]->{_mapq} * 2.0 / 3.0);
			$new_mapq2 = int($als->[$seli]->[1]->{_mapq} * 2.0 / 3.0) if $paired;
		} else {
			$new_mapq1 = 0;
			$new_mapq2 = 0 if $paired;
		}
	} elsif($n > 1) {
		$new_mapq1 = 0;
		$new_mapq2 = 0 if $paired;
	}
	(defined($new_mapq1) && $new_mapq1 == int($new_mapq1)) || die;
	if($paired) {
		(defined($new_mapq2) && $new_mapq2 == int($new_mapq2)) || die;
	}
	return ($max_mapq_l[$seli], $new_mapq1, $new_mapq2);
}

##
# Given a sam alignment string, compose a new string to be printed to the sam/bam output.  Includes key extra information.
#
sub prepare_sam($$$) {
	my ($line, $aligned, $is_watson) = @_;
	$line =~ s/\s+$//;       # Remove trailing whitespace
	my @ts = split(/\t/, $line);
	$ts[0] =~ /:([^:\t]*)$/; # Grab stuff we added to read name
	my $flag = undef;
	my $add = $1;
	defined($add) || croak("Couldn't get post-colon material from read name:\n$line\n");
	my $orig = substr($add, 0, length($add)-2);
	my $r = substr($add, -1);
	if($aligned) {
		$flag = "XB:Z:" . ($is_watson ? "W" : "C");
		$flag .= "R" if $r eq "-";
	}
	$ts[0] =~ s/:[^:]*$//; # Remove crud we added to read name
	my $ret = join("\t", @ts);
	$ret .= "\tYO:Z:$orig" if length($orig) > 0;
	$ret .= "\t$flag" if $aligned;
	return "$ret\n";
}

print STDERR "Parsing reference...\n";
for(@ref_fns) { read_fasta_to_hash($_, \%ref, { truncate_names => 1 }); }
print STDERR "  Parsed ".scalar(keys %ref)." sequences\n";

print STDERR "Processing SAM output from aligners...\n";
my ($last_al_wat, $last_al_cri) = (undef, undef);
my $iters = 0;
my ($naltot, $nmaxtot, $nunaltot) = (0, 0, 0);
my %mapq_hist = ();
my $nttot = 0;
my $t0 = time();
while(1) {
	$iters++;
	my $nnt = 0;
	# Output is in SAM format
	my ($line_wat, $line_cri) = next_lines($al_wat_fh, $al_cri_fh);
	last unless defined($line_wat);
	my $firstc_wat = substr($line_wat, 0, 1);
	my $firstc_cri = substr($line_cri, 0, 1);
	if($firstc_wat eq "\@") {
		# Header line
		$nheads++;
		$firstc_cri eq "\@" || die;
		# Send header lines to sam/bam
		print {$sam_fh_wat} $line_wat if defined($sam_fh_wat);
		print {$bam_fh_wat} $line_wat if defined($bam_fh_wat);
		print {$sam_fh_cri} $line_cri if defined($sam_fh_cri);
		print {$bam_fh_cri} $line_cri if defined($bam_fh_cri);
		print {$sam_fh_un}  $line_wat if defined($sam_fh_un);
		print {$bam_fh_un}  $line_wat if defined($bam_fh_un);
		next; # skip
	}
	if(!$stop_after_aln) {
		my ($nal, $nunal) = (0, 0); # reads with/without alignments
		my $al_paired = 0;  # true iff a paired alignment was found
		my $al = undef;     # the alignment/pair selected
		my @als = ();       # alignments found
		my @unals = ();     # records for unaligned
		my @lines = ();     # SAM lines
		my @unlines = ();   # SAM lines for unaligned
		my @is_watson = (); # 0/1 = Crick/Watson
		if(!$strand4 && !$paired) {
			my $al_wat = parse_sam_record($line_wat);
			my $al_cri = parse_sam_record($line_cri);
			$nnt = length($al_wat->{_rdseq});
			$al_wat->{_rdname} eq $al_cri->{_rdname} || die;
			if($al_wat->aligned) {
				push @als, [ $al_wat ];
				push @lines, [ $line_wat ];
				push @is_watson, 1;
				$nal++;
			}
			if($al_cri->aligned) {
				push @als, [ $al_cri ];
				push @lines, [ $line_cri ];
				push @is_watson, 0;
				$nal++;
			}
			if($nal == 0) {
				push @unals, $al_wat;
				push @unlines, $line_wat;
				$nunal++;
			}
			if($nal > 0) {
				# Select alignment with greatest MAPQ
				my ($seli, $nmapq1, $nmapq2) = select_al(\@als, \@is_watson, $strand4);
				$al = $als[$seli];
				$nmapq1 <= $al->[0]{_mapq} || die;
				$al->[0]{_mapq} = $nmapq1;
				if($bsc) {
					bsc_analyze_alignment($al->[0], $is_watson[$seli], \%ref, $ev_sink, undef, undef);
				} else {
					bscpg_analyze_alignment($al->[0], $is_watson[$seli], \%ref, $ev_sink, undef, undef);
				}
				my $rec = prepare_sam($lines[$seli]->[0], 1, $is_watson[$seli]);
				if($is_watson[$seli]) {
					print {$sam_fh_wat} $rec if defined($sam_fh_wat);
					print {$bam_fh_wat} $rec if defined($bam_fh_wat);
				} else {
					print {$sam_fh_cri} $rec if defined($sam_fh_cri);
					print {$bam_fh_cri} $rec if defined($bam_fh_cri);
				}
			}
			if($nunal > 0) {
				my $unline = $unlines[0];
				my $rec = prepare_sam($unline, 0, undef);
				print {$sam_fh_un} $rec if defined($sam_fh_un);
				print {$bam_fh_un} $rec if defined($bam_fh_un);
			}
		}
		elsif($strand4 && !$paired) {
			my $al_wat_fw = parse_sam_record($line_wat);
			my $al_cri_fw = parse_sam_record($line_cri);
			my ($line_wat_rc, $line_cri_rc) = next_lines($al_wat_fh, $al_cri_fh);
			my $al_wat_rc = parse_sam_record($line_wat_rc);
			my $al_cri_rc = parse_sam_record($line_cri_rc);
			$nnt = (length($al_wat_fw->{_rdseq}) + length($al_wat_rc->{_rdseq}));
			$al_wat_fw->{_rdname} eq $al_cri_fw->{_rdname} || die;
			$al_wat_rc->{_rdname} eq $al_cri_rc->{_rdname} || die;
			same_name($al_wat_fw, $al_wat_rc) || die;
			same_name($al_cri_fw, $al_cri_rc) || die;
			if($al_wat_fw->aligned) {
				push @als, [ $al_wat_fw ];
				push @lines, [ $line_wat ];
				push @is_watson, 1;
				$nal++;
			}
			if($al_wat_rc->aligned) {
				push @als, [ $al_wat_rc ];
				push @lines, [ $line_wat_rc ];
				push @is_watson, 1;
				$nal++;
			}
			if($al_cri_fw->aligned) {
				push @als, [ $al_cri_fw ];
				push @lines, [ $line_cri ];
				push @is_watson, 0;
				$nal++;
			}
			if($al_cri_rc->aligned) {
				push @als, [ $al_cri_rc ];
				push @lines, [ $line_cri_rc ];
				push @is_watson, 0;
				$nal++;
			}
			if($nal == 0) {
				push @unals, $al_wat_fw;
				push @unlines, $line_wat;
				$nunal++;
			}
			if($nal > 0) {
				# Select alignment with greatest MAPQ
				my ($seli, $nmapq1, $nmapq2) = select_al(\@als, \@is_watson, $strand4);
				$al = $als[$seli];
				# TODO: Be careful here: when the read aligns both to W and CR,
				# or both to C and WR, and the MAPQs are high for both
				# alignments, we may be looking at a highly methylated read
				# that doesn't deserve a low MAPQ.
				$nmapq1 <= $al->[0]{_mapq} || die;
				$al->[0]{_mapq} = $nmapq1;
				if($bsc) {
					bsc_analyze_alignment($al->[0], $is_watson[$seli], \%ref, $ev_sink, undef, undef);
				} else {
					bscpg_analyze_alignment($al->[0], $is_watson[$seli], \%ref, $ev_sink, undef, undef);
				}
				my $rec = prepare_sam($lines[$seli]->[0], 1, $is_watson[$seli]);
				if($is_watson[$seli]) {
					print {$sam_fh_wat} $rec if defined($sam_fh_wat);
					print {$bam_fh_wat} $rec if defined($bam_fh_wat);
				} else {
					print {$sam_fh_cri} $rec if defined($sam_fh_cri);
					print {$bam_fh_cri} $rec if defined($bam_fh_cri);
				}
			}
			if($nunal > 0) {
				my $unline = $unlines[0];
				my $rec = prepare_sam($unline, 0, undef);
				print {$sam_fh_un} $rec if defined($sam_fh_un);
				print {$bam_fh_un} $rec if defined($bam_fh_un);
			}
		}
		elsif(!$strand4 && $paired) {
			my $al_wat1 = parse_sam_record($line_wat);
			my $al_cri1 = parse_sam_record($line_cri);
			my ($line_wat2, $line_cri2) = next_lines($al_wat_fh, $al_cri_fh);
			my $al_wat2 = parse_sam_record($line_wat2);
			my $al_cri2 = parse_sam_record($line_cri2);
			$nnt = (length($al_wat1->{_rdseq}) + length($al_wat2->{_rdseq}));
			($al_wat1, $al_wat2) = order_mates($al_wat1, $al_wat2, $m1fw, $m2fw, 0);
			($al_cri1, $al_cri2) = order_mates($al_cri1, $al_cri2, $m1fw, $m2fw, 0);
			# Try committing paired-end evidence
			my ($nal1, $nal2) = (0, 0);
			if($al_wat1->aligned && $al_wat2->aligned) {
				push @als, [ $al_wat1, $al_wat2 ];
				push @lines, [ $line_wat, $line_wat2 ];
				push @is_watson, 1;
				$nal++; $nal1++; $nal2++;
			}
			if($al_cri1->aligned && $al_cri2->aligned) {
				push @als, [ $al_cri1, $al_cri2 ];
				push @lines, [ $line_cri, $line_cri2 ];
				push @is_watson, 0;
				$nal++; $nal1++; $nal2++;
			}
			if($nal == 0 && !$no_unpaired) {
				if($al_wat1->aligned || $al_wat2->aligned) {
					my $alx = $al_wat1->aligned ? $al_wat1 : $al_wat2;
					my $linex = $al_wat1->aligned ? $line_wat : $line_wat2;
					!$alx->mate_aligned || croak ("Mate aligned in unpaired alignment:\n$linex");
					push @als, [ $alx ];
					push @lines, [ $linex ];
					push @is_watson, 1;
					$nal++;
				}
				if($al_cri1->aligned || $al_cri2->aligned) {
					my $alx = $al_cri1->aligned ? $al_cri1 : $al_cri2;
					my $linex = $al_cri1->aligned ? $line_cri : $line_cri2;
					!$alx->mate_aligned || croak ("Mate aligned in unpaired alignment:\n$linex");
					push @als, [ $alx ];
					push @lines, [ $linex ];
					push @is_watson, 0;
					$nal++;
				}
				# See if only one or the other aligned
				if($al_wat1->aligned || $al_cri1->aligned) {
					$nal1++;
				}
				if($al_wat2->aligned || $al_cri2->aligned) {
					$nal2++;
				}
			}
			$al_paired = ($nal == 0 || scalar(@{$als[0]}) == 2);
			if($nal > 0) {
				# Select alignment with greatest MAPQ
				my ($seli, $nmapq1, $nmapq2) = select_al(\@als, \@is_watson, $strand4);
				$al = $als[$seli];
				$nmapq1 <= $al->[0]{_mapq} || die;
				if($al_paired) { $nmapq2 <= $al->[1]{_mapq} || die; }
				$al->[0]{_mapq} = $nmapq1;
				$al->[1]{_mapq} = $nmapq2 if $al_paired;
				if($bsc) {
					bsc_analyze_alignment($al->[0], $is_watson[$seli], \%ref, $ev_sink, undef, undef);
					bsc_analyze_alignment($al->[1], $is_watson[$seli], \%ref, $ev_sink, undef, undef) if $al_paired;
				} else {
					bscpg_analyze_alignment($al->[0], $is_watson[$seli], \%ref, $ev_sink, undef, undef);
					bscpg_analyze_alignment($al->[1], $is_watson[$seli], \%ref, $ev_sink, undef, undef) if $al_paired;
				}
				my $rec = prepare_sam($lines[$seli]->[0], 1, $is_watson[$seli]);
				my $rec2 = undef;
				$rec2 = prepare_sam($lines[$seli]->[1], 1, $is_watson[$seli]) if $al_paired;
				if($is_watson[$seli]) {
					print {$sam_fh_wat} $rec if defined($sam_fh_wat);
					print {$bam_fh_wat} $rec if defined($bam_fh_wat);
					if($al_paired) {
						print {$sam_fh_wat} $rec2 if defined($sam_fh_wat);
						print {$bam_fh_wat} $rec2 if defined($bam_fh_wat);
					}
				} else {
					print {$sam_fh_cri} $rec if defined($sam_fh_cri);
					print {$bam_fh_cri} $rec if defined($bam_fh_cri);
					if($al_paired) {
						print {$sam_fh_cri} $rec2 if defined($sam_fh_cri);
						print {$bam_fh_cri} $rec2 if defined($bam_fh_cri);
					}
				}
			}
			if($nal1 == 0) {
				my $unline = $line_wat;
				my $rec = prepare_sam($unline, 0, undef);
				print {$sam_fh_un} $rec if defined($sam_fh_un);
				print {$bam_fh_un} $rec if defined($bam_fh_un);
			}
			if($nal2 == 0) {
				my $unline = $line_wat2;
				my $rec = prepare_sam($unline, 0, undef);
				print {$sam_fh_un} $rec if defined($sam_fh_un);
				print {$bam_fh_un} $rec if defined($bam_fh_un);
			}
		}
		else {
			$paired || die;
			$strand4 || die;
			my $al_wat_fw1 = parse_sam_record($line_wat);
			my $al_cri_fw1 = parse_sam_record($line_cri);
			my ($line_wat_fw2, $line_cri_fw2) = next_lines($al_wat_fh, $al_cri_fh);
			my $al_wat_fw2 = parse_sam_record($line_wat_fw2);
			my $al_cri_fw2 = parse_sam_record($line_cri_fw2);
			my ($line_wat_rc1, $line_cri_rc1) = next_lines($al_wat_fh, $al_cri_fh);
			my $al_wat_rc1 = parse_sam_record($line_wat_rc1);
			my $al_cri_rc1 = parse_sam_record($line_cri_rc1);
			my ($line_wat_rc2, $line_cri_rc2) = next_lines($al_wat_fh, $al_cri_fh);
			my $al_wat_rc2 = parse_sam_record($line_wat_rc2);
			my $al_cri_rc2 = parse_sam_record($line_cri_rc2);
			$nnt = (length($al_wat_fw1->{_rdseq}) + length($al_wat_fw2->{_rdseq}));
			($al_wat_fw1, $al_wat_fw2) = order_mates($al_wat_fw1, $al_wat_fw2, $m1fw, $m2fw, 0);
			($al_wat_rc1, $al_wat_rc2) = order_mates($al_wat_rc1, $al_wat_rc2, $m1fw, $m2fw, 1);
			($al_cri_fw1, $al_cri_fw2) = order_mates($al_cri_fw1, $al_cri_fw2, $m1fw, $m2fw, 0);
			($al_cri_rc1, $al_cri_rc2) = order_mates($al_cri_rc1, $al_cri_rc2, $m1fw, $m2fw, 1);
			# Look for paired-end evidence
			my ($nal1, $nal2) = (0, 0);
			if($al_wat_fw1->aligned && $al_wat_fw2->aligned) {
				push @als, [ $al_wat_fw1, $al_wat_fw2 ];
				push @lines, [ $line_wat, $line_wat_fw2 ];
				push @is_watson, 1;
				$nal++; $nal1++; $nal2++;
			}
			if($al_cri_fw1->aligned && $al_cri_fw2->aligned) {
				push @als, [ $al_cri_fw1, $al_cri_fw2 ];
				push @lines, [ $line_cri, $line_cri_fw2 ];
				push @is_watson, 0;
				$nal++; $nal1++; $nal2++;
			}
			if($al_wat_rc1->aligned && $al_wat_rc2->aligned) {
				push @als, [ $al_wat_rc1, $al_wat_rc2 ];
				push @lines, [ $line_wat_rc1, $line_wat_rc2 ];
				push @is_watson, 1;
				$nal++; $nal1++; $nal2++;
			}
			if($al_cri_rc1->aligned && $al_cri_rc2->aligned) {
				push @als, [ $al_cri_rc1, $al_cri_rc2 ];
				push @lines, [ $line_cri_rc1, $line_cri_rc2 ];
				push @is_watson, 0;
				$nal++; $nal1++; $nal2++;
			}
			if($nal == 0 && !$no_unpaired) {
				if($al_wat_fw1->aligned || $al_wat_fw2->aligned) {
					push @als,   [ $al_wat_fw1->aligned ? $al_wat_fw1 : $al_wat_fw2 ];
					push @lines, [ $al_wat_fw1->aligned ? $line_wat : $line_wat_fw2 ];
					push @is_watson, 1;
					$nal++;
				}
				if($al_cri_fw1->aligned || $al_cri_fw2->aligned) {
					push @als,   [ $al_cri_fw1->aligned ? $al_cri_fw1 : $al_cri_fw2 ];
					push @lines, [ $al_cri_fw1->aligned ? $line_cri : $line_cri_fw2 ];
					push @is_watson, 0;
					$nal++;
				}
				if($al_wat_rc1->aligned || $al_wat_rc2->aligned) {
					push @als,   [ $al_wat_rc1->aligned ? $al_wat_rc1 : $al_wat_rc2 ];
					push @lines, [ $al_wat_rc1->aligned ? $line_wat_rc1 : $line_wat_rc2 ];
					push @is_watson, 1;
					$nal++;
				}
				if($al_cri_rc1->aligned || $al_cri_rc2->aligned) {
					push @als,   [ $al_cri_rc1->aligned ? $al_cri_rc1 : $al_cri_rc2 ];
					push @lines, [ $al_cri_rc1->aligned ? $line_cri_rc1 : $line_cri_rc2 ];
					push @is_watson, 0;
					$nal++;
				}
				# See if only one or the other aligned
				if($al_wat_fw1->aligned || $al_cri_fw1->aligned &&
				   $al_wat_rc1->aligned || $al_cri_rc1->aligned)
				{
					$nal1++;
				}
				if($al_wat_fw2->aligned || $al_cri_fw2->aligned &&
				   $al_wat_rc2->aligned || $al_cri_rc2->aligned)
				{
					$nal2++;
				}
			}
			$al_paired = ($nal == 0 || scalar(@{$als[0]}) == 2);
			if($nal > 0) {
				# Select alignment with greatest MAPQ
				my ($seli, $nmapq1, $nmapq2) = select_al(\@als, \@is_watson, $strand4);
				$al = $als[$seli];
				$nmapq1 <= $al->[0]{_mapq} || die;
				if($al_paired) { $nmapq2 <= $al->[1]{_mapq} || die; }
				$al->[0]{_mapq} = $nmapq1;
				$al->[1]{_mapq} = $nmapq2 if $al_paired;
				# TODO: Be careful here: when the read aligns both to W and CR,
				# or both to C and WR, and the MAPQs are high for both
				# alignments, we may be looking at a highly methylated read
				# that doesn't deserve a low MAPQ.
				if($bsc) {
					bsc_analyze_alignment($al->[0], $is_watson[$seli], \%ref, $ev_sink, undef, undef);
					bsc_analyze_alignment($al->[1], $is_watson[$seli], \%ref, $ev_sink, undef, undef) if $al_paired;
				} else {
					bscpg_analyze_alignment($al->[0], $is_watson[$seli], \%ref, $ev_sink, undef, undef);
					bscpg_analyze_alignment($al->[1], $is_watson[$seli], \%ref, $ev_sink, undef, undef) if $al_paired;
				}
				my $rec = prepare_sam($lines[$seli]->[0], 1, $is_watson[$seli]);
				my $rec2 = undef;
				$rec2 = prepare_sam($lines[$seli]->[1], 1, $is_watson[$seli]) if $al_paired;
				if($is_watson[$seli]) {
					print {$sam_fh_wat} $rec if defined($sam_fh_wat);
					print {$bam_fh_wat} $rec if defined($bam_fh_wat);
					if($al_paired) {
						print {$sam_fh_wat} $rec2 if defined($sam_fh_wat);
						print {$bam_fh_wat} $rec2 if defined($bam_fh_wat);
					}
				} else {
					print {$sam_fh_cri} $rec if defined($sam_fh_cri);
					print {$bam_fh_cri} $rec if defined($bam_fh_cri);
					if($al_paired) {
						print {$sam_fh_cri} $rec2 if defined($sam_fh_cri);
						print {$bam_fh_cri} $rec2 if defined($bam_fh_cri);
					}
				}
			}
			if($nal1 == 0) {
				my $unline = $line_wat;
				my $rec = prepare_sam($unline, 0, undef);
				print {$sam_fh_un} $rec if defined($sam_fh_un);
				print {$bam_fh_un} $rec if defined($bam_fh_un);
			}
			if($nal2 == 0) {
				my $unline = $line_wat_fw2;
				my $rec = prepare_sam($unline, 0, undef);
				print {$sam_fh_un} $rec if defined($sam_fh_un);
				print {$bam_fh_un} $rec if defined($bam_fh_un);
			}
		}
		if($nal == 0) {
			# None of the options aligned
			$nskip_unal += ($paired ? 2 : 1);
			$nunaltot += ($paired ? 2 : 1);
		} elsif($nal == 1) {
			defined($al) || die;
			# Got exactly one alignment from all the strand combos; could have
			# been paired or unpaired
			$naltot += ($al_paired ? 2 : 1);
			if($paired && !$al_paired) {
				# One mate aligned but the other didn't
				$nskip_unal++;
				$nunaltot++;
			}
			$mapq_hist{$al->[0]{_mapq}}++;
			$mapq_hist{$al->[1]{_mapq}}++ if $al_paired;
		} else {
			defined($al) || die;
			# More than one alignment from the strand combos
			$nskip_mult += ($al_paired ? 2 : 1);
			$nmaxtot += ($al_paired ? 2 : 1);
			if($paired && !$al_paired) {
				# One mate aligned but the other didn't
				$nskip_unal++;
				$nunaltot++;
			}
			$mapq_hist{$al->[0]{_mapq}}++;
			$mapq_hist{$al->[1]{_mapq}}++ if $al_paired;
		}
		$nrecs_p += ($paired ? 2 : 1);
		$nttot += $nnt;
		if((++$nrecs % $nrec_ival) == 0) {
			my $pct_al   = (100.00 * $naltot   / (1.0 * $nrecs_p || 1));
			my $pct_max  = (100.00 * $nmaxtot  / (1.0 * $nrecs_p || 1));
			my $pct_unal = (100.00 * $nunaltot / (1.0 * $nrecs_p || 1));
			my $t1 = time();
			my $ntps = (1.0 * $nttot) / (1000.0 * ($t1 - $t0));
			printf STDERR "  Processed $nrecs reads (%.3f knt/s, %.3f%% align in 1 orientation, %.3f%% in >1, %.3f%% in 0)...\n", $ntps, $pct_al, $pct_max, $pct_unal;
		}
	}
}
print STDERR "WARNING: Saw no SAM records\n" if $nrecs == 0;

print STDERR "Waiting for children if necessary\n";
waitpid($pid_feed, 0) if $pid_feed > 0;

print STDERR "Closing .sam/.bam filehandles\n";
for($sam_fh_wat, $sam_fh_cri, $sam_fh_un, $bam_fh_wat, $bam_fh_cri, $bam_fh_un) {
	close($_) if defined($_);
}

for my $ofh (%ev_out_fhs) { close($ofh); }
rmtree([$temp_dir]) unless defined($keep_dir);

# TODO: Print to a file as well as to STDERR

print STDERR "Summary:\n";
print STDERR "  SAM header lines parsed: $nheads\n";
print STDERR "  SAM record sets parsed: $nrecs\n";
print STDERR "  Loop iterations: $iters\n";
print STDERR "  Read-level measurements extracted: $nev\n" if !$stop_after_aln;
my $pct_al   = (100.00 * $naltot   / (1.0 * $nrecs_p || 1));
my $pct_max  = (100.00 * $nmaxtot  / (1.0 * $nrecs_p || 1));
my $pct_unal = (100.00 * $nunaltot / (1.0 * $nrecs_p || 1));
printf STDERR "  Of $nrecs input reads:\n";
printf STDERR "    $nunaltot (%.3f%%) reads/ends failed to align\n", $pct_unal;
printf STDERR "    $nmaxtot (%.3f%%) reads/ends aligned with bisulfite-strand ambiguity\n", $pct_max;
printf STDERR "    $naltot (%.3f%%) reads/ends aligned without bisulfite-strand ambiguity\n", $pct_al;
for my $k (sort {$a <=> $b} keys %mapq_hist) {
	print STDERR "      $mapq_hist{$k} at MAPQ=$k\n";
}
