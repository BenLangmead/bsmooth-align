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
# bs_merman_align.pl
#
# Given a set of reads, align the reads to the reference using merman in a way
# that does not penalize evidence for or against methylation.
#
# User input:
#
# There are two things that we would like the user to be able to input rather
# flexibly to this tool: a list of arguments for the aligner and a list of read
# files.  The aligner argument list is optional, but there must be at least one
# read file.  If possible, we would like to avoid the user having to specify
# options or read files one-by-one, each with its own argument.  One option is
# to use ARGV for both: the user specifies options for this script, then a
# separator, --, then a list of read files, then another separator, --, then a
# list of aligner arguments.
#
# Intermediate reads files versus named pipes:
#
# We'd like to avoid making intermediate files for the in-silico bisulfite
# converted reads.  The same set of converted reads need to be provided to both
# bowtie processes, regardless of whether the bowtie processes run in series or
# in parallel.  One option is to convert all the reads once and store them in
# intermediate files.  But another option is to convert them twice and pipe the
# output from the two conversion processes directly to the bowtie processes.
#
# Other issues:
#
# - Information about the original read is "threaded through" the alignment
#   step by appending the original read sequence to the read name.  This can
#   make the read name quite long.  In fact, it can make the read name longer
#   than the maximum permitted by the SAM spec (255 characters).  For this
#   reason, we have to be careful to ask the aligner to ignore the QNAME length
#   requirement when printing alignments.  Otherwise it could truncate the
#   saved sequnce.  Bowtie 1 and merman ignore it by default, but Bowtie 2
#   requires a special option (--sam-no-qname-trunc) to ignore it.
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
use BtlBio::Read::Qseq;
use BtlBio::Read::Fastq;
use BtlBio::Read::Fasta;
use BtlBio::Util::File;
use BtlBio::Util::Tool;
use BtlBio::Util::Temp;
use BtlBio::Format::SamBam::SamBam;
use POSIX;
use Carp;

my $ALIGNER = "merman";         # aligner binary name
my $aligner_exe = $ALIGNER;     # path to aligner bin
my $samtools_exe = "samtools";  # path to samtools bin; only for --keep-bam
my @aligner_args = ();          # arguments to pass to aligner
my $bscpg = 1;                  # 1 -> do CpG->YpG conversion
my $bsc = 0;                    # 1 -> do C->Y conversion
my $temp_dir = $ENV{TMPDIR};    # temp dir for intermediates; reads, alns
$temp_dir = $temp_dir || "/tmp";# temp dir for intermediates; reads, alns
my $sam_fn = undef;             # put SAM results here
my $bam_fn = undef;             # put BAM results here
my $format = "fastq";           # format of input files
my $use_named_pipe = 1;         # use named pipes, not intermediate aln files
my $echo_sam = 0;               # print out SAM as it comes in
my $ev_out_dir = ".";           # directory to write evidence
my $strand4 = 0;                # using 4-strand protocol?
my $keep_dir = undef;           # where to keep intermediate reads/alns
my $stop_after_aln = 0;         # stop after alignment step, leaving .sam/.bam
my $stop_after = 0;             # stop after this many reads
my @name_map_fn = ();           # list of files with long/short name maps
my %name_map = ();              # map from long to short names
my $compress = "";              # manner of compression ev output
my $skip_reads = 0;             # skip this many reads off the bat
my @mate1_fns = ();             # read files with unpaired reads/mate #1s
my @ref_fns = ();               # fasta files with references

my $usage = qq!
Usage:
  bs_merman_align.pl [options*] -- <refs> -- <merman-args> -- <reads>

Align reads to reference using bisulfite-aware aligner Merman.  CpG Cs and/or
non-CpG Cs in the genome are treated as Ys, matching either Cs or Ts with no
penalty.  Output is a directory of read-level measurement files, which can be
sorted with with the bsev_sort.pl script; the sorted directory can then be
tabulated with the bsev_tabulate.pl script.  Can also output alignments.

Required arguments:
  
  Required arguments are separated from options list and from each other by
  double dashes (--).  This allows more flexible use of wildcards (e.g. *.fa).
  
  <refs>
    Fasta files containing all original reference sequences
  <merman-args>
    Arguments to pass to Merman
  <reads>
    Files with unpaired reads to align.  Merman does not support paired reads.

Options (defaults in parentheses):
 
 Alignment:
  --merman=<path>          Path to merman executable (from PATH)
  --four-strand            Align as though read/fragment could have originated
                           from BSW, BSWR, BSC, BSCR (def: just BSW, BSC)
  
  --bscpg                  Extract read-level measurements from CpG Cs (on)
   OR:
  --bsc                    Extract read-level measurements from all Cs (off)
 
 Input:
  --fastq                  Reads files are FASTQ format (default)
  --qseq                   Reads files are in Illumina _qseq.txt format
  --fasta                  Reads files are in FASTA format
  --csfasta                Reads files are in CSFASTA format (if specified file
                           is XYZ, XYZ.csfasta/XYZ_qual.QV contain reads/quals)
  --csfastq                Reads files are in CSFASTQ format
  -s/--skip-reads <int>    Skip first <int> input reads
  -u/--stop-after <int>    Only align first <int> input reads
 
 Output:
  --name-map=<path>        Many-to-one map from raw reference names to bin
                           names.  Output will be binned by bin name.
  --output-dir=<path>      Write read-level measurements to <path> (.)
  --keep-sam=<prefix>      Store Merman output .sam (BS-SAM) in <prefix>.sam
  --keep-bam=<prefix>      Like --keep-sam but stores as .bam instead
  --samtools=<path>        samtools binary for converting SAM to BAM when
                           --keep-bam or --keep-merged-bam are used
  --output=<path>          Directory to store read-level measurements output
                           or --keep-fasta are specified (.)
  --gzip                   gzip compress read-level measurement output
  --bzip2                  bzip2 compress read-level measurement output
 
 Misc:
  --stop-after-alignment   Don't output read-level measurements
  --temp=<path>            Set path used for temporary FASTA files (TMPDIR)
  --keep=<path>            Keep intermediate/temporary files at <path>
  --help/--usage           Print this message and quit

!;

##
# Print and run given command.  Die if it returns non-zero.
#
sub run($) {
	my $cmd = shift;
	print STDERR "Running '$cmd'...\n";
	system($cmd);
	$? == 0 || confess("Command '$cmd' failed with exitlevel $?\n");
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

GetOptions (
	"merman=s"        => \$aligner_exe,
	"samtools=s"      => \$samtools_exe,
	"echo-sam"        => \$echo_sam,
	"merman-args=s"   => \@aligner_args,
	"aligner-args=s"  => \@aligner_args,
	"no-named-pipes"  => sub { $use_named_pipe = 0; },
	"keep-sam=s"      => \$sam_fn,
	"keep-bam=s"      => \$bam_fn,
	"bscpg"           => sub { $bscpg = 1; $bsc = 0; },
	"bsc"             => sub { $bscpg = 0; $bsc = 1; },
	"read=s"          => \@mate1_fns,
	"four-strand"     => \$strand4,
	"fastq"           => sub { $format = "fastq" },
	"qseq"            => sub { $format = "qseq"  },
	"fasta"           => sub { $format = "fasta" },
	"csfasta"         => sub { $format = "csfasta" },
	"csfastq"         => sub { $format = "csfastq" },
	"temp-directory|temporary-directory=s" => \$temp_dir,
	"output-directory|evidence-directory=s" => \$ev_out_dir,
	"keep=s"          => \$keep_dir,
	"stop-after-alignment" => \$stop_after_aln,
	"s|skip-reads=i"  => \$skip_reads,
	"u|stop-after=i"  => \$stop_after,
	"name-map=s"      => \@name_map_fn,
	"compress=s"      => \$compress,
	"gzip"            => sub { $compress = "gzip";  },
	"bzip2"           => sub { $compress = "bzip2"; },
	"help|usage"      => sub { print $usage; exit 0; }
) || dieusage("Bad option", 0);

mkpath($ev_out_dir);

# Parse out arguments from among the double-dashes
my $dds = 0;
for my $a (@ARGV) {
	if($a eq "--") {
		$dds++; next;
	}
	if($dds == 0) {
		push @ref_fns, $a;
	} elsif($dds == 1) {
		push @aligner_args, $a;
	} elsif($dds == 2) {
		push @mate1_fns, $a;
	}
}

# Inject default Bowtie arguments if none were specified
if(scalar(@aligner_args) == 0) {
	push @aligner_args, "-M 1";
	push @aligner_args, "--best";
}

mkpath($keep_dir) if defined($keep_dir);
my $inter1 = temp_name(".reads", $keep_dir);
my $inter2 = temp_name(".reads", $keep_dir);

print STDERR "Reference files:\n";
for(0..$#ref_fns) {print STDERR "  $ref_fns[$_]\n"}
print STDERR "$ALIGNER executable: $aligner_exe\n";
print STDERR "$ALIGNER arguments: ".join(" ", @aligner_args)."\n";
print STDERR "Methylation output directory: $ev_out_dir\n";
print STDERR "Keep .sam: ".(defined($sam_fn) ? $sam_fn : "no")."\n";
print STDERR "Keep .bam: ".(defined($bam_fn) ? $bam_fn : "no")."\n";
print STDERR "samtools executable override: $samtools_exe\n" if $bam_fn;
print STDERR "Read-level measurements: ".($bscpg ? "just CpGs" : "all Cs")."\n";
print STDERR "Four-strand?: ".($strand4 ? "yes":"no")."\n";
print STDERR "Unpaired inputs:\n";
for(0..$#mate1_fns) {print STDERR "  $mate1_fns[$_]\n"}
print STDERR "Input format: $format\n";
print STDERR "Skip first # reads: ".($skip_reads > 0 ? $skip_reads : "(no skipping)")."\n";
print STDERR "Stop after # reads: ".($stop_after > 0 ? $stop_after : "(no limit)")."\n";
print STDERR "Temporary directory: $temp_dir\n";
print STDERR "Intermediate files:\n";
print STDERR "  $inter1\n";
print STDERR "  $inter2\n";

print STDERR "Locating tools...\n";
$aligner_exe = find_tool($ALIGNER, $aligner_exe);
$samtools_exe = find_tool("samtools", $samtools_exe) if $bam_fn;
print STDERR "  Found $ALIGNER: $aligner_exe\n";
print STDERR "  Found samtools: $samtools_exe\n" if $bam_fn;

print STDERR "Checking reads...\n";
# Do all the files exist?
for(@mate1_fns) { $_ eq "-" || -f $_ || die "Could not find reads file '$_'"; }

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
# Given a filename and a format, open the file, process each read through
# in-silico bisulfite conversion, and write the converted read (in FASTQ
# format) to the given filehandle.
#
sub cat_file($$) {
	my ($ifn, $ofh) = @_;
	return if $stop_after > 0 && $nreads >= $stop_after;
	print STDERR "  Writing reads from file '$ifn' ...\n";
	my $ifh;
	if($ifn eq "-") {
		$ifh = *STDIN;
	} else {
		$ifh = openex($ifn);
	}
	if($format =~ /fastq/i) {
		while(my $read = parse_fastq_read($ifh)) {
			return if $stop_after > 0 && $nreads >= $stop_after;
			print {$ofh} $read->to_fastq;
			$nreads++;
		}
	} elsif($format =~ /qseq/i) {
		while(my $read = parse_qseq_read($ifh)) {
			return if $stop_after > 0 && $nreads >= $stop_after;
			print {$ofh} $read->to_fastq;
			$nreads++;
		}
	} else {
		$format =~ /fasta/i || die "Bad format: '$format'";
		while(my $read = parse_fasta_read($ifh)) {
			return if $stop_after > 0 && $nreads >= $stop_after;
			print {$ofh} $read->to_fastq;
			$nreads++;
		}
	}
}

# We can run two instances of Bowtie simultaneously, or we can run them one
# after the other.
if($use_named_pipe) {
	print STDERR "Creating named pipes...\n";
	mkfifo($inter1, 0700);
}

# fork off a process whose job it is to feed Bowtie
my ($pid1, $pid2) = (0, 0);
scalar(@mate1_fns) > 0 || die;
if($use_named_pipe) {
	$pid1 = fork();
	if($pid1 == 0) {
		# Open named pipe 1 for writing
		open(my $ofh, ">$inter1") ||
			die "Could not open '$inter1' for writing";
		for my $ifn (@mate1_fns) { cat_file($ifn, $ofh); }
		close($ofh);
		exit 0;
	}
	print STDERR "Forked off feeder process 1 (pid=$pid1)...\n";
} else {
	# Instead of using named pipes and forking off processes to feed them, we
	# simply dump all the reads to intermediate files within this process.  We
	# need only process the reads once - the intermediate files will then be
	# read twice.
	print STDERR "Processing reads for mate 1 into file '$inter1'...\n";
	open(my $ofh1, ">$inter1") ||
		die "Could not open '$inter1' for writing";
	for my $ifn (@mate1_fns) { cat_file($ifn, $ofh1); }
	close($ofh1);
}

my ($al_fh) = (undef, undef);

# Open a bowtie process that reads from the intermediate files/pipes and
# outputs alignments in .sam format.

my $aligner_req_args = "--sam --quiet";
if($bscpg) {
	if($strand4) {
		$aligner_req_args .= " --bisCpG-4";
	} else {
		$aligner_req_args .= " --bisCpG-2";
	}
} else {
	if($strand4) {
		$aligner_req_args .= " --bisC-4";
	} else {
		$aligner_req_args .= " --bisC-2";
	}
}

my $cmd = "$aligner_exe ".join(" ", @aligner_args)." $aligner_req_args";
# Add the references
$cmd .= " ";
$cmd .= join(",", @ref_fns);
$cmd .= " $inter1";

# We're going to run two aligners simultaneously
$cmd .= " |"; # Will read your output directly
print STDERR "Opening merman pipe:\n";
print STDERR "  $cmd\n";
open($al_fh, $cmd) || die "Could not open pipe '$cmd' for reading";

my ($sam_fh, $bam_fh) = (undef, undef);
if(defined($sam_fn)) {
	$sam_fn .= ".sam";
	open($sam_fh, ">$sam_fn") || die "Could not open '$sam_fn' for writing";
}
if(defined($bam_fn)) {
	$bam_fn .= ".bam";
	open($bam_fh, "| $samtools_exe view -Sb - > $bam_fn")
		|| die "Could not open '$bam_fh' for writing";
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
# Mark the evidence directory with a file indicating what type of evidence is
# there (either CpG-level or C-level)
open(MARK, ">$ev_out_dir/.".($bscpg ? "cpg" : "c").".ev") || die;
print MARK "\n"; close(MARK);
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
print STDERR "Processing SAM output from aligners...\n";
my $iters = 0;
my ($naltot, $nunaltot) = (0, 0, 0);
my %mapq_hist = ();
while(1) {
	$iters++;
	# Output is in SAM format
	my $line = readline $al_fh;
	last unless defined($line);
	# Pass SAM line through to .sam/.bam files if --keep-sam/--keep-bam are
	# specified
	print {$sam_fh} $line if defined($sam_fh);
	print {$bam_fh} $line if defined($bam_fh);
	print $line if $echo_sam;
	my $firstc = substr($line, 0, 1);
	if($firstc eq "\@") {
		# Header line
		$nheads++;
		next; # skip
	}
	if(!$stop_after_aln) {
		# Create Alignment objects for both
		my $al = parse_sam_record($line);
		if($al->aligned) { $naltot++;   }
		else             { $nunaltot++; }
		if($al->aligned) {
			bs_analyze_alignment($al, \&evidenceSink) ;
			defined($al->{_mapq}) || die "Couldn't parse mapping quality from line:\n$line";
			$al->{_mapq} == int($al->{_mapq}) || die "Couldn't parse mapping quality from line:\n$line";
			$mapq_hist{$al->{_mapq}}++;
		}
	}
	if((++$nrecs % $nrec_ival) == 0) {
		print STDERR "  Processed $nrecs SAM records...\n";
	}
}
close($al_fh);
$? == 0 || die "Non-zero exitlevel from aligner: $?";

print STDERR "Waiting for children if necessary\n";
waitpid($pid1, 0) if $pid1 > 0;
waitpid($pid2, 0) if $pid2 > 0;

if(defined($sam_fh)) {
	print STDERR "Closing .sam output filehandles\n";
	close($sam_fh);
	-f $sam_fn || die "Failed to create .sam file '$sam_fn'";
}
if(defined($bam_fh)) {
	print STDERR "Closing .bam output filehandles\n";
	close($bam_fh);
	-f $bam_fn || die "Failed to create .bam file '$bam_fn'";
}

##
# Remove a file and print a message
#
sub rm($) {
	defined($_[0]) || croak();
	print STDERR "  Removing $_[0]...\n";
	unlink($_[0]);
}

if(!defined($keep_dir)) {
	print STDERR "Removing named pipes/in-silico converted reads...\n";
	rm($inter1);
}

print STDERR "Summary:\n";
print STDERR "  SAM header lines parsed: $nheads\n";
print STDERR "  SAM record sets parsed: $nrecs\n";
print STDERR "  Loop iterations: $iters\n";
print STDERR "  Read-level measurements extracted: $nev\n" if !$stop_after_aln;

my $pct_al   = (100.00 * $naltot   / (1.0 * $nrecs || 1));
my $pct_unal = (100.00 * $nunaltot / (1.0 * $nrecs || 1));
if(!$stop_after_aln) {
	printf STDERR "  Of $nrecs input reads:\n";
	printf STDERR "    $naltot (%.3f%%) reads aligned without bisulfite-strand ambiguity\n", $pct_al;
	printf STDERR "    $nunaltot (%.3f%%) reads failed to align\n", $pct_unal;
}
for my $k (sort {$a <=> $b} keys %mapq_hist) {
	print STDERR "      $mapq_hist{$k} at MAPQ=$k\n";
}
for my $ofh (%ev_out_fhs) { close($ofh); }
