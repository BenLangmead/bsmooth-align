#!/usr/bin/env perl

#
# Copyright (C) 2012, Ben Langmead <blangmea@jhsph.edu>
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

use strict;
use warnings;
use Carp;
use FindBin qw($Bin);
use lib "$Bin/../../lib";
use BtlBio::Alphabet::DNA;
use BtlBio::Util::Temp;
use BtlBio::Align::Bisulfite::BsEvidence;
use BtlBio::Align::Bisulfite::BsSummary;
use File::Path;
use Data::Dumper;
use File::Which;

which("samtools") || die "Could not execute 'samtools'";
which("bowtie2") || die "Could not execute 'bowtie2'";
which("bowtie2-build") || die "Could not execute 'bowtie2-build'";

my $build_cmd_base    = "perl $Bin/../../bin/bswc_bowtie2_index.pl";
my $align_cmd_base    = "perl $Bin/../../bin/bswc_bowtie2_align.pl";
my $extract_cmd_base  = "perl $Bin/../../bin/bswc_sam_extract.pl";
my $sort_cmd_base     = "perl $Bin/../../bin/bsev_sort.pl";
my $tabulate_cmd_base = "perl $Bin/../../bin/bsev_tabulate.pl";

# Ref:      CATGCGTATTTTTTGCATCATCACGGTAGCGGGGCGCCGCGTAGTCATTCGACTAATTACGGTGCTGTGACTCTAACTG
#           0123456789012345678901234567890123456789012345678901234567890123456789012345678
#           0         1         2         3         4         5         6         7        
# BSW:      TATGTGTATTTTTTGTATTATTATGGTAGTGGGGTGTTGTGTAGTTATTTGATTAATTATGGTGTTGTGATTTTAATTG
# Read FW:     GCGTATTTTTTGTATCATTATGGT   GGGGTGCCGTGTAGTTATTTGATTAATTACGGTGCTGTG
# Read RC:     ACCATAATGATACAAAAAATACGC   CACAGCACCGTAATTAATCAAATAACTACACGGCACCCC

# Ref RC:   CAGTTAGAGTCACAGCACCGTAATTAGTCGAATGACTACGCGGCGCCCCGCTACCGTGATGATGCAAAAAATACGCATG
#           0123456789012345678901234567890123456789012345678901234567890123456789012345678
#           0         1         2         3         4         5         6         7        
# BSC:      TAGTTAGAGTTATAGTATTGTAATTAGTTGAATGATTATGTGGTGTTTTGTTATTGTGATGATGTAAAAAATATGTATG
# Read FW:      TAGAGTCACAGCATTGTAATTAGTTGAAT      GTGGTGCCCCGTTATTGTGATGATGCAAAAAATACGTATG
# Read RC:      ATTCAACTAATTACAATGCTGTGACTCTA      CATACGTATTTTTTGCATCATCACAATAACGGGGCACCAC

my $ref = "CATGCGTATTTTTTGCATCATCACGGTAGCGGGGCGCCGCGTAGTCATTCGACTAATTACGGTGCTGTGACTCTAACTG";

my $read1 = "ACCATAATGATACAAAAAATACGC";                 # from BSWR
my $qual1 = "AABBCCDDEEFFGGHHIIJJKKLL";
my $read2 = "GGGGTGCCGTGTAGTTATTTGATTAATTACGGTGCTGTG";  # from BSW
my $qual2 = "IJIJKMLASDHJAUYSGIAUDBKJSHFBAUYFGAUYNER";
my $read3 = "TAGAGTCACAGCATTGTAATTAGTTGAAT";            # from BSC
my $qual3 = "AERFUVGIBHNJMNBHIGVYUFCTRYGBH";
my $read4 = "CATACGTATTTTTTGCATCATCACAATAACGGGGCACCAC"; # from BSCR
my $qual4 = "PISDUAOSWBISFOUASYFGASFEHBoiufhgwoeifsad";
#             ^  ^                              ^  ^

# Read 1 fwbs: ATTATAATGATATAAAAAATATGT
#              ATTGTGATGATGTAAAAAATATGT
# Read 1 rcbs: GTGTATTTTTTGTATTATTATGGT
#              GTGTATTTTTTGTATTATTATGGT

# Read 3 fwbs: TAGAGTTATAGTATTGTAATTAGTTGAAT
#              TAGAGTTATAGTATTGTAATTAGTTGAAT
# Read 3 rcbs: ATTTGATTAATTATGGTGTTGTGATTTTA
#              ATTTGATTAATTACGGTGTTGTGATTTTA

# Read 1 (BSWR):
#                1         2
# Refoff: 345678901234567890123456
# Ref:    GCGTATTTTTTGCATCATCACGGT
# BsC:    |Y||||||||||Y||Y||Y|Y|||
# BsCpG:  |Y||||||||||||||||||Y|||
# Read:   GCGTATTTTTTGTATCATTATGGT
# Cyc:    321098765432109876543210
# Qual:   LLKKJJIIHHGGFFEEDDCCBBAA
#
# C Ev:   ref1,  4, "C", 1, 0, "L", -1, 22, 0, 255
# C Ev:   ref1, 15, "T", 1, 0, "F", -1, 11, 0, 255
# C Ev:   ref1, 18, "C", 1, 0, "E", -1,  8, 0, 255
# C Ev:   ref1, 21, "T", 1, 0, "C", -1,  5, 0, 255
# C Ev:   ref1, 23, "T", 1, 0, "B", -1,  3, 0, 255
#
# CpG Ev: ref1,  4, "C", 1, 0, "L", -1, 22, 0, 255
# CpG Ev: ref1, 23, "T", 1, 0, "B", -1,  3, 0, 255
#
# Read 2 (BSW):
#         3         4         5         6        
# Refoff: 012345678901234567890123456789012345678
# Ref:    GGGGCGCCGCGTAGTCATTCGACTAATTACGGTGCTGTG
# BsC:    ||||Y|YY|Y|||||Y|||Y||Y||||||Y||||Y||||
# BsCpG:  ||||Y||Y|Y|||||||||Y|||||||||Y|||||||||
# Read:   GGGGTGCCGTGTAGTTATTTGATTAATTACGGTGCTGTG
# Cyc:    012345678901234567890123456789012345678
# Qual:   IJIJKMLASDHJAUYSGIAUDBKJSHFBAUYFGAUYNER
#
# C Ev:   ref1, 34, "T", 1, 1, "K", -1,  4, 0, 255
# C Ev:   ref1, 36, "C", 1, 1, "L", -1,  6, 0, 255
# C Ev:   ref1, 37, "C", 1, 1, "A", -1,  7, 0, 255
# C Ev:   ref1, 39, "T", 1, 1, "D", -1,  9, 0, 255
# C Ev:   ref1, 45, "T", 1, 1, "S", -1, 15, 0, 255
# C Ev:   ref1, 49, "T", 1, 1, "U", -1, 19, 0, 255
# C Ev:   ref1, 52, "T", 1, 1, "K", -1, 22, 0, 255
# C Ev:   ref1, 59, "C", 1, 1, "U", -1, 29, 0, 255
# C Ev:   ref1, 64, "C", 1, 1, "U", -1, 34, 0, 255
#
# CpG Ev: ref1, 34, "T", 1, 1, "K", -1,  4, 0, 255
# CpG Ev: ref1, 37, "C", 1, 1, "A", -1,  7, 0, 255
# CpG Ev: ref1, 39, "T", 1, 1, "D", -1,  9, 0, 255
# CpG Ev: ref1, 49, "T", 1, 1, "U", -1, 19, 0, 255
# CpG Ev: ref1, 59, "C", 1, 1, "U", -1, 29, 0, 255
#
# Read 3 (BSC):
#             5         6         7
# Refoff: 67890123456789012345678901234
# Ref:    ATTCGACTAATTACGGTGCTGTGACTCTA
# BsC:    ||||R|||||||||RR|R||R|R||||||
# BsCpG:  ||||R|||||||||R||||||||||||||
# Read:   ATTCAACTAATTACAATGCTGTGACTCTA
# Cyc:    87654321098765432109876543210
# Qual:   HBGYRTCFUYVGIHBNMJNHBIGVUFREA
#
# C Ev:   ref1, 50, "A", 0, 1, "R", -1, 24, 0, 255
# C Ev:   ref1, 60, "A", 0, 1, "B", -1, 14, 0, 255
# C Ev:   ref1, 61, "A", 0, 1, "N", -1, 13, 0, 255
# C Ev:   ref1, 63, "G", 0, 1, "J", -1, 11, 0, 255
# C Ev:   ref1, 66, "G", 0, 1, "B", -1,  8, 0, 255
# C Ev:   ref1, 68, "G", 0, 1, "G", -1,  6, 0, 255
#
# CpG Ev: ref1, 50, "A", 0, 1, "R", -1, 24, 0, 255
# CpG Ev: ref1, 60, "A", 0, 1, "B", -1, 14, 0, 255
#
# Read 4 (BSCR):
#                   1         2         3
# Refoff: 0123456789012345678901234567890123456789
# Ref:    CATGCGTATTTTTTGCATCATCACGGTAGCGGGGCGCCGC
# BsC:    |||R|R||||||||R|||||||||RR||R|RRRR|R||R|
# BsCpG:  |||||R||||||||||||||||||R|||||R||||R||R|
# Read:   CATACGTATTTTTTGCATCATCACAATAACGGGGCACCAC
# Cyc:    0123456789012345678901234567890123456789
# Qual:   PISDUAOSWBISFOUASYFGASFEHBoiufhgwoeifsad
#
# C Ev:   ref1,  3, "A", 0, 0, "D", -1, 24, 0, 255
# C Ev:   ref1,  5, "G", 0, 0, "A", -1, 14, 0, 255
# C Ev:   ref1, 14, "G", 0, 0, "U", -1, 13, 0, 255
# C Ev:   ref1, 24, "A", 0, 0, "H", -1, 13, 0, 255
# C Ev:   ref1, 25, "A", 0, 0, "B", -1, 13, 0, 255
# C Ev:   ref1, 28, "A", 0, 0, "u", -1, 13, 0, 255
# C Ev:   ref1, 30, "G", 0, 0, "h", -1, 13, 0, 255
# C Ev:   ref1, 31, "G", 0, 0, "g", -1, 13, 0, 255
# C Ev:   ref1, 32, "G", 0, 0, "w", -1, 13, 0, 255
# C Ev:   ref1, 33, "G", 0, 0, "o", -1, 13, 0, 255
# C Ev:   ref1, 35, "A", 0, 0, "i", -1, 13, 0, 255
# C Ev:   ref1, 38, "A", 0, 0, "a", -1, 13, 0, 255

# CpG Ev: ref1,  5, "G", 0, 0, "A", -1, 14, 0, 255
# CpG Ev: ref1, 24, "A", 0, 0, "H", -1, 13, 0, 255
# CpG Ev: ref1, 30, "G", 0, 0, "h", -1, 13, 0, 255
# CpG Ev: ref1, 35, "A", 0, 0, "i", -1, 13, 0, 255
# CpG Ev: ref1, 38, "A", 0, 0, "a", -1, 13, 0, 255

my $fa1 = temp_fasta_filename([$ref], ["ref1"]);

my $fq1 = temp_fastq_filename(
	[ $read1,  $read2,  $read3,  $read4  ],
	[ $qual1,  $qual2,  $qual3,  $qual4  ],
	[ "read1", "read2", "read3", "read4" ]
);

my $ev_dir = ".tmp.ev";
rmtree($ev_dir);
my $sort_dir = ".tmpsort.ev";
rmtree($sort_dir);
my $tab_dir = ".tmptab.ev";
rmtree($tab_dir);

sub _test_evidence($$$$$$$$$$$$) {
	my ($ev, $tname, $off, $al, $watson, $fw, $flags, $qua1, $qua2, $cyc,
	    $score, $mapq) = @_;
	$off++; # change 0-based ref offset to 1-based
	$ev->ref_name eq $tname  || confess("Expected refname=$tname, got " .$ev->ref_name     . Dumper($ev));
	$ev->ref_off == $off     || confess("Expected offset=$off, got "    .$ev->ref_off      . Dumper($ev));
	$ev->ev eq $al           || confess("Expected al=$al, got "         .$ev->ev           . Dumper($ev));
	$ev->watson == $watson   || confess("Expected watson=$watson, got " .$ev->watson       . Dumper($ev));
	$ev->fw == $fw           || confess("Expected fw=$fw, got "         .$ev->fw           . Dumper($ev));
	$ev->flags == $flags     || confess("Expected flags=$flags, got "   .$ev->flags        . Dumper($ev));
	$ev->qual1 eq $qua1      || confess("Expected qual1=$qua1, got "    .$ev->qual1        . Dumper($ev));
	$ev->qual2 eq $qua2      || confess("Expected qual2=$qua2, got "    .$ev->qual2        . Dumper($ev));
	$ev->cycle == $cyc       || confess("Expected cycle=$cyc, got "     .$ev->cycle        . Dumper($ev));
	$ev->aln_score == $score || confess("Expected score=$score, got "   .$ev->aln_score    . Dumper($ev));
	if($mapq eq "hi") {
		$ev->aln_mapq >= 20  || confess("Expected mapq >= 20, got "     .$ev->aln_mapq     . Dumper($ev));
	} elsif($mapq eq "lo") {
		$ev->aln_mapq < 20   || confess("Expected mapq < 20, got "      .$ev->aln_mapq     . Dumper($ev));
	} else {
		$ev->aln_mapq  == $mapq  || confess("Expected mapq=$mapq, got "     .$ev->aln_mapq . Dumper($ev));
	}
	return 1;
}

for my $do_extract (undef, "sam", "bam") {
	for my $type ("bscpg", "bsc") {
		my $name = "test1";
		my $args = "--$type";
		my $sam_temp = temp_name();
		
		my $temp = temp_name();
		mkpath($temp);
		my $build_cmd = "$build_cmd_base --name=$name --temp=$temp $fa1";
		print STDERR "$build_cmd\n";
		system($build_cmd);
		$? == 0 || die "Bad exitlevel: $?";
		
		if($do_extract) {
			my $align_cmd =
				"$align_cmd_base ".
				"--echo-sam ".
				"--$do_extract=$sam_temp ".
				"--four-strand ".
				"$args ".
				"-- $name ".  # index
				"-- $fa1 ".   # refs
				"-- --score-min C,-10,0 ".
				"-- $fq1";    # reads
			print STDERR "$align_cmd\n";
			system($align_cmd);
			$? == 0 || die "Bad exitlevel: $?";
			
			rmtree($ev_dir);
			my $extract_cmd =
				"$extract_cmd_base ".
				"--echo-sam ".
				"--output=$ev_dir ".
				"$args ".
				"-- $fa1 ".
				"-- $sam_temp.watson.$do_extract ".
				"-- $sam_temp.crick.$do_extract";
			print STDERR "$extract_cmd\n";
			system($extract_cmd);
			$? == 0 || die "Bad exitlevel: $?";
		} else {
			rmtree($ev_dir);
			my $align_cmd =
				"$align_cmd_base ".
				"--echo-sam ".
				"--output=$ev_dir ".
				"--four-strand ".
				"$args ".
				"-- $name ".  # index
				"-- $fa1 ".   # refs
				"-- --score-min C,-10,0 ".
				"-- $fq1";    # reads
			print STDERR "$align_cmd\n";
			system($align_cmd);
			$? == 0 || die "Bad exitlevel: $?";
		}
		
		for(my $i = 0; $i <= 1; $i++) {
			my @ev = ();
			if($i == 0) {
				parse_meth_ev_records_dirname($ev_dir, \@ev);
			} else {
				# Sort it
				rmtree($sort_dir);
				my $sort_cmd =
					"$sort_cmd_base ".
					"--ev $ev_dir ".
					"--out $sort_dir ";
				print STDERR "$sort_cmd\n";
				system($sort_cmd);
				$? == 0 || die "Bad exitlevel: $?";
				parse_meth_ev_records_dirname($sort_dir, \@ev);
			}
			if($type eq "bsc") {
				scalar(@ev) == 32 || die "Expected 32 pieces of evidence, got ".scalar(@ev).":".Dumper(\@ev);
				@ev = sort {$a->rdid cmp $b->rdid or $a->ref_off <=> $b->ref_off} @ev;
				# BSWR
				_test_evidence($ev[0],  "ref1",  4, "C", 1, 0, 0, "L", -1, 22, 0, "hi");
				_test_evidence($ev[1],  "ref1", 15, "T", 1, 0, 0, "F", -1, 11, 0, "hi");
				_test_evidence($ev[2],  "ref1", 18, "C", 1, 0, 0, "E", -1,  8, 0, "hi");
				_test_evidence($ev[3],  "ref1", 21, "T", 1, 0, 0, "C", -1,  5, 0, "hi");
				_test_evidence($ev[4],  "ref1", 23, "T", 1, 0, 0, "B", -1,  3, 0, "hi");
				# BSW
				_test_evidence($ev[5],  "ref1", 34, "T", 1, 1, 0, "K", -1,  4, 0, "hi");
				_test_evidence($ev[6],  "ref1", 36, "C", 1, 1, 0, "L", -1,  6, 0, "hi");
				_test_evidence($ev[7],  "ref1", 37, "C", 1, 1, 0, "A", -1,  7, 0, "hi");
				_test_evidence($ev[8],  "ref1", 39, "T", 1, 1, 0, "D", -1,  9, 0, "hi");
				_test_evidence($ev[9],  "ref1", 45, "T", 1, 1, 0, "S", -1, 15, 0, "hi");
				_test_evidence($ev[10], "ref1", 49, "T", 1, 1, 0, "U", -1, 19, 0, "hi");
				_test_evidence($ev[11], "ref1", 52, "T", 1, 1, 0, "K", -1, 22, 0, "hi");
				_test_evidence($ev[12], "ref1", 59, "C", 1, 1, 0, "U", -1, 29, 0, "hi");
				_test_evidence($ev[13], "ref1", 64, "C", 1, 1, 0, "U", -1, 34, 0, "hi");
				# BSC
				_test_evidence($ev[14], "ref1", 50, "A", 0, 1, 0, "R", -1, 24, 0, "hi");
				_test_evidence($ev[15], "ref1", 60, "A", 0, 1, 0, "B", -1, 14, 0, "hi");
				_test_evidence($ev[16], "ref1", 61, "A", 0, 1, 0, "N", -1, 13, 0, "hi");
				_test_evidence($ev[17], "ref1", 63, "G", 0, 1, 0, "J", -1, 11, 0, "hi");
				_test_evidence($ev[18], "ref1", 66, "G", 0, 1, 0, "B", -1,  8, 0, "hi");
				_test_evidence($ev[19], "ref1", 68, "G", 0, 1, 0, "G", -1,  6, 0, "hi");
				# BSCR
				_test_evidence($ev[20], "ref1",  3, "A", 0, 0, 0, "D", -1,  3, 0, "hi");
				_test_evidence($ev[21], "ref1",  5, "G", 0, 0, 0, "A", -1,  5, 0, "hi");
				_test_evidence($ev[22], "ref1", 14, "G", 0, 0, 0, "U", -1, 14, 0, "hi");
				_test_evidence($ev[23], "ref1", 24, "A", 0, 0, 0, "H", -1, 24, 0, "hi");
				_test_evidence($ev[24], "ref1", 25, "A", 0, 0, 0, "B", -1, 25, 0, "hi");
				_test_evidence($ev[25], "ref1", 28, "A", 0, 0, 0, "u", -1, 28, 0, "hi");
				_test_evidence($ev[26], "ref1", 30, "G", 0, 0, 0, "h", -1, 30, 0, "hi");
				_test_evidence($ev[27], "ref1", 31, "G", 0, 0, 0, "g", -1, 31, 0, "hi");
				_test_evidence($ev[28], "ref1", 32, "G", 0, 0, 0, "w", -1, 32, 0, "hi");
				_test_evidence($ev[29], "ref1", 33, "G", 0, 0, 0, "o", -1, 33, 0, "hi");
				_test_evidence($ev[30], "ref1", 35, "A", 0, 0, 0, "i", -1, 35, 0, "hi");
				_test_evidence($ev[31], "ref1", 38, "A", 0, 0, 0, "a", -1, 38, 0, "hi");
			} else {
				scalar(@ev) == 14 || die "Expected 14 pieces of evidence, got ".scalar(@ev).":".Dumper(\@ev);
				@ev = sort {$a->rdid cmp $b->rdid or $a->ref_off <=> $b->ref_off} @ev;
				# BSWR
				_test_evidence($ev[0],  "ref1",  4, "C", 1, 0, 0, "L", -1, 22, 0, "hi");
				_test_evidence($ev[1],  "ref1", 23, "T", 1, 0, 0, "B", -1,  3, 0, "hi");
				# BSW
				_test_evidence($ev[2],  "ref1", 34, "T", 1, 1, 0, "K", -1,  4, 0, "hi");
				_test_evidence($ev[3],  "ref1", 37, "C", 1, 1, 0, "A", -1,  7, 0, "hi");
				_test_evidence($ev[4],  "ref1", 39, "T", 1, 1, 0, "D", -1,  9, 0, "hi");
				_test_evidence($ev[5],  "ref1", 49, "T", 1, 1, 0, "U", -1, 19, 0, "hi");
				_test_evidence($ev[6],  "ref1", 59, "C", 1, 1, 0, "U", -1, 29, 0, "hi");
				# BSC
				_test_evidence($ev[7],  "ref1", 50, "A", 0, 1, 0, "R", -1, 24, 0, "hi");
				_test_evidence($ev[8],  "ref1", 60, "A", 0, 1, 0, "B", -1, 14, 0, "hi");
				# BSCR
				_test_evidence($ev[9],  "ref1",  5, "G", 0, 0, 0, "A", -1,  5, 0, "hi");
				_test_evidence($ev[10], "ref1", 24, "A", 0, 0, 0, "H", -1, 24, 0, "hi");
				_test_evidence($ev[11], "ref1", 30, "G", 0, 0, 0, "h", -1, 30, 0, "hi");
				_test_evidence($ev[12], "ref1", 35, "A", 0, 0, 0, "i", -1, 35, 0, "hi");
				_test_evidence($ev[13], "ref1", 38, "A", 0, 0, 0, "a", -1, 38, 0, "hi");
			}
			if($i > 0) {
				# Tabulate
				rmtree($tab_dir);
				my $tab_cmd =
					"$tabulate_cmd_base ".
					"--cpg-out=$tab_dir ".
					"-- $sort_dir ".
					"-- $fa1 ";
				print STDERR "$tab_cmd\n";
				system($tab_cmd);
				$? == 0 || die "Bad exitlevel: $?";
				#my %mtab = ();
				#parse_meth_table_dirname($tab_dir, \%mtab);
				#scalar(keys %mtab) == 2 || die "Bad number of keys: ".Dumper(%mtab);
				#defined($mtab{"ref1"}{12}) || die;
				#defined($mtab{"ref1"}{15}) || die;
				#defined($mtab{"ref2"}{6})  || die;
				#defined($mtab{"ref2"}{13}) || die;
			}
		}
		
		print STDERR "PASSED $type\n";
		for (<$name*.bt2>) { unlink($_) };
		rmtree($ev_dir);
		if($do_extract) {
			for my $f (<$sam_temp*>) { unlink($f); }
		}
	}
	print STDERR "PASSED\n";
}
unlink($fa1);
unlink($fq1);
print STDERR "PASSED ALL\n";
