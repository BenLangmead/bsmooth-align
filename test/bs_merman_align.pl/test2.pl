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
use File::Path;
use Data::Dumper;

my $align_cmd_base = "perl $Bin/../../bin/bs_merman_align.pl";
my $extract_cmd_base = "perl $Bin/../../bin/bs_sam_extract.pl";

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

sub _test_evidence($$$$$$$$$$$) {
	my ($ev, $tname, $off, $al, $watson,
	    $fw, $qua1, $qua2, $cyc, $score,
	    $mapq) = @_;
	$off++; # change 0-based ref offset to 1-based
	$ev->ref_name eq $tname  || croak("Expected refname=$tname, got " .$ev->ref_name);
	$ev->ref_off == $off     || croak("Expected offset=$off, got "    .$ev->ref_off);
	$ev->ev eq $al           || croak("Expected al=$al, got "         .$ev->ev);
	$ev->watson == $watson   || croak("Expected watson=$watson, got " .$ev->watson);
	$ev->fw == $fw           || croak("Expected fw=$fw, got "         .$ev->fw);
	$ev->qual1 eq $qua1      || croak("Expected qual1=$qua1, got "    .$ev->qual1);
	$ev->qual2 eq $qua2      || croak("Expected qual2=$qua2, got "    .$ev->qual2);
	$ev->cycle == $cyc       || croak("Expected cycle=$cyc, got "     .$ev->cycle);
	$ev->aln_score == $score || croak("Expected score=$score, got "   .$ev->aln_score);
	$ev->aln_mapq  == $mapq  || croak("Expected mapq=$mapq, got "     .$ev->aln_mapq);
	return 1;
}

for my $do_extract (undef, "sam", "bam") {
	for my $type ("bscpg", "bsc") {
		my $name = "test1";
		my $args = "--$type";
		my $sam_temp = temp_name();
		
		if($do_extract) {
			my $align_cmd =
				"$align_cmd_base ".
				"--echo-sam ".
				"--keep-$do_extract=$sam_temp ".
				"--stop-after-alignment ".
				"--four-strand ".
				"$args ".
				"-- $fa1 ".
				"-- -v 4 -l 10 -L 10 -w 10 -M 1 ".
				"-- $fq1";
			print STDERR "$align_cmd\n";
			system($align_cmd);
			$? == 0 || die "Bad exitlevel: $?";
			
			my $extract_cmd =
				"$extract_cmd_base ".
				"--echo-sam ".
				"--output=$ev_dir ".
				"$args ".
				"$sam_temp.$do_extract";
			print STDERR "$extract_cmd\n";
			system($extract_cmd);
			$? == 0 || die "Bad exitlevel: $?";
		} else {
			my $align_cmd =
				"$align_cmd_base ".
				"--echo-sam ".
				"--output=$ev_dir ".
				"--four-strand ".
				"$args ".
				"-- $fa1 ".
				"-- -v 4 -l 10 -L 10 -w 10 -M 1 ".
				"-- $fq1";
			print STDERR "$align_cmd\n";
			system($align_cmd);
			$? == 0 || die "Bad exitlevel: $?";
		}
		
		my @ev = ();
		parse_meth_ev_records_dirname($ev_dir, \@ev);
		if($type eq "bsc") {
			scalar(@ev) == 32 || die "Expected 32 pieces of evidence, got ".scalar(@ev).":".Dumper(\@ev);
			# C Ev:   ref1,  4, "C", 1, 0, "L", -1, 22, 0, 255
			# C Ev:   ref1, 15, "T", 1, 0, "F", -1, 11, 0, 255
			# C Ev:   ref1, 18, "C", 1, 0, "E", -1,  8, 0, 255
			# C Ev:   ref1, 21, "T", 1, 0, "C", -1,  5, 0, 255
			# C Ev:   ref1, 23, "T", 1, 0, "B", -1,  3, 0, 255
			# C Ev:   ref1, 34, "T", 1, 1, "K", -1,  4, 0, 255
			# C Ev:   ref1, 36, "C", 1, 1, "L", -1,  6, 0, 255
			# C Ev:   ref1, 37, "C", 1, 1, "A", -1,  7, 0, 255
			# C Ev:   ref1, 39, "T", 1, 1, "D", -1,  9, 0, 255
			# C Ev:   ref1, 45, "T", 1, 1, "S", -1, 15, 0, 255
			# C Ev:   ref1, 49, "T", 1, 1, "U", -1, 19, 0, 255
			# C Ev:   ref1, 52, "T", 1, 1, "K", -1, 22, 0, 255
			# C Ev:   ref1, 59, "C", 1, 1, "U", -1, 29, 0, 255
			# C Ev:   ref1, 64, "C", 1, 1, "U", -1, 34, 0, 255
			# C Ev:   ref1, 50, "A", 0, 1, "R", -1, 24, 0, 255
			# C Ev:   ref1, 60, "A", 0, 1, "B", -1, 14, 0, 255
			# C Ev:   ref1, 61, "G", 0, 1, "N", -1, 13, 0, 255
			# C Ev:   ref1, 63, "G", 0, 1, "J", -1, 11, 0, 255
			# C Ev:   ref1, 66, "G", 0, 1, "B", -1,  8, 0, 255
			# C Ev:   ref1, 68, "G", 0, 1, "G", -1,  6, 0, 255
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
			#_test_evidence($ev[0], "ref1",  8, "T", 1, 1, "E", -1,  4, 0, 255);
		} else {
			scalar(@ev) == 14 || die "Expected 14 pieces of evidence, got ".scalar(@ev).":".Dumper(\@ev);
			# CpG Ev: ref1,  4, "C", 1, 0, "L", -1, 22, 0, 255
			# CpG Ev: ref1, 23, "T", 1, 0, "B", -1,  3, 0, 255
			# CpG Ev: ref1, 34, "T", 1, 1, "K", -1,  4, 0, 255
			# CpG Ev: ref1, 37, "C", 1, 1, "A", -1,  7, 0, 255
			# CpG Ev: ref1, 39, "T", 1, 1, "D", -1,  9, 0, 255
			# CpG Ev: ref1, 49, "T", 1, 1, "U", -1, 19, 0, 255
			# CpG Ev: ref1, 59, "C", 1, 1, "U", -1, 29, 0, 255
			# CpG Ev: ref1, 50, "A", 0, 1, "R", -1, 24, 0, 255
			# CpG Ev: ref1, 60, "A", 0, 1, "B", -1, 14, 0, 255
			# CpG Ev: ref1,  5, "G", 0, 0, "A", -1, 14, 0, 255
			# CpG Ev: ref1, 24, "A", 0, 0, "H", -1, 13, 0, 255
			# CpG Ev: ref1, 30, "G", 0, 0, "h", -1, 13, 0, 255
			# CpG Ev: ref1, 35, "A", 0, 0, "i", -1, 13, 0, 255
			# CpG Ev: ref1, 38, "A", 0, 0, "a", -1, 13, 0, 255
			#_test_evidence($ev[0], "ref1", 12, "T", 1, 1, "I", -1,  8, 0, 255);
		}
		
		print STDERR "PASSED $type\n";
		for (<$name*.ebwt>) { unlink($_) };
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
