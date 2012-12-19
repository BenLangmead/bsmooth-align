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
which("bowtie2") == 0 || die "Could not execute 'bowtie2'";
which("bowtie2-build") == 0 || die "Could not execute 'bowtie2-build'";

my $build_cmd_base    = "perl $Bin/../../bin/bswc_bowtie2_index.pl";
my $align_cmd_base    = "perl $Bin/../../bin/bswc_bowtie2_align.pl";
my $extract_cmd_base  = "perl $Bin/../../bin/bswc_sam_extract.pl";
my $sort_cmd_base     = "perl $Bin/../../bin/bsev_sort.pl";
my $tabulate_cmd_base = "perl $Bin/../../bin/bsev_tabulate.pl";

#                   *   *+ *+
#                     1         2         3         4         5         6
#            1234567890123456789012345678901234567890123456789012345678901234
#            ATCGTGTCATACGGCGATGCTAGCTNNNNNNNNNNNNCGCATCGCACCCCGAAAAACCCCCCCC
my $ref1  = "ATCGTGTCATACGGCGATGCTAGCTNNNNNNNNNNNNCGCATCGCACCCCGAAAAACCCCCCCC";
my $read1 =    "GTGTTATATGGCGA";
#               01234567890123
#                   U   U  M
my $qual1 =    "ABCDEFGHI12345";
# C evidence idx:   0   1  2
# CpG evidence idx:     0  1

#                                                      *+* ****+
#                     1         2         3         4         5         6
#            1234567890123456789012345678901234567890123456789012345678901234
#            ATCGTGTCATACGGCGATGCTAGCTNNNNNNNNNNNNCGCATCGCACCCCGAAAAACCCCCCCC
my $read2 =                                         "ATTGCACTTCGAAAAA";
#                                                    0123456789012345
#                                                      U M MUUM
my $qual2 =                                         "ABCDEFGHIJKLMNOP";
# C evidence idx:                                      3 4 5678
# CpG evidence idx:                                    2      3

my $fa1 = temp_fasta_filename(
	[ $ref1  ],
	[ "ref1" ]
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

for my $m1fw (1, 0) {
for my $m2fw (1, 0) {

my ($rd1, $qu1) = ($read1, $qual1);
if(!$m1fw) {
	$rd1 = revcomp($rd1);
	$qu1 = reverse $qu1;
}
my $fq1 = temp_fastq_filename([ $rd1 ], [ $qu1 ], [ "read1"]);

my ($rd2, $qu2) = ($read2, $qual2);
if(!$m2fw) {
	$rd2 = revcomp($rd2);
	$qu2 = reverse $qu2;
}
my $fq2 = temp_fastq_filename([ $rd2 ], [ $qu2 ], [ "read2"]);

my $mate_arg = "--";
$mate_arg .= ($m1fw ? "f" : "r");
$mate_arg .= ($m2fw ? "f" : "r");

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
				"$args ".
				"$mate_arg ".
				"-- $name ".
				"-- $fa1 ".
				"-- ".
				"-- $fq1 ".
				"-- $fq2";
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
				"$args ".
				"$mate_arg ".
				"-- $name ".
				"-- $fa1 ".
				"-- ".
				"-- $fq1 ".
				"-- $fq2";
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
			my $ref2plus = 37;
			my $thisRcFlag = 0;
			my $oppRcFlag = 0;
			if($type eq "bsc") {
				scalar(@ev) == 9 || die "Expected 9 pieces of evidence, got:".Dumper(\@ev);
				@ev = sort {$a->rdid cmp $b->rdid or $a->ref_off <=> $b->ref_off} @ev;
				#                                                al wat     fw ma  qu 
				_test_evidence($ev[0], "ref1",  8,             "T",  1, $m1fw, 1 | 2 | ($m1fw ? 0 : $thisRcFlag) | ($m2fw ? 0 : $oppRcFlag) | 64,  "E", -1, $m1fw ?  4 : 13- 4, 0, "hi");
				_test_evidence($ev[1], "ref1", 12,             "T",  1, $m1fw, 1 | 2 | ($m1fw ? 0 : $thisRcFlag) | ($m2fw ? 0 : $oppRcFlag) | 64,  "I", -1, $m1fw ?  8 : 13- 8, 0, "hi");
				_test_evidence($ev[2], "ref1", 15,             "C",  1, $m1fw, 1 | 2 | ($m1fw ? 0 : $thisRcFlag) | ($m2fw ? 0 : $oppRcFlag) | 64,  "3", -1, $m1fw ? 11 : 13-11, 0, "hi");
				
				_test_evidence($ev[3], "ref1",  6 + $ref2plus, "T",  1, $m2fw, 1 | 2 | ($m2fw ? 0 : $thisRcFlag) | ($m1fw ? 0 : $oppRcFlag) | 128, "C", -1, $m2fw ?  2 : 15- 2, 0, "hi");
				_test_evidence($ev[4], "ref1",  8 + $ref2plus, "C",  1, $m2fw, 1 | 2 | ($m2fw ? 0 : $thisRcFlag) | ($m1fw ? 0 : $oppRcFlag) | 128, "E", -1, $m2fw ?  4 : 15- 4, 0, "hi");
				_test_evidence($ev[5], "ref1", 10 + $ref2plus, "C",  1, $m2fw, 1 | 2 | ($m2fw ? 0 : $thisRcFlag) | ($m1fw ? 0 : $oppRcFlag) | 128, "G", -1, $m2fw ?  6 : 15- 6, 0, "hi");
				_test_evidence($ev[6], "ref1", 11 + $ref2plus, "T",  1, $m2fw, 1 | 2 | ($m2fw ? 0 : $thisRcFlag) | ($m1fw ? 0 : $oppRcFlag) | 128, "H", -1, $m2fw ?  7 : 15- 7, 0, "hi");
				_test_evidence($ev[7], "ref1", 12 + $ref2plus, "T",  1, $m2fw, 1 | 2 | ($m2fw ? 0 : $thisRcFlag) | ($m1fw ? 0 : $oppRcFlag) | 128, "I", -1, $m2fw ?  8 : 15- 8, 0, "hi");
				_test_evidence($ev[8], "ref1", 13 + $ref2plus, "C",  1, $m2fw, 1 | 2 | ($m2fw ? 0 : $thisRcFlag) | ($m1fw ? 0 : $oppRcFlag) | 128, "J", -1, $m2fw ?  9 : 15- 9, 0, "hi");
			} else {
				scalar(@ev) == 4 || die "Expected 4 pieces of evidence, got:".Dumper(\@ev);
				@ev = sort {$a->rdid cmp $b->rdid or $a->ref_off <=> $b->ref_off} @ev;
				#                                                al wat     fw ma  qu 
				
				_test_evidence($ev[0], "ref1", 12,             "T",  1, $m1fw, 1 | 2 | ($m1fw ? 0 : $thisRcFlag) | ($m2fw ? 0 : $oppRcFlag) | 64,  "I", -1, $m1fw ?  8 : 13- 8, 0, "hi");
				_test_evidence($ev[1], "ref1", 15,             "C",  1, $m1fw, 1 | 2 | ($m1fw ? 0 : $thisRcFlag) | ($m2fw ? 0 : $oppRcFlag) | 64,  "3", -1, $m1fw ? 11 : 13-11, 0, "hi");
				_test_evidence($ev[2], "ref1",  6 + $ref2plus, "T",  1, $m2fw, 1 | 2 | ($m2fw ? 0 : $thisRcFlag) | ($m1fw ? 0 : $oppRcFlag) | 128, "C", -1, $m2fw ?  2 : 15- 2, 0, "hi");
				_test_evidence($ev[3], "ref1", 13 + $ref2plus, "C",  1, $m2fw, 1 | 2 | ($m2fw ? 0 : $thisRcFlag) | ($m1fw ? 0 : $oppRcFlag) | 128, "J", -1, $m2fw ?  9 : 15- 9, 0, "hi");
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
unlink($fq1);
unlink($fq2);
}}
unlink($fa1);
print STDERR "PASSED ALL\n";
