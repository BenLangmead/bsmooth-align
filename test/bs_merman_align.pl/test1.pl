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

my $align_cmd_base    = "perl $Bin/../../bin/bs_merman_align.pl";
my $extract_cmd_base  = "perl $Bin/../../bin/bs_sam_extract.pl";
my $sort_cmd_base     = "perl $Bin/../../bin/bsev_sort.pl";
my $tabulate_cmd_base = "perl $Bin/../../bin/bsev_tabulate.pl";

#                   *   *+ *+
#            1234567890123456789012345
my $ref1  = "ATCGTGTCATACGGCGATGCTAGCT";
my $read1 =    "GTGTTATATGGCGA";
#               01234567890123
#                   U   U  M
my $qual1 =    "ABCDEFGHI12345";
# C evidence idx:   0   1  2
# CpG evidence idx:     0  1

#                     *+* ****+
#                123456789012345678901234567
my $ref2  =     "CGCATCGCACCCCGAAAAACCCCCCCC";
my $read2 =        "ATTGCACTTCGAAAAA";
#                   0123456789012345
#                     U M MUUM
my $qual2 =        "ABCDEFGHIJKLMNOP";
# C evidence idx:     3 4 5678
# CpG evidence idx:   2      3

my $fa1 = temp_fasta_filename(
	[ $ref1,  $ref2  ],
	[ "ref1", "ref2" ]
);

my $fq1 = temp_fastq_filename(
	[ $read1,  $read2 ],
	[ $qual1,  $qual2 ],
	[ "read1", "read2"]
);

my $ev_dir = ".tmp.ev";
rmtree($ev_dir);
my $sort_dir = ".tmp.sorted.tab";
rmtree($sort_dir);
my $tab_dir = ".tmp.tab";
rmtree($tab_dir);

sub _test_evidence($$$$$$$$$$$) {
	my ($ev, $tname, $off, $al, $watson,
	    $fw, $qua1, $qua2, $cyc, $score,
	    $mapq) = @_;
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
	for my $type ("bsc", "bscpg") {
		my $name = "test1";
		my $args = "--$type";
		my $sam_temp = temp_name();
		
		if($do_extract) {
			my $align_cmd =
				"$align_cmd_base ".
				"--echo-sam ".
				"--keep-$do_extract=$sam_temp ".
				"--stop-after-alignment ".
				"$args ".
				"-- $fa1 ".
				"-- -v 2 -M 1 -l 10 -L 10 -w 10 ".
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
				"$args ".
				"-- $fa1 ".
				"-- -v 2 -M 1 -l 10 -L 10 -w 10 ".
				"-- $fq1";
			print STDERR "$align_cmd\n";
			system($align_cmd);
			$? == 0 || die "Bad exitlevel: $?";
		}
		
		my @ev = ();
		parse_meth_ev_records_dirname($ev_dir, \@ev);
		if($type eq "bsc") {
			scalar(@ev) == 9 || die "Expected 9 pieces of evidence, got:".Dumper(\@ev);
			_test_evidence($ev[0], "ref1",  8, "T", 1, 1, "E", -1,  4, 0, 40);
			_test_evidence($ev[1], "ref1", 12, "T", 1, 1, "I", -1,  8, 0, 40);
			_test_evidence($ev[2], "ref1", 15, "C", 1, 1, "3", -1, 11, 0, 40);
			_test_evidence($ev[3], "ref2",  6, "T", 1, 1, "C", -1,  2, 0, 40);
			_test_evidence($ev[4], "ref2",  8, "C", 1, 1, "E", -1,  4, 0, 40);
			_test_evidence($ev[5], "ref2", 10, "C", 1, 1, "G", -1,  6, 0, 40);
			_test_evidence($ev[6], "ref2", 11, "T", 1, 1, "H", -1,  7, 0, 40);
			_test_evidence($ev[7], "ref2", 12, "T", 1, 1, "I", -1,  8, 0, 40);
			_test_evidence($ev[8], "ref2", 13, "C", 1, 1, "J", -1,  9, 0, 40);
		} else {
			scalar(@ev) == 4 || die "Expected 4 pieces of evidence, got:".Dumper(\@ev);
			_test_evidence($ev[0], "ref1", 12, "T", 1, 1, "I", -1,  8, 0, 40);
			_test_evidence($ev[1], "ref1", 15, "C", 1, 1, "3", -1, 11, 0, 40);
			_test_evidence($ev[2], "ref2",  6, "T", 1, 1, "C", -1,  2, -60, 40);
			_test_evidence($ev[3], "ref2", 13, "C", 1, 1, "J", -1,  9, -60, 40);
		}
		
		print STDERR "PASSED $type\n";
		
		# Sort the evidence
		my $sort_cmd =
			"$sort_cmd_base ".
			"--ev $ev_dir ".
			"--out $sort_dir";
		system($sort_cmd);
		$? == 0 || die "Bad exitlevel: $?";
		
		# Tabulate the evidence
		my $tab_cmd =
			"$tabulate_cmd_base ".
			"--cpg-out=$tab_dir ".
			"-- $sort_dir ".
			"-- $fa1";
		system($tab_cmd);
		$? == 0 || die "Bad exitlevel: $?";
		
		for (<$name*.ebwt>) { unlink($_) };
		rmtree($ev_dir);
		rmtree($sort_dir);
		rmtree($tab_dir);
		if($do_extract) {
			for my $f (<$sam_temp*>) { unlink($f); }
		}
	}
	
	print STDERR "PASSED\n";
}
unlink($fa1);
unlink($fq1);
print STDERR "PASSED ALL\n";
