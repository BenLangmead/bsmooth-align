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
my $mbias_cmd_base    = "perl $Bin/../../bin/bsev_mbias.pl";
my $tabulate_cmd_base = "perl $Bin/../../bin/bsev_tabulate.pl";

#                   *   *+ *+
#            1234567890123456789012345
my $ref1  = "ATCGTGTCATACGGCGATGCTAGCG";
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
my $read3 =        "TTTTTCGAAGTGCAAT";
#                   0123456789012345
#                     U M MUUM
my $qual2 =        "ABCDEFGHIJKLMNOP";
my $qual3 =        "PONMLKJIHGFEDCBA";
# C evidence idx:     3 4 5678
# CpG evidence idx:   2      3

my $fa1 = temp_fasta_filename(
	[ $ref1,  $ref2  ],
	[ "ref1", "ref2" ]
);

my $fq1 = temp_fastq_filename(
	[ $read1,  $read2,  $read3 ],
	[ $qual1,  $qual2,  $qual3 ],
	[ "read1", "read2", "read3"]
);

my $ev_dir = ".tmp.ev";
rmtree($ev_dir);
my $sort_dir = ".tmpsort.ev";
rmtree($sort_dir);
my $tab_dir = ".tmptab.ev";
rmtree($tab_dir);
my $mbias_dir = ".tmpmbias.ev";
rmtree($mbias_dir);

sub _test_evidence($$$$$$$$$$$$) {
	my ($ev, $tname, $off, $al, $watson, $fw, $flags, $qua1, $qua2, $cyc,
	    $score, $mapq) = @_;
	$ev->ref_name eq $tname  || croak("Expected refname=$tname, got " .$ev->ref_name);
	$ev->ref_off == $off     || croak("Expected offset=$off, got "    .$ev->ref_off);
	$ev->ev eq $al           || croak("Expected al=$al, got "         .$ev->ev);
	$ev->watson == $watson   || croak("Expected watson=$watson, got " .$ev->watson);
	$ev->fw == $fw           || croak("Expected fw=$fw, got "         .$ev->fw);
	$ev->flags == $flags     || croak("Expected flags=$flags, got "   .$ev->flags);
	$ev->qual1 eq $qua1      || croak("Expected qual1=$qua1, got "    .$ev->qual1);
	$ev->qual2 eq $qua2      || croak("Expected qual2=$qua2, got "    .$ev->qual2);
	$ev->cycle == $cyc       || croak("Expected cycle=$cyc, got "     .$ev->cycle);
	$ev->aln_score == $score || croak("Expected score=$score, got "   .$ev->aln_score);
	if($mapq eq "hi") {
		$ev->aln_mapq >= 20  || croak("Expected mapq >= 20, got "     .$ev->aln_mapq);
	} elsif($mapq eq "lo") {
		$ev->aln_mapq < 20   || croak("Expected mapq < 20, got "      .$ev->aln_mapq);
	} else {
		$ev->aln_mapq  == $mapq  || croak("Expected mapq=$mapq, got "     .$ev->aln_mapq);
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
				"--four-strand ".
				"--$do_extract=$sam_temp ".
				"$args ".
				"-- $name ".  # index
				"-- $fa1 ".   # refs
				"-- ".        # Bowtie 2 args
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
				"--four-strand ".
				"--output=$ev_dir ".
				"$args ".
				"-- $name ".  # index
				"-- $fa1 ".   # refs
				"-- ".        # Bowtie 2 args
				"-- $fq1";    # reads
			print STDERR "$align_cmd\n";
			system($align_cmd);
			$? == 0 || die "Bad exitlevel: $?";
		}
		
		# Generate M-bias table
		rmtree($mbias_dir);
		my $mbias_cmd =
			"$mbias_cmd_base ".
			"--ev=$ev_dir ".
			"--out=$mbias_dir";
		print STDERR "$mbias_cmd\n";
		system($mbias_cmd);
		$? == 0 || die "Bad exitlevel: $?";
		
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
				scalar(@ev) == 15 || die "Expected 15 pieces of evidence, got:".Dumper(\@ev);
				@ev = sort {$a->rdid cmp $b->rdid or $a->ref_off <=> $b->ref_off} @ev;
	 
				_test_evidence($ev[0],  "ref1",  8, "T", 1, 1, 0, "E", -1,  4, 0, "hi");
				_test_evidence($ev[1],  "ref1", 12, "T", 1, 1, 0, "I", -1,  8, 0, "hi");
				_test_evidence($ev[2],  "ref1", 15, "C", 1, 1, 0, "3", -1, 11, 0, "hi");
				
				_test_evidence($ev[3],  "ref2",  6, "T", 1, 1, 0, "C", -1,  2, 0, "hi");
				_test_evidence($ev[4],  "ref2",  8, "C", 1, 1, 0, "E", -1,  4, 0, "hi");
				_test_evidence($ev[5],  "ref2", 10, "C", 1, 1, 0, "G", -1,  6, 0, "hi");
				_test_evidence($ev[6],  "ref2", 11, "T", 1, 1, 0, "H", -1,  7, 0, "hi");
				_test_evidence($ev[7],  "ref2", 12, "T", 1, 1, 0, "I", -1,  8, 0, "hi");
				_test_evidence($ev[8],  "ref2", 13, "C", 1, 1, 0, "J", -1,  9, 0, "hi");
				_test_evidence($ev[9],  "ref2",  6, "T", 1, 0, 0, "C", -1, length($read3)-2-1, 0, "hi");
				_test_evidence($ev[10], "ref2",  8, "C", 1, 0, 0, "E", -1, length($read3)-4-1, 0, "hi");
				_test_evidence($ev[11], "ref2", 10, "C", 1, 0, 0, "G", -1, length($read3)-6-1, 0, "hi");
				_test_evidence($ev[12], "ref2", 11, "T", 1, 0, 0, "H", -1, length($read3)-7-1, 0, "hi");
				_test_evidence($ev[13], "ref2", 12, "T", 1, 0, 0, "I", -1, length($read3)-8-1, 0, "hi");
				_test_evidence($ev[14], "ref2", 13, "C", 1, 0, 0, "J", -1, length($read3)-9-1, 0, "hi");
			} else {
				scalar(@ev) == 6 || die "Expected 6 pieces of evidence, got:".Dumper(\@ev);
				@ev = sort {$a->rdid cmp $b->rdid or $a->ref_off <=> $b->ref_off} @ev;
				_test_evidence($ev[0], "ref1", 12, "T", 1, 1, 0, "I", -1,  8, 0, "hi");
				_test_evidence($ev[1], "ref1", 15, "C", 1, 1, 0, "3", -1, 11, 0, "hi");
				_test_evidence($ev[2], "ref2",  6, "T", 1, 1, 0, "C", -1,  2, 0, "hi");
				_test_evidence($ev[3], "ref2", 13, "C", 1, 1, 0, "J", -1,  9, 0, "hi");
				_test_evidence($ev[4], "ref2",  6, "T", 1, 0, 0, "C", -1, length($read3)-2-1, 0, "hi");
				_test_evidence($ev[5], "ref2", 13, "C", 1, 0, 0, "J", -1, length($read3)-9-1, 0, "hi");
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
				my %mtab = ();
				parse_meth_table_dirname($tab_dir, \%mtab);
				scalar(keys %mtab) == 2 || die "Bad number of keys: ".Dumper(%mtab);
				defined($mtab{"ref1"}{12}) || die;
				defined($mtab{"ref1"}{15}) || die;
				defined($mtab{"ref2"}{6})  || die;
				defined($mtab{"ref2"}{13}) || die;
				
				for my $ii (3, 24) {
					$mtab{"ref1"}{$ii}->[0] eq "M" || confess(Dumper(\%mtab));
					$mtab{"ref1"}{$ii}->[1] eq ""  || confess(Dumper(\%mtab));
					$mtab{"ref1"}{$ii}->[3] eq ""  || confess(Dumper(\%mtab));
					$mtab{"ref1"}{$ii}->[2] == 0   || confess(Dumper(\%mtab));
					$mtab{"ref1"}{$ii}->[4] == 0   || confess(Dumper(\%mtab));
					$mtab{"ref1"}{$ii}->[5] == 0   || confess(Dumper(\%mtab));
					$mtab{"ref1"}{$ii}->[6] == 0   || confess(Dumper(\%mtab));
					$mtab{"ref1"}{$ii}->[7] == 0   || confess(Dumper(\%mtab));
					$mtab{"ref1"}{$ii}->[8] == 0   || confess(Dumper(\%mtab));
				}
				
				$mtab{"ref1"}{12}->[0] eq "M" || confess(Dumper(\%mtab));
				$mtab{"ref1"}{12}->[1] eq ""  || confess(Dumper(\%mtab));
				$mtab{"ref1"}{12}->[3] eq "I" || confess(Dumper(\%mtab));
				$mtab{"ref1"}{12}->[2] == 0   || confess(Dumper(\%mtab));
				$mtab{"ref1"}{12}->[4] == 1   || confess(Dumper(\%mtab));
				$mtab{"ref1"}{12}->[5] == 0   || confess(Dumper(\%mtab));
				$mtab{"ref1"}{12}->[6] == 0   || confess(Dumper(\%mtab));
				$mtab{"ref1"}{12}->[7] == 0   || confess(Dumper(\%mtab));
				$mtab{"ref1"}{12}->[8] == 0   || confess(Dumper(\%mtab));

				$mtab{"ref1"}{15}->[0] eq "M" || confess(Dumper(\%mtab));
				$mtab{"ref1"}{15}->[1] eq "3" || confess(Dumper(\%mtab));
				$mtab{"ref1"}{15}->[3] eq ""  || confess(Dumper(\%mtab));
				$mtab{"ref1"}{15}->[2] == 1   || confess(Dumper(\%mtab));
				$mtab{"ref1"}{15}->[4] == 0   || confess(Dumper(\%mtab));
				$mtab{"ref1"}{15}->[5] == 0   || confess(Dumper(\%mtab));
				$mtab{"ref1"}{15}->[6] == 0   || confess(Dumper(\%mtab));
				$mtab{"ref1"}{15}->[7] == 0   || confess(Dumper(\%mtab));
				$mtab{"ref1"}{15}->[8] == 0   || confess(Dumper(\%mtab));
				
				for my $ii (1) {
					$mtab{"ref2"}{$ii}->[0] eq "M" || confess(Dumper(\%mtab));
					$mtab{"ref2"}{$ii}->[1] eq ""  || confess(Dumper(\%mtab));
					$mtab{"ref2"}{$ii}->[3] eq ""  || confess(Dumper(\%mtab));
					$mtab{"ref2"}{$ii}->[2] == 0   || confess(Dumper(\%mtab));
					$mtab{"ref2"}{$ii}->[4] == 0   || confess(Dumper(\%mtab));
					$mtab{"ref2"}{$ii}->[5] == 0   || confess(Dumper(\%mtab));
					$mtab{"ref2"}{$ii}->[6] == 0   || confess(Dumper(\%mtab));
					$mtab{"ref2"}{$ii}->[7] == 0   || confess(Dumper(\%mtab));
					$mtab{"ref2"}{$ii}->[8] == 0   || confess(Dumper(\%mtab));
				}
				
				$mtab{"ref2"}{6}->[0] eq "M"  || confess(Dumper(\%mtab));
				$mtab{"ref2"}{6}->[1] eq ""   || confess(Dumper(\%mtab));
				$mtab{"ref2"}{6}->[3] eq "CC" || confess(Dumper(\%mtab));
				$mtab{"ref2"}{6}->[2] == 0    || confess(Dumper(\%mtab));
				$mtab{"ref2"}{6}->[4] == 2    || confess(Dumper(\%mtab));
				$mtab{"ref2"}{6}->[5] == 0    || confess(Dumper(\%mtab));
				$mtab{"ref2"}{6}->[6] == 0    || confess(Dumper(\%mtab));
				$mtab{"ref2"}{6}->[7] == 0    || confess(Dumper(\%mtab));
				$mtab{"ref2"}{6}->[8] == 0    || confess(Dumper(\%mtab));

				$mtab{"ref2"}{13}->[0] eq "M"  || confess(Dumper(\%mtab));
				$mtab{"ref2"}{13}->[1] eq "JJ" || confess(Dumper(\%mtab));
				$mtab{"ref2"}{13}->[3] eq ""   || confess(Dumper(\%mtab));
				$mtab{"ref2"}{13}->[2] == 2    || confess(Dumper(\%mtab));
				$mtab{"ref2"}{13}->[4] == 0    || confess(Dumper(\%mtab));
				$mtab{"ref2"}{13}->[5] == 0    || confess(Dumper(\%mtab));
				$mtab{"ref2"}{13}->[6] == 0    || confess(Dumper(\%mtab));
				$mtab{"ref2"}{13}->[7] == 0    || confess(Dumper(\%mtab));
				$mtab{"ref2"}{13}->[8] == 0    || confess(Dumper(\%mtab));
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
