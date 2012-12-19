#!/usr/bin/env perl

##
# sim.pl
#
# Very simple script to simulate a small genome and some paired-end bisulfite
# reads.
#
#  Author: Ben Langmead
#    Date: 9/21/2011
# Contact: langmea@cs.jhu.edu
#

use strict;
use warnings;
use Math::Random qw(random_normal);

my %revcompMap = (
	"A" => "T",
	"T" => "A",
	"C" => "G",
	"G" => "C",
	"R" => "Y",
	"Y" => "R",
	"M" => "K",
	"K" => "M",
	"S" => "S",
	"W" => "W",
	"B" => "V",
	"V" => "B",
	"H" => "D",
	"D" => "H",
	"N" => "N",
	"." => ".",
	"a" => "t",
	"t" => "a",
	"c" => "g",
	"g" => "c",
	"r" => "y",
	"y" => "r",
	"m" => "k",
	"k" => "m",
	"s" => "s",
	"w" => "w",
	"b" => "v",
	"v" => "b",
	"h" => "d",
	"d" => "h",
	"n" => "n",
);

sub comp($) {
	my $ret = $revcompMap{$_[0]} || die "Can't reverse-complement '$_[0]'";
	return $ret;
}

sub revcomp {
	my ($ret, $color) = @_;
	$ret = reverse $ret;
	unless($color) {
		for(my $i = 0; $i < length($ret); $i++) {
			substr($ret, $i, 1) = comp(substr($ret, $i, 1));
		}
	}
	return $ret;
}

my $out_dir = "sim_output";
my $nreads = 200;
my $rdlen = 100;
my $frag_av = 300;
my $frag_sd = 25;
my @fraglens = random_normal($nreads, $frag_av, $frag_sd);

open(RD1, ">reads_1.fq") || die;
open(RD2, ">reads_2.fq") || die;

sub rand_dna($) {
	my $ret = "";
	for(1..$_[0]) { $ret .= substr("ACGT", int(rand(4)), 1); }
	return $ret;
}

sub rand_quals($) {
	my $ret = "";
	for(1..$_[0]) { $ret .= chr(33+int(rand(40))); }
	return $ret;
}

my %gn = ();
my $nchrs = int(rand(8));
my $concat = "";
for(0..$nchrs) {
	my $name = "sim".($_+1);
	$gn{$name} = rand_dna(int(rand(800))+1);
	my $fn = "ref_".($_+1).".fa";
	open(RF, ">$fn") || die;
	print RF ">$name\n$gn{$name}\n";
	close(RF);
	print STDERR "Made reference: $fn\n";
	$concat .= $gn{$name};
	$concat .= ("N" x int(rand(100)));
}

my $bsref_wat = $concat;
my $bsref_cri = revcomp($concat);

# Turn CpG Cs to Ys and remaining Cs to Ts
$bsref_wat =~ s/CG/YG/g; $bsref_wat =~ s/C/T/g;
$bsref_cri =~ s/CG/YG/g; $bsref_cri =~ s/C/T/g;

$concat = $bsref_wat.("N" x int(rand(100))).$bsref_cri;

my $rflen = length($concat);
for(0..$#fraglens) {
	my $flen = $fraglens[$_];
	my $off = int(rand($rflen - ($flen-1)));
	my $fstr = substr($concat, $off, $flen);
	# Select methylation state for Ys
	for(0..length($fstr)-1) {
		my $c = substr($fstr, $_, 1);
		substr($fstr, $_, 1) = ((rand() > 0.5) ? "C" : "T") if $c eq "Y";
	}
	$fstr = revcomp($fstr) if int(rand(2)) == 0;
	my $rd1 = substr($fstr, 0, $rdlen);
	my $rd2 = substr($fstr, length($fstr)-$rdlen);
	$rd2 = revcomp($rd2);
	my $qu1 = rand_quals($rdlen);
	my $qu2 = rand_quals($rdlen);
	print RD1 "\@r".($_+1)."\n$rd1\n+\n$qu1\n";
	print RD2 "\@r".($_+1)."\n$rd2\n+\n$qu2\n";
}

close(RD1);
close(RD2);
print STDERR "Made reads: reads_1.fq/reads_2.fq\n";

print STDERR "DONE\n";
