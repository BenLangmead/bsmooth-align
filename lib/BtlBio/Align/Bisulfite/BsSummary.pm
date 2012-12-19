#!/usr/bin/perl -w

#
# Copyright (C) 2012, Ben Langmead <langmea@cs.jhu.edu>
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

##
# BtlBio::Align::Bisulfite::BsSummary
#
# Author: Ben Langmead
# Email: langmea@cs.jhu.edu
#

package BtlBio::Align::Bisulfite::BsSummary;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(bs_summary_blank_tsv
             bs_summary_blank_tsv_nostrata
             parse_meth_table
             parse_meth_table_filename
             parse_meth_table_filehandle
             parse_meth_table_dirname);

use strict;
use warnings;
use Carp;
use FindBin qw($Bin);
use lib "$Bin/../../..";
use BtlBio::Align::Bisulfite::BsEvidence;
use BtlBio::Align::Bisulfite::BsFilter;
use BtlBio::Util::File;

##
# Create a new bisulfite evidence object.  An evidence object is just an array
# ref with fields for:
#
sub new {
	my ($class) = @_;
	my $ret = bless {
		_count     => {},   # allele => count
		_qual      => {},   # allele & quality strat => count
		_cycs      => {},   # allele & cyc => count
		_qualstr   => {},
		_mapqstr   => {},
		_filt_mapq => 0,
		_filt_qual => 0,
		_filt_cyc  => 0,
		_filt_rdl  => 0,
		_filt_nuc  => 0,
		_empty     => 1
	}, $class;
	$ret->rep_ok || die;
	return $ret;
}

sub rep_ok($$$) {
	my ($self, $qual_strata, $max_strat) = @_;
	for my $k (keys %{$self->{_qual}}) {
		defined($self->{_count}{$k}) || confess;
		defined($self->{_cycs}{$k}) || confess;
		scalar keys %{$self->{_cycs}{$k}} <= $self->{_count}{$k} || confess;
		my $tot = 0;
		my $idx = 0;
		for my $v (@{$self->{_qual}{$k}}) {
			$tot += ($v || 0);
		}
		$tot == $self->{_count}{$k} ||
			confess "tot=$tot, selfcount=$self->{_count}{$k}, max_strat=$max_strat, qual_strata=$qual_strata";
	}
	return 1;
}

##
# Add a new piece of evidence to this summary.
#
sub add($$$$$) {
	my ($self, $ev, $filter, $qual_strata, $max_strat) = @_;
	$qual_strata = 1 unless defined($qual_strata);
	$max_strat = 50 unless defined($max_strat);
	my ($filt_mapq, $filt_qual, $filt_cyc, $filt_rdl, $filt_nuc) = $filter->filter($ev);
	defined($filt_mapq) || die;
	defined($filt_qual) || die;
	defined($filt_cyc) || die;
	defined($filt_rdl) || die;
	defined($filt_nuc) || die;
	#print STDERR "Filter: $filt_mapq, $filt_qual, $filt_cyc, $filt_rdl, $filt_nuc\n";
	my $filt = $filt_mapq || $filt_qual || $filt_cyc || $filt_rdl || $filt_nuc;
	#print STDERR "Filter: $filt\n";
	$self->{_filt_mapq}++ if $filt_mapq;
	$self->{_filt_qual}++ if $filt_qual;
	$self->{_filt_cyc}++  if $filt_cyc;
	$self->{_filt_rdl}++  if $filt_rdl;
	$self->{_filt_nuc}++  if $filt_nuc;
	return 0 if $filt;
	# Evidence not filtered
	my $strat = ($ev->{_qual_summ}/$qual_strata);
	$strat = $max_strat if $strat > $max_strat;
	my $al = $ev->{_ev};
	$self->{_count}{$al}++;
	# Tally the base quality stratum
	$self->{_qual}{$al}[$strat]++;
	# Tally the cycle
	$self->{_cycs}{$al}{$ev->{_cy}}++;
	$self->{_empty} = 0;
	return 1;
}

##
# Merge this summary with another, placing the result in this summary.
#
sub merge_inplace {
	my ($self, $summ, $different_strand) = @_;
	# For each allele
	for my $k (keys %{$summ->{_qual}}) {
		# For each stratum
		for(my $i = 0; $i < scalar(@{$summ->{_qual}{$k}}); $i++) {
			if(defined($summ->{_qual}{$k}[$i])) {
				$self->{_qual}{$k}[$i] += $summ->{_qual}{$k}[$i];
			}
		}
		defined($summ->{_count}{$k}) || confess("Had _qual for $k, but not _count");
		$self->{_count}{$k} += $summ->{_count}{$k};
	}
	# For each allele
	for my $k (keys %{$summ->{_cycs}}) {
		# For each cycle
		for my $cy (keys %{$summ->{_cycs}{$k}}) {
			my $self_cy = $different_strand ? "+$cy" : "$cy";
			$self->{_cycs}{$k}{$self_cy} += $summ->{_cycs}{$k}{$cy};
		}
	}
	for my $k (keys %{$summ->{_qualstr}}) {
		$self->{_qualstr}{$k} = "" unless defined($self->{_qualstr}{$k});
		$self->{_qualstr}{$k} .= $summ->{_qualstr}{$k};
		$self->{_mapqstr}{$k} = "" unless defined($self->{_mapqstr}{$k});
		$self->{_mapqstr}{$k} .= $summ->{_mapqstr}{$k};
	}
	#$self->{_ncycs} += $summ->{_ncycs};
	$self->{_filt_mapq} += $summ->{_filt_mapq};
	$self->{_filt_qual} += $summ->{_filt_qual};
	$self->{_filt_cyc} += $summ->{_filt_cyc};
	$self->{_filt_rdl} += $summ->{_filt_rdl};
	$self->{_filt_nuc} += $summ->{_filt_nuc};
	$self->{_empty} = $self->{_empty} && $summ->{_empty};
}

##
# Return a tab-separated text summary of this summary.  Prints 16 fields.
#
sub summary_tsv($$$) {
	my ($self, $qual_strata, $max_strat) = @_;
	my @fs = ();
	my $tot = 0;
	my $nstrata = 0;
	for(my $i = 0; $i <= $max_strat; $i += $qual_strata) { $nstrata++; }
	if(defined($self->{_qual}{C})) {
		for(my $i = 0; $i < $nstrata; $i++) {
			push @fs, ($self->{_qual}{C}[$i] || 0); # stratified by quality
			$tot += $fs[$#fs];
		}
	} else {
		for(1..$nstrata) { push @fs, 0; }
	}
	if(defined($self->{_cycs}{C})) {
		my $ncyc = 0;
		for my $k (keys %{$self->{_cycs}{C}}) {
			$ncyc++ if $self->{_cycs}{C}{$k} > 0;
		}
		push @fs, $ncyc; # # distinct cycles
	} else {
		push @fs, 0;
	}
	if(defined($self->{_qual}{T})) {
		for(my $i = 0; $i < $nstrata; $i++) {
			push @fs, ($self->{_qual}{T}[$i] || 0); # stratified by quality
			$tot += $fs[$#fs];
		}
	} else {
		for(1..$nstrata) { push @fs, 0; }
	}
	if(defined($self->{_cycs}{T})) {
		my $ncyc = 0;
		for my $k (keys %{$self->{_cycs}{T}}) {
			$ncyc++ if $self->{_cycs}{T}{$k} > 0;
		}
		push @fs, $ncyc; # # distinct cycles
	} else {
		push @fs, 0;
	}
	push @fs, ($self->{_filt_cyc}  || 0); # filtered because of # cycles
	$tot += $fs[$#fs];
	push @fs, ($self->{_filt_rdl}  || 0); # filtered because of read length
	$tot += $fs[$#fs];
	push @fs, ($self->{_filt_nuc}  || 0); # filtered because nuc was not C/T
	$tot += $fs[$#fs];
	push @fs, ($self->{_filt_mapq} || 0); # filtered because of mapping quality
	$tot += $fs[$#fs];
	push @fs, ($self->{_filt_qual} || 0); # filtered because of alignment scoore
	$tot += $fs[$#fs];
	return (\@fs, $tot);
}

##
# Return a tab-separated text summary of this summary.  Prints 16 fields.
#
sub summary_tsv_nostrata($) {
	my ($self) = @_;
	my @fs = ();
	my $tot = 0;
	# Append a string encoding all the base qualities
	my $mstr = "";
	$mstr = $self->{_qualstr}{"C"} if defined($self->{_qualstr}{"C"});
	$tot += length($mstr);
	push @fs, $mstr;
	# Append Mcy
	if(defined($self->{_cycs}{C})) {
		my $ncyc = 0;
		for my $k (keys %{$self->{_cycs}{C}}) {
			$ncyc++ if $self->{_cycs}{C}{$k} > 0;
		}
		$ncyc <= length($mstr) || die "ncyc=$ncyc, len(mstr)=".length($mstr);
		push @fs, $ncyc; # # distinct cycles
	} else {
		push @fs, 0;
	}
	# Append a string encoding all the base qualities
	my $ustr = "";
	$ustr = $self->{_qualstr}{"T"} if defined($self->{_qualstr}{"T"});
	$tot += length($ustr);
	push @fs, $ustr;
	# Append Ucy
	if(defined($self->{_cycs}{T})) {
		my $ncyc = 0;
		for my $k (keys %{$self->{_cycs}{T}}) {
			$ncyc++ if $self->{_cycs}{T}{$k} > 0;
		}
		$ncyc <= length($ustr) || die "ncyc=$ncyc, len(ustr)=".length($ustr);
		push @fs, $ncyc; # # distinct cycles
	} else {
		push @fs, 0;
	}
	push @fs, ($self->{_filt_cyc}  || 0); # filtered because of # cycles
	$tot += $fs[$#fs];
	push @fs, ($self->{_filt_rdl}  || 0); # filtered because of read length
	$tot += $fs[$#fs];
	push @fs, ($self->{_filt_nuc}  || 0); # filtered because nuc was not C/T
	$tot += $fs[$#fs];
	push @fs, ($self->{_filt_mapq} || 0); # filtered because of mapping quality
	$tot += $fs[$#fs];
	push @fs, ($self->{_filt_qual} || 0); # filtered because of alignment scoore
	$tot += $fs[$#fs];
	return (\@fs, $tot);
}

sub num_measurements($) {
	my $self = shift;
	my $tot = 0;
	for my $k (keys %{$self->{_count}}) {
		$tot += $self->{_count}{$k};
	}
	$tot += $self->{_filt_mapq};
	$tot += $self->{_filt_qual};
	$tot += $self->{_filt_cyc};
	$tot += $self->{_filt_rdl};
	$tot += $self->{_filt_nuc};
}

sub bs_summary_blank_tsv($$) {
	my ($qual_strata, $max_strat) = @_;
	my $sz = 0;
	my $nstrata = 0;
	for(my $i = 0; $i <= $max_strat; $i += $qual_strata) { $nstrata++; }
	
	$sz += $nstrata; # Ms
	$sz++;           # Mcyc
	$sz += $nstrata; # Us
	$sz++;           # Ucyc
	
	$sz++;           # filt_cyc
	$sz++;           # filt_rdl
	$sz++;           # filt_nuc
	$sz++;           # filt_mapq
	$sz++;           # filt_qual
	return join("\t", (0) x $sz);
}

sub bs_summary_blank_tsv_nostrata() {
	# Mstr
	# Mcyc
	# Ustr
	# Ucyc
	# filt_cyc
	# filt_rdl
	# filt_nuc
	# filt_mapq
	# filt_qual
	return "\t0\t\t0\t0\t0\t0\t0\t0";
}

sub parse_meth_table($$) {
	my ($line, $h) = @_;
	chomp($line);
	my @ls = split(/\t/, $line);
	my %ret = ();
	if(scalar(@ls) == 12) {
		# No strata
		my ($ref, $off, $type, $mstr, $mcyc, $ustr, $ucyc,
		    $filt_cyc, $filt_rdl, $filt_nuc, $filt_mapq, $filt_baseq) = @ls;
		return if $off eq "off";
		$h->{$ref}{$off} = [$type, $mstr, $mcyc, $ustr, $ucyc, $filt_cyc, $filt_rdl, $filt_nuc, $filt_mapq, $filt_baseq];
	} else {
		
	}
}

##
# Given a filehandle for a methylation evidence file, parse it into a list of
# BsEvidence objects.
#
sub parse_meth_table_filehandle($$) {
	my ($fh, $h) = @_;
	while(readline $fh) {
		parse_meth_table($_, $h);
	}
}

##
# Given a filename for a methylation evidence file, parse it into a list of
# BsEvidence objects.
#
sub parse_meth_table_filename($$) {
	my ($fn, $h) = @_;
	return parse_meth_table_filehandle(openex($fn), $h);
}

##
# Given a directory name for a directory of methylation evidence, parse it into
# a list of BsEvidence objects.
#
sub parse_meth_table_dirname($$) {
	my ($dirn, $h) = @_;
	for my $fn (<$dirn/*.tsv>) { parse_meth_table_filename($fn, $h); }
}

#
# Simple tests
#

sub _test_1() {
	print STDERR "Testing BsSummary::add 1 ... ";
	my $filt = BtlBio::Align::Bisulfite::BsFilter->new(
		10,       # MAPQ min
		10,       # qual min
		0,        # cyc min
		9999999,  # cyc max
		3,        # cyc trim begin
		3,        # cyc trim end
		100);     # min read len
	my $ev = BtlBio::Align::Bisulfite::BsEvidence->new(
		"read1",  # read ID
		"ref1",   # ref name
		10,       # ref offset
		"C",      # allele
		1,        # watson
		1,        # fw
		1,        # flags
		chr(30+33),# qual 1
		chr(40+33),# qual 2
		9,        # cyc
		100,      # alignment len
		-55,      # alignment score
		40);      # MAPQ
	my $summ = BtlBio::Align::Bisulfite::BsSummary->new();
	$summ->add($ev, $filt, 10, 40);
	0 == ($summ->{_filt_mapq}  || 0) || die;
	0 == ($summ->{_filt_qual}  || 0) || die;
	0 == ($summ->{_filt_cyc}   || 0) || die;
	0 == ($summ->{_filt_rdl}   || 0) || die;
	0 == ($summ->{_filt_nuc}   || 0) || die;
	1 == ($summ->{_count}{C}   || 0) || die;
	1 == ($summ->{_qual}{C}[3] || 0) || die Dumper($summ);
	1 == ($summ->{_cycs}{C}{9} || 0) || die;
	print STDERR "PASSED\n";
}

sub _test_2() {
	print STDERR "Testing BsSummary::add 2 ... ";
	my $filt = BtlBio::Align::Bisulfite::BsFilter->new(
		10,       # MAPQ min
		10,       # qual min
		0,        # cyc min
		9999999,  # cyc max
		3,        # cyc trim begin
		3,        # cyc trim end
		101);     # min read length
	my $ev = BtlBio::Align::Bisulfite::BsEvidence->new(
		"read1",  # read ID
		"ref1",   # ref name
		10,       # ref offset
		"T",      # allele
		0,        # watson
		1,        # fw
		0,        # flags
		chr(7+33),# qual 1
		chr(4+33),# qual 2
		98,       # cyc
		99,       # alignment len
		-55,      # alignment score
		9);       # MAPQ
	my $summ = BtlBio::Align::Bisulfite::BsSummary->new();
	$summ->add($ev, $filt);
	1 == ($summ->{_filt_mapq}  || 0) || die $summ->{_filt_mapq};
	1 == ($summ->{_filt_qual}  || 0) || die $summ->{_filt_qual};
	1 == ($summ->{_filt_cyc}   || 0) || die $summ->{_filt_cyc};
	1 == ($summ->{_filt_rdl}   || 0) || die $summ->{_filt_rdl};
	1 == ($summ->{_filt_nuc}   || 0) || die $summ->{_filt_nuc};
	print STDERR "PASSED\n";
}

sub _test() {
	require Data::Dumper;
	_test_1();
	_test_2();
}

# Call _test() if the module is being run directly.
_test() unless caller();

1;
