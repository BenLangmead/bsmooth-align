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
# BtlBio::Alignment::Alignment
#
# Encapsulates an alignment.  The alignment is represented simply as the
# parallel query/reference strings with gap characters inserted.
#
# Author: Ben Langmead
# Email: langmea@cs.jhu.edu
#

package BtlBio::Alignment::Alignment;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(parse_generic_aln_record
             aln_opt_hash_to_string
             aln_string_to_opt_hash
             empty_alignment);

use strict;
use warnings;
use Carp;
use FindBin qw($Bin);
use lib "$Bin/../..";
use BtlBio::Alignment::Util;
use BtlBio::Read::Read;

##
# Create a new alignment object.  The object might represent an alignment to
# the genome, a record indicating that the read failed to align, or a record
# indicating that the read aligned repetitively.  If the object represents an
# alignment (repetitive or otherwise), then $aln->aligned() returns true.
#
# One subtlety is that SAM doesn't explicitly label alignments as being
# colorspace or not.  There are some tool-specific ways of detecting colorspace
# alignemnts; e.g. by checking for presence of the CS:Z and CQ:Z optional
# flags.  It's up to the user of these objects to check whether this is being
# deduced correctly for their app.
#
# BIG BIG PROBLEM: Doesn't save enough info to recreate a SAM record (RNEXT,
# PNEXT)
#
sub new {
	my (
		$class,
		$rdname,  # read name
		$rdseq,   # read sequence from alignment record
		$rdqual,  # read qualities from alignment record
		$color,   # whether read is colorspace
		$tname,   # name of reference sequence
		$toff,    # offset into reference
		$fw,      # true -> aligned to Watson strand
		$score,   # alignment score
		$mapq,    # mapping quality
		$md,      # MD:Z string
		$cigar,   # CIGAR string
		$mate1,   # is mate #1 in a pair?
		$mate2,   # is mate #2 in a pair?
		$flags,   # SAM flags (or undef if unavailable)
		$opts,    # options hash ref
		$orig,    # original string
		$unaligned # force to unaligned status
	) = @_;
	if($unaligned || (defined($cigar) && $cigar eq "*")) {
		$tname = undef;
		$toff  = undef;
		$md    = undef;
		$cigar = undef;
		$mapq  = undef;
	}
	defined($tname) == defined($toff)  || croak();
	defined($tname) == defined($cigar) || croak();
	defined($mate1) || croak();
	defined($mate2) || croak();
	my $ref_known = defined($md);
	my ($alread, $alref, $rdex, $rfex, $htriml, $striml, $htrimr, $strimr);
	my $al = bless {
		_rdname => $rdname,    # read name
		_rdseq  => $rdseq,     # read bases
		_rdqual => $rdqual,    # qualities
		_color  => $color,     # colorspace?
		_tname  => $tname,     # reference sequence name
		_toff   => $toff,      # reference sequence offset
		_fw     => $fw,        # true -> forward read aligned strand
		_score  => $score,     # alignment score
		_mapq   => $mapq,      # mapping quality
		_md     => $md,        # MD:Z string
		_cigar  => $cigar,     # CIGAR string
		_mate1  => $mate1,     # Mate #1
		_mate2  => $mate2,     # Mate #2
		_alread => $alread,    # read side of the stacked alignment
		_alref  => $alref,     # reference side of the stacked alignment
		_rdex   => $rdex,      # # read chars involved in alignment
		_rfex   => $rfex,      # # ref chars involved in alignment
		_rfknown=> $ref_known, # MD:Z string known?
		_htriml => $htriml,    # hard trimming on LHS
		_striml => $striml,    # soft trimming on LHS
		_htrimr => $htrimr,    # hard trimming on RHS
		_strimr => $strimr,    # soft trimming on RHS
		_flags  => $flags,     # sam flags (or undef)
		_opts   => $opts,      # optional flags
		_orig   => $orig
	}, $class;
	if(defined($tname) && $cigar ne "*") {
		($alread, $alref, $rdex, $rfex, $htriml, $striml, $htrimr, $strimr) =
			$al->md_cigar_to_al();
		$mapq  == int($mapq)  || die "MAPQ must be int, was $mapq";
		$score == int($score) || die "Score must be int, was $score";
		$toff  == int($toff)  || die "POS must be int, was $toff";
	}
	$al->{_alread} = $alread;
	$al->{_alref}  = $alref;
	$al->{_rdex}   = $rdex;
	$al->{_rfex}   = $rfex;
	
	$al->{_htriml} = $htriml;
	$al->{_striml} = $striml;
	$al->{_htrimr} = $htrimr;
	$al->{_strimr} = $strimr;
	
	return $al;
}

sub read_name   { return $_[0]->{_rdname} }
sub read_seq    { return $_[0]->{_rdseq}  }
sub read_qual   { return $_[0]->{_rdqual} }
sub read_extent { return $_[0]->{_rdex}   }
sub color       { return $_[0]->{_color}  }
sub ref_extent  { return $_[0]->{_rfex}   }
sub ref_known   { return $_[0]->{_rfknown}}
sub tname       { return $_[0]->{_tname}  }
sub toff        { return $_[0]->{_toff}   }
sub fw          { return $_[0]->{_fw}     }
sub score       { return $_[0]->{_score}  }
sub mapq        { return $_[0]->{_mapq}   }
sub md          { return $_[0]->{_md}     }
sub cigar       { return $_[0]->{_cigar}  }
sub mate1       { return $_[0]->{_mate1}  }
sub mate2       { return $_[0]->{_mate2}  }
sub flags       { return $_[0]->{_flags}  }
sub aligned     { return defined($_[0]->{_cigar}) }
sub htriml      { return $_[0]->{_htriml} }
sub striml      { return $_[0]->{_striml} }
sub htrimr      { return $_[0]->{_htrimr} }
sub strimr      { return $_[0]->{_strimr} }
sub triml       { return $_[0]->{_htriml} + $_[0]->{_striml} }
sub trimr       { return $_[0]->{_htrimr} + $_[0]->{_strimr} }

sub mate_aligned { return (($_[0]->{_flags} & 8) == 0) }
sub is_paired    { return (($_[0]->{_flags} & 1) != 0) }

sub mateid {
	return 1 if $_[0]->{_mate1};
	return 2 if $_[0]->{_mate2};
	return 0;
}

sub alread    { return $_[0]->{_alread} }
sub alref     { return $_[0]->{_alref}  }

sub empty_alignment() {
	return BtlBio::Alignment::Alignment->new(
		undef,  # read name
		undef,  # read sequence from alignment record
		undef,  # read qualities from alignment record
		undef,  # whether read is colorspace
		undef,  # name of reference sequence
		undef,  # offset into reference
		undef,  # true -> aligned to Watson strand
		undef,  # alignment score
		undef,  # mapping quality
		undef,  # MD:Z string
		undef,  # CIGAR string
		0,      # is mate #1 in a pair?
		0,      # is mate #2 in a pair?
		undef,  # SAM flags (if available)
		{},     # options hash ref
		undef); # original string
}

##
# Convert an options hash to a string.
#
sub aln_opt_hash_to_string($) {
	my $h = shift;
	my @toks = ();
	for my $k (sort keys %$h) {
		push @toks, join(":", ($k, $h->{$k}{type}, $h->{$k}{value}));
	}
	return join("\t", @toks);
}

##
# Convert a string-ized options hash back into a hash.
#
sub aln_string_to_opt_hash($) {
	my $str = shift;
	chomp($str);
	my @ts = split(/\t/, $str);
	my %h = ();
	for my $t (@ts) {
		my @cs = split(/:/, $t);
		scalar(@cs) >= 3 || croak("Malformed optional flag: '$t'");
		my ($key, $type, $value) = ($cs[0], $cs[1], join(":", @cs[2..$#cs]));
		defined($type) || die "Bad format for optional field: '$str'";
		defined($value) || die "Bad format for optional field: '$str'";
		$h{$key}{type} = $type;
		$h{$key}{value} = $value;
	}
	return \%h;
}

##
# Return the value of the given optional flag.  Optionally check that the
# flag's type is as expected.
#
sub opt {
	my ($self, $key, $typ) = @_;
	defined($self->{_opts}{$key}) ||
		croak("No such option: $key");
	!defined($typ) || $typ eq $self->{_opts}{$key}{type} ||
		croak("Expected type '$typ', but got '".$self->{_opts}{$key}{type}."'");
	defined($self->{_opts}{$key}{value}) || die;
	return $self->{_opts}{$key}{value};
}

##
# Remove trailing /1 or /2 if it's there.
#
sub fix_mate_name($) {
	$_[0]->{_rdname} =~ s/\/[12]$//;
}

##
# Given a read sequence (with characters in upstream-to-downstream order with
# respect to the reference - NOT necessarily 5'-to-3') and a CIGAR string and
# an MD:Z string, build the alignment strings.  The alignment strings will only
# contain the portion of the read that aligned.  Any portions that were either
# hard-trimmed or soft-trimmed are trimmed from this function's result.
#
# For now, I'm assuming that the MD:Z string only describes aligned characters,
# i.e. *after* trimming.
#
sub md_cigar_to_al($$) {
	my ($self) = @_;
	my ($seq, $md, $cigar) = ($self->{_rdseq}, $self->{_md}, $self->{_cigar});
	my $alread = "";
	my $alref = "";
	defined($seq) || die;
	defined($cigar) || die;
	$cigar ne "*" ||
		croak("CIGAR string should not be '*' - only ".
		      "aligned reads should call md_cigar_to_al");
	my ($parsed_cig, $nms, $nds, $nis) = parse_cigar($cigar);
	my $parsed_md;
	if(defined($md)) {
		$parsed_md = parse_md($md);
	} else {
		$parsed_md = " "x($nms+$nds);
	}
	my ($rdoff, $mdoff) = (0, 0);
	my ($htriml, $striml, $htrimr, $strimr) = (0, 0, 0, 0);
	my $nonsh = 0; # have I seen a non-S, non-H CIGAR op?
	my $nonh = 0;  # have I seen a non-H CIGAR op?
	my ($rdex, $rfex) = (0, 0);
	for(my $i = 0; $i < length($parsed_cig); $i++) {
		my $cigop = substr($parsed_cig, $i, 1);
		$rdoff < length($seq) || die "Bad rdoff:$rdoff for seq '$seq' cigop=$cigop";
		$nonh++ if $cigop ne "H";
		$nonsh++ if ($cigop ne "H" && $cigop ne "S");
		if($cigop eq "S") {
			if($nonsh) {
				$strimr++;
			} else {
				$striml++;
			}
			$rdoff++;
			next;
		}
		if($cigop eq "H") {
			if($nonh) {
				$htrimr++;
			} else {
				$htriml++;
			}
			next;
		}
		$cigop = "M" if $cigop eq "=" || $cigop eq "X";
		if($cigop eq "P" || $cigop eq "N") {
			# Padding
			$alread .= "-";
			$alref .= "-";
		} elsif($cigop eq "M") {
			my $rdc = substr($seq, $rdoff, 1);
			my $rfc = substr($parsed_md, $mdoff, 1);
			$rfc = $rdc if $rfc eq " ";
			$alread .= $rdc;
			$alref .= $rfc;
			$rdoff++;
			$mdoff++;
			$rdex++;
			$rfex++;
		} elsif($cigop eq "D") {
			# Read gap
			#  Read: AAA-AAA
			#   Ref: AAAAAAA
			my $rfc = substr($parsed_md, $mdoff, 1);
			if($self->{_rfknown}) {
				$rfc ne " " ||
					die "Must have a ref character opposite a gap in the read:\n".
					    "cig: $parsed_cig ($i)\nmd:  $parsed_md ($mdoff)\n";
			} else {
				$rfc = "?";
			}
			$alref .= $rfc;
			$alread .= "-";
			$mdoff++;
			$rfex++;
		} else {
			# Reference gap
			#  Read: AAAAAAA
			#   Ref: AAA-AAA
			$cigop eq "I" || die "Unsupported cigop: $cigop";
			my $rdc = substr($seq, $rdoff, 1);
			$alread .= $rdc;
			$alref .= "-";
			$rdoff++;
			$rdex++;
			# $mdoff SHOULD NOT be incremented here.  A simple test using BWA
			# confirms this
		}
	}
	return (
		$alread,
		$alref,
		$rdex,
		$rfex,
		$htriml,
		$striml,
		$htrimr,
		$strimr);
}

sub to_record {
	return join("\t", (
		$_[0]->{_rdname},
		$_[0]->{_rdseq},
		$_[0]->{_rdqual},
		$_[0]->{_color},
		$_[0]->{_tname},
		$_[0]->{_toff},
		$_[0]->{_fw},
		$_[0]->{_score},
		$_[0]->{_mapq},
		$_[0]->{_md},
		$_[0]->{_cigar},
		$_[0]->{_mate1},
		$_[0]->{_mate2},
		$_[0]->{_flags}))."\n";
}

sub parse_generic_aln_record($) {
	my $line = shift;
	chomp($line);
	my @ts = split(/\t/, $line, -1);
	return BtlBio::Alignment::Alignment->new(
		$ts[0],  # read name
		$ts[1],  # read sequence from alignment record
		$ts[2],  # read qualities from alignment record
		$ts[3],  # colorspace?
		$ts[4],  # name of ref sequence
		$ts[5],  # offset into ref
		$ts[6],  # true -> aligned to Watson
		$ts[7],  # alignment score
		$ts[8],  # mapping quality
		$ts[9],  # MD:Z string
		$ts[10], # CIGAR string
		$ts[11], # mate1?
		$ts[12], # mate2?
		$ts[13], # SAM flags (if available)
		{ },     # options
		undef);  # orig
}

##
# Return non-zero iff this Alignment and the given Alignment share at least one
# cell in common in the overall dynamic programming table.
#
sub overlaps($$) {
	my ($self, $al) = @_;
	return 0 if !$self->aligned || !$al->aligned;
	defined($self->{_alread}) || die;
	defined($al->{_alread}) || die;
	return 0 if $self->fw != $al->fw;
	return 0 if $self->tname ne $al->tname;
	# Now put all the cells traversed by this alignment into a hash
	my ($row, $col) = ($self->striml, $self->toff);
	my %cells = ();
	for(my $i = 0; $i < length($self->{_alread}); $i++) {
		$cells{$row}{$col} = 1;
		my $rdc = substr($self->{_alread}, $i, 1);
		my $rfc = substr($self->{_alref},  $i, 1);
		$row++ if $rdc ne "-";
		$col++ if $rfc ne "-";
	}
	# Now iterate through cells in the other alignment
	($row, $col) = ($al->striml, $al->toff);
	for(my $i = 0; $i < length($al->{_alread}); $i++) {
		return 1 if $cells{$row}{$col};
		my $rdc = substr($al->{_alread}, $i, 1);
		my $rfc = substr($al->{_alref},  $i, 1);
		$row++ if $rdc ne "-";
		$col++ if $rfc ne "-";
	}
	return 0;
}

sub _test_read_md_cigar_to_al($$$) {
	my ($seq, $md, $cig) = @_;
	my $al = BtlBio::Alignment::Alignment->new(
		"tmp",
		$seq,
		"I" x length($seq),
		0,     # colorspace?
		"ref",
		0,
		1,
		10,
		10,
		$md,
		$cig,
		0, 0,    # mate 1/2
		undef,   # SAM flags (if available)
		{ },     # options
		"");
	return $al->md_cigar_to_al();
}

sub _test_read_md_cigar_to_al_1() {
	print STDERR "Testing _read_md_cigar_to_al 1 ... ";
	my ($alread, $alref, $rdex, $rfex, $htriml, $striml, $htrimr, $strimr) =
		_test_read_md_cigar_to_al(
			"TCGATCCATCTAGCTGAGTCGCGACCATATATCGACTACGATCGATCGCTACG",
			"25^T28", "25M1D28M");
	#                     1         2               1         2       
	#           0123456789012345678901234 0123456789012345678901234567
	$alread eq "TCGATCCATCTAGCTGAGTCGCGAC-CATATATCGACTACGATCGATCGCTACG" || die;
	$alref  eq "TCGATCCATCTAGCTGAGTCGCGACTCATATATCGACTACGATCGATCGCTACG" || die;
	$rdex == (25+28) || die;
	$rfex == (25+28+1) || die;
	$htriml == 0 || die;
	$striml == 0 || die;
	$htrimr == 0 || die;
	$strimr == 0 || die;
	print STDERR "PASSED\n";
}

sub _test_read_md_cigar_to_al_2() {
	print STDERR "Testing _read_md_cigar_to_al 2 ... ";
	my ($alread, $alref, $rdex, $rfex, $htriml, $striml, $htrimr, $strimr) =
		_test_read_md_cigar_to_al(
			"ATCCATCTAGCTGAGTCACGTCTCATATATCGACTACGATCGATC",
			"17G2A24", "45M");
	$alread eq "ATCCATCTAGCTGAGTCACGTCTCATATATCGACTACGATCGATC" || die;
	$alref  eq "ATCCATCTAGCTGAGTCGCGACTCATATATCGACTACGATCGATC" || die;
	$rdex == 45 || die;
	$rfex == 45 || die;
	$htriml == 0 || die;
	$striml == 0 || die;
	$htrimr == 0 || die;
	$strimr == 0 || die;
	print STDERR "PASSED\n";
}

sub _test_read_md_cigar_to_al_3() {
	print STDERR "Testing _read_md_cigar_to_al 3 ... ";
	my ($alread, $alref, $rdex, $rfex, $htriml, $striml, $htrimr, $strimr) =
		_test_read_md_cigar_to_al(
			"ATCCATCTAGCTGAGTCGCGACTTTATATCGACTACGATCGATCGCTA",
			"23^C0A24", "23M1D25M");
	my $ex_alread = "ATCCATCTAGCTGAGTCGCGACT-TTATATCGACTACGATCGATCGCTA";
	my $ex_alref  = "ATCCATCTAGCTGAGTCGCGACTCATATATCGACTACGATCGATCGCTA";
	#        CATGATCGATCCATCTAGCTGAGTCGCGACTCATATATCGACTACGATCGATCGCTACGTACGT
	$alread eq $ex_alread || die "Expected/got:\n$ex_alread\n$alread";
	$alref  eq $ex_alref  || die "Expected/got:\n$ex_alref\n$alref";
	$rdex == 48 || die;
	$rfex == 49 || die;
	$htriml == 0 || die;
	$striml == 0 || die;
	$htrimr == 0 || die;
	$strimr == 0 || die;
	print STDERR "PASSED\n";
}

sub _test_read_md_cigar_to_al_4() {
	print STDERR "Testing _read_md_cigar_to_al 4 ... ";
	my ($alread, $alref, $rdex, $rfex, $htriml, $striml, $htrimr, $strimr) =
		_test_read_md_cigar_to_al(
			"ATCCATCTAGCTGAGTCGCGACTTTATATCGACTACGATCGATCGCTA",
			"21^C0A24", "2S21M1D25M");
	#                012345678901234567890  012345678901234567890123
	my $ex_alread = "CCATCTAGCTGAGTCGCGACT-TTATATCGACTACGATCGATCGCTA";
	my $ex_alref  = "CCATCTAGCTGAGTCGCGACTCATATATCGACTACGATCGATCGCTA";
	#        CATGATCGATCCATCTAGCTGAGTCGCGACTCATATATCGACTACGATCGATCGCTACGTACGT
	$alread eq $ex_alread || die "Expected/got:\n$ex_alread\n$alread";
	$alref  eq $ex_alref  || die "Expected/got:\n$ex_alref\n$alref";
	$rdex == 46 || die;
	$rfex == 47 || die;
	$htriml == 0 || die;
	$striml == 2 || die;
	$htrimr == 0 || die;
	$strimr == 0 || die;
	print STDERR "PASSED\n";
}

sub _test_read_md_cigar_to_al_5() {
	print STDERR "Testing _read_md_cigar_to_al 5 ... ";
	my ($alread, $alref, $rdex, $rfex, $htriml, $striml, $htrimr, $strimr) =
		_test_read_md_cigar_to_al(
			"AACCATCTAGCTGAGTCGCGACTCTTATATCGACTACGATCGATCGCTAAAAA",
	#        ^^ trimmed                                       ^^^^ trimmed
			"47", "2S21M1I25M4S");
	#                0123456789012345678901 012345678901234567890123
	my $ex_alread = "CCATCTAGCTGAGTCGCGACTCTTATATCGACTACGATCGATCGCTA";
	my $ex_alref  = "CCATCTAGCTGAGTCGCGACT-TTATATCGACTACGATCGATCGCTA";
	#        CATGATCGATCCATCTAGCTGAGTCGCGACTCATATATCGACTACGATCGATCGCTACGTACGT
	$alread eq $ex_alread || die "Expected/got:\n$ex_alread\n$alread";
	$alref  eq $ex_alref  || die "Expected/got:\n$ex_alref\n$alref";
	$rdex == 47 || die;
	$rfex == 46 || die;
	$htriml == 0 || die;
	$striml == 2 || die;
	$htrimr == 0 || die;
	$strimr == 4 || die;
	print STDERR "PASSED\n";
}

sub _test_read_md_cigar_to_al_6() {
	print STDERR "Testing _read_md_cigar_to_al 6 ... ";
	my ($alread, $alref, $rdex, $rfex, $htriml, $striml, $htrimr, $strimr) =
		_test_read_md_cigar_to_al(
			"AACCATCTAGCTGAGTCGCGACTCTTATATCGACTACGATCGATCGCTAAAAA",
	#        ^^ trimmed                                       ^^^^ trimmed
			"21A24", "2S21M1I25M4S");
	#                0123456789012345678901 012345678901234567890123
	my $ex_alread = "CCATCTAGCTGAGTCGCGACTCTTATATCGACTACGATCGATCGCTA";
	my $ex_alref  = "CCATCTAGCTGAGTCGCGACT-ATATATCGACTACGATCGATCGCTA";
	#        CATGATCGATCCATCTAGCTGAGTCGCGACTCATATATCGACTACGATCGATCGCTACGTACGT
	$alread eq $ex_alread || die "Expected/got:\n$ex_alread\n$alread";
	$alref  eq $ex_alref  || die "Expected/got:\n$ex_alref\n$alref";
	$rdex == 47 || die;
	$rfex == 46 || die;
	$htriml == 0 || die;
	$striml == 2 || die;
	$htrimr == 0 || die;
	$strimr == 4 || die;
	print STDERR "PASSED\n";
}

sub _test_read_md_cigar_to_al_7() {
	print STDERR "Testing _read_md_cigar_to_al 7 ... ";
	my ($alread, $alref, $rdex, $rfex, $htriml, $striml, $htrimr, $strimr) =
		_test_read_md_cigar_to_al(
			"ATCCATCTAGCTGAGTCGCGACTTTATATCGACTACGATCGATCGCTA",
			"21^C0A19", "2S21M1D20M5S");
	my $ex_alread = "CCATCTAGCTGAGTCGCGACT-TTATATCGACTACGATCGAT";
	my $ex_alref  = "CCATCTAGCTGAGTCGCGACTCATATATCGACTACGATCGAT";
	#        CATGATCGATCCATCTAGCTGAGTCGCGACTCATATATCGACTACGATCGATCGCTACGTACGT
	$alread eq $ex_alread || die "Expected/got:\n$ex_alread\n$alread";
	$alref  eq $ex_alref  || die "Expected/got:\n$ex_alref\n$alref";
	$htriml == 0 || die;
	$striml == 2 || die;
	$htrimr == 0 || die;
	$strimr == 5 || die;
	print STDERR "PASSED\n";
}

sub _test_read_md_cigar_to_al_8() {
	print STDERR "Testing _read_md_cigar_to_al 8 ... ";
	my ($alread, $alref, $rdex, $rfex, $htriml, $striml, $htrimr, $strimr) =
		_test_read_md_cigar_to_al("AAAACCAAAA", "5^A4", "3M1I2M1D4M");
	my $ex_alread = "AAAACC-AAAA";
	my $ex_alref  = "AAA-CCAAAAA";
	$alread eq $ex_alread || die "Expected/got:\n$ex_alread\n$alread";
	$alref  eq $ex_alref  || die "Expected/got:\n$ex_alref\n$alref";
	print STDERR "PASSED\n";
}

sub _test_parse_generic_aln_record_1() {
	print STDERR "Testing parse_generic_aln_record 1 ... ";
	my $line = "rd1	ACGTACGTACG	IIIIIIIIIII	0	rf1	20	1	-20	40	5^T6	5M1D6M	0	0	0\n";
	my $al = parse_generic_aln_record($line);
	$al->read_name eq "rd1" || die;
	my $rec = $al->to_record;
	$rec eq $line || die "Should have matched:\n$rec\n$line\n";
	print STDERR "PASSED\n";
}

sub _test_overlaps_1() {
	print STDERR "Testing overlaps 1 ... ";
	# length=11, toff=20, fw=true, score=-20, mapq=40
	my $al1 = parse_generic_aln_record(
		"read1	ACGTACGTACG	IIIIIIIIIII	1	ref1	20	1	-20	40	5T5	11M	0	0	0\n");
	my $al2 = parse_generic_aln_record(
		"read2	ACGTACGTACG	IIIIIIIIIII	1	ref1	20	0	-20	40	5T5	11M	0	0	0\n");
	!$al1->overlaps($al2) || die "Should not have overlapped:\n$al1\n$al2";
	print STDERR "PASSED\n";
}

sub _test_overlaps_2() {
	print STDERR "Testing overlaps 2 ... ";
	# length=11, toff=20, fw=true, score=-20, mapq=40
	my $al1 = parse_generic_aln_record(
		"read1	ACGTACGTACG	IIIIIIIIIII	1	ref1	20	0	-20	40	5T5	11M	0	0	0\n");
	my $al2 = parse_generic_aln_record(
		"read2	ACGTACGTACG	IIIIIIIIIII	1	ref1	21	0	-20	40	5T5	11M	0	0	0\n");
	!$al1->overlaps($al2) || die "Should not have overlapped:\n$al1\n$al2";
	print STDERR "PASSED\n";
}

sub _test_overlaps_3() {
	print STDERR "Testing overlaps 3 ... ";
	# length=11, toff=20, fw=true, score=-20, mapq=40
	my $al1 = parse_generic_aln_record(
		"read1	ACGTACGTACG	IIIIIIIIIII	0	ref1	21	0	-20	40	5T5	11M	0	0	0\n");
	my $al2 = parse_generic_aln_record(
		"read2	ACGTACGTACG	IIIIIIIIIII	0	ref2	21	0	-20	40	5T5	11M	0	0	0\n");
	!$al1->overlaps($al2) || die "Should not have overlapped:\n$al1\n$al2";
	print STDERR "PASSED\n";
}

sub _test_overlaps_4() {
	print STDERR "Testing overlaps 4 ... ";
	# length=11, toff=20, fw=true, score=-20, mapq=40
	my $al1 = parse_generic_aln_record(
		"read1	ACGTACGTACG	IIIIIIIIIII	0	ref1	21	0	-20	40	5T5	11M	0	0	0\n");
	my $al2 = parse_generic_aln_record(
		"read2	ACGTAGGTACG	IIIIIIIIIII	0	ref1	21	0	-20	40	5T5	11M	0	0	0\n");
	$al1->overlaps($al2) || die "Should have overlapped:\n$al1\n$al2";
	print STDERR "PASSED\n";
}

sub _test_overlaps_5() {
	print STDERR "Testing overlaps 5 ... ";
	# length=11, toff=20, fw=true, score=-20, mapq=40
	my $al1 = parse_generic_aln_record(
		"read1	ACGTACGTACG	IIIIIIIIIII	0	ref1	20	0	-20	40	1^T10	1M1D10M	0	0	0\n");
	my $al2 = parse_generic_aln_record(
		"read2	ACGTAGGTACG	IIIIIIIIIII	0	ref1	21	0	-20	40	5T5	11M	0	0	0\n");
	$al1->overlaps($al2) || die "Should have overlapped:\n$al1\n$al2";
	print STDERR "PASSED\n";
}

sub _test_overlaps_6() {
	print STDERR "Testing overlaps 6 ... ";
	# length=11, toff=20, fw=true, score=-20, mapq=40
	my $al1 = parse_generic_aln_record(
		"read1	ACGTACGTACG	IIIIIIIIIII	1	ref1	20	0	-20	40	10	1M1I9M	0	0	0\n");
	my $al2 = parse_generic_aln_record(
		"read2	ACGTAGGTACG	IIIIIIIIIII	1	ref1	21	0	-20	40	5T5	11M	0	0	0\n");
	!$al1->overlaps($al2) || die "Should not have overlapped:\n$al1\n$al2";
	print STDERR "PASSED\n";
}

sub _test_overlaps_7() {
	print STDERR "Testing overlaps 7 ... ";
	# length=11, toff=20, fw=true, score=-20, mapq=40
	my $al1 = parse_generic_aln_record(
		"read1	ACGTACGTACG	IIIIIIIIIII	1	ref1	21	0	-20	40	10	1M1I9M	0	0	0\n");
	my $al2 = parse_generic_aln_record(
		"read2	ACGTAGGTACG	IIIIIIIIIII	1	ref1	20	0	-20	40	5T5	11M	0	0	0\n");
	$al1->overlaps($al2) || die "Should have overlapped:\n$al1\n$al2";
	print STDERR "PASSED\n";
}

sub _test_overlaps_8() {
	print STDERR "Testing overlaps 8 ... ";
	# Test that soft clipping is taken into account when determining overlaps
	my $al1 = parse_generic_aln_record(
		"read1	ACGTACGTACG	IIIIIIIIIII	0	ref1	25	0	-20	40	6	5S6M	0	0	0\n");
	my $al2 = parse_generic_aln_record(
		"read2	ACGTAGGTACG	IIIIIIIIIII	0	ref1	20	0	-20	40	5	5M6S	0	0	0\n");
	!$al1->overlaps($al2) || die "Should not have overlapped:\n$al1\n$al2";
	print STDERR "PASSED\n";
}

sub _test_overlaps_9() {
	print STDERR "Testing overlaps 9 ... ";
	# Test that soft clipping is taken into account when determining overlaps
	my $al1 = parse_generic_aln_record(
		"read1	ACGTACGTACG	IIIIIIIIIII	1	ref1	22	0	-20	40	7	2S7M	0	0	0\n");
	my $al2 = parse_generic_aln_record(
		"read2	ACGTAGGTACG	IIIIIIIIIII	1	ref1	20	0	-20	40	5	5M6S	0	0	0\n");
	$al1->overlaps($al2) || die "Should have overlapped:\n$al1\n$al2";
	print STDERR "PASSED\n";
}

sub _test_overlaps_10() {
	print STDERR "Testing overlaps 10 ... ";
	# Test that soft clipping is taken into account when determining overlaps
	my $al1 = parse_generic_aln_record(
		"read1	ACGTACGTACG	IIIIIIIIIII	0	ref1	24	0	-20	40	7	4S7M	0	0	0\n");
	my $al2 = parse_generic_aln_record(
		"read2	ACGTAGGTACG	IIIIIIIIIII	0	ref1	20	0	-20	40	5	5M6S	0	0	0\n");
	$al1->overlaps($al2) || die "Should have overlapped:\n$al1\n$al2";
	print STDERR "PASSED\n";
}

sub _test_overlaps_11() {
	print STDERR "Testing overlaps 11 ... ";
	# Test that soft clipping is taken into account when determining overlaps
	my $al1 = parse_generic_aln_record(
		"read1	ACGTACGTACG	IIIIIIIIIII	1	ref1	21	0	-20	40	7	4S7M	0	0	0\n");
	my $al2 = parse_generic_aln_record(
		"read2	ACGTAGGTACG	IIIIIIIIIII	1	ref1	20	0	-20	40	8	3S8M	0	0	0\n");
	$al1->overlaps($al2) || die "Should have overlapped:\n$al1\n$al2";
	print STDERR "PASSED\n";
}

sub _test_opts_1() {
	print STDERR "Testing opts 1 ... ";
	my %h = ();
	$h{"MD"}{type}  = "Z";
	$h{"MD"}{value} = 16;
	$h{"CS"}{type}  = "Z";
	$h{"CS"}{value} = "A022201321013022";
	$h{"NM"}{type}  = "i";
	$h{"NM"}{value} = "40";
	my $str = aln_opt_hash_to_string(\%h);
	my $ex_str = "CS:Z:A022201321013022\tMD:Z:16\tNM:i:40";
	$str eq $ex_str || die "Expected:\n$ex_str\ngot:\n$str";
	my $h2 = aln_string_to_opt_hash($str);
	$h2->{"MD"}{type}  eq "Z"                || die;
	$h2->{"MD"}{value} eq "16"               || die;
	$h2->{"CS"}{type}  eq "Z"                || die;
	$h2->{"CS"}{value} eq "A022201321013022" || die;
	$h2->{"NM"}{type}  eq "i"                || die;
	$h2->{"NM"}{value} eq "40"               || die;
	print STDERR "PASSED\n";
}

sub _test_opts_2() {
	print STDERR "Testing opts 2 ... ";
	my %h = ();
	$h{"MD"}{type}  = "Z";
	$h{"MD"}{value} = "4";
	$h{"CS"}{type}  = "Z";
	$h{"CS"}{value} = "";
	$h{"NM"}{type}  = "i";
	$h{"NM"}{value} = "-10";
	my $al = BtlBio::Alignment::Alignment->new(
		"tmp",
		"ACGT",
		"IIII",
		1,      # colorspace?
		"ref",
		0,
		1,
		10,
		10,
		"4",
		"4M",
		0, 0,    # mate 1/2
		undef,   # SAM flags (if available)
		\%h,
		"");
	$al->opt("MD", "Z") eq "4" || die;
	$al->opt("CS", "Z") eq ""  || die;
	$al->opt("NM", "i") == -10 || die;
	print STDERR "PASSED\n";
}

sub _test_empty_1() {
	print STDERR "Testing empty_alignment 1 ... ";
	my $al = empty_alignment();
	!$al->aligned || die;
	print STDERR "PASSED\n";
}

sub _test() {
	_test_read_md_cigar_to_al_1();
	_test_read_md_cigar_to_al_2();
	_test_read_md_cigar_to_al_3();
	_test_read_md_cigar_to_al_4();
	_test_read_md_cigar_to_al_5();
	_test_read_md_cigar_to_al_6();
	_test_read_md_cigar_to_al_7();
	_test_read_md_cigar_to_al_8();
	_test_parse_generic_aln_record_1();
	_test_overlaps_1();
	_test_overlaps_2();
	_test_overlaps_3();
	_test_overlaps_4();
	_test_overlaps_5();
	_test_overlaps_6();
	_test_overlaps_7();
	_test_overlaps_8();
	_test_overlaps_9();
	_test_overlaps_10();
	_test_overlaps_11();
	_test_opts_1();
	_test_opts_2();
	_test_empty_1();
}

# Call _test() if the module is being run directly.
_test() unless caller();

1;
