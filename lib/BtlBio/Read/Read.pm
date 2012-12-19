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
# BtlBio::Read::Read
#
# Encapsulates a read, regardless of input format.  A read has the following
# properties:
#
# 1. A name
# 2. A mate id (1 or 2)
# 3. A sequence
# 4. A quality string
#
# Author: Ben Langmead
# Email: langmea@cs.jhu.edu
#

package BtlBio::Read::Read;

@ISA = qw(Exporter);
@EXPORT = qw(reads_to_tab5 reads_to_tab6);

use strict;
use warnings;
use Carp;

##
# Create a new read object
#
sub new {
	my (
		$class,
		$name,   # name of read
		$seq,    # sequence string
		$qual,   # quality string
		$mateid, # mate id (default: 0, indicating unpaired)
		$color,  # colorspace?
		$orig    # original record
	) = @_;
	$mateid = $mateid || 0;
	$orig = $orig || "(unspecified)";
	my $name_canonical = $name;
	$name_canonical =~ s/\/[12]$//;
	my $primer = '?';
	if(defined($seq) && $color) {
		# See if we can parse out the primer base
		my $ci = substr($seq, 0, 1);
		my $cf = substr($seq, -1, 1);
		if($ci =~ /[ACGT]/ && $cf =~ /[^ACGTN]/) {
			$primer = $ci;
			$seq = substr($seq, 1);
		}
		# Encode numeric colors as nucleotides
		$seq =~ tr/01234\./ACGTNN/;
	}
	$seq = uc $seq;
	length($seq) == 0 || $seq =~ /[ACGTN.]/ || croak("Bad characters in read sequence: '$seq'");
	return bless {
		_name           => $name,
		_name_canonical => $name_canonical,
		_seq            => $seq,
		_qual           => $qual,
		_mateid         => $mateid,
		_color          => $color,
		_primer         => $primer,
		_orig           => $orig
	}, $class;
}
sub name       { return $_[0]->{_name}   }
sub seq        { return $_[0]->{_seq}    }
sub qual       { return $_[0]->{_qual}   }
sub mateid     { return $_[0]->{_mateid} }
sub colorspace { return $_[0]->{_color}  }
sub primer     { return $_[0]->{_primer} }
sub orig       { return $_[0]->{_orig}   }

sub to_fastq {
	return "\@".$_[0]->{_name}."\n".
	       $_[0]->{_seq}."\n".
	       "+\n".
	       $_[0]->{_qual}."\n";
}

##
# Return a read record in tabbed format.
#
sub to_tab5 {
	return join("\t", ($_[0]->{_name}, $_[0]->{_seq}, $_[0]->{_qual}))."\n";
}

##
# Return a read record in tabbed format.
#
sub to_tab6 {
	return join("\t", ($_[0]->{_name}, $_[0]->{_seq}, $_[0]->{_qual}))."\n";
}

##
# Given a read and possibly a second read (in which case they're interpreted
# as a pair), return a read record in tab5 format.
#
sub reads_to_tab5 {
	my ($m1, $m2) = @_;
	defined($m1) || die;
	if(defined($m2)) {
		return join("\t", (
			$m1->{_name},
			$m1->{_seq},
			$m1->{_qual},
			$m2->{_seq},
			$m2->{_qual}))."\n";
	} else {
		return join("\t", (
			$m1->{_name},
			$m1->{_seq},
			$m1->{_qual}))."\n";
	}
}

##
# Given a read and possibly a second read (in which case they're interpreted
# as a pair), return a read record in tab6 format.
#
sub reads_to_tab6 {
	my ($m1, $m2) = @_;
	defined($m1) || die;
	if(defined($m2)) {
		return join("\t", (
			$m1->{_name},
			$m1->{_seq},
			$m1->{_qual},
			$m2->{_name},
			$m2->{_seq},
			$m2->{_qual}))."\n";
	} else {
		return join("\t", (
			$m1->{_name},
			$m1->{_seq},
			$m1->{_qual}))."\n";
	}
}

sub name_canonical { return $_[0]->{_name_canonical} }

##
# Turn a quality string made up of whitespace-separated numbers into an
# ASCII-coded sequence with ASCII offset=33.
#
sub asciiize_quals($) {
	my $self = shift;
	my $quals = $self->{_qual};
	$quals =~ s/^\s*//;
	$quals =~ s/\s*$//;
	my @qual_toks = split(/\s+/, $quals);
	$self->{_qual} = "";
	for my $q (@qual_toks) {
		$q == int($q) ||
			croak("Expected quality value to be an integer, got '$q'");
		$q = 0 if $q < 0;
		$q = 127 if $q > 127;
		$self->{_qual} .= chr($q+33);
	}
}

1;
