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
# BtlBio::Align::Bisulfite::BsAlign
#
# Routines for in-silico bisulfite conversion of sequence data and analysis of
# alignment results.
#
# Information about the original read sequence has to be "threaded through" the
# aligner somehow - otherwise, we can't discern any pro-methylation evidence
# from the alignment.  One idea is to have wrapper threads on either either
# side of the aligner; one that writes converted reads to the aligner, and one
# that reads alignments from it.  However, unless the read itself is modified
# to include information about the original sequence, our threaded scheme must
# have the property that the thread reading alignments must know or be able to
# look up each read's original sequence.  THAT is hard to imagine how to do
# unless an appropriate version of perl threading (e.g. ithreads) is available.
# That's not a very portable way to do it.
#
# Instead, what we do is simply append the original sequence to the read name.
# This pollutes the intermediate results with extra info, but doesn't require
# pre- and post-alignment software to share any information about the reads.
#
# Author: Ben Langmead
# Email: langmea@cs.jhu.edu
#

package BtlBio::Align::Bisulfite::BsAlign;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(bs_analyze_alignment
             bsc_ize_read
             bsc_analyze_alignment
             bscpg_ize_read
             bscpg_analyze_alignment);

use strict;
use warnings;
use Carp;

use FindBin qw($Bin);
use lib "$Bin/../../..";
use BtlBio::Util::File;
use BtlBio::Alphabet::DNA;
use BtlBio::Alignment::Alignment;
use BtlBio::Alignment::Colorspace;
use BtlBio::Align::Bisulfite::BsEvidence;
use List::Util qw[min max];

##
# Given a read, convert all Cs to Ts and remember which Ts are converted versus
# original by storing the original sequence in the {_bsc_orig_seq} field.  If
# the $rc parameter is defined and non-zero, reverse complement the read first.
#
sub bsc_ize_read {
	my ($rd, $rc, $mate) = @_;
	$mate = $mate || 0;
	$rd->{_seq} = normalize_nucs($rd->{_seq});  # turn non-ACGTs to Ns
	$rd->{_seq} = revcomp($rd->{_seq}) if $rc;
	$rd->{_qual} = reverse $rd->{_qual} if $rc;
	# Save the original sequence
	$rd->{_bsc_orig_seq} = $rd->{_seq};
	# In-silico convert the read
	$rd->{_seq} =~ tr/Cc/Tt/;
	# Append the original sequence on to the end of the name
	$rd->{_name} =~ s/\s/_/g;
	$rd->{_name} .= ":".$rd->{_bsc_orig_seq}.$mate.($rc ? "-" : "+");
	return $rd;
}

##
# Analyze an alignment output by a tool (e.g. Merman) that takes raw reference
# & read and, bisulfite aspects internally, and reports alignments with respect
# to the Watson strand.  Such a tool must also convey what strand was aligned
# to (BSW/BSC/BSWR/BSCR); we assume it does so via the XB:Z optional flag.
# Also, if the read is colorspace, the tool must convey the original quality
# values; we assume it does so via the CQ:Z flag.
#
sub bs_analyze_alignment($$) {
	my ($al, $sink) = @_;
	return unless $al->aligned;
	$al->{_toff} > 0 || confess("Bad offset: $al->{_toff}");
	my ($triml, $trimr) = ($al->triml,     $al->trimr);
	my ($alrd,  $alrf ) = ($al->{_alread}, $al->{_alref});
	my ($rdlen, $rflen) = ($al->{_rdex},   $al->{_rfex});
	my $len = length($alrd);
	my ($rfoff, $rdoff) = (0, 0);
	# Get strand info (so we know whether Y or R are relevant)
	defined($al->{_opts}{"XB"}) ||
		croak("No strand info (XB:Z) for aligned read");
	my $strand = $al->opt("XB", "Z");
	my $watson    = substr($strand, 0, 1) eq "W" ? 1 : 0;
	my $fw        = substr($strand, -1)   ne "R" ? 1 : 0;
	my $fivePleft = ($watson == $fw);
	my $refoff = $al->{_toff}-1;
	my ($cs, $cq) = (undef, undef);
	if($al->color) {
		# Get contents of CS:Z and CQ:Z optional flags.  Though they're
		# optional, here we require them for colorspace bisulfite alignments.
		($cs, $cq) = get_trimmed_cs_cq($al);
		defined($cq) || confess();
		$cq = reverse $cq if !$fivePleft;
		# Trim from LHS
		$cq = substr($cq, $al->{_striml}) if $al->{_striml} > 0;
		# Trim from RHS
		$cq = substr($cq, 0, length($cq)-$al->{_strimr}) if $al->{_strimr} > 0;
		# Now we have a trimmed color quality string that should either be one
		# longer or one shorter than the decoded nucleotide alignment.
		length($cq) == $rdlen-1 || length($cq) == $rdlen+1 ||
			croak("Length of the CQ:Z string after trimming (".
			      length($cq).") cannot be reconciled with ".
			      "alignment length ($rdlen)");
		$cq = " $cq " if length($cq) == $rdlen-1;
		length($cq) == $rdlen+1 || confess();
	}
	# Iterate over alignment positions
	for(my $i = 0; $i < $len; $i++) {
		# Get next characters in read and ref
		my ($alrdc, $alrfc) = (uc substr($alrd, $i, 1), uc substr($alrf, $i, 1));
		if($alrdc ne "-") {
			# Is alrfc the proper ambiguous characer?  For evidence on the
			# watson strand, it must be Y; for evidence on the crick strand
			# it must be R.
			# TODO: make it an option whether we report all alleles (as we do
			# here), or just the alleles that constitute methylation evidence.
			# Reporting all alleles might be helpful because it gives a hint as
			# to which CpGs are interrupted by SNPs.
			my $evpos = ( $watson && $alrfc eq "Y" && (1 || $alrdc eq "C" || $alrdc eq "T"));
			$evpos ||=  (!$watson && $alrfc eq "R" && (1 || $alrdc eq "G" || $alrdc eq "A"));
			if($evpos) {
				# This position was converted
				my $cy = $rdoff;
				$cy = $rdlen - $rdoff - 1 if !$fivePleft;
				my ($qu1, $qu2);
				if($al->color) {
					($qu1, $qu2) =
						(substr($cq, $rdoff, 1), substr($cq, $rdoff+1, 1));
					$qu1 = -1 if $qu1 eq " ";
					$qu2 = -1 if $qu2 eq " ";
				} else {
					($qu1, $qu2) = (substr($al->{_rdqual}, $rdoff, 1), -1);
				}
				defined($al->{_flags}) || confess("Expected SAM flags in alignment");
				$sink->(BtlBio::Align::Bisulfite::BsEvidence->new(
					$al->{_rdname},       # read ID
					$al->{_tname},        # ref name
					$refoff + $rfoff + 1, # 1-based ref off
					$alrdc,               # evidence
					$watson,              # whether evidence is from Watson
					$fw,                  # whether evidence is from fw read
					$al->{_flags},        # SAM flags
					$qu1,                 # positional quality #1
					$qu2,                 # positional quality #2
					$cy,                  # read cycle
					$rdlen,               # alignment length
					$al->{_score},        # alignment score
					$al->{_mapq}          # mapping quality
				));
			}
		}
		$rfoff++ if $alrfc ne "-";
		$rdoff++ if $alrdc ne "-";
	}
}

##
# Given a converted read's alignment (BtlBio::Alignment::Alignment), extract
# and commit any methylation evidence encountered.  Evidence is committed by
# calling the function ref provided as '$sink' and passing it the evidence
# wrapped in a BtlBio::Align::Bisulfite::BsEvidence object.
#
# TODO: What reference offset is being printed when the alignment is to the BS
# Crick strand?  And what allele is recorded?  I would propose that all
# referece offsets be reported w/r/t the Watson strand and that alleles be
# reported w/r/t the Watson strand.  E.g. a piece of methylation evidence on
# the Crick strand gets allele 'G', not 'C'.
#
# in: al:        alignment
# in: watson:    whether alignment is from Watson reference
# in: sink_meth: function to call with new methylation evidence
# in: sink_all:  function to call with new methylation evidence
# in: ref:       reference sequences, for checking what's a methylation site
#
sub bsc_analyze_alignment($$$$$$) {
	my ($al, $watson, $ref, $sink_meth, $sink_all, $sink_asm) = @_;
	return unless $al->aligned;
	$al->{_toff} > 0 || croak("Bad offset: $al->{_toff}");
	$al->{_fw} || croak("Expected only forward-strand alignments");
	# Try to get forward/revcomp from read name
	my $fw = substr($al->{_rdname}, -1);
	if($fw eq "+" || $fw eq "-") {
		$fw = ($fw eq "+" ? 1 : 0);
	} elsif(defined($al->{_opts}{XB})) {
		# Didn't work; try to get it from extra field
		$al->{_opts}{XB}{type} eq "Z" || croak("Expected type Z for extra field XB");
		$fw = (substr($al->{_opts}{XB}{value}, -1) eq "R" ? 0 : 1);
	} else {
		croak("Could not determine strand information");
	}
	# Get original read sequence
	my $orig_seq = undef;
	if(defined($al->{_opts}{YO})) {
		$al->{_opts}{YO}{type} eq "Z" || croak("Expected type Z for extra field YO");
		$orig_seq = $al->{_opts}{YO}{value};
	} else {
		my $cidx = rindex($al->{_rdname}, ":");
		$cidx >= 0 ||
			confess("Could not find : in read name: '".$al->{_rdname}."';".
			" send read through bsc_ize_read() first\n");
		$orig_seq = substr($al->{_rdname}, $cidx+1);
		# Remove trailing mate number and +/-
		$orig_seq = substr($orig_seq, 0, length($orig_seq)-2);
		# Remove original read sequence from end of name
		$al->{_rdname} = substr($al->{_rdname}, 0, $cidx);
	}
	# TODO: Hard trimming in the aligner could violate this assumption
	length($orig_seq) == length($al->{_rdseq}) ||
		confess("length of orig: '$orig_seq' didn't match seq: '$al->{_rdseq}'");
	# Possibly trim our original read sequence, both soft and hard
	my ($triml, $trimr) = ($al->triml, $al->trimr);
	$orig_seq = substr($orig_seq, $triml, length($orig_seq)-$triml-$trimr);
	# Get the stacked alignment
	my ($alrd, $alrf) = ($al->{_alread}, $al->{_alref});
	length($alrd) >= length($al->{_rdseq})-$al->triml-$al->trimr || confess();
	length($alrd) == length($alrf) || confess();
	# Get the reference offset of the first aligned character, i.e.
	# upstream-most character on Watson strand involved in the alignment, after
	# all trimming
	my $refoff = $al->{_toff}-1;
	my ($rdlen, $rflen) = ($al->{_rdex}, $al->{_rfex});
	$rdlen == length($orig_seq) ||
		confess("Couldn't reconcile lengths of trimmined, original sequence: ".
		    "'$orig_seq', and aligned sequence:\n$alrf\n$alrd");
	my $tname = $al->{_tname};
	defined($ref->{$tname}) || confess("Bad reference name: '$tname'");
	my $reflen = length($ref->{$tname});
	if(!$watson) {
		# If alignment is to the Crick strand, transform to Watson coordinates
		# before obtaining the reference substring.
		$refoff = $reflen - $refoff - 1;
		$refoff -= ($rflen-1);
	}
	# Get an extra character on the end
	my $refstr = substr($ref->{$tname}, $refoff, $rflen);
	defined($refstr) || confess();
	length($refstr) == $rflen || confess();
	$refstr = revcomp($refstr) if !$watson;
	# Iterate over each alignment position
	my $len = length($alrd);
	my ($rfoff, $rdoff) = (0, 0);
	for(my $i = 0; $i < $len; $i++) {
		# Get next characters in read and ref
		my $alrdc = substr($alrd, $i, 1);
		my $alrfc = substr($alrf, $i, 1);
		if($alrdc ne "-" && $alrfc ne "-") {
			# A read character appears opposite a reference character in this
			# column (no gap)
			my $orig_rdc = substr($orig_seq, $rdoff, 1);
			my $refc = uc substr($refstr, $rfoff, 1);
			my $convrdc = $orig_rdc;
			$convrdc = "T" if $orig_rdc eq "C";
			$convrdc eq $alrdc ||
				croak("Character from read side of alignment at offset $rdoff ".
				      "($alrdc) didn't match character from original read ".
				      "($convrdc):\n$orig_seq\n$alrd\n$alrf\n$refstr");
			my $convrefc = $refc;
			$convrefc = "T" if $refc eq "C";
			$convrefc eq $alrfc ||
				croak("Character from ref side of alignment at offset $rfoff ".
				      "($alrfc) didn't match character from ref ($convrefc) ".
				      ":\n$orig_seq\n$alrd\n$alrf\n$refstr");
			defined($orig_rdc) || confess();
			defined($refc) || confess();
			my $evpos = 0;
			$evpos = 1 if $refc eq "C";
			my $do_sink_meth = defined($sink_meth) && $evpos;
			my $do_sink_all = defined($sink_all) && $alrdc ne $alrfc;
			if($do_sink_meth || $do_sink_all) {
				# This position was converted
				my $cy = $rdoff;
				$cy = $rdlen - $rdoff - 1 if !$fw;
				my $qu = substr($al->{_rdqual}, $rdoff, 1);
				my $allele = $orig_rdc;
				my $evoff = $al->{_toff}-1 + $rfoff; # 0-based evidence offset
				if(!$watson) {
					# Reverse-complement the evidence and transform its offset,
					# forcing it to be with respect to the Watson strand
					$allele = revcomp($allele);
					$evoff = $reflen - $evoff - 1;
				}
				defined($al->{_flags}) || confess("Expected SAM flags in alignment");
				if($do_sink_meth) {
					$sink_meth->(BtlBio::Align::Bisulfite::BsEvidence->new(
						$al->{_rdname},       # read ID
						$tname,               # ref name
						$evoff + 1,           # 1-based ref off
						$allele,              # evidence
						$watson,              # whether read aligns to Watson
						$fw,                  # whether ev is from fw/rc read
						$al->{_flags},        # mate id
						$qu,                  # positional quality #1
						-1,                   # positional quality #2
						$cy,                  # read cycle
						$rdlen,               # alignment length
						$al->{_score},        # alignment score
						$al->{_mapq}));       # mapping quality
				}
				if($do_sink_all) {
					$sink_all->(BtlBio::Align::Bisulfite::BsEvidence->new(
						$al->{_rdname},       # read ID
						$tname,               # ref name
						$evoff + 1,           # 1-based ref off
						$allele,              # evidence
						$watson,              # whether read aligns to Watson
						$fw,                  # whether ev is from fw/rc read
						$al->{_flags},        # mate id
						$qu,                  # positional quality #1
						-1,                   # positional quality #2
						$cy,                  # read cycle
						$rdlen,               # alignment length
						$al->{_score},        # alignment score
						$al->{_mapq}));       # mapping quality
				}
			}
		}
		$rfoff++ if $alrfc ne "-";
		$rdoff++ if $alrdc ne "-";
	}
}

##
# Given a converted read's alignment (BtlBio::Alignment::Alignment), extract
# and commit any methylation evidence encountered.  Evidence is committed by
# calling the function ref provided as '$sink' and passing it the evidence
# wrapped in a BtlBio::Align::Bisulfite::BsEvidence object.
#
# in: al:        alignment
# in: watson:    whether alignment is from Watson reference
# in: sink_meth: function to call with new methylation evidence
# in: sink_all:  function to call with new methylation evidence
# in: ref:       reference sequences, for checking what's a methylation site
#
sub bscpg_analyze_alignment($$$$$$) {
	my ($al, $watson, $ref, $sink_meth, $sink_all, $sink_asm) = @_;
	return unless $al->aligned;
	$al->{_toff} > 0 || croak("Bad offset: $al->{_toff}");
	$al->{_fw} || croak("Expected only forward-strand alignments");
	# Try to get forward/revcomp from read name
	my $fw = substr($al->{_rdname}, -1);
	if($fw eq "+" || $fw eq "-") {
		$fw = ($fw eq "+" ? 1 : 0);
	} elsif(defined($al->{_opts}{XB})) {
		# Didn't work; try to get it from extra field
		$al->{_opts}{XB}{type} eq "Z" || croak("Expected type Z for extra field XB");
		$fw = (substr($al->{_opts}{XB}{value}, -1) eq "R" ? 0 : 1);
	} else {
		croak("Could not determine strand information");
	}
	# Get original read sequence
	my $orig_seq = undef;
	if(defined($al->{_opts}{YO})) {
		$al->{_opts}{YO}{type} eq "Z" || croak("Expected type Z for extra field YO");
		$orig_seq = $al->{_opts}{YO}{value};
	} else {
		my $cidx = rindex($al->{_rdname}, ":");
		$cidx >= 0 ||
			confess("Could not find : in read name: '".$al->{_rdname}."';".
			" send read through bsc_ize_read() first\n");
		$orig_seq = substr($al->{_rdname}, $cidx+1);
		# Remove trailing mate number and +/-
		$orig_seq = substr($orig_seq, 0, length($orig_seq)-2);
		# Remove original read sequence from end of name
		$al->{_rdname} = substr($al->{_rdname}, 0, $cidx);
	}
	# TODO: Hard trimming in the aligner could violate this assumption
	length($orig_seq) == length($al->{_rdseq}) ||
		confess("length of orig: '$orig_seq' didn't match seq: '$al->{_rdseq}'");
	# Possibly trim our original read sequence, both soft and hard
	my ($triml, $trimr) = ($al->triml, $al->trimr);
	$orig_seq = substr($orig_seq, $triml, length($orig_seq)-$triml-$trimr);
	# Get the stacked alignment
	my ($alrd, $alrf) = ($al->{_alread}, $al->{_alref});
	length($alrd) >= length($al->{_rdseq})-$al->triml-$al->trimr || confess();
	length($alrd) == length($alrf) || confess();
	# Get the reference offset of the first aligned character, i.e.
	# upstream-most character on Watson strand involved in the alignment, after
	# all trimming
	my $refoff = $al->{_toff}-1;
	my ($rdlen, $rflen) = ($al->{_rdex}, $al->{_rfex});
	$rdlen == length($orig_seq) ||
		confess("Couldn't reconcile lengths of trimmined, original sequence: ".
		    "'$orig_seq', and aligned sequence:\n$alrf\n$alrd");
	my $tname = $al->{_tname};
	defined($ref->{$tname}) || confess("Bad reference name: '$tname'");
	my $reflen = length($ref->{$tname});
	if(!$watson) {
		# If alignment is to the Crick strand, transform to Watson coordinates
		# before obtaining the reference substring.
		$refoff = $reflen - $refoff - 1;
		$refoff -= ($rflen-1);
	}
	# Get an extra character on the end
	my $refstr = substr($ref->{$tname}, $refoff, $rflen);
	defined($refstr) || confess();
	length($refstr) == $rflen || length($refstr) == $rflen+1 || confess();
	$refstr = revcomp($refstr) if !$watson;
	# Now pad on the right-hand side
	if($watson) {
		if(length($ref->{$tname}) > $refoff + $rflen) {
			$refstr .= substr($ref->{$tname}, $refoff+$rflen, 1);
		} else {
			$refstr .= "N";
		}
	} else {
		if($refoff > 0) {
			$refstr .= comp(substr($ref->{$tname}, $refoff-1, 1));
		} else {
			$refstr .= "N";
		}
	}
	# Iterate over each alignment position
	my $len = length($alrd);
	my ($rfoff, $rdoff) = (0, 0);
	my @ev_cpg = ();
	my @ev_snp = ();
	for(my $i = 0; $i < $len; $i++) {
		# Get next characters in read and ref
		my $alrdc = substr($alrd, $i, 1);
		my $alrfc = substr($alrf, $i, 1);
		if($alrdc ne "-" && $alrfc ne "-") {
			# A read character appears in this column (not a gap)
			my $orig_rdc = substr($orig_seq, $rdoff, 1);
			my $refc = uc substr($refstr, $rfoff, 1);
			my $convrdc = $orig_rdc;
			$convrdc = "T" if $orig_rdc eq "C";
			$convrdc eq $alrdc ||
				croak("Character from read side of alignment at offset $rdoff ".
				      "($alrdc) didn't match character from original read ".
				      "($convrdc):\n$orig_seq\n$alrd\n$alrf\n$refstr");
			my $convrefc = $refc;
			$convrefc = "T" if $refc eq "C";
			$convrefc eq $alrfc ||
				croak("Character from ref side of alignment at offset $rfoff ".
				      "($alrfc) didn't match character from ref ($convrefc) ".
				      ":\n$orig_seq\n$alrd\n$alrf\n$refstr");
			my $refcs = uc substr($refstr, $rfoff, 2);
			defined($orig_rdc) || confess();
			defined($refcs) || confess();
			length($refcs) == 2 || confess();
			my $evpos = 0;
			$evpos = 1 if $refcs eq "CG";
			my $do_sink_meth = defined($sink_meth) && $evpos;
			my $do_sink_all = defined($sink_all) && $alrdc ne $alrfc && substr($refcs, 0, 1) ne "C";
			if($do_sink_meth || $do_sink_all) {
				# This position was converted
				my $cy = $rdoff;
				$cy = $rdlen - $rdoff - 1 if !$fw;
				my $qu = substr($al->{_rdqual}, $rdoff, 1);
				my $allele = $orig_rdc;
				my $evoff = $al->{_toff}-1 + $rfoff; # 0-based evidence offset
				if(!$watson) {
					# Reverse-complement the evidence and transform its offset,
					# forcing it to be with respect to the Watson strand
					$allele = revcomp($allele);
					$evoff = $reflen - $evoff - 1;
				}
				defined($al->{_flags}) || confess("Expected SAM flags in alignment");
				my $ev = BtlBio::Align::Bisulfite::BsEvidence->new(
					$al->{_rdname},       # read ID
					$tname,               # ref name
					$evoff + 1,           # 1-based ref off
					$allele,              # evidence
					$watson,              # whether read aligns to Watson
					$fw,                  # whether ev is from fw/rc read
					$al->{_flags},        # mate id
					$qu,                  # positional quality #1
					-1,                   # positional quality #2
					$cy,                  # read cycle
					$rdlen,               # alignment length
					$al->{_score},        # alignment score
					$al->{_mapq});        # mapping quality
				if(substr($refcs, 0, 1) ne "C") {
					push @ev_snp, $ev;
					$sink_all->($ev, 1) if defined($sink_all);
				} else {
					push @ev_cpg, $ev;
					$sink_meth->($ev) if defined($sink_meth);
					$sink_all->($ev, 0) if defined($sink_all);
				}
			}
		}
		$rfoff++ if $alrfc ne "-";
		$rdoff++ if $alrdc ne "-";
	}
	# Output all pairs of CpG and SNP evidence
	if(defined($sink_asm)) {
		for(my $i = 0; $i < scalar(@ev_cpg); $i++) {
			for(my $j = 0; $j < scalar(@ev_snp); $j++) {
				$sink_asm->($ev_cpg[$i], $ev_snp[$j]);
			}
		}
	}
}

sub _test_bs_analyze_alignment($$$$$$$$$$$) {
	my ($read, $color, $md, $cigar, $cs,
	    $cq, $watson, $fw, $score, $mapq, $flags) = @_;
	my $qual = "I" x length($read);
	my %opts = ();
	$opts{CS}{type}  = "Z" if defined($cs);
	$opts{CS}{value} = $cs if defined($cs);
	$opts{CQ}{type}  = "Z" if defined($cq);
	$opts{CQ}{value} = $cq if defined($cq);
	$opts{XB}{type}  = "Z";
	$opts{XB}{value} = ($watson ? "W" : "C").($fw ? "" : "R");
	my $al = BtlBio::Alignment::Alignment->new(
		"read1",    # read name
		$read,      # read sequence from alignment record
		$qual,      # qualities from alignment record
		$color,     # colorspace?
		"ref1",     # name of reference sequence ("text")
		1,          # offset into reference
		$fw,        # true -> aligned to Watson strand
		$score,     # alignment score
		$mapq,      # mapping quality
		$md,        # MD:Z string
		$cigar,     # CIGAR string
		(($flags & 0x40) != 0), # mate 1
		(($flags & 0x80) != 0), # mate 2
		$flags,     # SAM flags (if available)
		\%opts,     # optional flags
		undef);     # original string
	my @ev = ();
	my $sink = sub { push @ev, $_[0]; };
	bs_analyze_alignment($al, $sink);
	return @ev;
}

sub _test_bsc_analyze_alignment($$$$$$$$$) {
	my ($orig_read, $ref1, $md, $cigar, $watson, $fw, $score, $mapq, $flags) = @_;
	$orig_read = revcomp($orig_read) unless $fw;
	my $read = $orig_read; $read =~ tr/C/T/;
	my $qual      = "I" x length($read);
	my %refs = ( "ref1" => $ref1 );
	$fw = $fw ? "+" : "-";
	my $mate = 0;
	$mate = 1 if (($flags & 0x40) != 0);
	$mate = 2 if (($flags & 0x80) != 0);
	my $al = BtlBio::Alignment::Alignment->new(
		"read1:$orig_read$mate$fw", # read name
		$read,      # read sequence from alignment record
		$qual,      # qualities from alignment record
		0,          # colorspace?
		"ref1",     # name of reference sequence ("text")
		1,          # offset into reference
		1,          # true -> aligned to Watson strand
		$score,     # alignment score
		$mapq,      # mapping quality
		$md,        # MD:Z string
		$cigar,     # CIGAR string
		(($flags & 0x40) != 0), # mate 1
		(($flags & 0x80) != 0), # mate 2
		$flags,     # SAM flags (if available)
		{},         # optional flags
		undef);     # original string
	my @ev = ();
	my $sink_meth = sub { push @ev, $_[0]; };
	my @ev2 = ();
	my $sink_all = sub { push @ev2, $_[0]; };
	bsc_analyze_alignment($al, $watson, \%refs, $sink_meth, $sink_all, sub {});
	return @ev;
}

sub _test_bscpg_analyze_alignment($$$$$$$$$$) {
	my ($orig_read, $ref1, $refoff, $md, $cigar, $watson, $fw, $score, $mapq, $flags) = @_;
	$orig_read = revcomp($orig_read) unless $fw;
	my $read = $orig_read; $read =~ tr/C/T/;
	my $qual      = "I" x length($read);
	my %refs = ( "ref1" => $ref1 );
	$fw = $fw ? "+" : "-";
	my $mate = 0;
	$mate = 1 if (($flags & 0x40) != 0);
	$mate = 2 if (($flags & 0x80) != 0);
	my $al = BtlBio::Alignment::Alignment->new(
		"read1:$orig_read$mate$fw", # read name
		$read,      # read sequence from alignment record
		$qual,      # qualities from alignment record
		0,          # colorspace?
		"ref1",     # name of reference sequence ("text")
		$refoff,    # offset into reference
		1,          # true -> aligned to Watson strand
		$score,     # alignment score
		$mapq,      # mapping quality
		$md,        # MD:Z string
		$cigar,     # CIGAR string
		(($flags & 0x40) != 0), # mate 1
		(($flags & 0x80) != 0), # mate 2
		$flags,     # SAM flags (if available)
		{},         # optional flags
		undef);     # original string
	my @ev_meth = ();
	my @ev_all = ();
	my $sink_meth = sub { push @ev_meth, $_[0]; };
	my $sink_all = sub { push @ev_all, $_[0]; };
	bscpg_analyze_alignment($al, $watson, \%refs, $sink_meth, $sink_all, sub { });
	return @ev_meth;
}

sub _test_evidence($$$$$$$$$$$$$) {
	my ($ev, $tname, $off, $rdid, $al, $watson,
	    $fw, $qua1, $qua2, $cyc, $score,
	    $mapq, $flags) = @_;
	$ev->read_id eq $rdid    || croak("Expected read_id=$rdid, got "  .$ev->read_id);
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
	$ev->flags  == $flags    || croak("Expected flags=$flags, got "   .$ev->flags);
	return 1;
}

sub _test_bs_analyze_alignment_1() {
	print STDERR "Testing bs_analyze_alignment 1 ... ";
	require Data::Dumper;
	my @ev = _test_bs_analyze_alignment(
		"AAAAA", # read
		0,       # colorspace?
		"5",     # MD:Z
		"5M",    # CIGAR
		undef,   # CS:Z
		undef,   # CQ:Z
		1,       # watson
		1,       # fw
		20,      # score
		40,      # MAPQ
		16);     # SAM flags
	scalar(@ev) == 0 || die "Expected 0 records, got ".Dumper(\@ev);
	print STDERR "PASSED\n";
}

sub _test_bs_analyze_alignment_2() {
	# Same as _test_bs_analyze_alignment_1 but with all nucleotides
	# represented in the read
	print STDERR "Testing bs_analyze_alignment 2 ... ";
	require Data::Dumper;
	my @ev = _test_bs_analyze_alignment(
		"TACAG", # read
		0,       # colorspace?
		"5",     # MD:Z
		"5M",    # CIGAR
		undef,   # CS:Z
		undef,   # CQ:Z
		1,       # watson
		1,       # fw
		20,      # score
		40,      # MAPQ
		16);     # SAM flags
	scalar(@ev) == 0 || die "Expected 0 records, got ".Dumper(\@ev);
	print STDERR "PASSED\n";
}

sub _test_bs_analyze_alignment_3() {
	print STDERR "Testing bs_analyze_alignment 3 ... ";
	require Data::Dumper;
	my @ev = _test_bs_analyze_alignment(
		"TACAG", # read
		0,       # colorspace?
		"2Y2",   # MD:Z
		"5M",    # CIGAR
		undef,   # CS:Z
		undef,   # CQ:Z
		1,       # watson
		1,       # fw
		20,      # score
		40,      # MAPQ
		16);     # SAM flags
	scalar(@ev) == 1 || die "Expected 1 records, got ".Dumper(\@ev);
	_test_evidence($ev[0],
		"ref1", # tname
		3,      # toff (1-based)
		"read1",# read ID
		"C",    # allele
		1,      # watson
		1,      # fw
		"I",    # qual 1
		-1,     # qual 2
		2,      # cyc (0-based)
		20,     # score
		40,      # MAPQ
		16);     # SAM flags
	print STDERR "PASSED\n";
}

sub _test_bs_analyze_alignment_4() {
	print STDERR "Testing bs_analyze_alignment 4 ... ";
	require Data::Dumper;
	my @ev = _test_bs_analyze_alignment(
		"TACAG", # read
		0,       # colorspace?
		"2Y0R1", # MD:Z
		"5M",    # CIGAR
		undef,   # CS:Z
		undef,   # CQ:Z
		0,       # watson
		0,       # fw
		24,      # score
		40,      # MAPQ
		16);     # SAM flags
	scalar(@ev) == 1 || die "Expected 1 records, got ".Dumper(\@ev);
	_test_evidence($ev[0],
		"ref1", # tname
		4,      # toff (1-based)
		"read1",# read ID
		"A",    # allele
		0,      # watson
		0,      # fw
		"I",    # qual 1
		-1,     # qual 2
		3,      # cyc (0-based)
		24,     # score
		40,      # MAPQ
		16);     # SAM flags
	print STDERR "PASSED\n";
}

sub _test_bs_analyze_alignment_5() {
	print STDERR "Testing bs_analyze_alignment 5 ... ";
	my ($read, $color, $md, $cigar, $cs, $cq, $watson, $fw, $score, $mapq) = @_;
	require Data::Dumper;
	my @ev = _test_bs_analyze_alignment(
		"TACAGGGGGGGG", # decoded nucleotides
		#311200000001T
		#  \|  \|  \
		#TAYRGGGRGGGR
		#LKJIHGFEDCBA
		#109876543210
		1,       # colorspace?
		"2Y0R3R3R0", # MD:Z
		"12M",       # CIGAR
		"T100000002113", # CS:Z
		"ABCDEFGHIJKL",  # CQ:Z
		0,       # watson
		1,       # fw
		17,      # score
		18,      # MAPQ
		16);     # SAM flags
	scalar(@ev) == 3 || die "Expected 3 records, got ".Dumper(\@ev);
	_test_evidence($ev[0],
		"ref1", # tname
		4,      # toff (1-based)
		"read1",# read ID
		"A",    # allele
		0,      # watson
		1,      # fw
		"J",    # qual 1
		"I",    # qual 2
		8,      # cyc (0-based)
		17,     # score
		18,      # MAPQ
		16);     # SAM flags
	_test_evidence($ev[1],
		"ref1", # tname
		8,      # toff (1-based)
		"read1",# read ID
		"G",    # allele
		0,      # watson
		1,      # fw
		"F",    # qual 1
		"E",    # qual 2
		4,      # cyc (0-based)
		17,     # score
		18,      # MAPQ
		16);     # SAM flags
	_test_evidence($ev[2],
		"ref1", # tname
		12,     # toff (1-based)
		"read1",# read ID
		"G",    # allele
		0,      # watson
		1,      # fw
		"B",    # qual 1
		-1,     # qual 2
		0,      # cyc (0-based)
		17,     # score
		18,      # MAPQ
		16);     # SAM flags
	print STDERR "PASSED\n";
}

sub _test_bs_analyze_alignment_6() {
	print STDERR "Testing bs_analyze_alignment 6 ... ";
	require Data::Dumper;
	#       9876-543210
	# read: AAAA-TCAAAA
	# ref:  AAAAAYY-AAA
	#       0123456-789
	my @ev = _test_bs_analyze_alignment(
		"AAAATCAAAA",  # read
		0,             # colorspace
		"4^A0Y0Y4",    # MD:Z
		"4M1D2M1I3M",  # CIGAR
		undef,         # CS:Z
		undef,         # CQ:Z
		1,             # watson
		0,             # fw
		-7,            # score
		31,            # MAPQ
		16);           # SAM flags
	scalar(@ev) == 2 || die "Expected 2 records, got ".Dumper(\@ev);
	_test_evidence($ev[0], "ref1", 6, "read1", "T", 1, 0, "I", -1, 5, -7, 31, 16);
	_test_evidence($ev[1], "ref1", 7, "read1", "C", 1, 0, "I", -1, 4, -7, 31, 16);
	print STDERR "PASSED\n";
}

sub _test_bs_analyze_alignment_7() {
	print STDERR "Testing bs_analyze_alignment 7 ... ";
	require Data::Dumper;
	#  10987654  3210
	# aACAAAACT--ATAAa
	#  /|   /||  /|
	# 011000123  3300
	# 012345678  9012
	#  AYA---YYAAAYAA
	#  123   45678901
	my ($watson, $fw) = (1, 0);
	my @ev = _test_bs_analyze_alignment(
		"ACAAAACTATAA",   # read
		1,                # colorspace
		"1Y1Y0Y^AA1Y2",   # MD:Z
		"3M3I2M2D4M",     # CIGAR
		"0033321000110",  # CS:Z
		"2109876543210",  # CQ:Z
		$watson,          # watson
		$fw,              # fw
		99,               # score
		20,               # MAPQ
		128);             # SAM flags
	scalar(@ev) == 4 || die "Expected 4 records, got ".Dumper(\@ev);
	my $cyc = 10; $cyc = 12 - $cyc - 1 if $fw;
	_test_evidence($ev[0], "ref1", 2, "read1", "C", $watson, $fw, "1", "2", $cyc, 99, 20, 128);
	$cyc = 5; $cyc = 12 - $cyc - 1 if $fw;
	_test_evidence($ev[1], "ref1", 4, "read1", "C", $watson, $fw, "6", "7", $cyc, 99, 20, 128);
	$cyc = 4; $cyc = 12 - $cyc - 1 if $fw;
	_test_evidence($ev[2], "ref1", 5, "read1", "T", $watson, $fw, "7", "8", $cyc, 99, 20, 128);
	$cyc = 2; $cyc = 12 - $cyc - 1 if $fw;
	_test_evidence($ev[3], "ref1", 9, "read1", "T", $watson, $fw, "9", "0", $cyc, 99, 20, 128);
	print STDERR "PASSED\n";
}

sub _test_bsc_analyze_alignment_1() {
	print STDERR "Testing bsc_analyze_alignment 1 ... ";
	require Data::Dumper;
	my @ev = _test_bsc_analyze_alignment(
		"AAAAA", # orig read
		"AAAAA", # ref
		"5",     # MD:Z
		"5M",    # CIGAR
		1,       # watson
		1,       # fw
		20,      # score
		40,      # MAPQ
		16);     # SAM flags
	scalar(@ev) == 0 || die "Expected 0 records, got ".Dumper(\@ev);
	print STDERR "PASSED\n";
}

sub _test_bsc_analyze_alignment_2() {
	print STDERR "Testing bsc_analyze_alignment 2 ... ";
	require Data::Dumper;
	my @ev = _test_bsc_analyze_alignment(
		"AACAA", # orig read
		"AACAA", # ref
		"5",     # MD:Z
		"5M",    # CIGAR
		1,       # watson
		1,       # fw
		-40,     # score
		40,      # MAPQ
		16);     # SAM flags
	scalar(@ev) == 1 || die "Expected 1 records, got ".Dumper(\@ev);
	_test_evidence(
		$ev[0],
		"ref1",
		3,
		"read1",
		"C",
		1,      # watson
		1,      # fw
		"I",    # qual 1
		-1,     # qual 2
		2,
		-40,
		40,
		16);
	print STDERR "PASSED\n";
}

sub _test_bsc_analyze_alignment_3() {
	print STDERR "Testing bsc_analyze_alignment 3 ... ";
	require Data::Dumper;
	my @ev = _test_bsc_analyze_alignment(
		"AACTAA",           # orig read
		#123456
		#TTGGTT
		revcomp("AACCAA"),  # ref
		"6",       # MD:Z
		"6M",      # CIGAR
		0,         # watson
		1,         # fw
		0,         # score
		35,      # MAPQ
		16);     # SAM flags
	scalar(@ev) == 2 || die "Expected 2 records, got ".Dumper(\@ev);
	_test_evidence($ev[0], "ref1", 4, "read1", "G", 0, 1, "I", -1, 2, 0, 35, 16);
	_test_evidence($ev[1], "ref1", 3, "read1", "A", 0, 1, "I", -1, 3, 0, 35, 16);
	print STDERR "PASSED\n";
}

sub _test_bsc_analyze_alignment_4() {
	print STDERR "Testing bsc_analyze_alignment 4 ... ";
	require Data::Dumper;
	my @ev = _test_bsc_analyze_alignment(
		revcomp("AACAA"), # orig read
		"AACAA", # ref
		"5",     # MD:Z
		"5M",    # CIGAR
		1,       # watson
		0,       # fw
		-24,     # score
		40,      # MAPQ
		16);     # SAM flags
	scalar(@ev) == 1         || die "Expected 1 records, got ".Dumper(\@ev);
	$ev[0]->fw == 0          || die "Expected fw=0, got "     .$ev[0]->fw;
	$ev[0]->watson == 1      || die "Expected watson=1, got " .$ev[0]->watson;
	$ev[0]->aln_score == -24 || die "Expected score=-24, got ".$ev[0]->aln_score;
	print STDERR "PASSED\n";
}

sub _test_bsc_analyze_alignment_5() {
	print STDERR "Testing bsc_analyze_alignment 5 ... ";
	require Data::Dumper;
	# AAAACC-AAAA
	# AAA-CCAAAAA
	my @ev = _test_bsc_analyze_alignment(
		"AAAACCAAAA", # orig read
		revcomp("AAACCAAAAA"), # ref
		#TTTTTGGTTT
		#1234567890
		"5^A4",       # MD:Z
		"3M1I2M1D4M", # CIGAR
		0,            # watson
		1,            # fw
		4,            # MD:Z
		33,      # MAPQ
		16);     # SAM flags
	scalar(@ev) == 2 || die "Expected 2 records, got ".Dumper(\@ev);
	_test_evidence($ev[0], "ref1", 7, "read1", "G", 0, 1, "I", -1, 4, 4, 33, 16);
	_test_evidence($ev[1], "ref1", 6, "read1", "G", 0, 1, "I", -1, 5, 4, 33, 16);
	print STDERR "PASSED\n";
}

sub _test_bsc_analyze_alignment_6() {
	print STDERR "Testing bsc_analyze_alignment 6 ... ";
	require Data::Dumper;
	# 9876 543210
	# AAAA-TCAAAA
	# AAAAACC-AAA
	my @ev = _test_bsc_analyze_alignment(
		revcomp("AAAATCAAAA"),  # orig read
		revcomp("AAAAACCAAA"),  # ref
		#TTTTGA-TTTT
		#TTT-GGTTTTT
		#123 4567890
		"4^A6",        # MD:Z
		"4M1D2M1I3M",  # CIGAR
		0,             # watson
		0,             # fw
		-7,
		31,      # MAPQ
		16);     # SAM flags
	scalar(@ev) == 2 || die "Expected 2 records, got ".Dumper(\@ev);
	_test_evidence($ev[0], "ref1", 5, "read1", "A", 0, 0, "I", -1, 5, -7, 31, 16);
	_test_evidence($ev[1], "ref1", 4, "read1", "G", 0, 0, "I", -1, 4, -7, 31, 16);
	print STDERR "PASSED\n";
}

sub _test_bsc_analyze_alignment_7() {
	print STDERR "Testing bsc_analyze_alignment 7 ... ";
	require Data::Dumper;
	# 10987654  3210
	# ACAAAACT--ATAA
	# ACA---CCAAACAA
	# 012   34567890
	for(my $watson = 0; $watson < 2; $watson++) {
		for(my $fw = 0; $fw < 2; $fw++) {
			my $orig_read = "ACAAAACTATAA";
			my $ref       = "ACACCAAACAA";
			my @ev = _test_bsc_analyze_alignment(
				$fw ? $orig_read : revcomp($orig_read),
				$watson ? $ref : revcomp($ref),
				"5^AA4",
				"3M3I2M2D4M",
				$watson,
				$fw,
				99,
				20,      # MAPQ
				16);     # SAM flags
			scalar(@ev) == 4 || die "Expected 4 records, got ".Dumper(\@ev);
			my $cyc = 10; $cyc = 12 - $cyc - 1 if $fw;
			my $refoff = 2; $refoff = 11 - $refoff + 1 if !$watson;
			_test_evidence($ev[0], "ref1", $refoff, "read1", $watson ? "C" : "G", $watson, $fw, "I", -1, $cyc, 99, 20, 16);
			$cyc = 5; $cyc = 12 - $cyc - 1 if $fw;
			$refoff = 4; $refoff = 11 - $refoff + 1 if !$watson;
			_test_evidence($ev[1], "ref1", $refoff, "read1", $watson ? "C" : "G", $watson, $fw, "I", -1, $cyc, 99, 20, 16);
			$cyc = 4; $cyc = 12 - $cyc - 1 if $fw;
			$refoff = 5; $refoff = 11 - $refoff + 1 if !$watson;
			_test_evidence($ev[2], "ref1", $refoff, "read1", $watson ? "T" : "A", $watson, $fw, "I", -1, $cyc, 99, 20, 16);
			$cyc = 2; $cyc = 12 - $cyc - 1 if $fw;
			$refoff = 9; $refoff = 11 - $refoff + 1 if !$watson;
			_test_evidence($ev[3], "ref1", $refoff, "read1", $watson ? "T" : "A", $watson, $fw, "I", -1, $cyc, 99, 20, 16);
		}
	}
	print STDERR "PASSED\n";
}

sub _test_bsc_analyze_alignment_8() {
	print STDERR "Testing bsc_analyze_alignment 8 ... ";
	#                                        0         1         2         3
	#                                        012345678901234567890123456789012345
	my $orig_read = uc                      "AAGAGTcGGtTATATATGtGGGGGtGtGTATGcTAT";
	my %refs = (
		"ref1" =>      "GCGTACATCATCTGTACGGCGCTTGCTACGTACTATTAAAAAAAACGATCGATCTACGCTAGCGTAAAA",
	#                   012345678901234567890123456789012345678901234567890123456789012345678
	#                   0         1         2         3         4         5         6
		"ref2" =>      "AGCGTATCACATCGATGCTAGAAGAGTCGGCTATATATGCGGGGGCGCGTATGCTATGGTCTAGCT",
	#                   012345678901234567890123456789012345678901234567890123456789012345
	#                   0         1         2         3         4         5         6
	);
	my $mapq = 30;
	my $score = 77;
	my $md = "36";
	my $cigar = "36M";
	for(my $watson = 1; $watson < 2; $watson++) {
		for(my $color = 0; $color < 2; $color++) {
			my $fw = 1;
			my $read = $orig_read;
			$read =~ tr/C/T/;
			my $qual = "I" x length($read);
			my $fwstr = $fw ? "+" : "-";
			my %myrefs = %refs;
			if(!$watson) {
				for my $k (keys %myrefs) {
					$myrefs{$k} = revcomp($myrefs{$k});
				}
			}
			my $al = BtlBio::Alignment::Alignment->new(
				"read1:${orig_read}1$fwstr", # read name w/ original read sequence
				$read,      # read sequence from alignment record
				$qual,      # qualities from alignment record
				$color,     # colorspace?
				"ref2",     # name of reference sequence ("text")
				22,         # offset into reference
				1,          # true -> aligned to Watson strand
				$score,     # alignment score
				$mapq,      # mapping quality
				$md,        # MD:Z string
				$cigar,     # CIGAR string
				1, 0,       # mate 1/2
				16,         # SAM flags
				{},         # optional flags
				undef);     # original string
			my @ev_meth = ();
			my @ev_all = ();
			bsc_analyze_alignment($al, $watson, \%refs, sub { push @ev_meth, $_[0]; }, sub { push @ev_all, $_[0]; }, sub {});
			scalar(@ev_meth) == 6 || die;
			scalar(@ev_all) == 0 || die;
			my $rdlen = 36;
			_test_evidence($ev_meth[0], "ref2", $watson ? 28 : 66-28+1, "read1", $watson ? "C" : "G", $watson, $fw, "I", -1, $fw ?  6 : $rdlen -  6 - 1, $score, $mapq, 16);
			_test_evidence($ev_meth[1], "ref2", $watson ? 31 : 66-31+1, "read1", $watson ? "T" : "A", $watson, $fw, "I", -1, $fw ?  9 : $rdlen -  9 - 1, $score, $mapq, 16);
			_test_evidence($ev_meth[2], "ref2", $watson ? 40 : 66-40+1, "read1", $watson ? "T" : "A", $watson, $fw, "I", -1, $fw ? 18 : $rdlen - 18 - 1, $score, $mapq, 16);
			_test_evidence($ev_meth[3], "ref2", $watson ? 46 : 66-46+1, "read1", $watson ? "T" : "A", $watson, $fw, "I", -1, $fw ? 24 : $rdlen - 24 - 1, $score, $mapq, 16);
			_test_evidence($ev_meth[4], "ref2", $watson ? 48 : 66-48+1, "read1", $watson ? "T" : "A", $watson, $fw, "I", -1, $fw ? 26 : $rdlen - 26 - 1, $score, $mapq, 16);
			_test_evidence($ev_meth[5], "ref2", $watson ? 54 : 66-54+1, "read1", $watson ? "C" : "G", $watson, $fw, "I", -1, $fw ? 32 : $rdlen - 32 - 1, $score, $mapq, 16);
		}
	}
	print STDERR "PASSED\n";
}

sub _test_bscpg_analyze_alignment_1() {
	print STDERR "Testing bscpg_analyze_alignment 1 ... ";
	require Data::Dumper;
	my @ev = _test_bscpg_analyze_alignment(
		"AAAAA", # orig read
		revcomp("AAAAA"), # ref
		1,       # ref off
		"5",     # MD:Z
		"5M",    # CIGAR
		0,       # watson
		1,       # fw
		-2,      # score
		40,      # MAPQ
		16);     # SAM flags
	scalar(@ev) == 0 || die "Expected 0 records, got ".Dumper(\@ev);
	print STDERR "PASSED\n";
}

sub _test_bscpg_analyze_alignment_2() {
	print STDERR "Testing bscpg_analyze_alignment 2 ... ";
	require Data::Dumper;
	my @ev = _test_bscpg_analyze_alignment(
		"AACAA", # orig read
		"AACAA", # ref
		1,       # ref off
		"5",     # MD:Z
		"5M",    # CIGAR
		1,       # watson
		1,       # fw
		-3,      # score
		40,      # MAPQ
		16);     # SAM flags
	scalar(@ev) == 0 || die "Expected 0 records, got ".Dumper(\@ev);
	print STDERR "PASSED\n";
}

sub _test_bscpg_analyze_alignment_3() {
	print STDERR "Testing bscpg_analyze_alignment 3 ... ";
	require Data::Dumper;
	# 01234
	# AACGA
	# AACGA
	# 01234
	for(my $watson = 1; $watson >= 0; $watson--) {
		for(my $fw = 0; $fw < 2; $fw++) {
			#print STDERR "Watson=$watson, fw=$fw\n";
			my @ev = _test_bscpg_analyze_alignment(
				$fw     ? "AACGA" : "TCGTT", # orig read
				$watson ? "AACGA" : "TCGTT", # ref
				1,       # ref off
				"5",     # MD:Z
				"5M",    # CIGAR
				$watson, # watson
				$fw,     # fw
				32,      # score
				40,      # MAPQ
				16);     # SAM flags
			scalar(@ev) == 1 || die "Expected 1 record, got ".Dumper(\@ev);
			my $cyc = 2; $cyc = 5 - $cyc - 1 unless $fw;
			_test_evidence(
				$ev[0],
				"ref1",
				$watson ? 3 : 5-3+1,
				"read1", 
				$watson ? "C" : "G",
				$watson,
				$fw,
				"I",
				-1,
				$cyc,
				32,
				40,
				16);
		}
	}
	print STDERR "PASSED\n";
}

sub _test_bscpg_analyze_alignment_4() {
	print STDERR "Testing bscpg_analyze_alignment 4 ... ";
	require Data::Dumper;
	# 10987654  3210
	# ACGAAACG--ATGA
	# ACG---CGAAACGA
	# 012   34567890
	for(my $watson = 0; $watson < 2; $watson++) {
		for(my $fw = 0; $fw < 2; $fw++) {
			my @ev = _test_bscpg_analyze_alignment(
				$fw ? "ACGAAACGATGA" : revcomp("ACGAAACGATGA"), # orig read
				$watson ? "ACGCGAAACGA" : revcomp("ACGCGAAACGA"),  # ref
				1, "5^AA4", "3M3I2M2D4M",
				$watson, $fw, 23, 20, 16);
			scalar(@ev) == 3 || die "Expected 3 records, got ".Dumper(\@ev);
			my $cyc = 10; $cyc = 12 - $cyc - 1 if $fw;
			_test_evidence($ev[0], "ref1", $watson ? 2 : 11-2+1, "read1", $watson ? "C" : "G", $watson, $fw, "I", -1, $cyc, 23, 20, 16);
			$cyc = 5; $cyc = 12 - $cyc - 1 if $fw;
			_test_evidence($ev[1], "ref1", $watson ? 4 : 11-4+1, "read1", $watson ? "C" : "G", $watson, $fw, "I", -1, $cyc, 23, 20, 16);
			$cyc = 2; $cyc = 12 - $cyc - 1 if $fw;
			_test_evidence($ev[2], "ref1", $watson ? 9 : 11-9+1, "read1", $watson ? "T" : "A", $watson, $fw, "I", -1, $cyc, 23, 20, 16);
		}
	}
	print STDERR "PASSED\n";
}

sub _test_bscpg_analyze_alignment_5() {
	print STDERR "Testing bscpg_analyze_alignment 5 ... ";
	require Data::Dumper;
	# 10987654  3210
	# TTTAAAGC--TGGA
	# ATC---GCGACGTA
	# 012   34567890
	for(my $watson = 1; $watson >= 0; $watson--) {
		for(my $fw = 1; $fw >= 0; $fw--) {
			my @ev = _test_bscpg_analyze_alignment(
				$fw     ? "ATTAAAGCTGGA" : revcomp("ATTAAAGCTGGA"), # orig read
				$watson ? "ATCGCGACGTA"  : revcomp("ATCGCGACGTA"),  # ref
				1,
				"A4^GA2T1",
				"3M3I2M2D4M",
				$watson,
				$fw,
				11,
				20,
				8);
			scalar(@ev) == 3 || die "Expected 3 records, got ".Dumper(\@ev);
			my $cyc = 9; $cyc = 12 - $cyc - 1 if $fw;
			_test_evidence($ev[0], "ref1", $watson ? 3 : 11-3+1, "read1", $watson ? "T" : "A", $watson, $fw, "I", -1, $cyc, 11, 20, 8);
			$cyc = 4; $cyc = 12 - $cyc - 1 if $fw;
			_test_evidence($ev[1], "ref1", $watson ? 5 : 11-5+1, "read1", $watson ? "C" : "G", $watson, $fw, "I", -1, $cyc, 11, 20, 8);
			$cyc = 3; $cyc = 12 - $cyc - 1 if $fw;
			_test_evidence($ev[2], "ref1", $watson ? 8 : 11-8+1, "read1", $watson ? "T" : "A", $watson, $fw, "I", -1, $cyc, 11, 20, 8);
		}
	}
	print STDERR "PASSED\n";
}

sub _test_bscpg_analyze_alignment_6() {
	print STDERR "Testing bscpg_analyze_alignment 6 ... ";
	require Data::Dumper;
	#  10987654  3210
	#  ATTAAAGC  TGGC
	# AATC   GCGACGTCG
	# 0123   456789012
	for(my $watson = 1; $watson < 2; $watson++) {
		for(my $fw = 0; $fw < 2; $fw++) {
			my @ev = _test_bscpg_analyze_alignment(
				$fw     ? "ATTAAAGCTGGC" : revcomp("ATTAAAGCTGGC"),  # orig read
				$watson ? "AATCGCGACGTCG" : revcomp("AATCGCGACGTCG"), # ref
				2,
				"A4^GA2T1",
				"3M3I2M2D4M",
				$watson,
				$fw,
				11,
				20,
				131);
			scalar(@ev) == 4 || die "Expected 4 records, got ".Dumper(\@ev);
			my $cyc = 9; $cyc = 12 - $cyc - 1 if $fw;
			_test_evidence($ev[0], "ref1", $watson ? 4 : 13-4+1, "read1", $watson ? "T" : "A", $watson, $fw, "I", -1, $cyc, 11, 20, 131);
			$cyc = 4; $cyc = 12 - $cyc - 1 if $fw;
			_test_evidence($ev[1], "ref1", $watson ? 6 : 13-6+1, "read1", $watson ? "C" : "G", $watson, $fw, "I", -1, $cyc, 11, 20, 131);
			$cyc = 3; $cyc = 12 - $cyc - 1 if $fw;
			_test_evidence($ev[2], "ref1", $watson ? 9 : 13-9+1, "read1", $watson ? "T" : "A", $watson, $fw, "I", -1, $cyc, 11, 20, 131);
			$cyc = 0; $cyc = 12 - $cyc - 1 if $fw;
			_test_evidence($ev[3], "ref1", $watson ? 12 : 13-12+1, "read1", $watson ? "C" : "G", $watson, $fw, "I", -1, $cyc, 11, 20, 131);
		}
	}
	print STDERR "PASSED\n";
}

sub _test() {
	_test_bs_analyze_alignment_1();
	_test_bs_analyze_alignment_2();
	_test_bs_analyze_alignment_3();
	_test_bs_analyze_alignment_4();
	_test_bs_analyze_alignment_5();
	_test_bs_analyze_alignment_6();
	_test_bs_analyze_alignment_7();
	
	_test_bsc_analyze_alignment_1();
	_test_bsc_analyze_alignment_2();
	_test_bsc_analyze_alignment_3();
	_test_bsc_analyze_alignment_4();
	_test_bsc_analyze_alignment_5();
	_test_bsc_analyze_alignment_6();
	_test_bsc_analyze_alignment_7();
	_test_bsc_analyze_alignment_8();
	
	_test_bscpg_analyze_alignment_1();
	_test_bscpg_analyze_alignment_2();
	_test_bscpg_analyze_alignment_3();
	_test_bscpg_analyze_alignment_4();
	_test_bscpg_analyze_alignment_5();
	_test_bscpg_analyze_alignment_6();
}

# Call _test() if the module is being run directly.
_test() unless caller();

1;
