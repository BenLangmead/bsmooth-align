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
# BtlBio::Index::Bowtie1::Index
#
# For building and checking the existence of Bowtie 1 indexes.
#
# Author: Ben Langmead
# Email: langmea@cs.jhu.edu
#

package BtlBio::Index::Bowtie1::Index;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(bowtie1_index_exists);

##
# Given index basename, check that all the index files are there.
#
sub bowtie1_index_exists($) {
	for my $ext ("1", "2", "3", "4", "rev.1", "rev.2") {
		my $fn = "$_[0].$ext.ebwt";
		if(! -f $fn) {
			print STDERR "Index at base '$_[0]' lacked file '$fn'\n";
			return 0;
		}
	}
	return 1;
}

1;
