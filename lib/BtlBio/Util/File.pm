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
# BtlBio::Util::File
#
# Basic file manipulation, including automatic setup of pipes through
# appropriate decompression tools for compressed files.
#
# Author: Ben Langmead
# Email: langmea@cs.jhu.edu
#

package BtlBio::Util::File;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(openex sanitize_filename);

use strict;
use warnings;
use Carp;

##
# Return a filehandle for a file.  If the file is compressed, the filehandle
# should correspond to a pipeline that reads the file using, for instance,
# 'gzip -dc'.
#
sub openex($) {
	my $fn = shift;
	$fn =~ s/^~/$ENV{HOME}/;
	my $pipe = $fn;
	if($fn =~ /\.gz$/) {
		$pipe = "gzip -dc $fn |";
	} elsif($fn =~ /\.bz2$/) {
		$pipe = "bzip2 -dc $fn |";
	}
	my $fh;
	open($fh, $pipe) || croak("Could not open '$pipe' for reading");
	defined($fh) || croak("Filehandle for '$pipe' should be defined");
	return $fh;
}

##
# Turn a potential name with wacky characters into a filename with only
# acceptable characters.
#
sub sanitize_filename($) {
	my $fn = shift;
	$fn =~ s/[^0-9a-zA-Z()_-]/_/g;
	return $fn;
}

1;
