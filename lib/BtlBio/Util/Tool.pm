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
# BtlBio::Tool::Tool
#
# Routines for finding needed tools.
#
# Author: Ben Langmead
# Email: langmea@cs.jhu.edu
#

package BtlBio::Util::Tool;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(find_tool);

##
# Find a tool, perhaps given some information about (a) a path where it might
# be, and (b) some environment variables that might tell us where it is.
#
sub find_tool {
	my ($binname, $userspec, $conf) = @_;
	# Try user-specified path first
	return $userspec if (defined($userspec) && -x $userspec);
	# Try PATH next
	my $wh = `which $binname`;
	$wh =~ s/^\s+//;
	$wh =~ s/\s+$//;
	return $wh if (defined($wh) && -x $wh);
	return "./$binname" if -x "./$binname";
	die "Could not find tool '$binname'";
}

1;
