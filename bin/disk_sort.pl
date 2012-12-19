#!/usr/bin/env perl

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
# disk_sort.pl
#
# Sort input according to a collection of sort keys and bins provided by the
# user.  Ues the DiskSort module, which uses files on disk to accomplish the
# sorting.  It first bins tuples into temporary files, then sorts and merges
# the temporary files into a directory of final sorted files, one per bin.
#

# Special case: Input 

use strict;
use warnings;
use Carp;
use Getopt::Long;
use File::Path;
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use BtlBio::Sort::DiskSort;

my $usage = qq!
Usage:
  disk_sort.pl --name=<string> [options*] <fasta-files>
!;

my $temp_dir = $ENV{TMPDIR} || "/tmp";
my $idir = undef;      # intermediate directory
my $keep_inter = 0;    # whether to keep intermediate directory
my $odir = undef;      # output directory
my $delim = "\t";      # delimited separating fields in a record
my @sort_keys = ();    # list of fields (0-based offset from left) to sort by
my @bin_keys = ();     # list of fields (0-based offset from left) to bin by
my %numerics = ();     # map indicating which fields are numeric
my $reverse = 0;       # reverse order? (descending alphebetical or numeric)
my $force = 0;         # delete directories if they already exist?
my %bin_map = ();      # map from bin names to bin names
my $sort_sz = 2;       # max memory footprint (GB) of any given sort process
my $compress = "";     # type of compression to use on output
my $prefix = "";       # prepend this to each output filename
my $suffix = "";       # append this to each output filename
my $threads = 1;       # concurrent threads/processes to use for bin & sort
my $min_field = 0;     # minimum # of fields in valid record
my $max_field = 99999; # maximum # of fields in valid record
my $allow_emp = 0;     # empty input records OK, or cause error?
my @input_dirs = ();   # directories with input tuples
my @name_map_fn = ();  # map from bin names to bin names

sub msg { print STDERR $_[0]; }

GetOptions (
	"delim=s"           => \$delim,
	"sort-key=s"        => \@sort_keys,
	"bin-key=s"         => \@bin_keys,
	"numeric=s"         => sub { $numerics{$_[1]} = 1; },
	"reverse"           => \$reverse,
	"force"             => \$force,
	"temp=s"            => \$temp_dir,
	"num-threads=i"     => \$threads,
	"keep-intermediate" => \$keep_inter,
	"intermediate=s"    => \$idir,
	"prefix=s"          => \$prefix,
	"suffix=s"          => \$suffix,
	"input=s"           => \@input_dirs,
	"output=s"          => \$odir,
	"name-map=s"        => \@name_map_fn,
	"compress=s"        => \$compress,
	"gzip"              => sub { $compress = "gzip"; },
	"bzip2"             => sub { $compress = "bzip2"; }
) || dieusage("Bad option", 0);

my $bin_map_arg = undef;
my $renamer = undef;
if(scalar(@name_map_fn) > 0) {
	for my $fn (@name_map_fn) {
		open(NM, $fn) || die "Could not open name-map file '$fn' for reading";
		while(<NM>) {
			chomp;
			next if /^\s*$/;
			next if /^\s*#/;
			my @s = split(/\t/);
			$bin_map{$s[0]} = $s[1];
		}
		close(NM);
	}
	$bin_map_arg = \%bin_map;
}

print STDERR "Delimiter: ".ord($delim)."\n";
print STDERR "Sort keys:\n  ".join("\n  ", @sort_keys)."\n";
print STDERR "Bin keys:\n  ".join("\n  ", @sort_keys)."\n";
print STDERR "Numerics:\n  ".join("\n  ", sort {$a <=> $b} keys %numerics)."\n";
print STDERR "Temp dir: $temp_dir\n";
print STDERR "Sort mem sz: ${sort_sz}GB\n";
print STDERR "Compression type (if any): $compress\n";
print STDERR "Input dirs:\n";
for(0..$#input_dirs) { print STDERR "  $input_dirs[$_]\n"; }
print STDERR "Output dir: $odir\n";
print STDERR "Output prefix: $prefix\n";
print STDERR "Output suffix: $suffix\n";
print STDERR "Name maps:\n";
for(0..$#name_map_fn) { print STDERR "  $name_map_fn[$_]\n"; }
print STDERR "Min/max # fields: [$min_field, $max_field]\n";
print STDERR "# threads: $threads\n";
print STDERR "Allow empty records: $allow_emp\n";

my $dsort = BtlBio::Sort::DiskSort->new(
	\@sort_keys, # list of fields (0-based offset from left) to sort by
	\%numerics,  # map indicating which fields are numeric
	$reverse,    # reverse order? (descending alphebetical or numeric)
	$idir,       # intermediate directory to use
	$keep_inter, # whether to keep intermediate directory
	$force,      # whether to delete directories if they already exist
	$odir,       # output directory
	\@bin_keys,  # # list of fields (0-based offset from left) to bin by
	$bin_map_arg,# map from raw bin names to final bin names
	undef,       # default bin name (for when bin isn't in map)
	$renamer,    # function that renames bin tuples
	$delim,      # delimiter
	$sort_sz,    # max memory footprint of any given sort process
	$threads,    # concurrent threads/processes to use for sort
	$prefix,     # string to prepend to bin filenames
	$suffix,     # string to append to bin filenames
	$compress,   # type of compression to use on output
	$min_field,  # minimum # of fields in valid record
	$max_field,  # maximum # of fields in valid record
	$allow_emp,  # empty input records OK, or cause error?
	\&msg);      # function to send log messages to

$threads = 0;
my $pm = new Parallel::ForkManager($threads);
# Set callback for when a child finishes up so we can get its exit code
my $childFailed = 0;
my $childFailedPid = 0;
$pm->run_on_finish(sub {
	my ($pid, $exit_code, $ident) = @_;
	if($threads > 0 && $exit_code != 0) {
		$childFailed = $exit_code;
		$childFailedPid = $pid;
	}
});

msg("Initializing DiskSort object...\n");
$dsort->initialize();

if(scalar(@input_dirs) > 0) {
	# Feed it a collection of files
	for my $dir (@input_dirs) { $dsort->next_batch_dir($dir); }
} else {
	# Feed it one record at a time
	my @buf = ();
	my $nchunk = 32 * 1024;
	my $tid = 0;
	my $i = 0;
	msg("Reading input records...\n");
	while(<>) {
		next if /^\s*$/;
		chomp;
		$buf[$i] = $_;
		++$i <= $nchunk || die;
		if($i == $nchunk) {
			$tid++ if $threads > 0; # TODO: good way to make this multithreaded?
			$i = 0;
			if($pm->start) {
				next;
			}
			my %memo = ();
			for my $b (@buf) { $dsort->next_record_ex($b, $tid, \%memo); }
			$pm->finish; # end of fork
		}
	}
	$tid++;
	my %memo = ();
	for(my $j = 0; $j < $i; $j++) {
		$dsort->next_record_ex($buf[$j], $tid, \%memo);
	}
	
	$pm->wait_all_children;
	if($childFailed) {
		croak("Child with PID $childFailedPid exited abnormally\n");
	} else { msg("All children succeeded\n"); }
}

# Now call finalize
msg("Finalize...\n");
$dsort->finalize();

msg("DONE\n");
