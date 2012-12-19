#!/usr/bin/env perl -w

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
# BtlBio::Sort::DiskSort
#
# An object that, given information about a tuple of interest, and given a
# collection of tuples, will sort the tuples according to a user-defined order.
# It will use the disk to accomplish this, first binning tuples into temporary
# files, then mergining temporary files into a final sorted file.
#
# Uses ForkManager for parallelism.
#
# Parallelism
# -----------
#
# There's a question as to whether the first step - where each record is
# assigned to a bin and appended to the appropriate file - should be
# parallelized.  This depends on the timing of how the records are being
# generated.  If all the records are sitting on the disk to begin with, this
# can be profitably parallelized.  If the records are being generated in
# real-time but slowly, there's no benefit to parallelizing this step.  For
# this reason, there are two function for providing records to the object:
# next_batch() and next_record().  next_batch() is used to provide large
# numbers of tuples to the DiskSorter at once, and they will be binned in
# parallel.  next_record() is used to provide records one-at-a-time, and it
# bins them serially.
#
# There is also tension between binning in parallel and supporting a notion of
# stabilized sort.  This class simply doesn't support stable sort.
#
# Binning
# -------
#
# Conceptually, there are two types of binning to do.  There's the type of
# binning we do simply to split the records up into smaller groups, which we
# call convenience binning.  Then there's the type of binning we do so that the
# output files are organized as the user requested, which we call final
# binning.
#
# We choose not to do any more convenience binning than the minimum amount
# dictated by the number of final bins and the number of threads.  We trust
# the sort subroutine to do merge-sorting as appropriate.
#
# Note that initially split the tuples into bins in a round-robin fashion, or
# something similar, is not a good option.  The second step (i.e. merge and
# split by bin) cannot be effectively parallelized because the merge must be
# done across all bins.
#
# Another option is to intially split tuples into bins based first on final
# bin, then on thread number (i.e. minimize convenience bins).  Then let the
# second step be paralell by the number of final bins.  This is the option used
# here.  Note that in order to get any parallelism in the second step, this
# option requires that the user specify final bins.
#

#
# Author: Ben Langmead
# Email: langmea@cs.jhu.edu
#

package BtlBio::Sort::DiskSort;

use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use Carp;
use FindBin qw($Bin);
use lib "$Bin/../..";
use BtlBio::Contrib::ForkManager;
use BtlBio::Util::File;
use BtlBio::Util::Temp;
use IO::File;
use List::Util qw[min max];
use File::Find;
use File::Path;
use File::Basename;

##
# Substitute out problematic characters.
#
sub _bin_subchar($) {
	$_[0] =~ s/[^a-zA-Z01-9\.-]/_/g;
}

##
# Default function for taking a list of field values and transforming them into
# a bin name.  By default, we convert all non-alphanumeric characters into
# underscores, and use underscores to join bin names.
#
sub _bin_renamer($) {
	my $binids = shift;
	defined($binids) || confess("_bin_renamer called with undef");
	my $name = join("_", @$binids);
	_bin_subchar($name);
	return $name;
}

sub _default_verbose($) { }

##
# Create a new DiskSort object.
#
sub new {
	my (
		$class,
		
		$sort_by,    # list of fields (0-based offset from left) to sort by
		$numeric,    # map indicating which fields are numeric
		$reverse,    # reverse order? (descending alphebetical or numeric)
		
		$inter,      # intermediate directory to use
		$keep_inter, # whether to keep intermediate directory
		$force,      # whether to delete directories if they already exist
		$odir,       # output directory
		
		$bin_by,     # list of fields (0-based offset from left) to bin by
		$bin_map,    # map from raw bin names to final bin names
		$bin_def,    # default bin name (for when bin isn't in map)
		
		$bin_ren,    # function that renames bin tuples
		
		$delim,      # delimiter
		$sort_sz,    # max memory footprint of any given sort process
		
		$threads,    # concurrent threads/processes to use for sort
		
		$prefix,     # string to prepend to bin filenames
		$suffix,     # string to append to bin filenames
		
		$compress,   # type of compression to use on output
		
		$min_field,  # minimum # of fields in valid record
		$max_field,  # maximum # of fields in valid record
		$allow_emp,  # empty input records OK, or cause error?
		$logfunc     # function to send log messages to
		
	) = @_;
	$bin_ren  = \&_bin_renamer     unless defined($bin_ren);
	$logfunc  = \&_default_verbose unless defined($logfunc);
	my $odir_is_tmp = !defined($odir);
	my $inter_is_tmp = !defined($inter);
	$odir     = temp_name()        unless defined($odir);
	$inter    = "$odir.pre"        unless defined($inter);
	my $binodir  = "$inter/bins";
	my $tempdir  = "$inter/temp";
	my $binedir  = "$inter/bins.err";
	if(defined($bin_map)) {
		# Make a new bin map where all the keys are already renamed
		my %new_map = ();
		for my $k (%$bin_map) {
			my $v = $bin_map->{$k};
			_bin_subchar($k);
			$new_map{$k} = $v;
		}
		$bin_map = \%new_map;
	}
	return bless {
		_sort_by      => $sort_by,
		_numeric      => $numeric,
		_reverse      => $reverse,
		
		_inter        => $inter,
		_inter_is_tmp => $inter_is_tmp,
		_keep_inter   => $keep_inter,
		_force        => $force,
		_odir         => $odir,
		_odir_is_tmp  => $odir_is_tmp,
		
		_binodir      => $binodir,
		_tempdir      => $tempdir,
		_binedir      => $binedir,
		
		_bin_by       => $bin_by,
		_bin_map      => $bin_map,
		_bin_def      => $bin_def,
		
		_bin_ren      => $bin_ren,
		
		_delim        => $delim,
		_sort_sz      => $sort_sz,
		
		_threads      => $threads,
		
		_prefix       => $prefix,
		_suffix       => $suffix,
		
		_compress     => $compress,
		
		_min_field    => $min_field,
		_max_field    => $max_field,
		_allow_emp    => $allow_emp,
		_logfunc      => $logfunc,
		
		_memo         => {},
		_first_rec    => 0,           # saw first record yet?
		_finalized    => 0            # whether finalize() was called
	}, $class;
}
sub _log($$) { return $_[0]->{_logfunc}->($_[1]); }

sub output_dir($)       { return $_[0]->{_odir};  }
sub intermediate_dir($) { return $_[0]->{_inter}; }

##
# Given values for the _sort_by, _numeric, and _reverse fields, generate the
# command-line option string to be passed to the 'sort' utility.
#
sub _sort_options($$) {
	my ($self, $sz) = @_;
	defined($sz) || confess();
	my $opts = "sort ";
	my $first = 1;
	$opts .= "-r " if $self->{_reverse};
	scalar(@{$self->{_sort_by}}) > 0 || confess();
	for my $fid (reverse @{$self->{_sort_by}}) {
		$opts .= "| sort -s " unless $first;
		$opts .= "-T $self->{_tempdir} ";
		$opts .= "-S $sz ";
		if($self->{_numeric}{$fid}) {
			$opts .= "-n ";
		}
		$opts .= "-k$fid,$fid ";
		$first = 0;
	}
	return $opts;
}

##
# Given an array ref of record fields, assign it to a bin.
#
sub _binid($$) {
	my ($self, $fs) = @_;
	my @bin_by_list = ();
	# Assembly the sort-by tuple
	for my $bini (@{$self->{_bin_by}}) {
		int($bini) == $bini ||
			croak("Bad bin identifier '$bini'; must be integer >= 1");
		$bini > 0 ||
			croak("Bad bin identifier '$bini'; must be integer >= 1");
		my $i = $bini - 1;
		defined($fs->[$i]) || confess("Record had no field with offset $i");
		push @bin_by_list, $fs->[$i];
	}
	# Call bin renamer
	my $orig_name = join("\t", @bin_by_list);
	my $name = $self->{_bin_ren}->(\@bin_by_list);
	defined($name) || die;
	# See if we should map that name to something else
	if(defined($self->{_bin_map})) {
		my $newname = $self->{_bin_map}{$name};
		if(!defined($newname) && !defined($self->{_bin_def})) {
			print STDERR "Error: No mapping for bin with name '$name'.  Mappings exist for:\n";
			for my $f (keys %{$self->{_bin_map}}) { print STDERR "  $f\n"; }
			die;
		}
		$newname = $self->{_bin_def} unless defined($newname);
		$name = $newname;
		_bin_subchar($name);
	}
	return ($name, $orig_name);
}

##
# Make sure directory exists, possibly removing it first if force is
# authorized.
#
sub check_dir($$) {
	my ($self, $dir) = @_;
	if(-d $dir) {
		croak("Output directory $dir already exists") unless $self->{_force};
		if($self->{_force}) {
			$self->_log("Removing directory $dir, since force is enabled\n");
			rmtree($dir);
			-d $dir && croak("Could not remove directory $dir");
		}
	}
	mkpath($dir);
	-d $dir || croak("Could not create new directory $dir");
}

##
#
#
sub initialize($) {
	my $self = shift;
	$self->setup_dirs();
	$self->{_first_rec} = 1;
}

##
# Set up all the temporary directories that we'll use for sorting.
#
sub setup_dirs($) {
	my $self = shift;
	$self->check_dir($self->{_odir});
	$self->check_dir($self->{_inter});
	$self->check_dir($self->{_binodir});
	$self->check_dir($self->{_tempdir});
	$self->check_dir($self->{_binedir});
}

##
# Process a whole batch of records in one file already on disk.  Use given
# state ($tid/$memory).
#
sub next_batch_file_ex($$$$) {
	my ($self, $batchfile, $tid, $memo) = @_;
	my $ifh = openex($batchfile);
	my $n = 0;
	while(readline $ifh) { $self->next_record_ex($_, $tid, $memo); $n++; }
	close($ifh);
	# Close output handles
	for (keys %{$memo->{fhs}}) { close($memo->{fhs}{$_}); }
}

##
# Process a whole batch of records in one file already on disk.  Use given
# state ($tid/$memory).
#
sub next_batch_file($$) {
	my ($self, $batchfile) = @_;
	return $self->next_batch_file_ex(
		$batchfile,       # file
		0,                # thread id (0 = master thread)
		$self->{_memo});  # where filehandles, other state is kept
}

##
# Process a whole batch of records in one or more files already on disk.  This
# can be done in parallel.  Each parallel thread must have its own filehandle
# map, and its own memory of the previous
#
sub next_batch_dir($$) {
	my ($self, $batchdir) = @_;
	if(!$self->{_first_rec}) {
		$self->initialize();
	}
	$self->_log("Starting fork manager with $self->{_threads} threads\n");
	my $pm = new Parallel::ForkManager($self->{_threads});
	# Set callback for when a child finishes up so we can get its exit code
	my $childFailed = 0;
	my $childFailedPid = 0;
	$pm->run_on_finish(sub {
		my ($pid, $exit_code, $ident) = @_;
		if($self->{_threads} > 0 && $exit_code != 0) {
			$childFailed = $exit_code;
			$childFailedPid = $pid;
		}
	});
	our @fns = ();
	sub add_fn {
		my $fn = fileparse($File::Find::name);
		if(substr($fn, 0, 1) ne ".") {
			push @fns, $File::Find::name;
			#$self->_log("  Input file: $File::Find::name\n");
		}
	}
	find(\&add_fn, $batchdir);
	# Sort filename by size so that we're doing large files before small ones
	@fns = sort { -s $b <=> -s $a } @fns;
	# First, determine the number of input files
	my $ninputs = scalar(@fns);
	$self->_log("Found $ninputs input files in directory $batchdir\n");
	# For each input dir
	my %filesDone = ();
	$self->_log("--- Batch bin ---\n");
	my $fi = 0;
	for my $fn (@fns) {
		$fi++;
		my $base = fileparse($fn);
		if($childFailed) {
			$self->_log("Aborting master loop because child failed\n");
			last;
		}
		$pm->start and next; # fork off a mapper for this input file
		# Keep memory of: previous binid, previous bin, and all filehandles
		my %memo = ();
		$self->next_batch_file_ex($fn, $fi, \%memo);
		$pm->finish; # end of fork
	}
	$self->_log("Aborted master loop because child failed\n") if $childFailed;
	$pm->wait_all_children;
	if($childFailed) {
		croak("Child with PID $childFailedPid exited abnormally\n");
	} else { $self->_log("All children succeeded\n"); }
	$self->_log("--- Finished batch bin ---\n");
}

##
# Feed another record into the DiskSorter, with some state provided.
#
sub next_record_ex($$$$) {
	my ($self, $rec, $tid, $memo) = @_;
	defined($memo) || die;
	chomp($rec);
	if(!$self->{_allow_emp} && $rec =~ /^\s*$/) {
		croak("next_record called with empty record!");
	}
	my @fs = split($self->{_delim}, $rec, -1);
	my $nfs = scalar(@fs);
	my ($min_field, $max_field) = ($self->{_min_field}, $self->{_max_field});
	if($nfs < $min_field || $nfs > $max_field) {
		croak("Expected # fields in [$min_field, $max_field]; got:\n$rec\n");
	}
	# Put it in a bin
	my ($binid_xformed, $binid) = $self->_binid(\@fs);
	defined($binid_xformed) || die;
	if(!$self->{_first_rec}) {
		$self->initialize();
	}
	if(!defined($memo->{outfhs}{$binid_xformed})) {
		# Open a new file for a new bin
		# TODO: avoid name collisions somehow - like by adding strings of
		# numbers on the end
		my $odir = "$self->{_binodir}/$binid_xformed";
		mkpath($odir);
		my $ofn = ">>$odir/$tid"; # has to be append!
		my $ofn_name = ">$odir/.$tid.name"; # has to be append!
		open($memo->{outfhs}{$binid_xformed}, $ofn) ||
			die "Could not open $ofn for writing\n";
		open (NM, ">$ofn_name") || die "Could not open $ofn_name for writing\n";
		print NM "$binid\n";
		close(NM);
		$self->_log("Opened '$ofn'; ".scalar(keys %{$memo->{outfhs}})." files open in PID $$\n");
	}
	print {$memo->{outfhs}{$binid_xformed}} "$rec\n";
}

##
# Feed another record into the DiskSorter.
#
sub next_record($$) {
	my ($self, $rec) = @_;
	$self->next_record_ex(
		$rec,             # record
		0,                # thread id (0 = master thread)
		$self->{_memo});  # where filehandles, other state is kept
}

##
# No more records coming - we can sort the bins and emit the final, sorted
# output.
#
sub finalize($) {
	my ($self) = @_;
	# Figure out: what were the final bins we saw records for?
	my %binvals = ();
	my $bindir = $self->{_binodir};
	my $tmpdir = $self->{_tempdir};
	my $errdir = $self->{_binedir};
	for my $fn (<$bindir/*>) {
		if(-d $fn) {
			# Found a directory for a bin
			$fn = fileparse($fn); # just get basename
			my @name_fns = <$bindir/$fn/.*.name>;
			scalar(@name_fns) > 0 || die "No name files in '$bindir/$fn'";
			my $long = undef;
			open(NM, "$name_fns[0]") ||
				die "Could not open name file '$name_fns[0]' for reading";
			$long = <NM>;
			chomp($long);
			close(NM);
			$self->_log("  Found bin $fn, with long name '$long'\n");
			$binvals{$long} = $fn;
		}
	}
	my $nbin = scalar(keys %binvals);
	if(scalar(keys %binvals) < 1) {
		$self->_log("No non-empty bins!  Returning...\n");
		return;
	}
	$self->_log("Found $nbin non-empty bins...\n");
	my $nth = $self->{_threads};
	$nth = 1 if $nth < 1;
	my $sortSize = int((1 * 1024 * 1024)/min($nth, $nbin));
	my $bi = 0;
	my $sortCmd = $self->_sort_options($sortSize);
	$self->_log("--- Sort ---\n");
	$self->_log("Sort command base: $sortCmd\n");
	my ($prefix, $suffix) = ($self->{_prefix}, $self->{_suffix});
	$suffix .= ".gz"  if $self->{_compress} eq "gzip";
	$suffix .= ".bz2" if $self->{_compress} eq "bzip2";
	my $childFaild = 0;
	my $pm = new Parallel::ForkManager($self->{_threads});
	# Set callback for when a child finishes up so we can get its exit code
	my $childFailed = 0;
	my $childFailedPid = 0;
	$pm->run_on_finish(sub {
		my ($pid, $exit_code, $ident) = @_;
		if($self->{_threads} > 0 && $exit_code != 0) {
			$childFailed = $exit_code;
			$childFailedPid = $pid;
		}
	});
	my %bin_to_sz = ();
	for my $k (keys %binvals) {
		my $v = $binvals{$k};
		my $sz = 0;
		for my $fn (<$bindir/$v/*>) { $sz += -s $fn; }
		$bin_to_sz{$k} = $sz;
	}
	for my $binval (sort { $bin_to_sz{$b} <=> $bin_to_sz{$a} } keys %binvals) {
		my $short = $binvals{$binval}; # Get the short name
		$bi++;
		if($childFailed) {
			$self->_log("Aborting master loop because child failed\n");
			last;
		}
		$pm->start and next; # fork off a mapper for this input file
		$self->_log("Pid $$ processing bin $binval [$bi of $nbin]...\n");
		my $outfn = $self->{_odir}."/$prefix$short$suffix";
		my $out_long_fn = $self->{_odir}."/.$prefix$short$suffix.full_name";
		open(LN, ">$out_long_fn") || die;
		print LN "$binval\n";
		close(LN);
		my $outpipe = "| cat > $outfn";
		$outpipe    = "| gzip  -c > $outfn" if $self->{_compress} eq "gzip";
		$outpipe    = "| bzip2 -c > $outfn" if $self->{_compress} eq "bzip2";
		my $sort_cmd = "cat $bindir/$short/* | $sortCmd 2>$errdir/$short $outpipe";
		$self->_log("Sort command: $sort_cmd\n");
		my $ret = system($sort_cmd);
		$ret == 0 || die "Non-0 exitlevel from sort command: $?";
		$pm->finish; # end of fork
	}
	$pm->wait_all_children;
	if($childFailed) {
		croak("Sort child with PID $childFailedPid exited abnormally\n");
	} else { $self->_log("All children succeeded\n"); }
	# Delete all the files that were inputs to the sort
	if(!$self->{_keep_inter}) {
		$self->_log("Deleting intermediate and temporary directories:\n");
		$self->_log("  Deleting temporary bin directories...\n");
		rmtree($bindir) unless $self->{_keep_inter};
		$self->_log("  Deleting error message directories...\n");
		rmtree($errdir) unless $self->{_keep_inter};
		$self->_log("  Deleting intermediate directories...\n");
		rmtree($self->{_inter}) unless $self->{_keep_inter};
		$self->_log("  Deleting temporary sort directories...\n");
		rmtree($tmpdir) unless $self->{_keep_inter};
	}
	$self->{_finalized} = 1;
}

##
# Reset the DiskSort object in preparation for a new stream of records.
#
sub reset($) {
	my $self = shift;
	if($self->{_odir_is_tmp}) {
		$self->{_odir} = temp_name();
	}
	if($self->{_inter_is_tmp}) {
		$self->{_inter} = $self->{_odir}.".pre";
	}
	$self->{_binodir}   = $self->{_inter}."/bins";
	$self->{_tempdir}   = $self->{_inter}."/temp";
	$self->{_binedir}   = $self->{_inter}."/bins.err";
	$self->{_memo}      = {};
	$self->{_first_rec} = 0;
	$self->{_finalized} = 0;
}

##
# Remove any temporary files and directories that might remain.
#
sub clean_up($) {
	my $self = shift;
	rmtree($self->{_binodir});
	rmtree($self->{_tempdir});
	rmtree($self->{_binedir});
}

sub _test_1() {
	print STDERR "Testing DiskSort 1 ... ";
	my $dsort = BtlBio::Sort::DiskSort->new(
		[1, 2, 3],  # Sort on first 3 fields
		{ 1 => 1 }, # First field is numeric
		0,          # Don't reverse (ascending, not descending)
		undef,      # Intermediate dir (use temp)
		0,          # Don't keep intermedaite dir
		0,          # Don't force
		undef,      # output dir (use temp)
		[1],        # Bin on just first field
		undef,      # No bin map - use first field directly
		"default",  # Bin default name
		undef,      # No bin renaming function
		"\t",       # Normal delimiter
		1,          # Small sort memory ceiling
		0,          # Don't use threads
		"prefix",   # Dummy prefix
		"suffix",   # Dummy suffix
		"gzip",     # Gzip-compress output files
		3,          # All tuples must have exactly 3 fields
		3,          # All tuples must have exactly 3 fields
		0,          # Empty tuples not allowed
		0 ? sub {print STDERR $_[0]; } : undef); # Simple logging function
	$dsort->{_finalized} == 0       || die "Expected not finalized";
	$dsort->{_delim}    eq "\t"     || die "Got: $dsort->{_delim}";
	$dsort->{_prefix}   eq "prefix" || die "Got: $dsort->{_prefix}";
	$dsort->{_suffix}   eq "suffix" || die "Got: $dsort->{_suffix}";
	$dsort->{_compress} eq "gzip"   || die "Got: $dsort->{_compress}";
	$dsort->{_min_field} == 3       || die "Got: $dsort->{_min_field}";
	$dsort->{_max_field} == 3       || die "Got: $dsort->{_max_field}";
	my $res = $dsort->_sort_options(1000);
	$res =~ s/\s+/ /g; $res =~ s/^\s*//; $res =~ s/\s*$//;
	$res =~ s/ -T [^ ]*//g;
	my $ex = "sort -S 1000 -k3,3 | sort -s -S 1000 -k2,2 | sort -s -S 1000 -n -k1,1";
	$res eq $ex || die "Got\n$res\nexpected\n$ex";
	my ($short, $long) = $dsort->_binid([1, "cat", "dog"]);
	$short eq "1" || die "Got $short";
	# Give the DiskSort a series of tuples via next_record method
	$dsort->next_record(join("\t", (1, "cat", "dog4"))."\n");
	$dsort->next_record(join("\t", (1, "cat", "dog1"))."\n");
	$dsort->next_record(join("\t", (1, "cat", "dog2"))."\n");
	$dsort->next_record(join("\t", (1, "cat", "dog3"))."\n");
	$dsort->next_record(join("\t", (2, "cheetah", "llama"))."\n");
	$dsort->next_record(join("\t", (2, "llama", "cheetah"))."\n");
	$dsort->next_record(join("\t", (2, "cheetah", "llama"))."\n");
	my $ofn = temp_name();
	open(OF, ">$ofn") || die "Could not open $ofn for writing";
	print OF join("\t", (3, "alpaca", "gorilla"))."\n";
	print OF join("\t", (3, "gorilla", "alpaca"))."\n";
	print OF join("\t", (6, "lemur01", "lemur77"))."\n";
	print OF join("\t", (20, "lemur1", "lemur007"))."\n";
	close(OF);
	$dsort->next_batch_file($ofn);
	
	$dsort->finalize();
	$dsort->{_finalized} == 1       || die "Expected finalized";
	
	$dsort->clean_up();
	my $odir = $dsort->output_dir();
	$dsort->reset();
	
	-f "$odir/prefix1suffix.gz" || die;
	-f "$odir/prefix2suffix.gz" || die;
	-f "$odir/prefix3suffix.gz" || die;
	-f "$odir/prefix6suffix.gz" || die;
	-f "$odir/prefix20suffix.gz" || die;
	
	for my $bin (1, 2, 3, 6, 20) {
		my @lines = ();
		open(IN, "gzip -dc $odir/prefix${bin}suffix.gz |");
		while(<IN>) { chomp; push @lines, $_; }
		close(IN);
		my @s;
		if($bin == 1) {
			scalar(@lines) == 4 || die;
			@s = split("\t", $lines[0]);
			$s[0] == 1 || die;
			$s[1] eq "cat" || die;
			$s[2] eq "dog1" || die;
		} elsif($bin == 2) {
			scalar(@lines) == 3 || die;
			@s = split("\t", $lines[0]);
			$s[0] == 2 || die;
			$s[1] eq "cheetah" || die;
			$s[2] eq "llama" || die;
		} elsif($bin == 3) {
			scalar(@lines) == 2 || die;
		} elsif($bin == 6) {
			scalar(@lines) == 1 || die;
		} elsif($bin == 20) {
			scalar(@lines) == 1 || die;
		}
	}
	
	print STDERR "PASSED\n";
}

sub _test_2() {
	print STDERR "Testing DiskSort 2 ... ";
	my %bin_map = (
		"Agua_1"   => "AOne",
		"Water_30" => "WThirty",
		"Eau_5"    => "EFive",
	);
	rmtree("/tmp/DiskSort.pm.intermediate");
	rmtree("/tmp/DiskSort.pm.output");
	my $dsort = BtlBio::Sort::DiskSort->new(
		[1, 2],     # Sort on first 2 fields
		{ 2 => 1 }, # Second field is numeric
		1,          # Reverse (descending)
		"/tmp/DiskSort.pm.intermediate", # Intermediate dir
		1,          # Keep intermedaite dir
		0,          # Don't force
		"/tmp/DiskSort.pm.output", # output dir (use temp)
		[1, 2],     # Bin on first two fields
		\%bin_map,  # No bin map - use first field directly
		"default",  # Bin default name
		undef,      # No bin renaming function
		",",        # Normal delimiter
		1,          # Small sort memory ceiling
		2,          # Don't use threads
		"Chr",      # Dummy prefix
		".tab",     # Dummy suffix
		"",         # No compression
		2,          # All tuples must have at least 2 fields
		10,         # All tuples must have <= 10 fields
		0,          # Empty tuples not allowed
		0 ? sub {print STDERR $_[0]; } : undef); # Simple logging function
	$dsort->{_finalized} == 0       || die "Expected not finalized";
	$dsort->{_delim}    eq ","      || die "Got: $dsort->{_delim}";
	$dsort->{_prefix}   eq "Chr"    || die "Got: $dsort->{_prefix}";
	$dsort->{_suffix}   eq ".tab"   || die "Got: $dsort->{_suffix}";
	$dsort->{_compress} eq ""       || die "Got: $dsort->{_compress}";
	$dsort->{_min_field} == 2       || die "Got: $dsort->{_min_field}";
	$dsort->{_max_field} == 10      || die "Got: $dsort->{_max_field}";
	my $ex = "sort -r -S 1000 -n -k2,2 | sort -s -S 1000 -k1,1";
	my $res = $dsort->_sort_options(1000);
	$res =~ s/\s+/ /g; $res =~ s/^\s*//; $res =~ s/\s*$//;
	$res =~ s/ -T [^ ]*//g;
	$res eq $ex || die "Got\n$res\nexpected\n$ex";
	my ($short, $long) = $dsort->_binid([1, "cat", "dog"]);
	$short eq "default" || die "Got $short";
	($short, $long) = $dsort->_binid(["Water", "30", "abracadabra"]);
	$short eq "WThirty" || die "Got $short";
	# Give the DiskSort a series of tuples via next_record method
	srand(3456);
	my @lines = ();
	my $ex_nlines = 0;
	for my $i1 ("Agua", "Water", "Eau", "alpha", "bravo", "Bravo", "charile") {
		for my $i2 ("001", "002", "01", "02", "1", "2", "3", "4",
		            "10", "20", "30", "40", "100", "0100")
		{
			my @toks = ();
			push @toks, $i1;
			push @toks, $i2;
			my $n = int(rand(8));
			for(0..$n) {
				push @toks, int(rand(900000));
			}
			push @lines, join(",", @toks)."\n";
			$ex_nlines++;
		}
	}
	my $dir = temp_name();
	mkpath($dir);
	for(my $i = 0; $i < scalar(@lines); $i += 10) {
		open(OUT, ">$dir/part$i") || die;
		for(my $j = 0; $j < 10; $j++) {
			print OUT $lines[$i+$j] if $i+$j < scalar(@lines);
		}
		close(OUT);
	}
	$dsort->next_batch_dir($dir);
	
	$dsort->finalize();
	$dsort->{_finalized} == 1 || die "Expected finalized";
	
	$dsort->clean_up();
	my $odir = $dsort->output_dir();
	$dsort->reset();
	
	-f "$odir/ChrAOne.tab" || die;
	-f "$odir/ChrWThirty.tab" || die;
	-f "$odir/Chrdefault.tab" || die;

	my $nlines = 0;
	for my $bin ("AOne", "WThirty", "default") {
		my @lines = ();
		open(IN, "$odir/Chr$bin.tab") || die;
		while(<IN>) { chomp; push @lines, $_; $nlines++; }
		close(IN);
	}
	$nlines eq $ex_nlines || die "Expected $ex_nlines lines, got $nlines";
	
	print STDERR "PASSED\n";
}

sub _test() {
	_test_1();
	_test_2();
}

# Call _test() if the module is being run directly.
_test() unless caller();

1;
