# BSmooth
An alignment pipeline.

This is a fork of [BSmooth](https://github.com/BenLangmead/bsmooth-align) by
Kasper D. Hansen, Ben Langmead, and Rafael A. Irizarry. The fork fixes some
problems which prevented the merman aligner from compiling.

---

BSmooth is by Kasper D. Hansen, Ben Langmead, Rafael A. Irizarry

Statistics pipeline is by Kasper D. Hansen
Alignment pipeline is by Ben Langmead

The BSmooth alignment pipeline is licensed under the GPLv3 license.  See
`LICENSE' file for details.

## Installation

### BSmooth
 
 To install the BSmooth alignment pipeline, permanently set the BSMOOTH_HOME
 environment variable to the main BSmooth directory.  I.e., the directory with
 'bin', 'example' and 'lib' subdirectories inside.

### Building merman
 
 The merman sources are included with BSmooth.  To use merman with BSmooth, you
 must first build the merman binary.  This requires standard GNU build tools
 such as GNU make and g++.
 
 To build merman, cd to the $BSMOOTH_HOME/merman directory and run 'make'.
 Finally, add the $BSMOOTH_HOME/merman directory to your PATH.
 
### Installing Bowtie 2
 
 Follow the standard instructions for installing those tools; see their
 respective manuals for details.  Once the tool is installed, make sure that
 the directory containing the key binaries ('bowtie2', 'bowtie2-build') is
 included in your PATH.

## Pipeline flow

The following 4 steps (or possibly 3 steps, depending on aligner) are needed to
get from raw input reads to CpG-level measurement summaries.  These summaries
constitute the input to the R pipeline for smoothing and calculating DMRs.

  1. Build bisulfite genome index
   (not necessary if using Merman)
   bswc_bowtie2_index.pl

  2. Align
   bs_merman_align.pl
   bswc_bowtie2_align.pl

  3. Sort evidence directory
   bsev_sort.pl
  
  4. Tabulate sorted evidence directory
   bsev_tabulate.pl

For an example, see example/Makefile, which uses example/sim.pl to create a
small, simulated problem, then goes through the steps listed above.
