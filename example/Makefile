#!gmake

##
# Makefile demonstrating two runs of the pipeine: one using the Merman aligner
# and one using the Bowtie 2 aligner + in-silico bisulfite conversion.
#
# Ultimate output of the merman pipeline is in ev_merman_tab and ultimate
# output of the Bowtie 2 pipeline is in ev_bt2_tab.
#

all: ev_bt2_tab ev_merman_tab

ev_bt2_tab: ev_bt2_sorted
	perl $(BSMOOTH_HOME)/bin/bsev_tabulate.pl \
		--cpg-output-dir=$@ \
		--cpg-stranded-output-dir=$(@)_stranded \
		-- \
		$^ \
		-- \
		*.fa

ev_merman_tab: ev_merman_sorted
	perl $(BSMOOTH_HOME)/bin/bsev_tabulate.pl \
		--cpg-output-dir=$@ \
		--cpg-stranded-output-dir=$(@)_stranded \
		-- \
		$^ \
		-- \
		*.fa

ev_bt2_sorted: ev_bt2_raw
	perl $(BSMOOTH_HOME)/bin/bsev_sort.pl \
		--evidence-directory=$< \
		--output-dir=$@

ev_merman_sorted: ev_merman_raw
	perl $(BSMOOTH_HOME)/bin/bsev_sort.pl \
		--evidence-directory=$< \
		--output-dir=$@

ev_bt2_raw: sim.watson.bt2
	perl $(BSMOOTH_HOME)/bin/bswc_bowtie2_align.pl \
		--four-strand \
		--fastq \
		--fr \
		--evidence-directory=$@ \
		--bscpg \
		-- \
		sim \
		-- \
		*.fa \
		-- \
		--sensitive \
		-- \
		reads_1.fq \
		-- \
		reads_2.fq

ev_merman_raw: ref_1.fa
	perl $(BSMOOTH_HOME)/bin/bs_merman_align.pl \
		--four-strand \
		--fastq \
		--evidence-directory=$@ \
		--bscpg \
		-- \
		*.fa \
		-- \
		-v 6 -M 1 \
		-- \
		reads_1.fq

ref_1.fa:
	perl sim.pl

sim.watson.bt2: ref_1.fa
	perl $(BSMOOTH_HOME)/bin/bswc_bowtie2_index.pl \
		--name=sim \
		--output=. \
		*.fa

clean:
	rm -f *.fa *.fq *.bt2
	rm -rf ev_bt2_raw ev_bt2_sorted* ev_bt2_tab*
	rm -rf ev_merman_raw ev_merman_sorted* ev_merman_tab*
