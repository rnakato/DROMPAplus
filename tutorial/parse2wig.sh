#!/bin/bash


# Get a sample BAM file
wget -nc http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwHistone/wgEncodeUwHistoneK562H3k4me3StdAlnRep1.bam

### Basic usage
# the simplest command
gt=../data/genometable/genometable.hg19.txt
parse2wig+ -i wgEncodeUwHistoneK562H3k4me3StdAlnRep1.bam -o H3K4me3 --gt $gt

# parse2wig+ allows multiple input files (separated by “,”):
parse2wig+ -i wgEncodeUwHistoneK562H3k4me3StdAlnRep1.bam,wgEncodeUwHistoneK562H3k4me3StdAlnRep1.bam -o H3K4me3 --gt $gt

# Output as bedGraph format
parse2wig+ -i wgEncodeUwHistoneK562H3k4me3StdAlnRep1.bam -o H3K4me3 --gt $gt --outputformat 2

# By default, parse2wig+ omits to output bins in which the value is zero to reduce the file size.
# To output all bins, add --outputzero:
parse2wig+ -i wgEncodeUwHistoneK562H3k4me3StdAlnRep1.bam -o H3K4me3 --gt $gt --outputzero

# Output 100-kbp bigWig file:
parse2wig+ -i wgEncodeUwHistoneK562H3k4me3StdAlnRep1.bam -o H3K4me3 --gt $gt --binsize 100000

# Add the --pair option for paired-end files:
parse2wig+ -i ChIP.paired.bam -o ChIP --gt genometable.txt --pair

### Check the statistics file
cat parse2wigdir+/H3K4me3.100.tsv

# Summarize the stats in one line (useful for summarization in TSV file)
~/git/DROMPAplus/scripts/parsestats4DROMPAplus.pl parse2wigdir+/H3K4me3.100.tsv

### Total read normalization
# Normalize reads for the total number of (nonredundant) mapped reads on the genome
gt=../data/genometable/genometable.hg19.txt
parse2wig+ -i wgEncodeUwHistoneK562H3k4me3StdAlnRep1.bam -o H3K4me3 --gt $gt -n GR

# Normalize reads for the total number of (nonredundant) mapped reads on each chromosome
# (useful in case that the reads are mapped on a specific chromosome only)
gt=../data/genometable/genometable.hg19.txt
parse2wig+ -i wgEncodeUwHistoneK562H3k4me3StdAlnRep1.bam -o H3K4me3 --gt $gt -n CR

# In default, the read number after normalization is 20000000 (20M).
# If the raw read number is much smaller than 20M, it is better to decrease the number (--nrpm) not to increase the noisy peaks.
parse2wig+ -i sample.sam -o sample --gt genometable.txt -n GR --nrpm 10000000

### High resolution with central regions of fragments
# Consider only 50 bp around the center of each fragment:
parse2wig+ -i wgEncodeUwHistoneK562H3k4me3StdAlnRep1.bam -o H3K4me3 --gt $gt --rcenter 50

### Mappability normalization
# Normalize by "mappable" genonic length (--mptable option)
gt=../data/genometable/genometable.hg19.txt
mptable=../data/mptable/mptable.UCSC.hg19.36mer.flen150.txt
parse2wig+ -i wgEncodeUwHistoneK562H3k4me3StdAlnRep1.bam -o H3K4me3 --gt $gt --mptable $mptable
# Base-pair normalization based on mappability (--mpdir)
# Note: the mappability files for several species are available on our Google Drive account
# https://drive.google.com/drive/folders/1GfKZkq3HIcMLQt-pZ_4bfwh21NyS2O-5?usp=sharing
parse2wig+ -i wgEncodeUwHistoneK562H3k4me3StdAlnRep1.bam -o H3K4me3 --gt $gt --mpdir <mpdir>

### Estimate GC content bias
# This option needs FASTA files for reference chromosomes.
# --chrdir option specifies the directory that contains the reference chromosome FASTA files.
parse2wig+ -i wgEncodeUwHistoneK562H3k4me3StdAlnRep1.bam -o H3K4me3 --gt $gt  --chrdir <chromosomedir>
