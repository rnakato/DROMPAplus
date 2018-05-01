Parse2wig
============

parse2wig preprocesses an input mapfile into bin data (the number of mapped read per bin).
The length of each read is calculated automatically. For single-end mode, mapped reads are extended to the expected DNA-fragment length.

---------------
1. Main options
---------------


1.1 Example
-------------------------------

The command::

  $ parse2wig -i sample.sam -o sample -gt genometable.txt -binsize 100

generates bin files from sample.sam with a bin size of 100 bp. The parse2wigdir directory is created and the bin files are outputted into the directory.
A bin file is outputted for each chromosome (in the case of binary and WIG outputs) or for the whole-genome (in the case of bedGraph and bigWig outputs).

The command::

  $ parse2wig -f BAM -i sample.bam -o sample -gt genometable.txt -of 3

reads BAM file and outputs bin data in the bedGraph format.
Furthermore, multiple mapfiles can be given as one sample (separated by a ’,’)::
  
  $ parse2wig -i sample_rep1.sam,sample_rep2.sam,sample_rep3.sam \
  $ -o sample -gt genometable.txt -binsize 100

1.2 Paired-end file
-------------------------------
You must supply the “-pair” option for paired-end files (when parsing paired-end mapfiles with single-end mode, warning messages are outputted). In paired-end mode, each fragment length is calculated from the mapfile automatically. Inter-chromosomal read-pairs and read-pairs longer than the maximum fragment length (specified by the “-maxins” option) are ignored.

Note: For TagAlign format, paired-end data is not supported.

1.3 Multiple mapped reads
-------------------------------

parse2wig automatically recognizes the uniquely mapped and multiple mapped reads. For multiple mapped reads, each mapped locus is weighted equally. Thus, the total number of reads mapped into bin x is r x = k∈R 1/n k where n k is the number of times that read k is mapped onto the reference genome and R is the full set of reads mapped onto bin x.

Note: For SAM and BAM format, while parse2wig uses the “NH” flag to check multiple mapped reads, some mapping tools such as Bowtie and BWA do not output the “NH” column. In those cases, all reads are considered as ‘uniquely mapped’. Therefore, we recommend the Bowtie format when treating multiple mapped reads.

1.4 Higher resolution with central regions of fragments
-------------------------------------------------------------
When high resolution is required (e.g., nucleosome-seq), it is better to consider only central regions of each fragment.
For such purpose, parse2wig has the option “-rcenter”. The command::

  $ parse2wig -i sample.sam -o sample -gt genometable.txt -flen 200 -rcenter 50

assumes the averaged fragment length is 200 bp and consider only 50 bp around the center of each fragment.

1.5 Statistics file (for quality check)
-------------------------------------------------------------
In addition to the bin files, parse2wig also outputs the statistics of the input file into the output directory, which are useful to check the quality of the sample. The commands in subsection 5.1.1 produce three types of statistics files: “sample.100.xls”, “sample.binarray dist.xls” and “sample.readlength dist.xls”.

1.5.1 sample.100.xls
-------------------------------------------------------------
The contents of “sample.100.xls” are the following:
- Input file name;
- Redundancy threshold: threshold used for PCR bias;
- Library complexity [2]
(the number of reads tested are in parenthesis)
(if the number of reads tested is insufficient, the score is put in parentheses);
- GC summit: the summit of GC distribution (when -GC is supplied);
- Poisson and Negative binomial: estimated parameter;

- (both whole-genome and chromosomal stats as follows);
- Length: total length;
- Mappable base and mappability: mappability calculated from specified mappability file;
- Total reads: the number of reads in the input file;
- Non-redundant reads: the number of reads remaining after PCR-bias filtering;
- Redundant reads: the number of reads filtered by PCR-bias filtering (mapped on forward, reverse and both strands are outputted);
- Reads (GCnormed): read number after GC normalization (when -GC is supplied);
- Read depth;
- Scaling weight: total weight of the read normalization;
- Normalized read number: read number after the total read normalization; (this read number is used in drompa peakcall and drompa draw)
- FRiP (fraction of reads in peaks) score [2] (when -bed is supplied).
  
For the quality check, library complexity and the number of non-redundant reads are especially important.

1.5.2 sample.binarray dist.xls
-----------------------------------
“sample.binarray dist.xls” describes the distribution of read numbers contained in each bin. Distributions of simulated data are also shown.

1.5.3 sample.readlength dist.xls
------------------------------------
“sample.readlength dist.xls” describes the read length distribution of the input mapfile. Only lengths in which the read number is nonzero are considered.

1.6 Filtering reads
---------------------
parse2wig filters “redundant reads” (reads starting exactly at the same 5’ ends) as “PCR bias” [1].
This filtering step can be omitted by supplying “-nofilter” option.
By default, the threshold of filtering is defined as:
thre pcr = max(1, 10 ∗ E genome )
where E genome is the averaged read depth. This is because E genome can be greater than 1 for a small genome.
thre_pcr can be supplied manually through the “-thre pb” option.
The number of redundant/non-redundant reads and library complexity [2] can be checked using the generated statistics file (see section 5.2). Since the library complexity depends on the number of mapped reads, parse2wig uses the library complexity for 10 million mapped reads.
This default number can be changed through the “-num4cmp” option.

1.7 Total read normalization
-------------------------------
For the comparison of multiple ChIP samples, read number normalization is necessary.
parse2wig has the “-n” option to normalize the bin data with the number of total mapped reads (after PCR-bias filtering).
  -n {NONE|GR|GD|CR|CD} (default:NONE)
  NONE; not normalize
  GR; for whole genome, read number
  GD; for whole genome, read depth
  CR; for each chromosome, read number
  CD; for each chromosome, read depth
  -np <int> read number after normalization
  (default: 10000000 (10 million))
  -nd <double>
  depth after normalization (default: 0.1)

  The users can choose total reads or read depth for normalization. For example, the command::

    $ parse2wig -i sample.sam -o sample -gt genometable.txt -n GR -np 20000000

scales bin data so that the total number of mapped reads (after filtering) onto the whole genome is 20 million.
The normalization for each chromosome (CR or CD) is useful when the large difference in one chromosome affects to whole-genome (e.g., rDNA regions in chromosome XII for Saccharomyces cerevisiae).

Note: it is not recommended to scale a small number of reads up to a larger number because
that will result in plenty of background noise (e.g., 1 million → 10 million).

1.8 Mappability
-------------------------------
parse2wig can normalize reads based on the genome mappability [3] by supplying mappability
files as follows::

  $ parse2wig -i sample.sam -o sample -gt genometable.txt \
  $ -mp mappability/map_fragL150

When “-mp” is not supplied, all bases are considered as mappable. The low mappability regions
(“-mpthre” option, < 0.3 (30%) as default) are ignored after ChIP-seq analysis.
DROMPA adopts the mappability files generated through the scripts provided by MOSAiCS [4].
See section 9.1 for details.

1.9 GC content
-------------------------------
Sometimes the sequenced data has much GC bias.
In those cases, GC normalization is necessary.
parse2wig can adopt a GC normalization similar to BEADS [5].
This procedure requires the FASTA files of chromosomes and the binary mappability files.
The command::

  $ parse2wig -i sample.sam -o sample -gt genometable.txt \
  $ -GC <chromosomedir> -mpbin mappability/map -flen4gc 100

calculates the GC contents of the input file using the central 100 bp of each fragment.
<chromosomedir> is the directory that contains the FASTA files of all chromosomes described
in genometable.txt with corresponding filenames. For example, if “chr1” is in genometable.txt,
there should be “chr1.fa” in <chromosomedir>. “-mpbin” specifies the binary mappability text
files (see section 9.1 for details).

Note: Since this GC normalization scheme is under development, if a sample has a GC dis-
tribution quite different from other samples, it is better to consider re-preparing the sample rather
than using it with GC normalization.

1.9.1 GC distribution file
-------------------------------
parse2wig uses the longest chromosome described in genometable.txt for GC bias estimation.
When using GC normalization, the GC distribution file “sample.GCdist.xls” is also outputted
into the output directory. The contents are the following:
- GC: the GC content;
- Genome prop: the proportion of the mappable bases containing the GC contents, then
prop GC = n GC/G, where n GC
are the number of positions containing the GC contents and G is the total number of mappable bases;
- Read prop: the proportion of the reads (fragments) containing the GC contents, then prop GC
= n GC /N, where n GC are the number of reads containing the GC contents and N is the total number of mapped reads;
- Depth: the ratio of GC contents between reads and genome sequence, namely, depth GC = reads genome
n GC /n GC ;
- Weight: the ratio of the proportion between reads and genome sequence, namely, weight = reads prop GC /prop GC

Because GC contents with low depth (depth GC ) cause background noise, by default parse2wig 
sets a weight of 1 to the GC content with depth GC less than 0.001, and a weight of 0 to the GC
genome content having prop GC less than 0.00001.
When supplying the “-gcdepthoff” option, the former threshold is ignored.
Using the GC distribution file, the user can draw GC and weight distribution of the input file
and the genome sequence. Figure 2 shows an example.

1.9.2 Ignore peak regions
----------------------------
For ChIP samples, it is necessary to ignore fragments that overlap with peak regions and use back-
ground reads only, because ChIP reads can have different GC distribution from the background.
To do that, specify a peak list using the “-bed” option::

  $ parse2wig -i sample.sam -o sample -gt genometable.txt \
  $ -GC <chromosomedir> -mpbin mappability/map -flen4gc 100 -bed peaklist.bed

1.10 Cross-correlation analysis
--------------------------------
Optionally, parse2wig can output the cross-correlation profile [2] with a strategy similar to spp
[6] by specifying the “-ccp” option.
The command
$ parse2wig -i sample.sam -o sample -gt genometable.txt -ccp
generates “sample.ccp.xls” in the output directory, which describes the cross-correlation plot be-
tween the read number of forward and reverse strands from -500 to 1500 bp with a 5 bp step.
In version 3.0.0, the value of bins that hace above the 95 th percentile is reduced to 95 th per-
centile on the cross-correlation analysis.