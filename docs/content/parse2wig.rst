============
Parse2wig
============

parse2wig preprocesses an input mapfile into bin data (the number of mapped read per bin).


.. contents::


Example
-------------------------------

The command below generates a bigWig data ``ChIP.100.bw`` and a statistics file ``ChIP.100.tsv`` in ``parse2wigdir+`` directory::

  $ parse2wig+ -i ChIP.bam -o ChIP --gt genometable.txt
  (...)
  $ ls parse2wigdir+
  parse2wigdir+/ChIP.100.bw   parse2wigdir+/ChIP.100.tsv

The default bin size is 100 bp.
The input file format is automatically detected by postfix (.sam/.bam/.cram/.bowtie/.tagalign(.gz)).
If the detection does not work well, supply ``-f`` option (e.g., ``-f BAM``).

parse2wig+ allows multiple input files (separated by ",")::

  $ parse2wig+ -i ChIP1.bam,ChIP2.bam,ChIP3.bam -o ChIP --gt genometable.txt

If you want to generate it as bedGraph format, supply ``--outputformat 2``::

  $ parse2wig+ -i ChIP.bam -o ChIP --gt genometable.txt --outputformat 2

In default, parse2wig+ omits to output bins in which the value is zero to reduce the file size. Supply ``--outputzero`` to output all bins::

  $ parse2wig+ -i ChIP.bam -o ChIP --gt genometable.txt --outputzero

For bin size 100kbp::

  $ parse2wig+ -i ChIP.bam -o ChIP --gt genometable.txt --binsize 100000

To use multiple cpus, supply ``-p``::

  $ parse2wig+ -i ChIP.bam -o ChIP --gt genometable.txt -p 4

.. note::

    * Multithreading is activated only in strand-shift profile for estimating frangment length. When suppying ``--nomodel`` option, multithreading will make no differece.


Quality check
------------------------

The statistics file ``parse2wigdir+/ChIP.100.tsv`` describes various quality values and other information of the input mapfile as follows:

- **Redundancy threshold**: threshold for PCR bias
- **Library complexity**: non-redundant read fraction for 10 million mapped reads (supply ``--ncmp`` to change this mapped read number). This score is put in parentheses when the number of mapped reads is insufficient.
- **GC summit** (when supplying ``--GC`` option): the summit of GC distribution

- **Length**: total length of the genome and chromosomes
- **Mappable base and mappability**: mappability calculated from specified mappability file
- **Total reads**: the number of mapped reads
- **Non-redundant reads**: the number of reads remaining after PCR-bias filtering
- **Redundant reads**: the number of reads filtered by PCR-bias filtering (mapped on forward, reverse and both strands are outputted)
- **Reads (GCnormed)** (when supplying ``--GC`` option): read number after GC normalization
- **Read depth**: the expected number of mapped reads per base pair
- **Scaling weight**: scaling weight for the read normalization
- **Normalized read number**: read number after the normalization
- **FRiP score** (when supplying ``--bed`` option): fraction of reads in peaks. .

Fragment length estimation
+++++++++++++++++++++++++++++++++++

For single-end data, DROMPAplus internally uses `SSP <https://github.com/rnakato/SSP>`_ to estimate averaged fragment length and extends to the length.
When ``--nomodel`` option is supplied, DROMPAplus omits to use SSP and extends read to a predetermined length (150 bp in default. Supply ``--flen`` option to change this value). 

Paired-end file
+++++++++++++++++++++++++++++++++++

Supply ``--pair`` option for paired-end files::

  $ parse2wig+ --pair -i ChIP.paired.bam -o ChIP --gt genometable.txt

In ``--pair`` mode, the fragment length of each read pair is calculated automatically.
parse2wig+ discards read pairs that are mapped onto different chromosomes or the fragment length is longer than 500bp (in default, specified ``-maxins`` to change).

.. note::

   * When parsing paired-end mapfiles with single-end mode, warning messages will be outputted.
   * In TagAlign format, paired-end data is not supported.

PCR bias filtering
++++++++++++++++++++++

parse2wig+ filters "redundant reads" (reads starting exactly at the same 5' ends) as "PCR bias".
This filtering step can be omitted by supplying ``-nofilter`` option.

By default, the threshold of filtering is defined as::

	threshold = max(1, 10 * E_genome)

where E\ :sub:`genome`\  is the averaged read depth.
10 * E\ :sub:`genome`\  can be greater than 1 for a small genome (e.g., yeast).
Additionally, ``--thre_pb`` option can be used to fix this threshold.


Multiple mapped reads
++++++++++++++++++++++++++++++

parse2wig+ recognizes the uniquely mapped and multiple mapped reads using 'NH' flag in SAM/BAM/CRAM format. For multiple mapped reads, each mapped locus is weighted equally.

Some mapping tools (e.g., Bowtie and BWA) do not output the 'NH' column. In this case, all reads are considered as uniquely mapped reads.

Total read normalization
---------------------------------

parse2wig+ has the ``-n`` option to normalize the read distribution based on the number of nonredundant reads

* **-n NONE** (default); not normalize
* **-n GR**; for whole genome, read number
* **-n GD**; for whole genome, read depth
* **-n CR**; for each chromosome, read number
* **-n CD**; for each chromosome, read depth

``-n GR`` is a typical total read normalization.
If the mapped read number is quite different among chromosomes (e.g., mapfile contains chrX only), consider to use ``-n CR``. Also, use ``-np`` option to change read number after normalization (default: 10 million). 

For example, the command below scales bin data so that the total number of nonredundant reads is 20 million::

    $ parse2wig+ -i sample.sam -o sample --gt genometable.txt -n GR --np 20000000

.. note::

       it is not recommended to scale a small number of reads up to a larger number (e.g., 1 million → 10 million) because that will result in plenty of background noise.

Higher resolution with central regions of fragments
-------------------------------------------------------------

When high resolution is required (e.g., nucleosome-seq), consider ``--rcenter`` option that focus on the the central region of each fragment. 
For example, the command below considers only 50 bp around the center of each fragment::

  $ parse2wig+ -i ChIP.bam -o ChIP --gt genometable.txt --rcenter 50

Mappability
-------------------------------

parse2wig+ can normalize reads based on the genome mappability by supplying mappability files::

  $ parse2wig+ -i ChIP.bam -o ChIP --gt genometable.txt --mp mappability/map_fragL150

When ``--mp`` is not supplied, all bases are considered as mappable.
The low mappability regions (``--mpthre`` option, < 0.3 as default) are ignored in mappability normalization (and GC normalization).

The mappability files for several species are available in /DROMPAplus/data/mptable/ directory.

GC content
-------------------------------

Sometimes the sequenced data has much GC bias.
In those cases, GC normalization is necessary.
parse2wig+ can implement a GC normalization by supplying the chromosome FASTA files by ``--GC`` option and the binary mappability files by ``--mpbin`` option.
The command::

  $ parse2wig+ -i sample.sam -o sample --gt genometable.txt \
  $ --GC <chromosomedir> --mpbin mappability/map

calculates the GC contents of the input file using the central 100 bp of each fragment.
``<chromosomedir>`` is the directory that contains the FASTA files of all chromosomes described
in ``genometable.txt`` with corresponding filenames. For example, if ``chr1`` is in ``genometable.txt``,
there should be ``chr1.fa`` in ``<chromosomedir>``. ``-mpbin`` specifies the binary mappability text
files (see section 9.1 for details).

.. note:: 
    
    Since this GC normalization scheme is under development, if a sample has a GC distribution quite different from other samples, it is better to consider re-preparing the sample rather than using it with GC normalization.
