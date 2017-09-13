DROMPA+
============
drompa+ can visualize multiple ChIP samples with specified genome annotation, using modes for the implementation of various types of ChIP-seq analysis:
PC_SHARP peak-calling (for sharp mode)
PC_BROAD peak-calling (for broad mode)
PC_ENRICH peak-calling (enrichment ratio)
GV global-view visualization
PD peak density 
CI compare peak-intensity between two samples
CG output ChIP-reads in each gene body
GOVERLOOK genome-wide overlook of peak positions
PROFILE make R script of averaged read density
HEATMAP make heatmap of multiple samples
TR calculate the travelling ratio (pausing index) for each gene

7.1 Output of parameter used
-------------------------------
When excuting drompa peakcall and drompa draw, the summary of parameters specified is outputted to STDOUT. The users can check whether the command is specified as expected.

7.2 Read distribution visualization (PC SHARP)
--------------------------------------------------------------
drompa draw can take multiple ChIP-input pairs as input. Each pair should be specified with the option “-i”, as in the drompa peakcall mode. For example, the command::
  
  $ IP1=’parse2wigdir/ChIP1’
  $ IP2=’parse2wigdir/ChIP2’
  $ IP3=’parse2wigdir/ChIP3’
  $ IP4=’parse2wigdir/ChIP4’
  $ Input=’parse2wigdir/Input’
  $ drompa_draw PC_SHARP -p ChIPseq -gt genometable.txt \
  $ -i $ChIP1,$Input,ChIP1 \
  $ -i $ChIP2,$Input,ChIP2 \
  $ -i $ChIP3,$Input,ChIP3 \
  $ -i $ChIP4,$Input,ChIP4 \
  $ -gene refFlat.txt -ls 1000 -lpp 2 -show_itag 2 -scale_tag 30

generates the PDF files “ChIPseq*.pdf” 45 for four ChIP samples (ChIP1, 2, 3 and 4) and using the same Input sample (Input), as shown in Figure 3a.
By default, drompa draw visualizes ChIP-read lines only. The “-show itag 1” option displays input lines for all ChIP samples while the “-show itag 2” option displays only the line for first input (Figure 3a). The latter is recommended when the same input sample is used for all ChIP samples.
Similarly, to display the lines of p-value and enrichment line (Figure 3b), type::

  $ drompa_draw PC_SHARP -p ChIPseq -gt genometable.txt \
  $ -i $ChIP1,$Input,ChIP1 \
  $ -i $ChIP2,$Input,ChIP2 \
  $ -i $ChIP3,$Input,ChIP3 \
  $ -i $ChIP4,$Input,ChIP4 \
  $ -gene refFlat.txt -showratio 1 -showpinter 1 -showpenrich 1 \
  $ -scale_tag 30 -scale_ratio 3 -scale_pvalue 3

where the “-scale tag”, “-scale ratio” and “-scale pvalue” options change the maximum values for the y axis of the corresponding lines.

7.2.1 Specify different parameter for each sample pair
--------------------------------------------------------------
For drompa draw, the option “-i” can take the following comma-separated multiple fields:
1. ChIP sample (required);
2. Input control sample;
3. Sample name to be shown in figure;
4. peak list to be highlighted;
5. binsize;
6. scale tag;
7. scale ratio;
8. scale pvalue.

Except for the “ChIP sample”, all the other fields can be omitted. When the peak list is specified, drompa draw highlights the specified peak regions instead of using the internal peak-calling engine, which is useful when comparing among multiple peak-calling programs. The rest of the options can be used to specify different parameters for each sample pair. The command::

  $ drompa_draw PC_SHARP -p ChIPseq -gt genometable.txt \
  $ -i $ChIP1,$Input,ChIP1,ChIP1peak.bed,,50 \
  $ -i $ChIP2,$Input,ChIP2,ChIP2peak.bed,,80 \
  $ -i $ChIP3,$Input,ChIP3,ChIP3peak.bed,1000,60 \
  $ -i $ChIP4,$Input,ChIP4,ChIP4peak.bed,1000,60 \
  $ -gene refFlat.txt -lpp 1 -chr 4 -ls 10000 -rmchr -sm 2000

generates the results presented in Figure 4. The parameter for each sample is superior to the global parameters.

The previous command can be rewritten as::

  $ s1="-i $ChIP1,$Input,ChIP1,ChIP1peak.bed,,50"
  $ s2="-i $ChIP2,$Input,ChIP2,ChIP2peak.bed,,80"
  $ s3="-i $ChIP3,$Input,ChIP3,ChIP3peak.bed,1000,60"
  $ s4="-i $ChIP4,$Input,ChIP4,ChIP4peak.bed,1000,60"
  $ drompa_draw PC_SHARP -p ChIPseq -gt genometable.txt $s1 $s2 $s3 $s4 \
  $ -gene refFlat.txt -lpp 1 -chr 4 -ls 10000 -rmchr -sm 2000

where “ $s1 $s2 $s3 $s4” are used as ChIP-input sample pairs.

7.3 Read distribution visualization (PC BROAD)
--------------------------------------------------------------
To identify broadly enriched region, use the PC BROAD mode as follows::
  
  $ drompa_draw PC_BROAD -p ChIPseq_broad -gt genometable.txt $s1 $s2 $s3 $s4 \
  $ -gene refFlat.txt -showratio 1 -showpinter 1 -showpenrich 1

7.4 Enrichment visualization (PC ENRICH)
--------------------------------------------------------------
For a small genome, such as the yeast’s, the sequencing depth is generally enough (> 10 fold).
In such cases, the genome-wide ChIP/Input enrichment distribution is informative because the
technical and biological bias in high throughput sequencing can be minimized.
To make a PDF file of the enrichment distribution for S. cerevisiae (Figure 5a), type::

  $ drompa_draw PC_ENRICH -p ChIPseq_enrich $s1 $s2 $s3 $s4 $s5 $s6 \
  $ -gt genometable.txt -gene SGD_features.tab -gftype 3 -ars ARS-oriDB.txt \
  $ -lpp 1 -scale_ratio 4 -ls 100

Supply “-showratio 2” to use logratio. If you want to check the enrichment precisely, it is good to adjust the y-axis (Figure 5b) as follows::

  $ drompa_draw PC_ENRICH -p ChIPseq_enrich $s1 $s2 $s3 $s4 $s5 $s6 \
  $ -gt genometable.txt -gene SGD_features.tab -gftype 3 -ars ARS-oriDB.txt \
  $ -lpp 1 -scale_ratio 4 -ls 100 -bn 5 -ystep 10

7.5 Annotation data for drompa draw
--------------------------------------------------------------
DROMPA accepts annotation data from the publicly accessible websites listed below. These annotation files can also be downloaded from the DROMPA website.

7.5.1 Gene annotation data
--------------------------------------------------------------
DROMPA+ accepts the Gtf, refFlat and “SGD features.tab” obtained from the Saccharomyces Genome Database (SGD) format for gene annotation.

- RefSeq annotation (refFlat format) obtained from the UCSC Genome Browser website [8].
• Ensembl gene data. The data for several species can be downloaded from the DROMPA
website.
• For the gene annotation data of S. pombe, download a GFT-formatted file (e.g., “schizosac-
charomyces pombe.EF1.62.gtf”) from the Ensembl website.
Supply the option “-gene” to specify gene data.

7.5.2 Replication origin data
--------------------------------------------------------------
DROMPA can visualize DNA replication origin data (ARS) available for S. cerevisiae and S.pombe.
The annotation data can be obtained from OriDB 7 . Download the origin list and supply with the option “-ars”.

7.5.3 Mappability and Gap-region data
--------------------------------------------------------------
If the mappability file and/or gap regions (filled with “Ns”) are supplied through the “-mp” and “-gap” options, the low mappable regions and gap regions are shaded in purple and gray in the figure, respectively. See section 9.1 for details on how to generate these data.::

  $ drompa_draw PC_SHARP -p ChIPseq -gt genometable.txt $s1 $s2 $s3 $s4 \
  $ -gene refFlat.txt -mp mappability/map_fragL150 -gap mappability/N_fragL150

7.5.4 Showing limited regions
--------------------------------------------------------------
When the “-chr” option specified, only the specified chromosome is outputted::

  $ drompa_draw PC_SHARP -p ChIPseq -gt genometable.txt $s1 $s2 $s3 $s4 \
  $ -gene refFlat.txt -chr 12

This command outputs the result of chromosome 12 only 8 .
To focus on specific regions (in this example, the HOX A cluster region), supply a BED file describing the regions to be shown with the option “-r”as follows::
  $ echo "chr7 27100000 27280000" > HOXA.txt
  $ drompa_draw PC_SHARP -gene refFlat.txt $s1 $s2 $s3 $s4 -p HOXA \
  $ -gt genometable.txt -r HOXA.txt -ls 300

7.5.5 Repeat data (RepBase) and GC contents
--------------------------------------------------------------
DROMPA can incorporate the BED-formatted GC content files and RepBase files using the options “-repeat” and “-GC”, respectively. These data can be obtained from the Table Browser of the UCSC Genome Browser [8].
$ drompa_draw PC_SHARP -p ChIPseq -gt genometable.txt $s1 $s2 $s3 $s4 \
$ -repeat RepeatMasker_hg19.txt -GC GCcontents/ -gcsize 1000
where “-gcsize” specifies the window size of GC contents. GC content files should be chromosome-separated in the specified directory (chr*-bs*).

To supply an arbitrary window size, the DROMPA website provides the program GCcount.pl to generate these files from a FASTA-formatted file.

7.5.6 BED annotation and long-range interactions
--------------------------------------------------------------
drompa draw accepts annotation data in BED or BED6 format (e.g., ChromHMM results [9]) with the “-bed” option.
The long-range interactions file such as ChIA-PET results are also allowed
with the “-inter” option, which takes tab-separated files with six columns: head chr, head start, head end, tail chr, tail start, and tail end. The intra- and inter-chromosomal interactions are shown in red and green, respectively.

For example, the following command generates the PDF file shown in Figure 6::

  $ drompa_draw PC_SHARP -p ChIP-seq -gt genometable.txt $s1 $s2 $s3 $s4 \
  $ -gene refFlat.txt -bed chromhmm.bed,emission \
  $ -inter ChIA-PET.bed,interaction
