=====================================================
DROMPAplus: a pipeline tool for ChIP-seq analysis
=====================================================

.. image:: images/top.png

DROMPA (DRaw and Observe Multiple enrichment Profiles and Annotation) is a program for user-friendly and flexible ChIP-seq pipelining.
DROMPA can be used for quality check, PCR-bias filtering, normalization, peak calling, visualization and other multiple analyses of ChIP-seq data. DROMPA is specially designed so that it is easy to handle, and for users without a strong
bioinformatics background.

The main features of DROMPAplus are:
* Applicable to any species whose genomic sequence is available;
* Accepts multiple input/output file formats (SAM, BAM, Bowtie, WIG, BED, TagAlign(.gz), bigWig, and bedGraph);
* Normalization using mappability and GC content;
* Visualize ChIP/input enrichment and p-value of statistical tests along with the read distribution itself;
* Output PDF format that is suitable to share the visualization on a cloud;
* In addtion to typical peak calling, various types of ChIP-seq analysis are available.

Contents:
---------
.. toctree::
   :maxdepth: 2
	     content/Install


Citation:
---------
Nakato, R., Itoh, T. and Shirahige, K. (2013). DROMPA: easy-to-handle peak calling and visualization software for the computational analysis and validation of ChIP-seq data. Genes Cells, 18(7), 589-601.

  R
yuichiro Nakato <rnakato AT iam.u-tokyo.ac.jp>
