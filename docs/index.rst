=====================================================
DROMPAplus: a pipeline tool for ChIP-seq analysis
=====================================================

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

.. image:: img/workflow.jpg
   :scale: 100%
   :align: center

Contents:
---------
.. toctree::
   :numbered:
   :glob:
   :maxdepth: 1

   content/*


Citation:
---------
* Nakato R, Shirahige K, Statistical Analysis and Quality Assessment of ChIP-seq Data with DROMPA, Methods in Molecular Biology, 1672, 631-643, 2018.
* Nakato R, Shirahige K, Recent advances in ChIP-seq analysis: from quality management to whole-genome annotation, Briefings in Bioinformatics, vol. 18, issue 2, 279-290, 2017.

Contact:
---------

:Mail: rnakato AT iam.u-tokyo.ac.jp
:Twitter: @RyuichiroNakato
