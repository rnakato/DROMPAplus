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
   :width: 600px
   :align: center

Contents:
---------
.. toctree::
   :numbered:
   :glob:
   :maxdepth: 1

   content/*


Citations:
---------
Please cite DROMPAplus as follows:

* Nakato R, Shirahige K, Methods for ChIP-seq analysis: A practical workflow and advanced applications, *in prep*.

Contact:
---------

:Mail: rnakato AT iam.u-tokyo.ac.jp
:Twitter: @RyuichiroNakato
