=====================================================
DROMPAplus: a pipeline tool for ChIP-seq analysis
=====================================================

DROMPA (DRaw and Observe Multiple enrichment Profiles and Annotation) is a ChIP-seq pipeline tool that meets various needs, including quality check, analysis, and visualization of multiple ChIP samples.

DROMPAplus is an update of `DROMPA3 <https://github.com/rnakato/DROMPA3>`_. It is written in C++ and runs from a single launch command on conventional Linux systems.

* The main features of DROMPAplus:

    * Accepts multiple map file formats (SAM, BAM, CRAM, Bowtie, TagAlign(.gz)) and read distribution formats (WIG(.gz), bigWig,  bedGraph).
    * Support spike-in normalization as well as total read normalization. 
    * Output various quality metrics for ChIP-seq analysis.
    * Visualize read distribution in conventional PDF format, which is preferable for many users, especially when sharing results (e.g., on a cloud storage) with other collaborators who do not have a strong bioinformatics background, because no additional programs are required.
    * Automatic estimation of fragment length from single-end reads using `SSP <https://github.com/rnakato/SSP>`_.
    * Can visualize two samples in one line, which delineates the co-occurrence (e.g., H3K4me3 and H3K27ac) and exclusivity (e.g., H3K27me3 and H3K36me3) of read enrichment. Transparency (alpha) of read color can be specified.
    * Support chromatin loops from ChIA-PET (Mango format) and Hi-C (HICCUPS format) with colors corresponding to the p-values.
    * HEATMAP command is now based on Python3, which enables flexible custamization.

.. image:: img/workflow.jpg
   :width: 600px
   :align: center

Contents:
---------

.. toctree::
   :numbered:
   :glob:
   :maxdepth: 1

   content/Install
   content/parse2wig
   content/drompa
   content/Appendix


Citation:
---------

* Nakato R., Sakata T., Methods for ChIP-seq analysis: A practical workflow and advanced applications, *submitted*.

Contact:
---------

:Mail: rnakato AT iam.u-tokyo.ac.jp
:Twitter: @RyuichiroNakato
