DROMPAplus: read distritbution
====================================

drompa+ can visualize multiple ChIP samples with specified genome annotation, using modes for the implementation of various types of ChIP-seq analysis:

* **PC_SHARP** peak-calling (for sharp mode)
* **PC_BROAD** peak-calling (for broad mode)
* **PC_ENRICH** peak-calling (enrichment ratio)
* **GV** global-view visualization
* **PROFILE** make pdf file and corresponding R script of averaged read density
* **HEATMAP** make png file and corresponding Python script


Note that the algorithm of **PC_SHARP**, **PC_BROAD** and **PC_ENRICH** is identical. Just the default parameter set is different. 
A sample script file for tutorial is in "tutorial" derectory.
``sample.human.sh`` describes the tutorial of these commands of drompa+ for histone modification data of human K562 cells.




.. toctree::
   :numbered:
   :glob:
   :maxdepth: 1

   drompa/PCSHARP
   drompa/PCENRICH
   drompa/PEAKCALL
   drompa/GV
   drompa/Profile
   drompa/Heatmap
