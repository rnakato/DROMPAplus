====================================
DROMPAplus: read distritbution
====================================

drompa+ can visualize multiple ChIP samples with specified genome annotation, using modes for the implementation of various types of ChIP-seq analysis:

* **PC_SHARP** peak-calling (for sharp mode)
* **PC_BROAD** peak-calling (for broad mode)
* **PC_ENRICH** peak-calling (enrichment ratio)
* **GV** global-view visualization
* **PROFILE** make R script of averaged read density

Note that the algorithm of **PC_SHARP**, **PC_BROAD** and **PC_ENRICH** is identical. Just the default parameter set is different. 


.. toctree::
   :numbered:
   :glob:
   :maxdepth: 1

   drompa/PCSHARP
   drompa/PCENRICH
   drompa/GV
   drompa/Profile
   drompa/Heatmap
