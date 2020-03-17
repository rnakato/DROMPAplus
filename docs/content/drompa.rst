DROMPAplus: read distritbution
====================================

drompa+ can visualize multiple ChIP samples with specified genome annotation, using modes for the implementation of various types of ChIP-seq analysis:

* **PC_SHARP**: peak-calling (for sharp mode)
* **PC_BROAD**: peak-calling (for broad mode)
* **PC_ENRICH**: peak-calling (enrichment ratio)
* **GV**: global-view visualization
* **PROFILE**: makes a PDF file and corresponding R script of the averaged read density (also used to make a PNG file and corresponding Python script)

Note that the algorithms for **PC_SHARP**, **PC_BROAD** and **PC_ENRICH** are identical. Only the default parameter sets differ.
A sample script file can be found in the "tutorial" directory. 
``sample.human.sh`` describes the tutorial of these commands of drompa+ for the histone modification data of human K562 cells.


.. toctree::
   :glob:
   :maxdepth: 1

   drompa/PCSHARP
   drompa/PCENRICH
   drompa/PEAKCALL
   drompa/GV
   drompa/Profile
   drompa/Heatmap
