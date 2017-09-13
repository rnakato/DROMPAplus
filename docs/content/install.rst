Installation
============

---------------
1. Dependencies
---------------

1.1 Install required libraries
-------------------------------
DROMPAplus is written in C++11 and requires the following programs and libraries:

- Boost C++ library
- Cairo libraries
- GTK library
- GNU Scientific Library
- zlib
- SAMtools (for BAM format)
- UCSC utility (for BigWig format)
- R (for PROFILE command)

for Ubuntu::

  sudo apt-get install git build-essential libgtkmm-3.0-dev libboost-all-dev libgsl-dev libz-dev samtools r-base

for CentOS::

  sudo yum -y install git gcc-c++ boost-devel zlib-devel gsl-devel gtkmm30-devel

and install samtools from the website.

SAMtools and R can also be installed by Anaconda.

1.2 Install required software
-------------------------------

UCSC utility tools::

  wget ftp://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
  wget ftp://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph
  chmod +x bedGraphToBigWig bigWigToBedGraph

Coherent PDF (https://www.coherentpdf.com/)::

  git clone https://github.com/coherentgraphics/cpdf-binaries.git
  
2.4. Add the PATH for DROMPA, samtools, cpdf and UCSC tools
------------------------------------------------------------

For example, if you downloaded DROMPA and cpdf into the $HOME/my_chipseq_exp directory, type::

   export PATH = $PATH:$HOME/my_chipseq_exp/DROMPAplus/bin:$HOME/my_chipseq_exp/cpdf-binaries/Linux-Intel-64bit
