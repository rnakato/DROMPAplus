# DROMPAplus

# 1. Overview
DROMPA (DRaw and Observe Multiple enrichment Profiles and Annotation) is a ChIP-seq pipeline tool that meets various needs, including quality check, analysis, and visualization of multiple ChIP samples.
DROMPAplus is written in C++ and has many valuable features:
* accept various map file formats including CRAM and wigfile (Wig, bedGraph and bigWig)
* ability to visualize two samples in one line, which delineates the co-occurrence and exclusivity of read enrichment
* highly customizable for track heights, axis limits, and display features
* output the visualization in conventional PDF format, which is preferable for many users, especially when sharing results (e.g., on a cloud storage)


# 2. Install
DROMPAplus is written in C++ and requires the following programs and libraries:
* [Boost C++ library](http://www.boost.org/)
* [Cairo libraries](http://www.cairographics.org/)
* [GTK library](http://www.gtk.org/)
* [GNU Scientific Library](http://www.gnu.org/software/gsl/)
* [zlib](http://www.zlib.net/)
* [SAMtools](http://samtools.sourceforge.net/) (for BAM/CRAM formatted input)
* [R](http://www.r-project.org/) (for PROFILE command)

#### 2.1. Install required libraries
for Ubuntu:

    sudo apt install git build-essential libgtkmm-3.0-dev libboost-all-dev \
    libgsl-dev libz-dev samtools r-base

for CentOS:

    sudo yum -y install git gcc-c++ boost-devel zlib-devel gsl-devel gtkmm30-devel
and install samtools from [the website](http://samtools.sourceforge.net/).

#### 2.2. Install DROMPAplus
    git clone --recursive https://github.com/rnakato/DROMPAplus
    cd DROMPAplus
    make

If you get an installation error, make sure that all required libraries are installed.

#### 2.3 Submodules
DROMPAplus contains submodules in 'submodules' directory as below:

* SSP (https://github.com/rnakato/SSP)
* Coherent PDF (http://community.coherentpdf.com/)

#### 2.3. Add the PATH environment variable
For example, if you downloaded DROMPA+ into the $HOME/my_chipseq_exp directory, type:

    export PATH = $PATH:$HOME/my_chipseq_exp/DROMPAplus/bin
    export PATH = $PATH:$HOME/my_chipseq_exp/DROMPAplus/otherbins
    export PATH = $PATH:$HOME/my_chipseq_exp/DROMPAplus/submodules/cpdf/Linux-Intel-64bit

# 3. Usage
 See http://drompaplus.readthedocs.io/ for detail.

# 4. Reference
1. Nakato R., Shirahige K. Recent advances in ChIP-seq analysis: from quality management to whole-genome annotation, Briefings in Bioinformatics, 2016.

2. Nakato, R., Itoh T. and Shirahige K.: DROMPA: easy-to-handle peak calling and visualization software for the computational analysis and validation of ChIP-seq data, Genes to Cells, vol.18, issue 7, 2013.
