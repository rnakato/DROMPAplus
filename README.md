# DROMPAplus

# 1. Overview
DROMPA (DRaw and Observe Multiple enrichment Profiles and Annotation) is a ChIP-seq pipeline tool that meets various needs, including quality check, analysis, and visualization of multiple ChIP samples.

DROMPAplus is written in C++ and has many valuable features. DROMPAplus:
* Accepts multiple map file formats (SAM, BAM, CRAM, Bowtie, TagAlign(.gz)) and read distribution formats (WIG(.gz), bigWig, bedGraph).
* Supports spike-in normalization and total read normalization.
* Outputs various quality metrics for ChIP-seq analysis.
* Visualizes read distributions in conventional PDF format; therefore, no additional programs are required which is preferable for many users, especially when sharing results (e.g., on cloud storage) with collaborators who do not have a strong bioinformatics background.
* Automatically estimates the fragment length from single-end reads using SSP.
* Can visualize two samples in a single line, which delineates the co-occurrence (e.g., H3K4me3 and H3K27ac) and exclusivity (e.g., H3K27me3 and H3K36me3) of read enrichment. Transparency (alpha) of read color can be specified.
* Supports chromatin loops from ChIA-PET (Mango format) and Hi-C (HICCUPS format) with colors corresponding to the p-values.
* Bases the HEATMAP command on Python3, which enables flexible customization.

See http://drompaplus.readthedocs.io/ for the detailed Manual.

# 2. Requirements
DROMPAplus requires the following programs and libraries:
* [Boost C++ library](http://www.boost.org/)
* [Cairo libraries](http://www.cairographics.org/)
* [GTK library](http://www.gtk.org/)
* [GNU Scientific Library](http://www.gnu.org/software/gsl/)
* [zlib](http://www.zlib.net/)
* [SAMtools](http://samtools.sourceforge.net/) (for BAM/CRAM formatted input)

and also contains two submodules:

* [SSP](https://github.com/rnakato/SSP)
* [Coherent PDF](http://community.coherentpdf.com/)

# 3. Install

### 3.1. Docker image

We recommend to use the latest Docker image of DROMPAplus from [DockerHub](https://hub.docker.com/r/rnakato/ssp_drompa) using:

    docker pull rnakato/ssp_drompa
    docker run -it --rm rnakato/ssp_drompa drompa+

For Singularity:

    singularity pull ssp_drompa.img docker://rnakato/ssp_drompa
    singularity exec ssp_drompa.img drompa+

### 3.2. Building from source

#### 3.2.1. Install required libraries
On Ubuntu:

    sudo apt install git build-essential libgtkmm-3.0-dev libboost-all-dev \
    libgsl-dev libz-dev samtools r-base

On CentOS:

    sudo yum -y install git gcc-c++ boost-devel zlib-devel gsl-devel gtkmm30-devel

#### 3.2.2. Install DROMPAplus
    git clone --recursive https://github.com/rnakato/DROMPAplus
    cd DROMPAplus
    make

If you get an installation error, make sure that all required libraries are successfully installed.

#### 3.2.3. Add the PATH environment variable
Permanently set the path to the DROMPAplus download directory by updating your **~/.bashrc** file. For example, if you downloaded DROMPAplus into $HOME directory, add the following lines to **~/.bashrc**:

    export PATH = $PATH:$HOME/DROMPAplus/bin:$HOME/DROMPAplus/otherbins

# 4. Reference
For DROMPAplus:
* Nakato R., Sakata T., [Methods for ChIP-seq analysis: A practical workflow and advanced applications](https://www.sciencedirect.com/science/article/pii/S1046202320300591), *Methods*, 2020.

For DROMPA:
* Nakato R, Shirahige K., Statistical Analysis and Quality Assessment of ChIP-seq Data with DROMPA, *Methods in Molecular Biology*, 2018.
