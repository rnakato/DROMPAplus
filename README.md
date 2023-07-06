# DROMPAplus

## 0. Changelog
See [Changelog](https://github.com/rnakato/DROMPAplus/blob/master/ChangeLog.md)

## 1. Overview
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

## 2. Install

### 2.1. Docker image

We recommend to use the latest Docker image of DROMPAplus from [DockerHub](https://hub.docker.com/r/rnakato/ssp_drompa).

#### 2.1.1. Docker 
 To use docker command, type:

    docker pull rnakato/ssp_drompa
    docker run -it --rm rnakato/ssp_drompa drompa+

**Note**: When using the docker image, it is necessary to mount the directory by ``-v`` option to access the input files as follows:

    docker run -it --rm -v $(pwd):/mnt rnakato/ssp_drompa parse2wig+ \
             -i /mnt/ChIP.bam -o ChIP --odir /mnt/parse2wigdir+ --gt /mnt/genometable.txt

This command mounts the current directory to ``/mnt`` directory in the container. 
Please see also the [document of Docker](https://docs.docker.com/storage/volumes/).

#### 2.1.2. Singularity

Singularity can also be used to execute the docker image:

    singularity build ssp_drompa.sif docker://rnakato/ssp_drompa
    singularity exec ssp_drompa.sif drompa+

Singularity mounts the current directory automatically. If you access the files in the other directory, please mount by `--bind` option, for instance:

    singularity exec --bind /work ssp_drompa.sif drompa+
    
This command mounts `/work` directory.

### 2.2. Building from source

#### 2.2.1 Requirements
DROMPAplus requires the following programs and libraries:
* [Boost C++ library (>1.53.0)](http://www.boost.org/)
* [Cairo libraries](http://www.cairographics.org/)
* [GTK library](http://www.gtk.org/)
* [GNU Scientific Library (>1.15)](http://www.gnu.org/software/gsl/)
* [zlib (>1.2.7)](http://www.zlib.net/)
* [CMake (>2.8)](https://cmake.org/)
* [HTSlib (1.10.2)](https://github.com/samtools/htslib) (for SAM/BAM/CRAM formatted input)

and also contains two submodules:

* [SSP](https://github.com/rnakato/SSP)
* [Coherent PDF](http://community.coherentpdf.com/)

#### 2.2.2. Install required libraries
On Ubuntu:

    sudo apt install git build-essential libgtkmm-3.0-dev libboost-all-dev \
    libz-dev libbz2-dev libgsl-dev libcurl4-gnutls-dev cmake

On CentOS:

    sudo yum -y install git gcc-c++ clang boost-devel zlib-devel gsl-devel gtkmm30-devel bzip2-devel cmake

On Mac:

     brew install gsl gtk gtkmm cairo pkgconfig curl xz zlib boost cmake

#### 2.2.3. Install DROMPAplus

    git clone --recursive https://github.com/rnakato/DROMPAplus
    cd DROMPAplus
    make

If you get an installation error, make sure that all required libraries are successfully installed.

#### 2.2.4. Add the PATH environment variable
Permanently set the path to the DROMPAplus download directory by updating your **~/.bashrc** file. For example, if you downloaded DROMPAplus into $HOME directory, add the following lines to **~/.bashrc**:

    export PATH=$PATH:$HOME/DROMPAplus/bin:$HOME/DROMPAplus/otherbins

# 3. Reference
For DROMPAplus:
* Nakato R., Sakata T., [Methods for ChIP-seq analysis: A practical workflow and advanced applications](https://www.sciencedirect.com/science/article/pii/S1046202320300591), *Methods*, 2020.

For DROMPA:
* Nakato R, Shirahige K., Statistical Analysis and Quality Assessment of ChIP-seq Data with DROMPA, *Methods in Molecular Biology*, 2018.
