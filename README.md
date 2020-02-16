# DROMPAplus

# 1. Overview
DROMPA (DRaw and Observe Multiple enrichment Profiles and Annotation) is a ChIP-seq pipeline tool that meets various needs, including quality check, analysis, and visualization of multiple ChIP samples.

DROMPAplus is written in C++ and has many valuable features:
* accept various map file formats including CRAM and wigfile (Wig, bedGraph and bigWig)
* ability to visualize two samples in one line, which delineates the co-occurrence and exclusivity of read enrichment
* highly customizable for track heights, axis limits, and display features
* output the visualization in conventional PDF format, which is preferable for many users, especially when sharing results (e.g., on a cloud storage)

# 2. Requirements
DROMPAplus requires the following programs and libraries:
* [Boost C++ library](http://www.boost.org/)
* [Cairo libraries](http://www.cairographics.org/)
* [GTK library](http://www.gtk.org/)
* [GNU Scientific Library](http://www.gnu.org/software/gsl/)
* [zlib](http://www.zlib.net/)
* [SAMtools](http://samtools.sourceforge.net/) (for BAM/CRAM formatted input)

and also contains two submodules:

* SSP (https://github.com/rnakato/SSP)
* Coherent PDF (http://community.coherentpdf.com/)

# 3. Install

### 3.1. Docker image

We recommend to use the latest Docker image of DROMPAplus from DockerHub using:

    docker pull rnakato/ssp_drompa
    docker run -it --rm rnakato/ssp_drompa drompa+
    
For Singularity:

    singularity pull ssp_drompa.img docker://rnakato/ssp_drompa
    singularity exec ssp_drompa.img drompa+

### 3.2. Building from source

#### 3.2.1. Install required libraries
for Ubuntu:

    sudo apt install git build-essential libgtkmm-3.0-dev libboost-all-dev \
    libgsl-dev libz-dev samtools r-base

for CentOS:

    sudo yum -y install git gcc-c++ boost-devel zlib-devel gsl-devel gtkmm30-devel

#### 3.2.2. Install DROMPAplus
    git clone --recursive https://github.com/rnakato/DROMPAplus
    cd DROMPAplus
    make

If you get an installation error, make sure that all required libraries are successfully installed.

#### 3.2.3. Add the PATH environment variable
For example, if you downloaded DROMPA+ into $HOME directory, type:

    export PATH = $PATH:$HOME/DROMPAplus/bin
    export PATH = $PATH:$HOME/DROMPAplus/otherbins
    export PATH = $PATH:$HOME/DROMPAplus/submodules/cpdf/Linux-Intel-64bit

# 4. Usage
 See http://drompaplus.readthedocs.io/ for detail.

# 5. Reference
For DROMPAplus:
* Nakato R., Sakata T., Methods for ChIP-seq analysis: A practical workflow and advanced applications, ***submitted***.

For DROMPA:
* Nakato R, Shirahige K., Statistical Analysis and Quality Assessment of ChIP-seq Data with DROMPA, ***Methods in Molecular Biology***, 2018.
