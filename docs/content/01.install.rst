Installation
================

Requirements
-------------------------------

DROMPAplus requires the following programs and libraries:

- Boost C++ library
- Cairo libraries
- GTK library
- GNU Scientific Library
- zlib
- SAMtools

and also contains two submodules:

- SSP
- Coherent PDF


Docker image
---------------------------------

We recommend to use the latest Docker image of DROMPAplus from DockerHub using:

.. code-block:: bash

   docker pull rnakato/ssp_drompa
   docker run -it --rm rnakato/ssp_drompa drompa+

For Singularity:

.. code-block:: bash

   singularity pull ssp_drompa.img docker://rnakato/ssp_drompa
   singularity exec ssp_drompa.img drompa+

Building from source
---------------------------------

Install required libraries
++++++++++++++++++++++++++++++

For Ubuntu:

.. code-block:: bash

    sudo apt install git build-essential libgtkmm-3.0-dev libboost-all-dev libgsl-dev libz-dev samtools

For CentOS:

.. code-block:: bash

    sudo yum -y install git gcc-c++ boost-devel zlib-devel gsl-devel gtkmm30-devel

Install DROMPAplus
+++++++++++++++++++++++++

.. code-block:: bash

   git clone --recursive https://github.com/rnakato/DROMPAplus
   cd DROMPAplus
   make

Add the PATH
+++++++++++++++++++++++++

If you downloaded DROMPAplus into $HOME/ directory, type:

.. code-block:: bash

   export PATH = $PATH:$HOME/DROMPAplus/bin
   export PATH = $PATH:$HOME/DROMPAplus/otherbins
   export PATH = $PATH:$HOME/DROMPAplus/submodules/cpdf/Linux-Intel-64bit

Add these lines to ~/.bashrc to permanently set PATH.
