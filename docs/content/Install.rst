Installation
================

Docker image
---------------------------------

We recommend using the latest Docker image of DROMPAplus from DockerHub as follows:

.. code-block:: bash

   docker pull rnakato/ssp_drompa
   docker run -it --rm rnakato/ssp_drompa drompa+

For Singularity:

.. code-block:: bash

   singularity pull ssp_drompa.img docker://rnakato/ssp_drompa
   singularity exec ssp_drompa.img drompa+


Building from source
---------------------------------

Requirements
++++++++++++++++++++++++++++++

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

Permanently set the path to the DROMPAplus download directory by updating your **~/.bashrc** file. For example, if you downloaded DROMPAplus into the **$HOME** directory, add the following lines to **~/.bashrc**:

.. code-block:: bash

   export PATH = $PATH:$HOME/DROMPAplus/bin:$HOME/DROMPAplus/otherbins
