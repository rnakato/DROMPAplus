Installation
================

Docker image
---------------------------------

We recommend using the latest Docker image of DROMPAplus from DockerHub as follows:

.. code-block:: bash

   docker pull rnakato/ssp_drompa
   docker run -it --rm rnakato/ssp_drompa drompa+

.. note::

    If you are using the docker image to execute parse2wig+, it is necessary to mount the directory by ``-v`` option to access the input files as follows::

        docker run -it --rm -v $(pwd):/mnt rnakato/ssp_drompa parse2wig+ \
               -i /mnt/ChIP.bam -o /mnt/ChIP --gt /mnt/genometable.txt

    This command mounts the current directory to /mnt directory in the container. 
    Please see also `the document of Docker <https://docs.docker.com/storage/volumes/>`_.

For Singularity:

.. code-block:: bash

   singularity build ssp_drompa.img docker://rnakato/ssp_drompa
   singularity exec ssp_drompa.img drompa+

Singularity mounts the current directory automatically. If you access the files in the other directory,
please mount by ``--bind`` option, for instance::

    singularity exec --bind /work ssp_drompa.img drompa+

This command mounts ``/work`` directory.


Building from source
---------------------------------

Requirements
++++++++++++++++++++++++++++++

DROMPAplus requires the following programs and libraries:

- Boost C++ library (>1.53.0)
- Cairo libraries
- GTK library
- GNU Scientific Library (>1.15)
- zlib (>1.2.7)
- CMake (>2.8)
- HTSlib (1.10.2) (for SAM/BAM/CRAM formatted input)

and also contains two submodules:

- SSP
- Coherent PDF

Install required libraries
++++++++++++++++++++++++++++++

On Ubuntu:

.. code-block:: bash

    sudo apt install git build-essential libgtkmm-3.0-dev libboost-all-dev \
    libgsl-dev libz-dev libbz2-dev libgsl-dev libcurl4-gnutls-dev cmake

On CentOS:

.. code-block:: bash

    sudo yum -y install git gcc-c++ clang boost-devel zlib-devel gsl-devel gtkmm30-devel bzip2-devel cmake

On Mac:

.. code-block:: bash

    brew install gsl gtk gtkmm cairo pkgconfig curl xz zlib boost cmake

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
