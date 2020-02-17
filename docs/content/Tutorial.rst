Tutorial
=================


Sample scripts
--------------------

A sample script file for tutorial is in "tutorial" derectory.
``sample.human.sh`` includes the drompa+ **PC_SHARP**, **GV** and **PROFILE** commands for histone modification data of human K562 cells.

Data
----------------

The "data" directory contains several types of data:

- Genometable files for several species. To generate genometable files for other species (abd genome builds), use ``makegenometable.pl`` in the "scripts" directory as follows::

        makegenometable.pl genome.fa > genometable.txt

- Ideogram files for visualization of chromosomes in **GV** command.
- Mptable files that describe the mappable bases of each chromsome.
