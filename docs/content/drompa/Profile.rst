PROFILE: Aggregation plot (using R)
-----------------------------------------

The **PROFILE** mode makes an aggregation plot by supplying:

-  ``--ptype 0``: around transcription start sites (TSS)
-  ``--ptype 1``: around transcription termination sites (TTS)
-  ``--ptype 2``: divide genes into 100 subregions
-  ``--ptype 3``: around peaks

In **PROFILE** mode, short genes (< 1kbp) are ignored.

The following command outputs a PDF file (aroundtss.pdf) and a corresponding R script (aroundtss.R)::

    dir=parse2wigdir+
    drompa+ PROFILE --ptype 0 \
       -i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3 \
       -i $dir/H3K27me3.100.bw,$dir/Input.100.bw,H3K27me3 \
       -i $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3 \
       -o aroundtss -g refFlat.txt --gt genometable.txt

.. figure:: img/aroundtss.jpg
   :width: 400px
   :align: center
   :alt: Alternate

   The read density as a function of the distance from the TSS.


The following command outputs a PDF file (aroundgene.pdf) and a corresponding R script (aroundgene.R)::

    dir=parse2wigdir+
    drompa+ PROFILE --ptype 2 \
       -i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3 \
       -i $dir/H3K27me3.100.bw,$dir/Input.100.bw,H3K27me3 \
       -i $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3 \
       -o aroundgene -g refFlat.txt --gt genometable.txt

.. figure:: img/aroundgene.jpg
   :width: 400px
   :align: center
   :alt: Alternate

   The read density as a function of the percentage of the gene length from the TSS.

The shaded regions indicate the 95% confidence interval.


Modifying the plot parameters
++++++++++++++++++++++++++++++++

To modify the plot parameters (e.g., max value of y-axis), change the parameters in the generated R script and remake the PDF file as follows::

    $ R --vanilla < aroundtss.R
    $ R --vanilla < aroundgene.R

The generated TSV files (aroundtss.tsv and aroundgene.tsv) describe the read number of each site and can be used for further analysis, e.g., hierarchical clustering.


Normalizing the read number for the specified regions
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The averaged read number in the background regions sometimes varies among samples. The option ``--ntype 1`` normalizes the read number::

    dir=parse2wigdir+
    drompa+ PROFILE --ptype 2 --ntype 1 \
       -i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3 \
       -i $dir/H3K27me3.100.bw,$dir/Input.100.bw,H3K27me3 \
       -i $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3 \
       -o aroundgene2 -g refFlat.txt --gt genometable.txt

.. figure:: img/aroundgene2.jpg
   :width: 400px
   :align: center
   :alt: Alternate

   The read enrichment with ``--ntype 1`` option.

Using the ``--stype 1`` option, drompa+ plots the averaged ChIP/Input enrichment (note that this is the averaged enrichment for all sites, not the enrichment of the averaged read density.)::

    dir=parse2wigdir+
    drompa+ PROFILE --ptype 2 --stype 1 \
       -i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3 \
       -i $dir/H3K27me3.100.bw,$dir/Input.100.bw,H3K27me3 \
       -i $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3 \
       -o aroundgene.enrich -g refFlat.txt --gt genometable.txt

.. figure:: img/aroundgene.enrich.jpg
   :width: 400px
   :align: center
   :alt: Alternate

   The read enrichment as a function of the percentage of the gene length from the TSS.
