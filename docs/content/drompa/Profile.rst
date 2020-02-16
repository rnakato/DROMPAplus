PROFILE: Aggregation plot (using R)
====================================

**PROFILE** mode makes an aggregation plot by supplying:

-  ``--ptype 0``: around transcription start sites (TSS)
-  ``--ptype 1``: around transcription termination sites (TTS)
-  ``--ptype 2``: divide genes into 100 subregions 
-  ``--ptype 3``: around peaks

In **PROFILE** mode, short genes (< 1kbp) are ignored.

The following command outputs a PDF file (aroundtss.pdf) and a corresponding R script (aroundtss.R)::

    dir=parse2wigdir+
    drompa+ PROFILE --ptype 2 \
       -i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3 \
       -i $dir/H3K27me3.100.bw,$dir/Input.100.bw,H3K27me3 \
       -i $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3 \
       -o aroundtss -g refFlat.txt --gt genometable.txt

.. image:: img/aroundtss.jpg
   :width: 400px
   :align: center


The following command outputs a PDF file (aroundgene.pdf) and a corresponding R script (aroundgene.R)::

    dir=parse2wigdir+
    drompa+ PROFILE --ptype 2 \
       -i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3 \
       -i $dir/H3K27me3.100.bw,$dir/Input.100.bw,H3K27me3 \
       -i $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3 \
       -o aroundgene -g refFlat.txt --gt genometable.txt

.. image:: img/aroundgene.jpg
   :width: 400px
   :align: center

The shaded regions indicate the 95% confidence interval.

Modify the parameter for plot
++++++++++++++++++++++++++++++++

If you want to modify the parameters of the plot (e.g., max value of y-axis),
change the parameters in the generated R script and remake the PDF file as follows::

    $ R --vanilla < aroundtss.R
    $ R --vanilla < aroundgene.R


The generated tsv files (aroundtss.tsv and aroundgene.tsv) describe the read number of each site and can be used to further analysis such as hierarchical Clustering.

Normalize read number for the specified regions
++++++++++++++++++++++++++++++++++++++++++++++++++++++

The averaged read number in the background regions sometimes varies among samples. The option ``--ntype 1`` normalizes the read number so that the number of reads mapped within the specified regions::

    dir=parse2wigdir+
    drompa+ PROFILE --ptype 2 --ntype 1 \
       -i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3 \
       -i $dir/H3K27me3.100.bw,$dir/Input.100.bw,H3K27me3 \
       -i $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3 \
       -o aroundgene2 -g refFlat.txt --gt genometable.txt

.. image:: img/aroundgene2.jpg
   :width: 400px
   :align: center

``--stype 1`` option plot the averaged ChIP/Input enrichment::

    dir=parse2wigdir+
    drompa+ PROFILE --ptype 2 --stype 1 \
       -i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3 \
       -i $dir/H3K27me3.100.bw,$dir/Input.100.bw,H3K27me3 \
       -i $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3 \
       -o aroundgene.enrich -g refFlat.txt --gt genometable.txt

.. image:: img/aroundgene.enrich.jpg
   :width: 400px
   :align: center
