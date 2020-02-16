GV: Whole-chromosome view
====================================

The **GV** mode shows a chromosome-wide overview of the ChIP-seq data. Type::

    $ drompa+ GV $s1 $s2 $s3 $s4 -p ChIPseq-wholegenome -gt genometable.txt \verb|\|
    $ --GC GCcontents/   --gcsize 500000 \  # GC contents directory and window size (500 kbp)
    $ --GD gene_density/ --gdsize 500000    # Gene-density directory and window size (500 kbp)

In **GV** mode, ChIP/Input enrichment for 100k-bp bin is shown by default. 
**GV** mode does not perform the significance test but simply
highlights the bins containing ChIP/Input enrichments above the middle of
y axis (specified with ``-scale_ratio``) in red.
