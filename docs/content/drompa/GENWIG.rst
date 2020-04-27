GENWIG: Generate wig data
-----------------------------------------

The **GENWIG** mode generates wig data of ChIP/Input enrichment and p-value distributions.
This mode is useful to analyze these distributions with Python and R by the user.

The following command outputs a PDF file (aroundtss.pdf) and a corresponding R script (aroundtss.R)::

    dir=parse2wigdir+
    drompa+ GENWIG \
       -i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3 \
       -i $dir/H3K27me3.100.bw,$dir/Input.100.bw,H3K27me3 \
       -i $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3 \
       -o aroundtss -g refFlat.txt --gt genometable.txt


-  ``--outputvalue 0``: ChIP/Input enrichment
-  ``--outputvalue 1``: P-value (ChIP internal)
-  ``--outputvalue 2``: P-value (ChIP/Input enrichment)

