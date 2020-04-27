GENWIG: Generate wig data
-----------------------------------------

The **GENWIG** mode generates wig data of ChIP/Input enrichment and p-value distributions.
This mode is useful to analyze these distributions with Python and R by the user.

The following command outputs a PDF file (aroundtss.pdf) and a corresponding R script (aroundtss.R)::

    dir=parse2wigdir+
    drompa+ GENWIG \
	-i $dir/YST1019_Gal_60min.100.bw,$dir/YST1019_Gal_0min.100.bw,YST1019_Gal \
	-i $dir/YST1019_Raf_60min.100.bw,$dir/YST1019_Raf_0min.100.bw,YST1019_Raf \
	-i $dir/YST1053_Gal_60min.100.bw,$dir/YST1053_Gal_0min.100.bw,YST1053_Gal \
	-o drompa-yeast --gt genometable.sacCer3.txt --outputformat 3 --outputvalue 0

then "drompa-yeast.YST1019_Gal.enrich.100.bw", "drompa-yeast.YST1019_Raf.enrich.100.bw" and "drompa-yeast.YST1053_Gal.enrich.100.bw" are generated.

-  ``--outputvalue 0`` (default): ChIP/Input enrichment
-  ``--outputvalue 1``: P-value (ChIP internal)
-  ``--outputvalue 2``: P-value (ChIP/Input enrichment)

-  ``--outputformat 0``: compresed wig (.wig.gz)
-  ``--outputformat 1``: uncompresed wig (.wig)
-  ``--outputformat 2``: bedGraph (.bedGraph)
-  ``--outputformat 3`` (default): bigWig (.bw)
