MULTICI: Generate matrix of averaged read density
----------------------------------------------------

The **MULTICI** generates a matrix file that describes the averaged read density in specified BED site::

    dir=parse2wigdir+
    gt=../data/genometable/genometable.sacCer3.txt
    drompa+ MULTICI \
	-i $dir/YST1019_Gal_60min.100.bw,$dir/YST1019_Gal_0min.100.bw,YST1019_Gal \
	-i $dir/YST1019_Raf_60min.100.bw,$dir/YST1019_Raf_0min.100.bw,YST1019_Raf \
	-i $dir/YST1053_Gal_60min.100.bw,$dir/YST1053_Gal_0min.100.bw,YST1053_Gal \
	-o drompa-multici --gt $gt --bed test.bed

The output file (drompa-multici.ChIPread.tsv) contains the averaged ChIP-read intensity per bin size for each BED site.
Add ``--stype 1`` to generate averaged ChIP/Input enrichment matrix.
