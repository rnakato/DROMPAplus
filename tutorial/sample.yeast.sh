#!/bin/bash
# DROMPAplus sample scripts
# Replication analysis of S. cerevisiae

### Download cram files from the URLs below:
# YST1019_Gal_0min-n2-k1.sort.cram <https://drive.google.com/open?id=1-T5eq6ZgY5kUmRDlbIfFYNMBuwVCQKkc>
# YST1019_Gal_60min-n2-k1.sort.cram <https://drive.google.com/open?id=1-ZabbvJdAYxdaH5MD8s-sMKHkNiEsTQu>
# YST1019_Raf_0min-n2-k1.sort.cram <https://drive.google.com/open?id=1-aCu7DSfEG8EIfYCdwPwsmzjsLXfelY6>
# YST1019_Raf_60min-n2-k1.sort.cram <https://drive.google.com/open?id=1-bzJi5MNjVee5J3chQ7kbhlHTMgqqFut>
# YST1053_Gal_0min-n2-k1.sort.cram <https://drive.google.com/open?id=1-c-FHugX47qUHx50KWshhYS8V2Ns0n_c>
# YST1053_Gal_60min-n2-k1.sort.cram <https://drive.google.com/open?id=1-i0Vbpbd5nLfiTtieD3LLACCf5mVYUiG>

# parse2wig (make BigWig)
gt=../data/genometable/genometable.sacCer3.txt
mptable=../data/mptable/mptable.UCSC.sacCer3.50mer.flen150.txt
for cell in YST1019_Gal YST1019_Raf YST1053_Gal; do
for min in 0min 60min; do
cram=${cell}_${min}-n2-k1.sort.cram
parse2wig+ -i $cram  -o ${cell}_${min} --gt $gt --mptable $mptable -n GR
done
done

# Make pdf
dir=parse2wigdir+
gene=../data/S_cerevisiae/SGD_features.tab
gt=../data/genometable/genometable.sacCer3.txt
drompa+ PC_ENRICH \
	-i $dir/YST1019_Gal_60min.100.bw,$dir/YST1019_Gal_0min.100.bw,YST1019_Gal,,,200 \
	-i $dir/YST1019_Raf_60min.100.bw,$dir/YST1019_Raf_0min.100.bw,YST1019_Raf,,,200 \
	-i $dir/YST1053_Gal_60min.100.bw,$dir/YST1053_Gal_0min.100.bw,YST1053_Gal,,,200 \
	-o drompa-yeast --gt $gt -g $gene --gftype 2 \
	--scale_ratio 1 --ls 200 --sm 10 --lpp 3

ars=../data/S_cerevisiae/ARS-oriDB_scer.txt
drompa+ PC_ENRICH \
	-i $dir/YST1019_Gal_60min.100.bw,$dir/YST1019_Gal_0min.100.bw,YST1019_Gal,,,200 \
	-i $dir/YST1019_Raf_60min.100.bw,$dir/YST1019_Raf_0min.100.bw,YST1019_Raf,,,200 \
	-i $dir/YST1053_Gal_60min.100.bw,$dir/YST1053_Gal_0min.100.bw,YST1053_Gal,,,200 \
	-o drompa-yeast-ARS --gt $gt --ars $ars \
	--scale_ratio 1 --ls 200 --sm 10 --lpp 3

drompa+ PC_ENRICH \
	-i $dir/YST1019_Gal_60min.100.bw,$dir/YST1019_Gal_0min.100.bw,YST1019_Gal,,,200 \
	-i $dir/YST1019_Raf_60min.100.bw,$dir/YST1019_Raf_0min.100.bw,YST1019_Raf,,,200 \
	-i $dir/YST1053_Gal_60min.100.bw,$dir/YST1053_Gal_0min.100.bw,YST1053_Gal,,,200 \
	--callpeak \
	-o drompa-yeast-ARS-peak --gt $gt --ars $ars \
	--scale_ratio 1 --ls 200 --sm 10 --lpp 3

drompa+ PC_ENRICH \
	-i $dir/YST1019_Gal_60min.100.bw,$dir/YST1019_Gal_0min.100.bw,YST1019_Gal,,,200 \
	-i $dir/YST1019_Raf_60min.100.bw,$dir/YST1019_Raf_0min.100.bw,YST1019_Raf,,,200 \
	-i $dir/YST1053_Gal_60min.100.bw,$dir/YST1053_Gal_0min.100.bw,YST1053_Gal,,,200 \
	--callpeak --ethre 1.5 \
	-o drompa-yeast-ARS-peak2 --gt $gt --ars $ars \
	--scale_ratio 1 --ls 200 --sm 10 --lpp 3



exit

drompa+ PC_ENRICH \
	-i $dir/YST1019_Gal_60min.100.bw,$dir/YST1019_Gal_0min.100.bw,YST1019_Gal,,,200 \
	-i $dir/YST1019_Raf_60min.100.bw,$dir/YST1019_Raf_0min.100.bw,YST1019_Raf,,,200 \
	-i $dir/YST1053_Gal_60min.100.bw,$dir/YST1053_Gal_0min.100.bw,YST1053_Gal,,,200 \
	-o drompa-yeast-logratio --gt $gt --ars $ars \
	--showratio 2 --scale_ratio 1 --ls 200 --sm 10 --lpp 3

drompa+ PC_ENRICH \
	-i $dir/YST1019_Gal_60min.100.bw,$dir/YST1019_Gal_0min.100.bw,YST1019_Gal,,,200 \
	-i $dir/YST1019_Raf_60min.100.bw,$dir/YST1019_Raf_0min.100.bw,YST1019_Raf,,,200 \
	-i $dir/YST1053_Gal_60min.100.bw,$dir/YST1053_Gal_0min.100.bw,YST1053_Gal,,,200 \
	-o drompa-yeast-detail --gt $gt --ars $ars \
	--scale_ratio 1 --bn 5 --ystep 10 --ls 200 --sm 10 --lpp 3
