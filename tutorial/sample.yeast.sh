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

dir=parse2wigdir+
gene=../data/S_cerevisiae/SGD_features.tab
gt=../data/genometable/genometable.sacCer3.txt
drompa+ PC_ENRICH \
	-i $dir/YST1019_Gal_60min.100.bw,$dir/YST1019_Gal_0min.100.bw,YST1019_Gal,,,200 \
	-i $dir/YST1019_Raf_60min.100.bw,$dir/YST1019_Raf_0min.100.bw,YST1019_Raf,,,200 \
	-i $dir/YST1053_Gal_60min.100.bw,$dir/YST1053_Gal_0min.100.bw,YST1053_Gal,,,200 \
	-o drompa-yeast --gt $gt -g $gene --gftype 2 \
	--scale_ratio 1 --ls 200 --sm 10 --lpp 3

# parse2wig (make BigWig)
gt=../data/genometable/genometable.sacCer3.txt
mptable=../data/mptable/mptable.UCSC.sacCer3.50mer.flen150.txt
for cell in YST1019_Gal YST1019_Raf YST1053_Gal; do
    for min in 0min 60min; do
	cram=${cell}_${min}-n2-k1.sort.cram
	parse2wig+ -i $cram -o ${cell}_${min} --gt $gt --mptable $mptable -n GR
	# parse statsfile
#	parsestats4DROMPAplus.pl parse2wigdir+/${cell}_${min}.100.tsv
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
	-o drompa-yeast-detail --gt $gt --ars $ars \
	--scale_ratio 1 --bn 5 --ystep 10 --ls 200 --sm 10 --lpp 3

# Highlight peaks
dir=parse2wigdir+
gt=../data/genometable/genometable.sacCer3.txt
ars=../data/S_cerevisiae/ARS-oriDB_scer.txt
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

# logratio
dir=parse2wigdir+
gt=../data/genometable/genometable.sacCer3.txt
ars=../data/S_cerevisiae/ARS-oriDB_scer.txt
drompa+ PC_ENRICH \
	-i $dir/YST1019_Gal_60min.100.bw,$dir/YST1019_Gal_0min.100.bw,YST1019_Gal,,,200 \
	-i $dir/YST1019_Raf_60min.100.bw,$dir/YST1019_Raf_0min.100.bw,YST1019_Raf,,,200 \
	-i $dir/YST1053_Gal_60min.100.bw,$dir/YST1053_Gal_0min.100.bw,YST1053_Gal,,,200 \
	-o drompa-yeast-log2ratio --gt $gt --ars $ars \
	--showratio 2 --scale_ratio 2 --ls 200 --sm 10 --bn 4 --lpp 3 --chr I

drompa+ PC_ENRICH \
	-i $dir/YST1019_Gal_60min.100.bw,$dir/YST1019_Gal_0min.100.bw,YST1019_Gal,,,200 \
	-i $dir/YST1019_Raf_60min.100.bw,$dir/YST1019_Raf_0min.100.bw,YST1019_Raf,,,200 \
	-i $dir/YST1053_Gal_60min.100.bw,$dir/YST1053_Gal_0min.100.bw,YST1053_Gal,,,200 \
	-o drompa-yeast-log2ratio2 --gt $gt --ars $ars \
	--showratio 2 --scale_ratio 2 --ls 200 --sm 10 --bn 4 --lpp 3 --chr I --callpeak


drompa+ PC_ENRICH \
	-i $dir/YST1019_Gal_60min.100.bw,$dir/YST1019_Gal_0min.100.bw,YST1019_Gal,,,200 \
	-i $dir/YST1019_Raf_60min.100.bw,$dir/YST1019_Raf_0min.100.bw,YST1019_Raf,,,200 \
	-i $dir/YST1053_Gal_60min.100.bw,$dir/YST1053_Gal_0min.100.bw,YST1053_Gal,,,200 \
	-o drompa-yeast-log10ratio --gt $gt --ars $ars \
	--showratio 2 --scale_ratio 10 --ls 200 --sm 10 --bn 4 --lpp 3 --chr I

# GENWIG
dir=parse2wigdir+
gene=../data/S_cerevisiae/SGD_features.tab
gt=../data/genometable/genometable.sacCer3.txt
drompa+ PC_ENRICH \
	-i $dir/YST1019_Gal_60min.100.bw,$dir/YST1019_Gal_0min.100.bw,YST1019_Gal,,,200 \
	-i $dir/YST1019_Raf_60min.100.bw,$dir/YST1019_Raf_0min.100.bw,YST1019_Raf,,,200 \
	-i $dir/YST1053_Gal_60min.100.bw,$dir/YST1053_Gal_0min.100.bw,YST1053_Gal,,,200 \
	-o drompa-yeast --gt $gt -g $gene --gftype 2  --showpenrich 1 --showpinter 1 \
	--scale_ratio 5 --ls 200 --sm 10 --lpp 3

drompa+ GENWIG \
	-i $dir/YST1019_Gal_60min.100.bw,$dir/YST1019_Gal_0min.100.bw,YST1019_Gal \
	-i $dir/YST1019_Raf_60min.100.bw,$dir/YST1019_Raf_0min.100.bw,YST1019_Raf \
	-i $dir/YST1053_Gal_60min.100.bw,$dir/YST1053_Gal_0min.100.bw,YST1053_Gal \
	-o drompa-yeast --gt $gt --outputformat 3 --outputvalue 0
drompa+ GENWIG \
	-i $dir/YST1019_Gal_60min.100.bw,$dir/YST1019_Gal_0min.100.bw,YST1019_Gal \
	-i $dir/YST1019_Raf_60min.100.bw,$dir/YST1019_Raf_0min.100.bw,YST1019_Raf \
	-i $dir/YST1053_Gal_60min.100.bw,$dir/YST1053_Gal_0min.100.bw,YST1053_Gal \
	-o drompa-yeast --gt $gt --outputformat 3 --outputvalue 1
drompa+ GENWIG \
	-i $dir/YST1019_Gal_60min.100.bw,$dir/YST1019_Gal_0min.100.bw,YST1019_Gal \
	-i $dir/YST1019_Raf_60min.100.bw,$dir/YST1019_Raf_0min.100.bw,YST1019_Raf \
	-i $dir/YST1053_Gal_60min.100.bw,$dir/YST1053_Gal_0min.100.bw,YST1053_Gal \
	-o drompa-yeast --gt $gt --outputformat 3 --outputvalue 2

s=""
for sample in YST1019_Gal YST1019_Raf YST1053_Gal; do
    for str in pinter penrich enrich; do
	s="$s -i drompa-yeast.$sample.$str.100.bw"
    done
done

drompa+ PC_SHARP $s -o drompa-yeast_genwig --gt $gt -g $gene --gftype 2 \
	--scale_tag 5 --ls 200 --sm 10 --lpp 3
