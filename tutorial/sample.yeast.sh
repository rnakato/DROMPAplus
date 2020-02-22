#!/bin/bash
# DROMPAplus sample scripts
# Replication analysis of S. cerevisiae

### Download cram files from the URLs below:
# https://drive.google.com/open?id=1f5H-umPgjzmDYLVHBdlIMWQXtt91S6Fc
# https://drive.google.com/open?id=1f991hi-V9ITAWeUzr3AkqyRY8CbFTBXO
# https://drive.google.com/open?id=1fDUgeo8IvhI-fikQv_PDCiDmNc9cBfWT
# https://drive.google.com/open?id=1fGjY5nlhjePk_0TXKC5bLOzmdrer7IHh
# https://drive.google.com/open?id=1fJh-f3CDNc7wAkuN79rndjptVq_qYU6k
# https://drive.google.com/open?id=1fNRt1uvA1CQrIfb9NSJq3hYL9-GGJnXv

# parse2wig (make BigWig)
gt=../data/genometable/genometable.sacCer3.txt
mptable=../data/mptable/mptable.UCSC.sacCer3.50mer.flen150.txt
for cell in YST1019_Gal YST1019_Raf YST1053_Gal; do
for min in 0min 60min; do
cram=${cell}_${min}-n2-k1.sort.cram
#parse2wig+ -i $cram  -o ${cell}_${min} --gt $gt --mptable $mptable -n GR
done
done

# Make pdf
dir=parse2wigdir+
gene=../data/S_cerevisiae/SGD_features.tab

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

exit


# Overlayed pdf
dir=parse2wigdir+
drompa+ PC_SHARP \
	-i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3,,,100 \
	--ioverlay $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3,,,10 \
	-o drompa-K562-overlay -g refFlat.txt --gt $gt \
	--lpp 3 --chr 1 \
	--alpha 0.6 \

# GV
gt=../data/genometable/genometable.hg19.txt
mptable=../data/mptable/mptable.UCSC.hg19.36mer.flen150.txt
parse2wig+ -i wgEncodeUwHistoneK562H3k4me3StdAlnRep1.bam  -o H3K4me3  --gt $gt --mptable $mptable -n GR --binsize 100000
parse2wig+ -i wgEncodeUwHistoneK562H3k27me3StdAlnRep1.bam -o H3K27me3 --gt $gt --mptable $mptable -n GR --binsize 100000
parse2wig+ -i wgEncodeUwHistoneK562H3k36me3StdAlnRep1.bam -o H3K36me3 --gt $gt --mptable $mptable -n GR --binsize 100000
parse2wig+ -i wgEncodeUwHistoneK562InputStdAlnRep1.bam    -o Input    --gt $gt --mptable $mptable -n GR --binsize 100000



# PROFILE
drompa+ PROFILE \
	-i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3 \
	-i $dir/H3K27me3.100.bw,$dir/Input.100.bw,H3K27me3 \
	-i $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3 \
	-o profile-K562-aroundTSS -g refFlat.txt --gt $gt --ptype 0

drompa+ PROFILE \
	-i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3 \
	-i $dir/H3K27me3.100.bw,$dir/Input.100.bw,H3K27me3 \
	-i $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3 \
	-o profile-K562-aroundGene -g refFlat.txt --gt $gt --ptype 2
