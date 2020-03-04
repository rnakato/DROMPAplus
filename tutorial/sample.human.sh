#!/bin/bash
# DROMPAplus sample scripts
# Visualization of H3K4me3, H3K27me3, H3K9me3 and Input samples for K562 cells from ENCODE project

# Get BAM files
wget -nc http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwHistone/wgEncodeUwHistoneK562H3k4me3StdAlnRep1.bam
wget -nc http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwHistone/wgEncodeUwHistoneK562H3k27me3StdAlnRep1.bam
wget -nc http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwHistone/wgEncodeUwHistoneK562H3k36me3StdAlnRep1.bam
wget -nc http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwHistone/wgEncodeUwHistoneK562InputStdAlnRep1.bam

# parse2wig (make BigWig)
gt=../data/genometable/genometable.hg19.txt
mptable=../data/mptable/mptable.UCSC.hg19.36mer.flen150.txt
parse2wig+ -i wgEncodeUwHistoneK562H3k4me3StdAlnRep1.bam  -o H3K4me3  --gt $gt --mptable $mptable -n GR
parse2wig+ -i wgEncodeUwHistoneK562H3k27me3StdAlnRep1.bam -o H3K27me3 --gt $gt --mptable $mptable -n GR
parse2wig+ -i wgEncodeUwHistoneK562H3k36me3StdAlnRep1.bam -o H3K36me3 --gt $gt --mptable $mptable -n GR
parse2wig+ -i wgEncodeUwHistoneK562InputStdAlnRep1.bam    -o Input    --gt $gt --mptable $mptable -n GR

# Make pdf
gt=../data/genometable/genometable.hg19.txt
dir=parse2wigdir+
gene=../data/geneannotation/Homo_sapiens.GRCh37.87.chr.gene.name.refFlat
drompa+ PC_SHARP \
	-i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3 \
	-i $dir/H3K27me3.100.bw,$dir/Input.100.bw,H3K27me3 \
	-i $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3 \
	-o drompa1 -g $gene --gt $gt \
	--lpp 2 --showitag 2 \
	--chr 1

drompa+ PC_SHARP \
	-i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3 \
	-i $dir/H3K27me3.100.bw,$dir/Input.100.bw,H3K27me3 \
	-i $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3 \
	-o drompa2 -g $gene --gt $gt \
	--lpp 2 --showitag 2 --scale_tag 10 \
	--chr 1

drompa+ PC_SHARP \
	-i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3,,,100 \
	-i $dir/H3K27me3.100.bw,$dir/Input.100.bw,H3K27me3,,,10 \
	-i $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3,,,10 \
	-o drompa3 -g $gene --gt $gt \
	--lpp 2 --showitag 2 \
	--chr 1

drompa+ PC_SHARP \
	-i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3 \
	-i $dir/H3K27me3.100.bw,$dir/Input.100.bw,H3K27me3 \
	-i $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3 \
	-o drompa_pvalue -g $gene --gt $gt \
        --showratio 1 --showpinter 1 --showpenrich 1 \
        --scale_ratio 3 --scale_pvalue 3 \
        --chr 1

# Visualize specific regions
gt=../data/genometable/genometable.hg19.txt
dir=parse2wigdir+
gene=../data/geneannotation/Homo_sapiens.GRCh37.87.chr.gene.name.refFlat

### make BED file "HOXA.txt"
echo -e "chr7\t27100000\t27280000" > HOXA.txt
#cat HOXA.txt

### supply "HOXA.txt" with -r option
drompa+ PC_SHARP \
	-i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3,,,100 \
	-i $dir/H3K27me3.100.bw,$dir/Input.100.bw,H3K27me3,,,10 \
	-i $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3,,,10 \
	-o drompa_HOXA -g $gene --gt $gt \
	--lpp 2 --showitag 2 \
        -r HOXA.txt

# Overlayed pdf
gt=../data/genometable/genometable.hg19.txt
dir=parse2wigdir+
gene=Homo_sapiens.GRCh37.87.chr.gene.name.refFlat
drompa+ PC_SHARP \
	-i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3,,,100 \
	--ioverlay $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3,,,10 \
	-o drompa-overlay -g $gene --gt $gt \
	--lpp 3 --chr 1 \
	--alpha 0.6 \

# with chromatin loops
gt=../data/genometable/genometable.hg19.txt
dir=parse2wigdir+
gene=../data/geneannotation/Homo_sapiens.GRCh37.87.chr.gene.name.refFlat
drompa+ PC_SHARP \
	-i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3,,,200 \
	-i $dir/H3K27me3.100.bw,$dir/Input.100.bw,H3K27me3,,,10 \
	-i $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3,,,10 \
	-o drompa_loops -g $gene --gt $gt \
	--inter interactions.all.mango,ChIA-PET,mango \
	--inter HICCUPS_looplist.txt,Hi-C,hiccups \
	--lpp 2 --chr 20 --ls 5000

# GV
gt=../data/genometable/genometable.hg19.txt
mptable=../data/mptable/mptable.UCSC.hg19.36mer.flen150.txt
parse2wig+ -i wgEncodeUwHistoneK562H3k4me3StdAlnRep1.bam  -o H3K4me3  --gt $gt --mptable $mptable -n GR --binsize 100000
parse2wig+ -i wgEncodeUwHistoneK562H3k27me3StdAlnRep1.bam -o H3K27me3 --gt $gt --mptable $mptable -n GR --binsize 100000
parse2wig+ -i wgEncodeUwHistoneK562H3k36me3StdAlnRep1.bam -o H3K36me3 --gt $gt --mptable $mptable -n GR --binsize 100000
parse2wig+ -i wgEncodeUwHistoneK562InputStdAlnRep1.bam    -o Input    --gt $gt --mptable $mptable -n GR --binsize 100000

dir=parse2wigdir+
drompa+ GV \
	-i $dir/H3K4me3.100000.bw,$dir/Input.100000.bw,H3K4me3   \
	-i $dir/H3K27me3.100000.bw,$dir/Input.100000.bw,H3K27me3 \
	-i $dir/H3K36me3.100000.bw,$dir/Input.100000.bw,H3K36me3 \
	-o drompaGV-K562 --gt $gt

drompa+ GV \
	-i $dir/H3K4me3.100000.bw,$dir/Input.100000.bw,H3K4me3   \
	-i $dir/H3K27me3.100000.bw,$dir/Input.100000.bw,H3K27me3 \
	-i $dir/H3K36me3.100000.bw,$dir/Input.100000.bw,H3K36me3 \
	-o drompaGV-K562_2 --gt $gt \
	--GC GCcontents --gcsize 500000 \
	--GD genedensity --gdsize 500000 \
	--ideogram ../data/ideogram/hg19.tsv

# PROFILE
dir=parse2wigdir+
gt=../data/genometable/genometable.hg19.txt
gene=../data/geneannotation/Homo_sapiens.GRCh37.87.chr.gene.name.refFlat
drompa+ PROFILE \
	-i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3 \
	-i $dir/H3K27me3.100.bw,$dir/Input.100.bw,H3K27me3 \
	-i $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3 \
	-o profile-aroundTSS -g $gene --gt $gt --ptype 0

drompa+ PROFILE \
	-i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3 \
	-i $dir/H3K27me3.100.bw,$dir/Input.100.bw,H3K27me3 \
	-i $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3 \
	-o profile-aroundGene -g $gene --gt $gt --ptype 2

drompa+ PROFILE \
	-i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3 \
	-i $dir/H3K27me3.100.bw,$dir/Input.100.bw,H3K27me3 \
	-i $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3 \
	-o profile-aroundGene -g $gene --gt $gt --ptype 2 --stype 1

#HEATMAP
../otherbins/drompa.heatmap.py \
    -o heatmap-aroundTSS \
    profile-aroundTSS.ChIPread.H3K4me3.tsv \
    profile-aroundTSS.ChIPread.H3K27me3.tsv \
    profile-aroundTSS.ChIPread.H3K36me3.tsv

../otherbins/drompa.heatmap.py \
    -o heatmap-aroundGene \
    profile-aroundGene.ChIPread.H3K4me3.tsv \
    profile-aroundGene.ChIPread.H3K27me3.tsv \
    profile-aroundGene.ChIPread.H3K36me3.tsv
