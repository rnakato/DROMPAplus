# DROMPAplus sample scripts
# Visualization of H3K4me3, H3K27me3, H3K9me3 and Input samples for K562 cells from ENCODE project

# Get BAM files
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwHistone/wgEncodeUwHistoneK562H3k4me3StdAlnRep1.bam
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwHistone/wgEncodeUwHistoneK562H3k27me3StdAlnRep1.bam
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwHistone/wgEncodeUwHistoneK562H3k36me3StdAlnRep1.bam
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwHistone/wgEncodeUwHistoneK562InputStdAlnRep1.bam

# Get refFlat (gene annotation)
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz
gunzip refFlat.txt.gz

# parse2wig (make BigWig)
### When the current directory is DROMPAplus/scripts
gt=../src/SSP/data/genometable/genometable.hg19.txt
mptable=../src/SSP/data/mptable/mptable.UCSC.hg19.36mer.flen150.txt
parse2wig+ --nomodel -i wgEncodeUwHistoneK562H3k4me3StdAlnRep1.bam  -o H3K4me3  --gt $gt --mptable $mptable -n GR
parse2wig+ --nomodel -i wgEncodeUwHistoneK562H3k27me3StdAlnRep1.bam -o H3K27me3 --gt $gt --mptable $mptable -n GR
parse2wig+ --nomodel -i wgEncodeUwHistoneK562H3k36me3StdAlnRep1.bam -o H3K36me3 --gt $gt --mptable $mptable -n GR
parse2wig+ --nomodel -i wgEncodeUwHistoneK562InputStdAlnRep1.bam    -o Input    --gt $gt --mptable $mptable -n GR


# Make pdf
dir=parse2wigdir+
drompa+ PC_SHARP \
	-i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3,,,200 \
	-i $dir/H3K27me3.100.bw,$dir/Input.100.bw,H3K27me3,,,10 \
	-i $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3,,,10 \
	-o drompa-K562 -g refFlat.txt --gt $gt \
	--lpp 2 --showitag 1

# Overlayed pdf
dir=parse2wigdir+
drompa+ PC_SHARP \
	-i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3,,,100 \
	--ioverlay $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3,,,10 \
	-o drompa-K562-overlay -g refFlat.txt --gt $gt \
	--lpp 3 --chr 1 --alpha 0.6

# PROFILE
drompa+ PROFILE \
	-i $dir/H3K4me3.100.bw,$dir/Input.100.bw,H3K4me3 \
	-i $dir/H3K27me3.100.bw,$dir/Input.100.bw,H3K27me3 \
	-i $dir/H3K36me3.100.bw,$dir/Input.100.bw,H3K36me3 \
	-o profile-K562 -g refFlat.txt --gt $gt --ptype 2
