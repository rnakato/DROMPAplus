gene=/home/Database/Ensembl/homo_sapiens/Homo_sapiens.GRCh38.83.chr.refFlat
gt=/home/Database/Ensembl/homo_sapiens/genome_table

for prefix in drompa/2016_010A_D_KT55-n2-m1-hg38-raw-mpbl-sm0
do
    /share/git/DROMPAplus/compare_bed2tss -g $gene --refFlat -b $prefix.bed -m 2 --gt $gt > $prefix.comp.bed
    s="-b $prefix.bed"
    /share/git/DROMPAplus/peak_occurance -g $gene -r -b $prefix.bed --genelist Matrix.KD7.genes.count.homo_sapiens.edgeR.DEGs.csv > $prefix.prmt.bed
done

#./multibed2gene -g $gene -r $s
