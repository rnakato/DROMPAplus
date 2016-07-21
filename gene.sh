gene=/home/Database/Ensembl/homo_sapiens/*refFlat

for prefix in
do
    ./compare_bed2tss -g $gene --refFlat -b $prefix.bed -m 2 > $prefix.comp.bed
    s="-b $prefix.bed"
    ./peak_occurance -g $gene -r -b $prefix.bed --genelist genelist.txt > $prefix.prmt.bed
done

./multibed2gene -g $gene -r $s
