sam=~/RNA-seq/Cohesin-iPS/bam/CDLS1-1.homo_sapiens.Aligned.sortedByCoord.out.bam
#sam=/home/rnakato/ChIP-seq/MCF7/Schmidt/bam/MCF7_CTCF_treat_rep1-n2-m1-hg38.sort.bam
#sam=/home/rnakato/NIPBL_30min-pair-bowtie2-hg19.sort.bam
#sam=/home/rnakato/temp/MCF7_CTCF_treat_rep1-n2-m1-hg19.sort.bam
#sam=example.sam
make DEBUG=1 && ./parse2wig+ -i $sam -o output-n2-m1 --gt /home/rnakato/Dropbox/temp/genome_table -f SAM --bed test.bed --ccp --ntype GR --mp mappability --pair #--genome GCgenome --mpbin mpbin 
