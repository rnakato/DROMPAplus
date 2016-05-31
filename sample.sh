sam=/home/rnakato/temp/MCF7_CTCF_treat_rep1-n2-m1-hg19.sort.bam
#sam=example.sam
make DEBUG=1 && ./parse2wig+ -i $sam -o output-n2-m1 --gt /home/rnakato/Dropbox/temp/genome_table -f SAM --bed test.bed --ccp --ntype GR --mp mappability #--genome GCgenome --mpbin mpbin 
