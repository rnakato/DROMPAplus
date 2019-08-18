library(IdeoViz)
list = list("mm9", "mm10", "hg19", "hg38",)
for (build in list) {
	file <- paste(build, ".tsv", sep = "")
	print (file)
	write.table(getIdeo(build), file, quote=F, sep="\t", col.names=F, row.names = F)
}