log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("edgeR")

cts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene", check.names=FALSE)

coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample_name", check.names=FALSE)

group <- as.factor(coldata$condition)
dge <- DGEList(counts=cts, group=group)
dge <- calcNormFactors(dge, method='TMM')
tmm <- cpm(dge)

write.table(tmm, file=snakemake@output[["table"]], quote=FALSE, sep='\t', col.names = NA)
