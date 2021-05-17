log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(TCC)
library(dplyr)

LFC_THRESHOLD <- 1.5
FDR_CUTOFF <- 0.05

cts <- read.table(snakemake@input[["counts"]], sep = '\t', header = TRUE, row.names = 'gene')
samples <- read.table(snakemake@params[["samples"]], sep = '\t', header = TRUE)

# TCC pipeline
tcc <- new("TCC", cts, samples$condition)
tcc <- filterLowCountGenes(tcc)
tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger", iteration=3, FDR=FDR_CUTOFF, floorPDEG=0.05) # Iterative normalization and elimination of non-DEGs

# Get normalized counts
normalized.cts <- as.data.frame(getNormalizedData(tcc))
normalized.cts <- tibble::rownames_to_column(normalized.cts, "gene_id")
normalized.cts$gene_id <- as.factor(normalized.cts$gene_id)

# Get DGEs and filter by L2FC
tcc <- estimateDE(tcc, test.method = "edger", FDR=FDR_CUTOFF) # DEG identification
res <- getResult(tcc, sort = TRUE)
res$estimatedDEG <- as.numeric(res$estimatedDEG & (abs(res$m.value) > LFC_THRESHOLD))

# Merge DEG table to include normalized counts
dge.table <- inner_join(res, normalized.cts, by="gene_id")

# Store results
write.table(dge.table, file=snakemake@output[["dge"]], row.names=FALSE, sep='\t')