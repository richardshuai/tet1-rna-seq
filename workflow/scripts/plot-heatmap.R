library(pheatmap)
library(dplyr)

LFC_THRESHOLD <- 1
FDR_CUTOFF <- 0.05



# dge <- read.table('results/tables/ko-vs-wt.table.tsv', sep = '\t', header=TRUE, row.names = NULL)

library(TCC)
cts <- read.table('results/counts/all.tsv', sep = '\t', header = TRUE, row.names = 'gene')
samples <- read.table('config/samples.tsv', sep = '\t', header = TRUE)

tcc <- new("TCC", cts, samples$condition)
tcc <- filterLowCountGenes(tcc)

tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger", iteration=3, FDR=FDR_CUTOFF, floorPDEG=0.05) # Iterative normalization and elimination of non-DEGs
tcc <- estimateDE(tcc, test.method = "edger", FDR=FDR_CUTOFF) # DEG identification
res <- getResult(tcc, sort = TRUE)


# Reorder dge and samples
dge <- dge[, c(seq(1, 7), 11, 12, 13, 8, 9, 10)]
samples <- samples[c(4, 5, 6, 1, 2, 3), ]

# Filter DGE by padj and L2FC
dge <- dplyr::filter(dge, dge$padj < FDR_CUTOFF)
dge <- dplyr::filter(dge, abs(dge$log2FoldChange) > LFC_THRESHOLD)

rownames(dge) <- dge$gene_symbol

tmm <- dge[8:(8 + nrow(samples) - 1)]
z.scores <- t(scale(t(tmm)))
z.scores <- z.scores[order(-dge$padj),]
# z.scores <- z.scores[order(-z.scores[, 1]),]

comparison.groups <- as.data.frame(factor(samples$condition))
colnames(comparison.groups) <- "Condition"
rownames(comparison.groups) <- samples$sample_name

dge.negative.l2fc <- dge[dge$log2FoldChange < 0, ]
dge.positive.l2fc <- dge[dge$log2FoldChange > 0, ]

select.negative <- Filter(function(x) !is.na(x), match(rownames(dge.negative.l2fc), rownames(z.scores)))
select.positive <- Filter(function(x) !is.na(x), match(rownames(dge.positive.l2fc), rownames(z.scores)))
select.all <- c(select.negative, select.positive)

## Create heatmap and save to file
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
col_palette <- colorRampPalette(c("#4233AA", "#FCFFD2", "#F93203"))(n=200)
Condition <- c("#0CA402", "#1B4FAB") # WU vs KU
names(Condition) <- levels(as.factor(samples$condition))
anno_colors <- list(Condition = Condition)

heatmap <- pheatmap(z.scores[select.all,], 
                    cluster_rows=FALSE,
                    show_rownames=TRUE,
                    show_colnames=TRUE,
                    cluster_cols=FALSE,
                    cellwidth=35,
                    cellheight=9,
                    col=col_palette,
                    border_color="white",
                    annotation=comparison.groups,
                    annotation_color=anno_colors,
                    annotation_names_col=FALSE,
                    legend=TRUE,
                    fontsize=8,
                    gaps_row=0)

save_pheatmap_pdf(heatmap, paste0(g1.name, "_vs_", g2.name, "/", g1.name, "_vs_", g2.name, "_at_lfc_1.5_p_.05.pdf"))

