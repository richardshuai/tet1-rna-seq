log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(pheatmap)
library(edgeR)

dge <- read.table(snakemake@input[["dge"]], sep = '\t', header = TRUE, row.names = NULL)
samples <- read.table(snakemake@params[["samples"]], sep = '\t', header = TRUE)


# Filter for DGEs
dge <- dge[dge$estimatedDEG == 1, ]
rownames(dge) <- dge$gene_symbol

# Compute row-wise z-scores of log cpm normalized counts
normalized.cts <- dge[, -(1:7)]
z.scores <- t(scale(t(cpm(normalized.cts, log=TRUE))))
z.scores <- z.scores[order(-dge$q.value),]

comparison.groups <- as.data.frame(factor(samples$condition))
colnames(comparison.groups) <- "Condition"
rownames(comparison.groups) <- samples$sample_name

dge.negative.l2fc <- dge[dge$m.value < 0, ]
dge.positive.l2fc <- dge[dge$m.value > 0, ]

# Partition genes into up-regulated and down-regulated
select.negative <- Filter(function(x) !is.na(x), match(rownames(dge.negative.l2fc), rownames(z.scores)))
select.positive <- Filter(function(x) !is.na(x), match(rownames(dge.positive.l2fc), rownames(z.scores)))
select.all <- c(select.positive, select.negative)
z.scores <- z.scores[select.all, ]

# Create heatmap and save to file
col_palette <- colorRampPalette(c("#4233AA", "#FCFFD2", "#F93203"))(n=500)
# Condition <- c("#000000", "#804B09")
Condition <- c("#1B4FAB", "#0CA402")
names(Condition) <- levels(as.factor(samples$condition))
anno.colors <- list(Condition = Condition)

## Swap ordering of columns
z.scores <- z.scores[, c(4, 5, 6, 1, 2, 3)]

heatmap <- pheatmap(z.scores,
                    cluster_rows=FALSE,
                    show_rownames=TRUE,
                    show_colnames=TRUE,
                    cluster_cols=FALSE,
                    cellwidth=35,
                    cellheight=9,
                    col=col_palette,
                    border_color="white",
                    annotation=comparison.groups,
                    annotation_color=anno.colors,
                    annotation_names_col=FALSE,
                    legend=TRUE,
                    fontsize=8,
                    gaps_row=0,
                    filename=snakemake@output[['heatmap']])