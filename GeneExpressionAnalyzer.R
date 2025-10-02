# RNA-seq analysis for cancer types using normalized data
# load stuff
library(pheatmap)
library(RColorBrewer)
library(genefilter)

# read in the data
data <- read.csv("data.csv", row.names = 1)
head(data)

# read sample info
labels <- read.csv("labels.csv", row.names = 1)
colnames(labels) <- "condition"

# make sure samples match up
sample_names <- colnames(data)
labels <- labels[match(sample_names, rownames(labels)), , drop = FALSE]
if (!all(sample_names == rownames(labels))) {
    stop("Sample names don't match between files!")
}

# set up sample groups (5 cancer types)
groups <- factor(labels$condition, levels = c("LUAD", "PRAD", "BRCA", "KIRC", "COAD"))
sample_info <- data.frame(row.names = sample_names, condition = groups)

# demonstrate loop through cancer types
cat("\nSample counts by cancer type:\n")
for(cancer_type in levels(groups)) {
  count <- sum(groups == cancer_type, na.rm=TRUE)
  cat(cancer_type, ":", count, "samples\n")
}

# clean up any weird values
data[is.na(data)] <- 0
data[data == Inf | data == -Inf] <- 0

# heatmap 1 - top 50 variable genes
top_genes <- order(rowVars(as.matrix(data)), decreasing = TRUE)[1:50]
pheatmap(data[top_genes, ],
    cluster_rows = TRUE,
    show_rownames = TRUE,
    cluster_cols = TRUE,
    annotation_col = sample_info,
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
    main = "Top 50 Variable Genes",
    filename = "heatmap_top_genes.png",
    width = 10, height = 8
)

# heatmap 2 - high variance genes (top 5%)
var_threshold <- quantile(rowVars(as.matrix(data)), 0.95, na.rm = TRUE)
high_var_genes <- which(rowVars(as.matrix(data)) >= var_threshold)
if (length(high_var_genes) > 0) {
    pheatmap(data[high_var_genes, ],
        cluster_rows = TRUE,
        show_rownames = ifelse(length(high_var_genes) < 50, TRUE, FALSE),
        cluster_cols = TRUE,
        annotation_col = sample_info,
        color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
        main = "High Variance Genes",
        filename = "heatmap_high_var_genes.png",
        width = 10, height = 8
    )
} else {
    print("No high variance genes found")
}

# saving the gene list
high_var_gene_names <- rownames(data)[high_var_genes]
write.csv(data.frame(Gene = high_var_gene_names), "high_variance_genes.csv", row.names = FALSE)
