```r
# RNA-seq analysis for cancer types using normalized data
# load stuff
library(pheatmap)
library(RColorBrewer)
library(genefilter)

# read in the data
# Updated to skip empty first column header (learned from zero-length name error)
data <- read.csv("data.csv", row.names = 1, check.names = FALSE)
# Check for valid gene and sample names
if (any(is.na(rownames(data)) | rownames(data) == "")) {
    stop("Empty or NA gene names in data.csv!")
}
if (any(is.na(colnames(data)) | colnames(data) == "")) {
    stop("Empty or NA sample names in data.csv!")
}

# read sample info
# Updated to skip empty first column header
labels <- read.csv("labels.csv", row.names = 1, check.names = FALSE)
colnames(labels) <- "condition"
# Check for valid sample names
if (any(is.na(rownames(labels)) | rownames(labels) == "")) {
    stop("Empty or NA sample names in labels.csv!")
}

# make sure samples match up
sample_names <- colnames(data)
labels <- labels[match(sample_names, rownames(labels)), , drop = FALSE]
if (!all(sample_names == rownames(labels))) {
    stop("Sample names don't match between files!")
}

# set up sample groups (removed COAD since 0 samples)
groups <- factor(labels$condition, levels = c("LUAD", "PRAD", "BRCA", "KIRC"))
sample_info <- data.frame(row.names = sample_names, condition = groups)

# demonstrate loop through cancer types
cat("\nSample counts by cancer type:\n")
for (cancer_type in levels(groups)) {
    count <- sum(groups == cancer_type, na.rm = TRUE)
    cat(cancer_type, ":", count, "samples\n")
}

# clean up any weird values
data[is.na(data)] <- 0
data[data == Inf | data == -Inf] <- 0

# Added variance filter to avoid hclust errors from zero-variance rows (learned from crashes!)
vars <- rowVars(as.matrix(data), na.rm = TRUE)
vars <- vars[!is.na(vars) & vars > 0] # Remove NA or zero-variance rows
data_filt <- data[names(vars), ] # Keep only valid genes for both heatmaps

# Heatmap 1 - top 50 variable genes
if (length(vars) >= 2) { # Need at least 2 genes for clustering
    top_genes <- order(vars, decreasing = TRUE)[1:min(50, length(vars))]
    if (length(top_genes) >= 2) {
        pheatmap(data_filt[top_genes, ],
            cluster_rows = TRUE,
            show_rownames = TRUE,
            cluster_cols = length(top_genes) >= 2, # Disable col clustering if not enough
            annotation_col = sample_info,
            color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
            main = "Top 50 Variable Genes",
            filename = "heatmap_top_genes.png",
            width = 10, height = 8
        )
    } else {
        print("Not enough variable genes for heatmap 1")
    }
} else {
    print("Not enough valid genes for heatmap 1")
}

# Heatmap 2 - high variance genes (top 5%)
# Updated to use filtered data (learned to reuse vars to avoid recomputing)
var_threshold <- quantile(vars, 0.95, na.rm = TRUE)
high_var_genes <- which(vars >= var_threshold)
if (length(high_var_genes) >= 2) {
    pheatmap(data_filt[high_var_genes, ],
        cluster_rows = TRUE,
        show_rownames = ifelse(length(high_var_genes) < 50, TRUE, FALSE),
        cluster_cols = length(high_var_genes) >= 2,
        annotation_col = sample_info,
        color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
        main = "High Variance Genes",
        filename = "heatmap_high_var_genes.png",
        width = 10, height = 8
    )
} else {
    print("Not enough high variance genes for heatmap 2")
}

# saving the gene list
# Updated to use filtered data (learned to match heatmap genes)
if (length(high_var_genes) > 0) {
    high_var_gene_names <- rownames(data_filt)[high_var_genes]
    write.csv(data.frame(Gene = high_var_gene_names), "high_variance_genes.csv", row.names = FALSE)
} else {
    print("No high variance genes to save")
}
```