#!/usr/bin/env Rscript

# ----------------------------------------------
# DESeq2 differential expression analysis
# Compatible with R 3.4.1 / DESeq2 1.16
# ----------------------------------------------

suppressPackageStartupMessages({
    library("DESeq2")
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop("Usage: run_deseq2.R <samples.tsv> <counts_files...> <output.csv>")
}

# First argument = samples.tsv
samples_file <- args[1]
output_file <- args[length(args)]
count_files <- args[2:(length(args)-1)]

cat("Samples file:", samples_file, "\n")
cat("Count files:\n")
print(count_files)
cat("\nOutput:", output_file, "\n\n")

# ----------------------------------------------
# 1) Read and merge featureCounts tables
# ----------------------------------------------

read_fc <- function(file) {
    df <- read.table(file, header = TRUE, sep = "\t", comment.char = "#")
    # Keep only gene ID and counts
    df <- df[, c("Geneid", grep("^\\S+\\.bam$", names(df), value = TRUE))]
    names(df)[2] <- file  # label column by filename
    return(df)
}

list_df <- lapply(count_files, read_fc)

# Merge on Geneid
counts <- Reduce(function(x, y) merge(x, y, by="Geneid"), list_df)

# Remove Geneid column (store separately)
gene_ids <- counts$Geneid
count_matrix <- counts[, -1]

# Clean column names (should match SRA IDs)
colnames(count_matrix) <- gsub(".*/|\\.bam$", "", colnames(count_matrix))

# Convert to integer matrix
count_matrix <- as.matrix(count_matrix)
mode(count_matrix) <- "integer"

# ----------------------------------------------
# 2) Read sample metadata and match conditions
# ----------------------------------------------

samples <- read.table(samples_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(samples) <- samples$sra

# Ensure order matches count_matrix columns
coldata <- samples[colnames(count_matrix), , drop = FALSE]
coldata$condition <- factor(coldata$condition)

cat("Sample table:\n")
print(coldata)
cat("\n")

# ----------------------------------------------
# 3) DESeq2 analysis
# ----------------------------------------------

dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = coldata,
    design = ~ condition
)

dds <- DESeq(dds)
res <- results(dds)

# Add gene IDs back
res_df <- as.data.frame(res)
res_df$gene_id <- gene_ids

# ----------------------------------------------
# 4) Output
# ----------------------------------------------

write.csv(res_df, file = output_file, row.names = FALSE)
cat("DESeq2 results written to:", output_file, "\n")