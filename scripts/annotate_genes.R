#!/usr/bin/env Rscript
if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr", repos = "https://cloud.r-project.org/")
}

suppressPackageStartupMessages({
    library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop("Usage: annotate_genes.R <deseq2_results.csv> <reference.gff3> <output.csv>")
}

deseq_file <- args[1]
gff_file   <- args[2]
outfile    <- args[3]

cat("Lecture du fichier DESeq2 :", deseq_file, "\n")
deg <- read.csv(deseq_file, stringsAsFactors = FALSE)

if (!"gene_id" %in% colnames(deg)) {
    stop("Le fichier DESeq2 doit contenir une colonne 'gene_id'.")
}

# ----------------------------------------------
# Lecture et filtrage du GFF
# ----------------------------------------------

cat("Lecture du GFF :", gff_file, "\n")

gff <- read.delim(gff_file, comment.char="#", header=FALSE, sep="\t",
                  stringsAsFactors=FALSE)

if (ncol(gff) < 9) stop("Format GFF inattendu.")

colnames(gff)[1:9] <- c(
    "seqid","source","type","start","end",
    "score","strand","phase","attributes"
)

# Equivalent de filter(type %in% ...)
gff_genes <- gff[gff$type %in% c("gene", "CDS"), ]

# ----------------------------------------------
# Extraction des attributs
# ----------------------------------------------

extract_attr <- function(attr, key) {
    match <- str_match(attr, paste0(key, "=([^;]+)"))
    ifelse(is.na(match[,2]), NA, match[,2])
}

gff_annot <- data.frame(
    locus_tag = extract_attr(gff_genes$attributes, "locus_tag"),
    gene_name = extract_attr(gff_genes$attributes, "gene"),
    product   = extract_attr(gff_genes$attributes, "product"),
    stringsAsFactors = FALSE
)

# Equivalent de distinct() : suppression des doublons
gff_annot <- gff_annot[!duplicated(gff_annot$locus_tag), ]

cat("Annotations extraites :", nrow(gff_annot), "gènes.\n")

# ----------------------------------------------
# Merge DESeq2 + GFF
# ----------------------------------------------

annotated <- merge(
    deg, gff_annot,
    by.x = "gene_id",
    by.y = "locus_tag",
    all.x = TRUE
)

cat("Résultats annotés :", nrow(annotated), "entrées.\n")

# ----------------------------------------------
# Sauvegarde
# ----------------------------------------------

write.csv(annotated, outfile, row.names = FALSE)
cat("Fichier annoté écrit dans :", outfile, "\n")