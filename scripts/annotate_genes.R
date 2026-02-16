#!/usr/bin/env Rscript
if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr", repos = "https://cloud.r-project.org/")
}

suppressPackageStartupMessages({
    library(dplyr)
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
deg <- read.csv(deseq_file)

# Vérifications minimales
if (!"gene_id" %in% colnames(deg)) {
    stop("Le fichier DESeq2 doit contenir une colonne 'gene_id'.")
}

# ----------------------------------------------
# Lecture et parsing du GFF
# ----------------------------------------------

cat("Lecture du GFF :", gff_file, "\n")

gff <- read.delim(gff_file, comment.char="#", header=FALSE, sep="\t", stringsAsFactors=FALSE)

# Le GFF3 contient 9 colonnes
if (ncol(gff) < 9) stop("Format GFF inattendu.")

colnames(gff)[1:9] <- c("seqid","source","type","start","end","score","strand","phase","attributes")

# On garde uniquement les gènes ou CDS
gff_genes <- gff %>% filter(type %in% c("gene","CDS"))

# Extraction attributs sous forme clé=valeur
extract_attr <- function(attr, key) {
    match <- str_match(attr, paste0(key, "=([^;]+)"))
    return(ifelse(is.na(match[,2]), NA, match[,2]))
}

gff_annot <- gff_genes %>% 
    mutate(
        locus_tag = extract_attr(attributes, "locus_tag"),
        gene_name = extract_attr(attributes, "gene"),
        product   = extract_attr(attributes, "product")
    ) %>%
    select(locus_tag, gene_name, product) %>%
    distinct()

cat("Annotations extraites :", nrow(gff_annot), "gènes.\n")

# ----------------------------------------------
# Merge DESeq2 + annotation
# ----------------------------------------------

annotated <- deg %>%
    left_join(gff_annot, by = c("gene_id" = "locus_tag"))

cat("Résultats annotés :", nrow(annotated), "entrées.\n")

# ----------------------------------------------
# Sauvegarde
# ----------------------------------------------

write.csv(annotated, outfile, row.names = FALSE)
cat("Fichier annoté écrit dans :", outfile, "\n")
