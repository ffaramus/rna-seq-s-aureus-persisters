#!/usr/bin/env Rscript
library(KEGGREST)

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]

gene_pathway <- keggLink(target = "brite", source = "sao")

Gene_ID    <- sub("^sao:", "", names(gene_pathway))
Pathway_ID <- sub("^br:" , "", gene_pathway)

gene_pathway <- data.frame(Gene_ID, Pathway_ID)

gene_pathway <- aggregate(. ~ Gene_ID, data = gene_pathway,
                          FUN = function(x) paste(unique(x), collapse = ";"))

write.table(
  gene_pathway,
  file = file.path(out_dir, "gene_pathway_table.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

