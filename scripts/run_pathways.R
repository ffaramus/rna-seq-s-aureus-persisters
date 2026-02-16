#!/usr/bin/env Rscript

# ---------------------------------------
# KEGG ORA (Over-Representation Analysis)
# Reproduction of Fig. 2c
# ---------------------------------------

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(EnrichmentBrowser)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Usage: run_pathways.R <deseq2_results_annotated.csv> <output_dir>")
}

deseq_file <- args[1]
outdir     <- args[2]

dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

cat("Lecture du fichier DESeq2 annoté :", deseq_file, "\n")
deg <- read_csv(deseq_file, show_col_types = FALSE)

# ---------------------------------------
# 1) Filtrage DEGs UP / DOWN
# ---------------------------------------

deg <- deg %>% filter(!is.na(padj))

universe <- unique(deg$gene_id)

deg_up   <- deg %>% filter(padj < 0.05, log2FoldChange > 0)
deg_down <- deg %>% filter(padj < 0.05, log2FoldChange < 0)

genes_up   <- unique(deg_up$gene_id)
genes_down <- unique(deg_down$gene_id)

cat("DEGs UP   :", length(genes_up), "\n")
cat("DEGs DOWN :", length(genes_down), "\n")

# ---------------------------------------
# 2) Téléchargement KEGG gene sets
# ---------------------------------------

cat("Récupération des gene-sets KEGG (organisme 'sao')...\n")
gs_kegg <- get.kegg.gs("sao")

# Nettoyage des ID : enlever prefix "sao:"
clean_id <- function(x) gsub("^sao:", "", x)
gs_kegg  <- lapply(gs_kegg, clean_id)

# Récupérer noms des pathways
pw_names <- get.kegg.pathways("sao")

# ---------------------------------------
# 3) Fonction ORA (Fisher Test)
# ---------------------------------------

run_ora <- function(de_genes, universe) {

    results <- lapply(names(gs_kegg), function(pid) {
        
        gset <- intersect(gs_kegg[[pid]], universe)
        if (length(gset) == 0) return(NULL)

        overlap <- length(intersect(de_genes, gset))
        if (overlap == 0) return(NULL)

        k <- overlap
        m <- length(de_genes) - k
        n <- length(gset)     - k
        t <- length(universe) - k - m - n

        mat <- matrix(c(k, m, n, t), nrow=2, byrow=TRUE)
        ft  <- fisher.test(mat, alternative="greater")

        tibble(
            pathway_id   = pid,
            pathway_name = pw_names[pid],
            n_in_pathway = length(gset),
            n_de_in_path = k,
            odds_ratio   = unname(ft$estimate),
            pvalue       = ft$p.value
        )
    })

    df <- bind_rows(results)
    if (nrow(df) == 0) return(df)

    df <- df %>%
        mutate(padj = p.adjust(pvalue, method="BH")) %>%
        arrange(pvalue)

    return(df)
}

# ---------------------------------------
# 4) Exécution ORA pour UP / DOWN
# ---------------------------------------

cat("Analyse KEGG pour DEGs UP...\n")
res_up <- run_ora(genes_up, universe)

cat("Analyse KEGG pour DEGs DOWN...\n")
res_down <- run_ora(genes_down, universe)

# Sauvegarde CSV
write_csv(res_up,  file.path(outdir, "kegg_pathways_up.csv"))
write_csv(res_down,file.path(outdir, "kegg_pathways_down.csv"))

cat("Résultats écrits dans :", outdir, "\n")

# ---------------------------------------
# 5) Graphiques
# ---------------------------------------

make_plot <- function(df, outfile, title) {
    if (is.null(df) || nrow(df) == 0) {
        warning("Aucun pathway significatif pour :", title)
        return(NULL)
    }
    
    df_top <- df %>% slice_min(pvalue, n=min(20, n()))

    p <- ggplot(df_top, aes(x=reorder(pathway_name, -log10(pvalue)),
                            y=-log10(pvalue))) +
        geom_bar(stat="identity", fill="steelblue") +
        coord_flip() +
        xlab("KEGG Pathway") +
        ylab("-log10(p-value)") +
        ggtitle(title)

    ggsave(outfile, p, width=8, height=6)
}

make_plot(res_up,   file.path(outdir, "kegg_up_barplot.png"),   "KEGG ORA – DEGs UP")
make_plot(res_down, file.path(outdir, "kegg_down_barplot.png"), "KEGG ORA – DEGs DOWN")

cat("Plots générés.\n")
