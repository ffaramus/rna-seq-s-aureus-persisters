
#!/usr/bin/env Rscript

# ---------------------------------------
# KEGG ORA (Over-Representation Analysis) sans EnrichmentBrowser
# ---------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(httr)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: run_pathways_noEB.R <deseq2_results_annotated.csv> <output_dir>")
}

deseq_file <- args[1]
outdir     <- args[2]

dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

cat("Lecture du fichier DESeq2 annoté :", deseq_file, "\n")

# Lecture CSV (base R)
deg <- read.csv(deseq_file, stringsAsFactors = FALSE)

# Filtrage
deg <- deg[!is.na(deg$padj), ]

universe <- unique(deg$gene_id)

deg_up   <- deg[deg$padj < 0.05 & deg$log2FoldChange > 0, ]
deg_down <- deg[deg$padj < 0.05 & deg$log2FoldChange < 0, ]

genes_up   <- unique(deg_up$gene_id)
genes_down <- unique(deg_down$gene_id)

cat("DEGs UP   :", length(genes_up), "\n")
cat("DEGs DOWN :", length(genes_down), "\n")

# --- Fonction pour récupérer les pathways KEGG et leurs gènes (organisme 'sao') via REST API KEGG

get_kegg_pathways <- function(organism) {
  url <- paste0("http://rest.kegg.jp/list/pathway/", organism)
  res <- httr::GET(url)
  stopifnot(res$status_code == 200)
  txt <- httr::content(res, "text", encoding = "UTF-8")
  lines <- strsplit(txt, "\n")[[1]]
  
  pw <- sapply(lines, function(line) {
    parts <- strsplit(line, "\t")[[1]]
    pid <- gsub("path:", "", parts[1])
    pname <- parts[2]
    c(pid, pname)
  })
  
  pathway_ids <- pw[1,]
  pathway_names <- pw[2,]
  names(pathway_names) <- pathway_ids
  
  return(pathway_names)
}

get_kegg_genesets <- function(organism, pathway_ids) {
  gs <- list()
  for (pid in pathway_ids) {
    url <- paste0("http://rest.kegg.jp/link/genes/", pid)
    res <- httr::GET(url)
    if (res$status_code != 200) next
    txt <- httr::content(res, "text", encoding = "UTF-8")
    lines <- strsplit(txt, "\n")[[1]]
    genes <- sapply(lines, function(line) {
      parts <- strsplit(line, "\t")[[1]]
      # Format: pathway_id   gene_id (ex: path:sao00010  sao:gene_id)
      gene <- gsub("^sao:", "", parts[2])
      return(gene)
    })
    genes <- unique(genes)
    gs[[pid]] <- genes
  }
  return(gs)
}

cat("Téléchargement des pathways KEGG...\n")
pw_names <- get_kegg_pathways("sao")
gs_kegg <- get_kegg_genesets("sao", names(pw_names))

# --- ORA avec test de Fisher

run_ora <- function(de_genes, universe) {
  results_list <- list()
  
  for (pid in names(gs_kegg)) {
    gset <- intersect(gs_kegg[[pid]], universe)
    if (length(gset) == 0) next
    
    overlap <- length(intersect(de_genes, gset))
    if (overlap == 0) next
    
    k <- overlap
    m <- length(de_genes) - k
    n <- length(gset) - k
    t <- length(universe) - k - m - n
    
    mat <- matrix(c(k, m, n, t), nrow=2, byrow=TRUE)
    ft <- fisher.test(mat, alternative = "greater")
    
    results_list[[pid]] <- data.frame(
      pathway_id   = pid,
      pathway_name = pw_names[pid],
      n_in_pathway = length(gset),
      n_de_in_path   = k,
      odds_ratio   = unname(ft$estimate),
      pvalue       = ft$p.value,
      stringsAsFactors = FALSE
    )
  }
  
  if (length(results_list) == 0) {
    return(data.frame())
  }
  
  df <- do.call(rbind, results_list)
  
  df$padj <- p.adjust(df$pvalue, method = "BH")
  df <- df[order(df$pvalue), ]
  return(df)
}

cat("Analyse KEGG pour DEGs UP...\n")
res_up <- run_ora(genes_up, universe)

cat("Analyse KEGG pour DEGs DOWN...\n")
res_down <- run_ora(genes_down, universe)

write.csv(res_up,  file.path(outdir, "kegg_pathways_up.csv"), row.names = FALSE)
write.csv(res_down,file.path(outdir, "kegg_pathways_down.csv"), row.names = FALSE)

cat("Résultats écrits dans :", outdir, "\n")

# --- Graphiques

make_plot <- function(df, outfile, title) {
  if (is.null(df) || nrow(df) == 0) {
    warning("Aucun pathway significatif pour :", title)
    return(NULL)
  }
  
  n_top <- min(20, nrow(df))
  df_top <- df[order(df$pvalue), ][1:n_top, ]
  
  p <- ggplot(df_top, aes(x = reorder(pathway_name, -log10(pvalue)), y = -log10(pvalue))) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    xlab("KEGG Pathway") +
    ylab("-log10(p-value)") +
    ggtitle(title)
  
  ggsave(outfile, p, width = 8, height = 6)
}

make_plot(res_up,   file.path(outdir, "kegg_up_barplot.png"),   "KEGG ORA – DEGs UP")
make_plot(res_down, file.path(outdir, "kegg_down_barplot.png"), "KEGG ORA – DEGs DOWN")

cat("Plots générés.\n")