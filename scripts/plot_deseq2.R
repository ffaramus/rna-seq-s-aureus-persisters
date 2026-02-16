##!/usr/bin/env Rscript


# Chargement des packages nécessaires
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

# Arguments en ligne de commande (1er argument = fichier CSV, 2e argument = fichier image sortie)
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2){
  stop("Usage: Rscript plot_deseq2_results.R <input_csv> <output_plot>")
}

input_csv <- args[1]
output_plot <- args[2]

# Chargement du fichier CSV
res_df <- read.csv(input_csv, stringsAsFactors = FALSE)

# Vérification que les colonnes nécessaires existent
required_cols <- c("baseMean", "log2FoldChange", "padj", "gene_name")  # adapter 'gene_name' selon colonne réelle
missing_cols <- setdiff(required_cols, colnames(res_df))
if(length(missing_cols) > 0){
  stop(paste("Colonnes manquantes dans le CSV :", paste(missing_cols, collapse = ", ")))
}

# Ajouter colonne de signification statistique
res_df$significant <- ifelse(!is.na(res_df$padj) & res_df$padj < 0.05, "Significatif", "Non significatif")

# Construction du plot
p <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Significatif" = "red", "Non significatif" = "grey")) +
  theme_minimal() +
  labs(
    x = "Abondance moyenne du gène (niveau basal contrôle)",
    y = "Différence d'expression (Persisters vs Contrôles, log2FC)",
    color = "Significativité (p-adj < 0.05)",
    title = "Expression différentielle des gènes"
  ) + geom_text_repel(data = subset(res_df, significant == "Significatif"),
                  aes(label = gene_name),
                  size = 3,
                  max.overlaps = 10)

# Sauvegarde dans un fichier 
ggsave(filename = output_plot, plot = p, width = 10, height = 7, units = "in", dpi = 300)

cat("Plot enregistré dans :", output_plot, "\n")
