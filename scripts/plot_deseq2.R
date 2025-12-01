#!/usr/bin/env Rscript

# Chargement des packages nécessaires
suppressPackageStartupMessages({
  library(ggplot2)
})
#options(bitmapType='cairo')

# Arguments en ligne de commande (1er argument = fichier CSV, 2e argument = fichier image sortie)
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) {
  stop("Usage: Rscript plot_deseq2_results.R <input_csv> <output_plot>")
}
input_csv <- args[1]
output_plot <- args[2]

# Chargement du fichier CSV
res_df <- read.csv(input_csv, stringsAsFactors = FALSE)

# Vérification que les colonnes nécessaires existent
required_cols <- c("baseMean", "log2FoldChange", "padj", "gene_id")
missing_cols <- setdiff(required_cols, colnames(res_df))
if(length(missing_cols) > 0) {
  stop(paste("Colonnes manquantes dans le CSV :", paste(missing_cols, collapse = ", ")))
}

# Ajouter colonne de signification statistique
res_df$significant <- ifelse(!is.na(res_df$padj) & res_df$padj < 0.05, "Significatif", "Non significatif")

# Filtrer les gènes significatifs pour éviter le surchargement du plot
significant_genes <- subset(res_df, significant == "Significatif" & !is.na(gene_id))

# Limiter le nombre d'étiquettes pour éviter le chevauchement
if (nrow(significant_genes) > 20) {
  significant_genes <- significant_genes[order(significant_genes$log2FoldChange, decreasing = TRUE), ]
  significant_genes <- head(significant_genes, 20)
}

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
  ) +
  geom_text(
    data = significant_genes,
    aes(label = gene_id),
    vjust = -1,  # Ajuste la position verticale des étiquettes
    hjust = 0.5, # Centre les étiquettes horizontalement
    size = 3,
    check_overlap = TRUE  # Évite le chevauchement des étiquettes
  )+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = c(1, NA)  # commence à 1, limite supérieure libre
  )
# Sauvegarde dans un fichier
pdf(file = sub("\\.png$", ".pdf", output_plot), width = 10, height = 7)
print(p)
dev.off()

cat("Plot enregistré dans :", output_plot, "\n")