#!/usr/bin/env Rscript
library(ggplot2)
library(VennDiagram)
library(ggrepel)
library(grid)
library(DESeq2)

# ==== Arguments Nextflow ====
args <- commandArgs(trailingOnly = TRUE)
paper_path         <- args[1]
deseq_results_path <- args[2]
counts_file        <- args[3]
mapping_path       <- args[4]


# Lire la table du papier (attention appel robuste sinon ne marche pas)
res_paper <- read.delim(
  file = paper_path,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fill = TRUE,
  na.strings = c("", "NA", "N/A"),
  strip.white = TRUE
)


# Lire deseq results
res <- read.csv(deseq_results_path,row.names =1)

# Lire counts mat
cts <- read.table(
  counts_file,
  header       = TRUE,
  comment.char = "#",
  row.names    = 1,
  check.names  = FALSE
)

# Lire table de mapping 
mapping <- read.table(
  mapping_path,
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE,
  stringsAsFactors = FALSE,
  check.names = F
)
# On ajoute 'pth' à la main car il n'est pas présent dans la table d'origine
mapping[mapping$product == "peptidyl-tRNA hydrolase", ]$symbol="pth"


# === Pré-traitement ===
# Retirer les colonnes d'annot
annot_cols <- c("Chr", "Start", "End", "Strand", "Length")
cts <- cts[, !(colnames(cts) %in% annot_cols), drop = FALSE]

# Clean noms de colonne : à enlever quand on aura fix directement sur le main.nf
old_names <- colnames(cts)
new_names <- gsub("_trimmed.sorted.bam$", "", old_names)
colnames(cts) <- new_names
# =====================


# On merge les dataframes:
res_combined <- cbind(
  res,
  cts[rownames(res), ], 
  mapping[match(rownames(res), mapping$locus_tag), -1]
)

# Je rajoute "symbol" à la table du papier
res_paper <- merge(
  res_paper,
  mapping[, c("locus_tag", "symbol")],
  by.x = "Name",
  by.y = "locus_tag",
  all.x = TRUE
)


# ==== Compare Results ====
# Créer une colonne commune
res_combined$Name <- rownames(res_combined)
# Gènes communs 
common_genes <- intersect(res_paper$Name, res_combined$Name)
# Gènes présents dans le papier mais absents de nos résultats 
missing_in_repro <- setdiff(res_paper$Name, res_combined$Name)
# Gènes présents dans nos résultats mais absents du papier 
extra_in_repro <- setdiff(res_combined$Name, res_paper$Name)

# Tableau
compare_summary <- data.frame(
  Category = c("Total genes (paper)", 
               "Total genes (reproduction)", 
               "Common genes",
               "Missing in reproduction",
               "Extra in reproduction"),
  Count = c(
    length(res_paper$Name),
    nrow(res_combined),
    length(common_genes),
    length(missing_in_repro),
    length(extra_in_repro)
  )
)

write.table(compare_summary,
            file = "comparison_summary_genes.tsv",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
# =========================


# === Comparaison des valeurs numériques : ===
df_compare <- merge(
  res_paper[, c("Name", "baseMean", "log2FoldChange", "padj")],
  res_combined[, c("Name", "baseMean", "log2FoldChange", "padj")],
  by = "Name",
  suffixes = c("_paper", "_repro")
)

cor_baseMean <- cor(df_compare$baseMean_paper, df_compare$baseMean_repro, use = "complete.obs")
cor_log2FC   <- cor(df_compare$log2FoldChange_paper, df_compare$log2FoldChange_repro, use = "complete.obs")
# ===========================================


# === PCA : papier ===
# Mapping
sample_mapping <- c(
  "SRR10379721" = "IP1",
  "SRR10379722" = "IP2",
  "SRR10379723" = "IP3",
  "SRR10379724" = "ctrl4",
  "SRR10379725" = "ctrl5",
  "SRR10379726" = "ctrl6"
)

# Colonnes de comptage (ctrl1, IP2, ...)
rownames<-res_paper$Name
count_cols <- grep("^ctrl|^IP", colnames(res_paper), value = TRUE)
paper_counts <- res_paper[, count_cols]
rownames(paper_counts)=rownames

# Reverse mapping pour obtenir les noms SRR
reverse_mapping <- setNames(names(sample_mapping), sample_mapping)

# Renommer les colonnes
new_colnames <- reverse_mapping[colnames(paper_counts)]
colnames(paper_counts) <- new_colnames
paper_counts <- paper_counts[!apply(is.na(paper_counts), 1, any), ]

# Créer le coldata avec conditions et replicates
condition <- ifelse(grepl("SRR1037972[123]", new_colnames), "persister", "control")
coldata <- data.frame(
  row.names = new_colnames,
  condition = factor(condition)
)

# Création DESeqDataSet
dds_paper <- DESeqDataSetFromMatrix(
  countData = paper_counts, 
  colData = coldata,
  design = ~ condition
)

# Transformation VST
vst_paper <- vst(dds_paper, blind = FALSE)

pdf("PCA_paper.pdf")
plotPCA(vst_paper, intgroup="condition") +
  labs(color = "Sample",
       title = "PCA - Paper") +
  geom_label(
    aes(label = name),
    size = 3,
    label.padding = unit(0.15, "lines"),
    alpha = 0.7,
    nudge_y = 1.1
  ) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_cartesian(clip = "off")             # labels en dehors visibles
dev.off()
# ====================


# === Venn Diagram gènes DE entre le papier et la repro ===
# Gènes DE
de_paper <- unique(na.omit(df_compare$Name[df_compare$padj_paper < 0.05]))
de_repro <- unique(na.omit(df_compare$Name[df_compare$padj_repro < 0.05]))

pdf("venn_DEgenes_paper_repro.pdf")
grid.newpage()
pushViewport(viewport(width = unit(10, "cm"), height = unit(10, "cm")))
draw.pairwise.venn(
  area1 = length(de_paper),
  area2 = length(de_repro),
  cross.area = length(intersect(de_paper, de_repro)),
  
  category = c("DE genes in the paper", "DE genes in the repro"),
  
  lwd = 1,
  col = c("#440154ff", "#21908dff"),
  fill = c("#440154ff", "#21908dff"),
  alpha = c(0.3, 0.3),
  
  cex = 0.8,
  fontfamily = "sans",
  
  cat.pos = c(-30, 30),      # angle
  cat.dist = c(0.05, 0.05),  # distance au cercle
  cat.cex = 1,
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", "#21908dff")
)
dev.off()
# =========================================================


# === Quadrant Plot sur les DE genes ===
df_de <- df_compare[
  df_compare$padj_paper < 0.05 &
    df_compare$padj_repro < 0.05,
]

# Définir les catégories quadrant (UP/DOWN)
df_de$category <- NA

df_de$category[
  df_de$log2FoldChange_paper > 0 &
    df_de$log2FoldChange_repro > 0
] <- "UP-UP"

df_de$category[
  df_de$log2FoldChange_paper < 0 &
    df_de$log2FoldChange_repro < 0
] <- "DOWN-DOWN"

df_de$category[
  df_de$log2FoldChange_paper > 0 &
    df_de$log2FoldChange_repro < 0
] <- "UP-DOWN"

df_de$category[
  df_de$log2FoldChange_paper < 0 &
    df_de$log2FoldChange_repro > 0
] <- "DOWN-UP"


# Plot
pdf("scatter_quadrants_DE_genes.pdf")
ggplot(df_de, aes(x = log2FoldChange_paper, y = log2FoldChange_repro)) +
  geom_point(aes(color = category), alpha = 0.6, size = 1.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  theme_classic() +
  scale_color_manual(values = c(
    "UP-UP" = "#E41A1C",     # rouge
    "DOWN-DOWN" = "#377EB8", # bleu
    "UP-DOWN" = "#984EA3",   # violet
    "DOWN-UP" = "#FF7F00"    # orange
  )) +
  labs(
    title = "Quadrant scatter (DE-only)",
    subtitle = paste0("n = ", nrow(df_de), " DE genes"),
    x = "log2FC (paper)",
    y = "log2FC (repro)",
    color = "Category"
  ) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
# ==========================



# ==== VOLCANO PLOT DU PAPIER (bonus) ====
# Préparer les colonnes utiles
res_paper$minusLog10Padj <- -log10(res_paper$padj)

res_paper$diffexp <- "Not significant"
res_paper$diffexp[res_paper$log2FoldChange > 0 & res_paper$padj<0.05] <- "Up"
res_paper$diffexp[res_paper$log2FoldChange < 0 & res_paper$padj<0.05] <- "Down"

# Top 10 gènes les plus signifs
top10 <- res_paper[order(res_paper$padj), ][1:10, ]

pdf("volcano_paper.pdf")
# Plot
volcano <- ggplot(res_paper, aes(x = log2FoldChange, y = minusLog10Padj)) +
  geom_point(aes(color = diffexp), alpha = 0.7, size = 1.5) +
  
  # Seuils
  geom_vline(xintercept = c(-0.6, 0.6), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  
  # Labels des top gènes
  geom_text_repel(
    data = top10,
    aes(label = symbol),
    size = 2.2,
    segment.color = "black",
    segment.size = 0.8,
    min.segment.length = 0,
    box.padding = 0.8,
    point.padding = 0
  ) +
  
  # Couleurs
  scale_color_manual(
    values = c("Down" = "blue", "Not significant" = "grey70", "Up" = "red")
  ) +
  
  theme_classic() +
  labs(
    title = "Volcano plot : Paper results",   
    x = expression(log[2]~"Fold Change"),
    y = expression(-log[10]~"adjusted p-value"),
    color = "Differential expression"
  )+
  theme(plot.title = element_text(hjust = 0.5))


print(volcano)
dev.off()
#======================================


# === Scatter plots ===
pdf("scatterplot_log2FC_paper_repro.pdf")
ggplot(df_compare, aes(x = log2FoldChange_paper, y = log2FoldChange_repro)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  theme_classic() +
  labs(title = paste("log2FC correlation (r =", round(cor_log2FC, 3), ")"))+ 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("scatterplot_basemean_paper_repro.pdf")
ggplot(df_compare, aes(x = baseMean_paper, y = baseMean_repro)) +
  geom_point(alpha = 0.4) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  theme_classic() +
  labs(title = paste("BaseMean correlation (r =", round(cor_baseMean, 3), ")")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
# ====================
