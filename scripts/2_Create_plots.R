#!/usr/bin/env Rscript
library(ggplot2)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
deseq_results_path  <- args[1]
genes_pathways_path <- args[2]
mapping_path        <- args[3]

res <- read.csv(deseq_results_path, row.names = 1)

kegg_all <- read.table(
  genes_pathways_path, header = TRUE, sep = "\t",
  stringsAsFactors = FALSE
)

mapping <- read.table(
  mapping_path, header = TRUE, sep = "\t", quote = "",
  fill = TRUE, stringsAsFactors = FALSE, check.names = FALSE
)

mapping[mapping$product == "peptidyl-tRNA hydrolase", ]$symbol <- "pth"

res$log2BaseMean <- ifelse(res$baseMean > 0, log2(res$baseMean), NA)
res$is_sig <- !is.na(res$padj) & res$padj < 0.05
res$signif <- ifelse(res$is_sig, "Significant", "Non-Significant")

res$locus_tag <- rownames(res)

res_annot <- merge(res, mapping, by = "locus_tag", all.x = TRUE)

colnames(kegg_all) <- c("locus_tag", "Pathway_ID")
res_annot <- merge(res_annot, kegg_all, by = "locus_tag", all.x = TRUE)

translation_genes <- !is.na(res_annot$Pathway_ID) &
  (grepl("(^|;)sao03011(;|$)", res_annot$Pathway_ID) |
   grepl("(^|;)sao03009(;|$)", res_annot$Pathway_ID) |
   grepl("(^|;)sao03016(;|$)", res_annot$Pathway_ID) |
   grepl("(^|;)sao03012(;|$)", res_annot$Pathway_ID))

AA_tRNA_synthetases <- !is.na(res_annot$Pathway_ID) &
  grepl("(^|;)sao03016(;|$)", res_annot$Pathway_ID) &
  grepl("-tRNA synthetase", res_annot$product)

res_annot$is_AA_tRNA <- AA_tRNA_synthetases

typical_members <- subset(
  res_annot,
  symbol %in% c("pth", "infA", "infB", "infC", "frr", "tsf")
)

res_annot$minusLog10Padj <- -log10(res_annot$padj)
res_annot$diffexp <- "Not significant"
res_annot$diffexp[res_annot$log2FoldChange > 0 & res_annot$is_sig] <- "Up"
res_annot$diffexp[res_annot$log2FoldChange < 0 & res_annot$is_sig] <- "Down"

top10 <- res_annot[order(res_annot$padj), ][1:10, ]

pdf("volcano_repro.pdf")
volcano <- ggplot(res_annot, aes(x = log2FoldChange, y = minusLog10Padj)) +
  geom_point(aes(color = diffexp), alpha = 0.7, size = 1.5) +
  geom_vline(xintercept = c(-0.6, 0.6), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_text_repel(
    data = top10,
    aes(label = symbol),
    size = 2.2
  ) +
  scale_color_manual(values = c(
    "Down" = "blue", "Not significant" = "grey70", "Up" = "red")) +
  theme_classic()

print(volcano)
dev.off()

plot_df <- res_annot[translation_genes, ]

pdf("MA_plot_translationgenes.pdf")
p <- ggplot(plot_df, aes(x = log2BaseMean, y = log2FoldChange, color = signif)) +
  geom_point(alpha = 0.8, size = 1.2) +
  geom_point(
    data = subset(plot_df, is_AA_tRNA),
    aes(shape = "AA_tRNA_synthetases"),
    color = "black", size = 1.3, stroke = 1.3
  ) +
  scale_color_manual(values = c("Non-Significant" = "grey70",
                                "Significant" = "red")) +
  geom_text_repel(
    data = typical_members,
    aes(label = symbol),
    size = 4
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

print(p)
dev.off()

