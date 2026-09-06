library(tidyverse)
library(httr2)
library(fgsea)
library(dplyr)



##### settlement genes ####
setwd("/home/gospozha/haifa/hiba/pl_align/settlement/")

res <- read.csv("../res_annot_full.csv", check.names = FALSE, row.names = 1)

#gene_set_dir <- "./gene_sets/styl_genes"
gene_set_dir <- "./gene_sets/relaxed"


# styl_gene_sets <- list.files(gene_set_dir, pattern = "\\.txt$", full.names = TRUE) %>%
#   set_names(~ tools::file_path_sans_ext(basename(.x))) %>%
#   map(~ readLines(.x, warn = FALSE) %>%
#         str_trim() %>%
#         discard(~ is.na(.x) || .x == "") %>%
#         unique())

styl_gene_sets <- list.files(gene_set_dir, pattern = "\\genes.txt$", full.names = TRUE) %>%
  set_names(~ tools::file_path_sans_ext(basename(.x))) %>%
  map(~ readLines(.x, warn = FALSE) %>%
        str_trim() %>%
        discard(~ is.na(.x) || .x == "") %>%
        unique())

names(styl_gene_sets)
lengths(styl_gene_sets)

# run fgsea analysis
run_fgsea <- function(res, gene_sets, contrast_name = "HF_vs_NF") {
  
  res_df <- as.data.frame(res) %>%
    filter(!is.na(Symbol), is.finite(stat)) %>%
    group_by(Symbol) %>%
    slice_max(abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  ranks <- setNames(res_df$stat, res_df$Symbol) %>%
    sort(decreasing = TRUE)
  
  fgseaMultilevel(
    pathways = gene_sets,
    stats = ranks,
    minSize = 5,
    maxSize = 500,
    eps = 1e-10
  ) %>%
    as_tibble() %>%
    mutate(
      contrast = contrast_name,
      enriched_in = if_else(NES > 0, "HF", "NF")
    ) %>%
    arrange(padj)
}

styl_fgsea <- run_fgsea(
  res = res,
  gene_sets = styl_gene_sets,
  contrast_name = "HF_vs_NF"
)

signif_styl_fgsea <- styl_fgsea %>%
  filter(padj < 0.05) %>%
  mutate(leadingEdge = map_chr(leadingEdge, paste, collapse = ";")) %>%
  arrange(padj)

# no significant pathways
# write.csv(
#   signif_styl_fgsea,
#   "GSEA_styl_proteins_HF_vs_NF.csv",
#   row.names = FALSE
# )

# plot_styl_fgsea <- styl_fgsea %>%
#   filter(padj < 0.05) %>%
#   mutate(
#     pathway_plot = pathway %>%
#       str_replace_all("_", " ") %>%
#       str_wrap(width = 30)
#   ) %>%
#   arrange(enriched_in, NES) %>%
#   mutate(pathway_plot = factor(pathway_plot, levels = rev(unique(pathway_plot))))
# 
# fgsea.styl.proteins <- ggplot(
#   plot_styl_fgsea,
#   aes(
#     x = NES,
#     y = pathway_plot,
#     color = enriched_in,
#     size = -log10(padj)
#   )
# ) +
#   geom_point(alpha = 0.9) +
#   geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5) +
#   scale_color_manual(values = c("HF" = "#F393C3", "NF" = "#00A6ED")) +
#   scale_size_continuous(range = c(2.5, 7)) +
#   labs(
#     x = "Normalized enrichment score",
#     y = NULL,
#     color = "Enriched in",
#     size = expression(-log[10](adjusted~P))
#   ) +
#   theme_classic(base_size = 8) +
#   theme(
#     legend.position = "top",
#     axis.text.y = element_text(size = 7),
#     axis.ticks.y = element_blank()
#   )
# 
# ggsave(
#   "styl_proteins.gsea.jpg",
#   fgsea.styl.proteins,
#   width = 6,
#   height = 4,
#   dpi = 300
# )



# intersect with DE results

sig_de <- res %>%
  as.data.frame() %>%
  filter(
    !is.na(Symbol),
    !is.na(padj),
    !is.na(log2FoldChange),
    padj < 0.05,
    abs(log2FoldChange) > 1
  ) %>%
  group_by(Symbol) %>%
  slice_min(pvalue, n = 1, with_ties = FALSE) %>%
  ungroup()

gene_set_de_overlap <- imap_dfr(
  styl_gene_sets,
  ~ sig_de %>%
    filter(Symbol %in% .x) %>%
    mutate(
      gene_set = .y,
      direction = if_else(log2FoldChange > 0, "HF", "NF")
    )
) %>%
  select(
    gene_set,
    Symbol,
    log2FoldChange,
    pvalue,
    padj,
    direction,
    everything()
  ) %>%
  arrange(gene_set, pvalue)

gene_set_de_overlap

write.csv(
  gene_set_de_overlap,
  "Stylophora_gene_sets_DE_overlap.csv",
  row.names = FALSE
)

gene_set_overlap_summary <- imap_dfr(
  styl_gene_sets,
  ~ tibble(
    gene_set = .y,
    gene_set_size = length(unique(.x)),
    genes_tested = sum(unique(.x) %in% res$Symbol),
    significant_de_genes = sum(unique(.x) %in% sig_de$Symbol),
    significant_up_HF = sum(unique(.x) %in% sig_de$Symbol[sig_de$log2FoldChange > 0]),
    significant_up_NF = sum(unique(.x) %in% sig_de$Symbol[sig_de$log2FoldChange < 0])
  )
) %>%
  arrange(desc(significant_de_genes))

gene_set_overlap_summary

write.csv(
  gene_set_overlap_summary,
  "Stylophora_gene_sets_DE_overlap_summary.csv",
  row.names = FALSE
)