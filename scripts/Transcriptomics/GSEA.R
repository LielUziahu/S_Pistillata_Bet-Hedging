library(tidyverse)
library(httr2)
library(fgsea)
library(dplyr)

# setting working directory 
setwd("/home/gospozha/haifa/hiba/pl_align/GSEA")


#### REACTOME ####
anno <- read.csv("~/haifa/S.pistillata/rna/genomes/anno/gene_annotations_reactome.tsv", sep = "\t", check.names = FALSE)
res <- read.csv("../res_annot_full.csv", check.names = FALSE, row.names = 1)

names(anno)
names(res)

# term2gene
reactome_gene_ids <- anno %>%
  dplyr::select(gene, ipr_pathways) %>%
  filter(!is.na(ipr_pathways), ipr_pathways != "") %>%
  separate_rows(ipr_pathways, sep = "\\|") %>%
  mutate(ipr_pathways = str_trim(ipr_pathways)) %>%
  filter(str_starts(ipr_pathways, "Reactome:")) %>%
  transmute(gene = as.character(gene),
            reactome_id = str_remove(ipr_pathways, "^Reactome:")) %>%
  distinct()


# proper names
reactome_names <- read.delim(
  "https://reactome.org/download/current/ReactomePathways.txt",
  header = FALSE,
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE
)

colnames(reactome_names) <- c("reactome_id", "pathway_name", "species")

# add pathways only from invertebrates
selected_species <- c(
  "Drosophila melanogaster",
  "Caenorhabditis elegans",
  "Arenicola marina",
  "Homarus americanus"
)

reactome_annot <- reactome_gene_ids %>%
  left_join(reactome_names, by = "reactome_id") %>%
  filter(species %in% selected_species)

term2gene <- reactome_annot %>%
  filter(!is.na(pathway_name), pathway_name != "") %>%
  dplyr::select(pathway = pathway_name, gene) %>%
  distinct()

# split into sets
reactome_sets <- split(term2gene$gene, term2gene$pathway)

#### GSEA on Reactome ####
term2gene <- read_tsv("term2gene.reactome.tsv")
# gene ranking
rank_table <- res %>%
  transmute(gene = as.character(Symbol), stat = as.numeric(stat))

# remove NA and duplicates
rank_table <- rank_table %>%
  filter(!is.na(gene), is.finite(stat)) %>%
  group_by(gene) %>%
  slice_max(abs(stat), n = 1, with_ties = FALSE) %>%
  ungroup()

# named ranking
ranks <- setNames(rank_table$stat, rank_table$gene)
ranks <- sort(ranks, decreasing = TRUE)

# filter gene sets
reactome_sets <- reactome_sets[
  lengths(reactome_sets) >= 10 &
    lengths(reactome_sets) <= 500
]

# run multilevel fgsea 
fg.reactome <- fgseaMultilevel(
  pathways = reactome_sets,
  stats = ranks,
  minSize = 10,
  maxSize = 500,
  eps = 0
) 

fg.reactome.arranged <- fg.reactome %>%
  as_tibble() %>%
  mutate(
    enriched_in = if_else(NES > 0, "HF", "NF"),
    leadingEdge = map_chr(leadingEdge, paste, collapse = ";")
  ) %>%
  arrange(padj)

write.csv(
  fg.reactome.arranged,
  "GSEA_Reactome_HF_vs_NF.csv",
  row.names = FALSE
)

fg.reactome.sign <- fg.reactome %>%
  dplyr::filter(padj < 0.05) 

fg.reactome.sign.export <-fg.reactome.sign %>%
  as_tibble() %>%
  mutate(
    enriched_in = if_else(NES > 0, "HF", "NF"),
    leadingEdge = map_chr(leadingEdge, paste, collapse = ";")
  ) %>%
  arrange(enriched_in, padj)

write.csv(
  fg.reactome.sign.export,
  "GSEA_Reactome_HF_vs_NF_significant.csv",
  row.names = FALSE
)

# # collapse pathways
# fg.hf <- fg.reactome.sign %>% filter(NES > 0)
# fg.nf <- fg.reactome.sign %>% filter(NES < 0)
# 
# 
# 
# collapsed.hf <- collapsePathways(
#   fgseaRes = fg.hf,
#   pathways = reactome_sets,
#   stats = ranks,
#   pval.threshold = 0.05
# )
# 
# collapsed.nf <- collapsePathways(
#   fgseaRes = fg.nf,
#   pathways = reactome_sets,
#   stats = ranks,
#   pval.threshold = 0.05
# )
# 
# main.pathways <- c(collapsed.hf$mainPathways, collapsed.nf$mainPathways)
# 
# fg.reactome.main <- fg.reactome.sign %>%
#   as_tibble() %>%
#   filter(pathway %in% main.pathways) %>%
#   mutate(enriched_in = if_else(NES > 0, "HF", "NF")) %>%
#   arrange(enriched_in, padj)
# 
# fg.reactome.main.export <- fg.reactome.main %>%
#   mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1)))
# 
# write.csv(
#   as.data.frame(fg.reactome.main.export),
#   "GSEA_Reactome_HF_vs_NF_main_pathways.csv",
#   row.names = FALSE
# )
# 
# # plot 
# fg.plot <- fg.reactome.main %>%
#   mutate(enriched_in = if_else(NES > 0, "HF", "NF")) %>%
#   group_by(enriched_in) %>%
#   slice_min(padj, n = 10, with_ties = FALSE) %>%
#   ungroup() %>%
#   mutate(pathway_plot = str_wrap(pathway, width = 42)) %>%
#   arrange(enriched_in, NES) %>%
#   mutate(pathway_plot = factor(pathway_plot, levels = pathway_plot))
# 
# plot <- ggplot(fg.plot, aes(x = NES, y = pathway_plot, fill = enriched_in)) +
#   geom_col(width = 0.75) +
#   geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5) +
#   scale_fill_manual(values = c("HF" = "#F393C3", "NF" = "#00A6ED")) +
#   labs(
#     x = "Normalized enrichment score",
#     y = NULL,
#     fill = "Enriched in"
#   ) +
#   theme_classic(base_size = 8) +
#   theme(
#     legend.position = "top",
#     axis.text.y = element_text(size = 7),
#     axis.ticks.y = element_blank()
#   )
# 
# ggsave("reactome_gsea.jpg", plot)


#### GSEA on cell atlas ####

gene_id_col <- "gene.ID"
gene_set_dir <- "./gene_sets"  

xp_to_loc <- read.delim("~/haifa/S.pistillata/rna/genomes/xp_to_loc.tsv", header = FALSE, stringsAsFactors = FALSE,
                        col.names = c("XP", "Symbol")) %>%
  mutate(XP = as.character(XP), Symbol = as.character(Symbol)) %>%
  filter(!is.na(XP), !is.na(Symbol), XP != "", Symbol != "") %>%
  distinct(XP, Symbol)

auto_gene_sets <- list.files(gene_set_dir, pattern = "\\.tsv$", full.names = TRUE) %>%
  set_names(~ tools::file_path_sans_ext(basename(.x))) %>%
  map(~ read.delim(.x, stringsAsFactors = FALSE) %>%
        transmute(
          original_id = as.character(.data[[gene_id_col]]),
          XP = str_replace(original_id, "^Spis_(XP_[0-9]+)_([0-9]+)$", "\\1.\\2")
        ) %>%
        left_join(xp_to_loc, by = "XP") %>%
        pull(Symbol) %>%
        unique() %>%
        na.omit() %>%
        as.character())

names(auto_gene_sets)
lengths(auto_gene_sets)

biomin_gene_set <- readLines(
  file.path(gene_set_dir, "genes.biomin.accessions.txt"),
  warn = FALSE
) %>%
  str_trim() %>%
  discard(~ is.na(.x) || .x == "") %>%
  unique()

auto_gene_sets[["biomineralization"]] <- biomin_gene_set

names(auto_gene_sets)
lengths(auto_gene_sets)
# if we need to select certain sets
# sets_to_keep <- c(
#   "oocytes",
#   "calicoblasts",
#   "gastro_algae",
#   "biomin toolkit",
#   "biomin proteome"
# )
# 
# bio_gene_set <- bio_gene_set[sets_to_keep]


run_fgsea <- function(res, gene_sets, contrast_name = "HF_vs_NF") {
  res_df <- as.data.frame(res) %>%
    filter(!is.na(Symbol), is.finite(stat)) %>%
    group_by(Symbol) %>%
    slice_max(abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  ranks <- setNames(res_df$stat, res_df$Symbol)
  ranks <- sort(ranks, decreasing = TRUE)
  
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

celltype_fgsea <- run_fgsea(
  res = res,
  gene_sets = auto_gene_sets,
  contrast_name = "HF_vs_NF"
)


all_results_final <- celltype_fgsea %>%
  mutate(global_padj = p.adjust(pval, method = "BH"))

signif_fgsea <- all_results_final %>%
  filter(global_padj < 0.05)

signif_fgsea.arranged <- signif_fgsea %>%
  as_tibble() %>%
  mutate(
    enriched_in = if_else(NES > 0, "HF", "NF"),
    leadingEdge = map_chr(leadingEdge, paste, collapse = ";")
  ) %>%
  arrange(padj)

write.csv(
  signif_fgsea.arranged,
  "GSEA_cell.atlas_HF_vs_NF.csv",
  row.names = FALSE
)

# 2. Save as a standard CSV

plot_fgsea <- signif_fgsea %>%
  mutate(
    enriched_in = if_else(NES > 0, "HF", "NF"),
    pathway_plot = str_replace_all(pathway, "_", " "),
    pathway_plot = str_wrap(pathway_plot, width = 30)
  ) %>%
  arrange(enriched_in, NES) %>%
  mutate(pathway_plot = factor(pathway_plot, levels = rev(unique(pathway_plot))))

fgsea.cell.atlas <- ggplot(plot_fgsea, aes(x = NES, y = pathway_plot, color = enriched_in,
                       size = -log10(global_padj))) +
  geom_point(alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5) +
  scale_color_manual(values = c("HF" = "#F393C3", "NF" = "#00A6ED")) +
  scale_size_continuous(range = c(2.5, 7)) +
  labs(
    x = "Normalized enrichment score",
    y = NULL,
    color = "Enriched in",
    size = expression(-log[10](adjusted~P))
  ) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 7),
    axis.ticks.y = element_blank()
  )

ggsave("cell.atlas.gsea.jpg", fgsea.cell.atlas, width = 6, height = 4)



