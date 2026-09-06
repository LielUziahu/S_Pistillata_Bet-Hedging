library(goseq)
library(tidyverse)
library(GSEABase)               #BiocManager::install("GSEABase")
library(data.table)
library(ggplot2)
library(cowplot)                #install.packages("cowplot")
library(patchwork)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(ontologyIndex)
library(GOSemSim)
library(simplifyEnrichment)
library(org.Hs.eg.db)
library(rlang)
library(readr)


#save.image(file="220126.GO.planulae.RData") 
load("/home/gospozha/haifa/hiba/pl_align/clusterprofiler/220126.GO.planulae.RData")

# setting working directory 
setwd("/home/gospozha/haifa/hiba/pl_align/clusterprofiler/")

library(GO.db)
library(dplyr)

term2name <- AnnotationDbi::select(
  GO.db,
  keys = unique(term2gene$term),
  columns = c("TERM"),
  keytype = "GOID"
) %>%
  rename(term = GOID, name = TERM) %>%
  distinct()

write_tsv(x = term2name, file = "term2name_GO.tsv")


term2gene <- read_tsv("term2gene.tsv")
term2name <- read_tsv("term2name_GO.tsv")


term2gene <- term2gene %>%
  filter(term %in% term2name$term)

#### reading count matrix from a file ####
countData  <- read.csv2('../CountMatrix.csv', header=TRUE, row.names=1, sep=',', check.names = F)
# reading metadata file
MetaData <- read.csv2('../Metadata', header=TRUE, sep="\t")

MetaData$condition <- as.factor(MetaData$condition)
#MetaData$origin <- as.factor(MetaData$sample)

MetaData <- MetaData[match(colnames(countData), MetaData$id), ]

countData$geneID <- rownames(countData)
smallestGroupSize <- 3
keep <- rowSums(countData >= 5) >= smallestGroupSize
countData <- countData[keep,]

background_genes <- countData %>%
  dplyr::select("geneID") %>%
  unlist() %>%
  as.vector()

file_path <- ("../res_annot_stranded.csv")
interesting_set <- read_csv(file_path, show_col_types = FALSE) %>%               # if you want only significant genes
  pull(Symbol) %>%                        # <-- THIS is your LOC column
  unique() %>%
  na.omit()
  
  #### Run enrichment ####
interesting_set_FL <- read_csv(file_path, show_col_types = FALSE) %>%     
  filter(padj < 0.05, log2FoldChange > 0) %>%
  pull(Symbol) %>%                       
  unique() %>%
  na.omit()

interesting_set_NF <- read_csv(file_path, show_col_types = FALSE) %>%     
  filter(padj < 0.05, log2FoldChange < 0) %>%
  pull(Symbol) %>%                       
  unique() %>%
  na.omit()

enrichment_FL <- enricher(interesting_set_FL,
                         TERM2GENE = term2gene,
                         TERM2NAME = term2name,
                         pvalueCutoff = 0.05,
                         universe = background_genes,
                         pAdjustMethod = "fdr",
                         qvalueCutoff = 0.2)
  
  # Save enrichment results
write_csv(enrichment_FL@result, "enrichment_results.lfc1.FL.csv")

enrichment_NF <- enricher(interesting_set_NF,
                          TERM2GENE = term2gene,
                          TERM2NAME = term2name,
                          pvalueCutoff = 0.05,
                          universe = background_genes,
                          pAdjustMethod = "fdr",
                          qvalueCutoff = 0.2)

# Save enrichment results
write_csv(enrichment_NF@result, "enrichment_results.lfc1.NF.csv")
# not significant

p <- dotplot(enrichment,
                            x = "geneRatio",
                            color = "p.adjust",
                            orderBy = "x",
                            showCategory = 100,
                            font.size = 7) 

ggsave(paste0("enrichment_dotplot.pdf"), p, width = 4, height = 9)
             

#### bar plot ####
deg_results <- list(
  FLvNF   = read_csv("../res_annot.csv")
)

# Enrichment results
enrich_results <- list(
  FLvNF   = read_csv("enrichment_results.lfc1.FL.csv")
)


# Function to compute mean logFC per GO term
compute_mean_logfc <- function(enrich_df, deg_df) {
  
  deg_small <- deg_df %>%
    dplyr::select(Symbol, log2FoldChange)
  
  enrich_df %>%
    dplyr::select(ID, Description, geneID, p.adjust) %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    dplyr::left_join(deg_small, by = c("geneID" = "Symbol")) %>%
    dplyr::group_by(ID, Description, p.adjust) %>%
    dplyr::summarise(
      mean_logFC = mean(log2FoldChange, na.rm = TRUE),
      n_genes_matched = sum(!is.na(log2FoldChange)),
      .groups = "drop"
    )
}

# Apply to all comparisons
mean_logfc_list <- map2(enrich_results, deg_results, compute_mean_logfc)

# Add comparison name
mean_logfc_named <- map2_dfr(mean_logfc_list, names(mean_logfc_list), ~mutate(.x, contrast = .y))


# You can optionally filter top N GO terms across all contrasts
top_terms <- mean_logfc_named %>%
  mutate(z_logFC = scale(mean_logFC)[,1]) %>%
  dplyr::filter(p.adjust<0.05) %>%
  arrange(-abs(mean_logFC))

go_terms <- unique(top_terms$ID)

top_terms <- mean_logfc_named %>%
  mutate(z_logFC = scale(mean_logFC)[,1]) %>%
  dplyr::filter(p.adjust<0.05)%>%
  mutate(Description_ordered = factor((Description), levels = rev
                                      (unique(Description))))
write_csv(top_terms, "enrichment_results.lfc1.FL.top.csv")

# Plot barplot
# either z score of logfc or pure logfc

clust.barplot <- ggplot(
  top_terms,
  aes(x = Description_ordered,
      y = mean_logFC,
      fill = contrast,
      shape = contrast)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.6) +
  geom_point(size = 3, position = position_dodge(width = 0.8), color = "black") +
  coord_flip() +
  #facet_grid(cluster ~ ., scales = "free_y", space = "free_y") +
  scale_y_continuous(
    name = "Mean log2FC\n(Down in S  ←  0  →  Up in S)",
    breaks = scales::breaks_width(5)
  ) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 35), name = "GO term")+
  #scale_fill_manual(values = contrast_colors) +
  #scale_shape_manual(values = contrast_shapes) +
  theme_minimal() +
  guides(fill = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE)) +
  labs(title = "GO term enrichment") +
  theme(
    plot.title = element_text(size = 13, margin = margin(t = 9, b = 6)),
    plot.margin = margin(10, 10, 8, 8),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    legend.position = "right"
  )

clust.barplot

#### go slim ####
library(GSEABase)
# 1. Read enrichment results
enrich_results <- read_csv("enrichment_results.lfc1.FL.csv")

# 2. Combine into one dataframe

all_enrich <- enrich_results %>%
  select(ID, Description, Count, p.adjust) %>%
  filter(!is.na(p.adjust), p.adjust < 0.05) %>%   # <-- filter here
  mutate(
    logp = -log10(p.adjust),
    Description = str_trunc(Description, 60)
  ) %>%
  left_join(mean_logfc_named %>% select(ID, mean_logFC, n_genes_matched),
              by = "ID")%>%
  mutate(
    mean_logFC = as.numeric(mean_logFC),
    n_genes_matched = as.numeric(n_genes_matched)
  )

# Load GO ontology
go <- get_ontology("http://purl.obolibrary.org/obo/go.obo",
                   extract_tags = "everything")

# Load generic GO slim
goslim <- get_ontology("http://current.geneontology.org/ontology/subsets/goslim_generic.obo",
                       extract_tags = "everything")

# 3. Map enriched GO terms to GO slim parents
# Each GO term maps to one or more slim terms via ancestors
slim_map <- lapply(all_enrich$ID, function(go_id) {
  if (!go_id %in% names(go$ancestors)) return(NA_character_)
  ancestors <- unique(c(go_id, go$ancestors[[go_id]]))  # include self
  slim_hits <- intersect(ancestors, names(goslim$name))
  if (length(slim_hits) == 0) return(NA_character_)
  slim_hits[1]
})

all_enrich$GOSLIM <- unlist(slim_map)

# 4. Aggregate by slim term
slim_agg <- all_enrich %>%
  mutate(
    GOSLIM_name = goslim$name[GOSLIM],
    plot_term   = if_else(is.na(GOSLIM_name), Description, GOSLIM_name),
  ) %>%
  group_by(plot_term) %>%
  summarise(
    mean_logp = mean(logp, na.rm = TRUE),
    n_terms      = n(),
    mean_logFC   = mean(mean_logFC, na.rm = TRUE),
    #mean_logFC_w = weighted.mean(mean_logFC, w = w, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_terms))

write_csv(slim_agg, "GOslim.lfc1.FL.top.csv")


slim_agg %>%
  slice_max(n_terms, n = 10) %>%
  mutate(plot_term = factor(plot_term, levels = rev(plot_term))) %>%
  ggplot(aes(x = n_terms, y = plot_term)) +
  geom_col() +
  labs(x = "# enriched GO terms in slim bucket", y = NULL) +
  theme_classic() -> GOplot
ggsave("GOslim.top10.jpg", GOplot, width = 10, height = 10)  
