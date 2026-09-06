library(tidyverse)
library(dplyr)
library(tidyr)
library(purrr)
library(GO.db)

setwd("/home/gospozha/haifa/hiba/pl_align/")

res <- read.csv("res_annot.csv", check.names = FALSE, row.names = 1)

# GO annotation table
term2gene <- read_tsv("clusterprofiler/term2gene.tsv")
term2name <- read_tsv("clusterprofiler/term2name_GO.tsv")


#### GO slims ####
category_roots <- list(
  
  ECM_adhesion = c(
    "GO:0031012", # extracellular matrix
    "GO:0030198", # extracellular matrix organization
    "GO:0007155", # cell adhesion
    "GO:0005201"  # extracellular matrix structural constituent
  ),
  
  cytoskeleton = c(
    "GO:0007010", # cytoskeleton organization
    "GO:0005856", # cytoskeleton
    "GO:0003774"  # cytoskeletal motor activity
  ),
  
  metabolism = c(
    "GO:0006091", # generation of precursor metabolites and energy
    "GO:0005975", # carbohydrate metabolic process
    "GO:0006629", # lipid metabolic process
    "GO:0006520",  # cellular amino acid metabolic process
    "GO:0008028"  # monocarboxylic acid transmembrane transporter activity
  ),
  
  morphogenesis = c(
    "GO:0009653", # anatomical structure morphogenesis
    "GO:0032502", # developmental process
    "GO:0048856", # anatomical structure development
    "GO:0007275", # multicellular organism development
    "GO:0007389"  # pattern specification process
  ),
  
  redox = c(
    "GO:0045454", # cell redox homeostasis
    "GO:0006979", # response to oxidative stress
    "GO:0098869", # cellular oxidant detoxification
    "GO:0016209", # antioxidant activity
    "GO:0004497", # monooxygenase activity
    "GO:0016705",  # oxidoreductase activity involving molecular oxygen
    "GO:0004601", # peroxidase activity
    "GO:0016684", # oxidoreductase activity, acting on peroxide as acceptor
    "GO:0055114",  # oxidation-reduction process
    "GO:0020037" # heme binding
  ),
  
  immunity = c(
    "GO:0002376", # immune system process
    "GO:0006952", # defense response
    "GO:0004866"  # endopeptidase inhibitor activity
  ),
  
  cell_cycle_apoptosis = c(
    "GO:0007049", # cell cycle
    "GO:0051276", # chromosome organization
    "GO:0006260", # DNA replication
    "GO:0006281", # DNA repair
    "GO:0006915", # apoptotic process
    "GO:0042981", # regulation of apoptotic process
    "GO:0097190", # apoptotic signaling pathway
    "GO:0012501"  # programmed cell death
  ),
  
  neuro_sensing = c(
    "GO:0050877", # nervous system process
    "GO:0007600",  # sensory perception
    "GO:0022834", # ligand-gated channel activity / gated channel branch
    "GO:0006816", # calcium ion transport
    "GO:0019722",  # calcium-mediated signaling
    "GO:0004930", # G protein-coupled receptor activity
    "GO:0007186"  # G protein-coupled receptor signaling pathway
  )
)


# get descendants of parent terms
get_descendants <- function(go_id) {
  
  term <- GOTERM[[go_id]]
  
  if (is.null(term))
    return(go_id)
  
  ont <- Ontology(term)
  
  descendants <- switch(
    ont,
    BP = GOBPOFFSPRING[[go_id]],
    CC = GOCCOFFSPRING[[go_id]],
    MF = GOMFOFFSPRING[[go_id]]
  )
  
  descendants <- descendants[!is.na(descendants)]
  
  unique(c(go_id, descendants))
}

category_GO <- lapply(
  category_roots,
  function(roots) {
    unique(unlist(lapply(roots, get_descendants)))
  }
)




#### classify ####


genes_to_classify <- unique(res$Symbol)

gene_go <- term2gene %>%
  filter(
    term %in% term2name$term,
    gene %in% genes_to_classify
  ) %>%
  distinct(term, gene)


str(category_GO)


go_categories <- imap_dfr(
  category_GO,
  ~ tibble(
    term = .x,
    category = .y
  )
)

# category overlap
go_categories <- go_categories %>%
  distinct(term, category)
go_categories %>%
  count(term) %>%
  filter(n > 1)

go_categories %>%
  group_by(term) %>%
  filter(n_distinct(category) > 1) %>%
  arrange(term)

go_categories %>%
  group_by(term) %>%
  summarise(
    categories = paste(sort(unique(category)), collapse = "; "),
    n_categories = n_distinct(category),
    .groups = "drop"
  ) %>%
  filter(n_categories > 1) %>%
  count(categories, sort = TRUE)

# combining category and gene
gene_categories_long <- gene_go %>%
  inner_join(go_categories, by = "term")

gene_categories <- gene_categories_long %>%
  group_by(gene) %>%
  summarise(
    category = paste(sort(unique(category)), collapse = "; "),
    .groups = "drop"
  )

# define primary category
category_scores <- gene_categories_long %>%
  distinct(gene, term, category) %>%
  count(gene, category, name = "n_GO")

primary_categories <- category_scores %>%
  group_by(gene) %>%
  filter(n_GO == max(n_GO)) %>%
  mutate(
    ambiguous = n() > 1
  ) %>%
  ungroup()

primary_categories %>%
  count(ambiguous)

# no match to other 
gene_categories <- tibble(gene = genes_to_classify) %>%
  left_join(gene_categories, by = "gene") %>%
  mutate(
    category = replace_na(category, "other")
  )

# add to res
res.categories <- res %>%
  left_join(
    gene_categories,
    by = c("Symbol" = "gene")
  )

# add primary category

res.categories <- res.categories %>%
  left_join(
    primary_categories %>%
      dplyr::select(gene, category) %>%
      dplyr::rename(primary_category = category),
    by = c("Symbol" = "gene")
  )%>%
  dplyr::mutate(
    primary_category = tidyr::replace_na(primary_category, "other")
  )

write.csv(
  res.categories,
  "res.annot.categorized.csv",
  row.names = FALSE
)