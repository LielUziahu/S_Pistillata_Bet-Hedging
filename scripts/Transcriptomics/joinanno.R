library(dplyr)
library(tidyr)
library(readr)
library(stringr)


setwd("~/haifa/S.pistillata/rna/genomes/anno/")

eggnog <- read_tsv("Spis_eggnog.emapper.annotations", comment = "#")
length(unique(eggnog$query)) 
ipr <- read_tsv(
  "GCF_002571385.2_Stylophora_pistillata_v1.1_protein.faa.tsv",col_names = c(
    "query","md5","length","analysis","signature_acc","signature_desc",
    "start","end","score","status","date",
    "ipr_acc","ipr_desc","go_terms","pathways"
  )
)
length(unique(ipr$query)) 

# FILTER 
ipr_filt <- ipr %>%
  filter(
    # Pfam / SMART / TIGRFAMs: keep confident HMM hits
    (analysis %in% c("Pfam", "SMART", "TIGRFAMs") & score <= 1e-5) |
      
      # SUPERFAMILY: very permissive → tighten
      (analysis == "SUPERFAMILY" & score <= 1e-10) |
      
      # PANTHER: score often 0, keep as-is
      (analysis == "PANTHER") |
      
      # CDD: curated, usually fine
      (analysis == "CDD")
  )%>%
  filter((end - start + 1) >= 30)

# How much did filtering remove?
nrow(ipr) - nrow(ipr_filt)

# Any proteins lost entirely?
sum(!ipr$query %in% ipr_filt$query)


ipr_collapsed <- ipr_filt %>%
  group_by(query) %>%
  summarise(
    ipr_accessions = paste(unique(na.omit(ipr_acc)), collapse = ";"),
    ipr_descriptions = paste(unique(na.omit(ipr_desc)), collapse = ";"),
    domains = paste(unique(paste0(analysis, ":", signature_acc)), collapse = ";"),
    domain_names = paste(unique(na.omit(signature_desc)), collapse = ";"),
    ipr_go = paste(unique(unlist(strsplit(na.omit(go_terms), "\\|"))), collapse = ","),
    ipr_pathways = paste(unique(unlist(strsplit(na.omit(pathways), "\\|"))), collapse = "|"),
    .groups = "drop"
  )

# clean GO
clean_go <- function(x) {
  if (is.na(x)) return(NA_character_)
  
  x %>%
    strsplit(",") %>%
    unlist() %>%
    # remove "-" and empty
    .[. != "-" & . != ""] %>%
    # remove "(PANTHER)" "(InterPro)" etc
    gsub("\\(.*?\\)", "", .) %>%
    # keep only valid GO IDs
    .[grepl("^GO:\\d{7}$", .)] %>%
    unique() %>%
    paste(collapse = ",")
}

clean_semicolon <- function(x) {
  if (is.na(x)) return(NA_character_)
  
  x %>%
    strsplit(";") %>%
    unlist() %>%
    .[. != "-" & . != ""] %>%
    unique() %>%
    paste(collapse = ";")
}

clean_pathways <- function(x) {
  if (is.na(x)) return(NA_character_)
  
  x %>%
    strsplit("\\|") %>%
    unlist() %>%
    .[. != "-" & . != ""] %>%
    unique() %>%
    paste(collapse = "|")
}


ipr_collapsed2 <- ipr_collapsed %>%
  mutate(ipr_go = vapply(ipr_go, clean_go, character(1)),
         domains = vapply(domains, clean_semicolon, character(1)),
         domain_names = vapply(domain_names, clean_semicolon, character(1)),
         ipr_accessions = vapply(ipr_accessions, clean_semicolon, character(1)),
         ipr_pathways = vapply(ipr_pathways, clean_pathways, character(1)))


# merge annotations
annotation_merged <- eggnog %>%
  left_join(ipr_collapsed2, by = "query")%>%
  dplyr::select("query", "Description", "GOs", "ipr_descriptions", "ipr_go")%>%
  mutate(
    GOs    = na_if(GOs, "-"),
    ipr_go = na_if(ipr_go, "-"),
    
    # remove any " -;" fragments and standalone "-" fragments
    ipr_go = str_replace_all(ipr_go, "(^|,)-;?", "\\1"),
    ipr_go = str_replace_all(ipr_go, ",{2,}", ","),   # compress duplicate commas
    ipr_go = str_replace_all(ipr_go, "^,|,$", ""),    # trim leading/trailing commas
    
    ipr_descriptions = str_replace_all(ipr_descriptions, "(^|,)-;?", "\\1"),
    ipr_descriptions = str_replace_all(ipr_descriptions, ",{2,}", ","),   # compress duplicate commas
    ipr_descriptions = str_replace_all(ipr_descriptions, "^,|,$", ""),    # trim leading/trailing commas
    
    # remove origin tags after GO term
    GOs    = str_replace_all(GOs, "\\([^\\)]*\\)", ""),
    ipr_go = str_replace_all(ipr_go, "\\([^\\)]*\\)", "")
  )

annotation_merged %>% summarise(
  n_GOs_na = sum(is.na(GOs)),
  n_ipr_go_na = sum(is.na(ipr_go)),
  example_ipr = first(na.omit(ipr_go))
)


annotation_combined <- annotation_merged %>%
  rowwise() %>%
  mutate(
    GO_merged = {
      combined <- paste(c(GOs, ipr_go), collapse = ",")
      combined <- str_replace_all(combined, "\\s+", "")  # remove spaces
      parts <- unlist(strsplit(combined, ",", fixed = TRUE), use.names = FALSE)
      parts <- parts[parts != "" & parts != "-"]
      parts <- unique(parts)
      paste(parts, collapse = ",")
    }
  ) %>%
  ungroup()


annotation_combined$GO_merged <- gsub("(^NA,|,NA$|,NA,|^NA$)", "", annotation_combined$GO_merged)

annotation_combined <- annotation_combined %>%
  dplyr::select("query", "GO_merged")


write_tsv(annotation_combined, "protein_annotations_eggnog_interpro_merged.tsv")

anno <- read_tsv("protein_annotations_eggnog_interpro_merged.tsv")
View(anno)
colnames(anno) <- c("protein", "GO")
# convert to LOC
loc <- read_tsv("../xp_to_loc.tsv",  col_names = FALSE)
View(loc)
colnames(loc) <- c("protein", "gene")

anno_loc <- anno %>%
  left_join(loc, by = "protein") %>%        # add LOC
  mutate(GO = str_replace_all(GO, "\\s+", "")) %>%
  separate_rows(GO, sep = "[,;/|]") %>%      # one GO per row
  distinct(gene, GO) %>%                     # avoid duplicates
  group_by(gene) %>%
  summarise(
    GO = paste(sort(unique(GO)), collapse = ","),
    .groups = "drop"
  )

write_tsv(anno_loc, "gene_annotations_eggnog_interpro_merged.tsv")

length(unique(loc$protein))
length(unique(loc$gene))
length(unique(anno$protein)) 
sum(unique(anno$protein) %in% unique(loc$protein))

term2gene <- anno %>%
  mutate(GO = str_replace_all(GO, "\\s+", "")) %>%
  separate_rows(GO, sep = "[,;/|]") %>%
  filter(str_detect(GO, "^GO:\\d{7}$")) %>%
  inner_join(loc, by = "protein") %>%
  distinct(term = GO, gene)
View(term2gene)

write_tsv(term2gene, "term2gene.tsv")
