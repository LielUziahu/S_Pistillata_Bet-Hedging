library(tidyverse)
library(readr)

setwd("~/haifa/hiba/pl_align/symbio/")
# Path to quant.sf files
quant_dir <- "./salmon"  

# Read files
files <- list.files(quant_dir, pattern = "\\.quant\\.sf$", full.names = TRUE)
stopifnot(length(files) > 0)

samples <- basename(files) |>
  sub("\\.quant\\.sf$", "", x = _)
names(files) <- samples

read_quant <- function(f, sample_name) {
  read_tsv(f, show_col_types = FALSE) |>
    mutate(sample = sample_name)
}

quant_long <- purrr::map2_dfr(files, names(files), read_quant)

# Parse "species" from transcript name: take text before first underscore
quant_long <- quant_long |>
  mutate(
    species = sub("_.*", "", Name),
    # Flag names that dont contain underscore 
    has_underscore = grepl("_", Name)
  )

# Sanity checks
cat("Total transcripts loaded:", nrow(quant_long), "\n")
cat("Fraction with underscore in Name:", mean(quant_long$has_underscore), "\n\n")
cat("Species tags found:\n")
print(sort(unique(quant_long$species)))

# # Summarize per sample x species
# species_summary <- quant_long |>
#   group_by(sample, species) |>
#   summarise(
#     TPM = sum(TPM, na.rm = TRUE),
#     NumReads = sum(NumReads, na.rm = TRUE),
#     n_tx = n(),
#     .groups = "drop"
#   ) |>
#   group_by(sample) |>
#   mutate(
#     prop_TPM = TPM / sum(TPM),
#     prop_reads = NumReads / sum(NumReads)
#   ) |>
#   ungroup()
# 
# # OPTIONAL: exclude some categories (edit this list if you added HOST, etc.)
# exclude_species <- c("HOST")  # add decoys if you want, e.g. c("HOST","DecoyX")
# plot_df <- species_summary |>
#   filter(!species %in% exclude_species)
# 
# # Plot: proportions by TPM
# p1 <- ggplot(plot_df, aes(x = sample, y = prop_TPM, fill = species)) +
#   geom_col() +
#   coord_flip() +
#   labs(x = NULL, y = "Proportion of TPM", title = "Symbiont composition (TPM)") +
#   theme_bw()
# 
# # Plot: absolute TPM (log10) to see weak tails
# p2 <- ggplot(plot_df, aes(x = sample, y = TPM, fill = species)) +
#   geom_col() +
#   coord_flip() +
#   scale_y_continuous(trans = "log10") +
#   labs(x = NULL, y = "TPM (log10)", title = "Symbiont signal (TPM, log10)") +
#   theme_bw()
# 
# print(p1)
# print(p2)
# 
# # Wide tables for export
# species_wide_prop <- species_summary |>
#   select(sample, species, prop_TPM) |>
#   pivot_wider(names_from = species, values_from = prop_TPM, values_fill = 0)
# 
# species_wide_TPM <- species_summary |>
#   select(sample, species, TPM) |>
#   pivot_wider(names_from = species, values_from = TPM, values_fill = 0)
# 
# # Save outputs
# write_csv(species_summary, file.path(quant_dir, "symbiont_species_summary_long.csv"))
# write_csv(species_wide_prop, file.path(quant_dir, "symbiont_species_propTPM_wide.csv"))
# write_csv(species_wide_TPM, file.path(quant_dir, "symbiont_species_TPM_wide.csv"))

# Summarise READ COUNTS per species per sample
species_reads <- quant_long |>
  group_by(sample, species) |>
  summarise(
    reads = sum(NumReads, na.rm = TRUE),
    .groups = "drop"
  )

# Convert to proportions (within each sample)
species_prop <- species_reads |>
  group_by(sample) |>
  mutate(
    total_reads = sum(reads),
    prop_reads = reads / total_reads
  ) |>
  ungroup()

# Optional: remove unnecessary samples
# here i dont remove anything
exclude_species <- c("HOST")  
species_prop_plot <- species_prop |>
  filter(!species %in% exclude_species)

# Plot: proportions
ggplot(species_prop_plot,
       aes(x = sample, y = prop_reads, fill = species)) +
  geom_col() +
  coord_flip() +
  labs(
    x = NULL,
    y = "Proportion of mapped reads",
    title = "Symbiont species composition"
  ) +
  theme_bw()

# Wide table 
species_prop_wide <- species_prop_plot |>
  select(sample, species, prop_reads) |>
  pivot_wider(
    names_from = species,
    values_from = prop_reads,
    values_fill = 0
  )

# Save outputs
write_csv(species_prop_plot,
          file.path(quant_dir, "symbiont_species_read_proportions_long.csv"))

write_csv(species_prop_wide,
          file.path(quant_dir, "symbiont_species_read_proportions_wide.csv"))

mlt2 <- species_prop_plot %>%
  mutate(
    #  order samples 
    sample = factor(sample, levels = unique(sample)),
    # order species by overall abundance
    species = factor(
      species,
      levels = species_prop_plot %>%
        group_by(species) %>%
        summarise(mean_prop = mean(prop_reads)) %>%
        #arrange(desc(mean_prop)) %>%
        pull(species)
    )
  )
hm1 <- ggplot(
  data = mlt2,
  mapping = aes(
    x = sample,
    y = species,
    fill = prop_reads
  )
) +
  geom_raster() +
  labs(
    #title = "Proportion of reads mapped to symbiont databases",
    x = "Sample",
    y = NULL,
    fill = "Proportion"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, size = 8, vjust = 0.8, hjust = 0.8),
    panel.grid = element_blank()
  ) +
  scale_fill_gradientn(
    colours = c("white", "darkorange1", "darkred"),
    limits = c(0, 0.5),
    oob = scales::squish
  ) +
  scale_y_discrete(
    labels = c(
      "Dtrenchii" = "Durusdinium trenchii",
      "symb" = "Breviolum minutum",
      "syma" = "Symbiodinium clade A3",
      "Stri" = "Symbiodinium tridacnidorum",
      "SPilosum" = "Symbiodinium pilosum",
      "Snec" = "Symbiodinium necroappetens",
      "Snat" = "Symbiodinium natans",
      "SMicUQSCI" = "Symbiodinium microadriaticum (UQ)",
      "SMicUQCass" = "Symbiodinium microadriaticum (UQ Cass)",
      "SMicReefgen" = "Symbiodinium microadriaticum (ReefGenomics)",
      "Slin" = "Symbiodinium linuacheae",
      "SKawaF" = "Fugacium kawagutii",
      "SGoreauC" = "Cladocopium goreaui",
      "CladC15" = "Cladocopium C15",
      "BrevMin" = "Breviolum minutum",
      "SymbA3" = "Symbiodinium A3"
    )
  )

hm1
ggsave("symb.prop.jpg", hm1, width = 10, height = 6)  
ggsave("symb.prop.pdf", hm1, width = 10, height = 6)  
