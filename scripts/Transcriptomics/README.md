### Filtering, mapping, quality control, annotation, settlement genes identification: see [bash_scripts](https://github.com/LielUziahu/S_Pistillata_Bet-Hedging/tree/main/scripts/Transcriptomics/bash_scripts) folder  
### R analysis:

- *DESeq2.planulae.R*: DE analysis results, DE heatmap, GFP genes, biomineralization genes (uses [*S. pistillata* candidate biomineralization genes](genes.biomin.accessions.txt) and [*S. pistillata* NCBI Annotation](sp_genes.tsv))
- *joinanno.R*: combining eggnog and interpro annotations
- *clusterprof.planulae.R*: GO and GOslim analyses (uses newly obtained [GO annotations](gene_annotations_eggnog_interpro_merged.tsv))
- *symbio.salmon.R*: symbiont proportion calculation
- *GSEA.R*: GSEA enrichment of Reactome pathways, biomineralization set and cell atlas gene sets  
- *settlement.R*: GSEA enrichment of settlement-related gene sets and overlap with DE results
- *classification.R*: GO-based classification of DE results
