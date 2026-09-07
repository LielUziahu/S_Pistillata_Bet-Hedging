### QC and filtering 
- FastQC check: *fastqc.sh* and *Sequencing_adaptors.fasta*
- Remove adaptors and low-quality bases: *trimmomatic.multi.sh*
- Combine QC metrics into reports: *multiqc.sh*

### rRNA filtering and QC
- Detect and filter rRNA: *ribodetector.sh*
- Calculate rRNA statistics: *calc.perc.sh*

### Mapping to the host genome
- Change meaningless ids to sample ids: *ids.csv* and *rename.sh*
- STAR mapping: *mapping.basic.sh*

### Mapping to the symbiont transcriptome database
- Salmon indexing and mapping: *mapping.salmon.sh*

### *S. pistillata* genome annotation 
- Annotation through eggnog: *eggnog.sh*
- Annotation through interproscan: *iprscan.sh*

### Settlement genes identification
- *S. pistillata* protein to gene conversion: *xp_to_loc.tsv*
- Gene sets: *transdec.sh*, *nucleotide_proteinort.sh*
- Protein sets: *protein_proteinort.sh*
