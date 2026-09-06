#!/bin/bash
#################################################################################################################
#SBATCH --job-name=sbatchTemplate ## Name of your job
#SBATCH --ntasks=16 ## number of cpu's to allocate for a job
#SBATCH --ntasks-per-node=16 ## number of cpu's to allocate per each node
#SBATCH --nodes=1 ## number of nodes to allocate for a job
#SBATCH --mem=64G ## memory to allocate for your job in MB
#SBATCH --time=1-00:00:00 ## time to allocate for your job in format: DD-HH:MM:SS
#SBATCH --error=%J.errors ## stderr file name(The %J will print job ID number)
#SBATCH --output=%J.output ## stdout file name(The %J will print job ID number)
#SBATCH --mail-type=NONE ## Send your job status via e-mail: Valid type values are NONE, BEGIN, END, FAIL, REQUEUE, ALL
########### Job information #############
echo "================================"
echo "Start at `date`"
echo "Job id is $SLURM_JOBID"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NTASKS processors."
echo "================================"
#########################################

######## Load required modules ##########
#. /etc/profile.d/modules.sh # Required line for modules environment to work
#module load openmpi/1.8.4 python/2.7 # Load modules that are required by your program
source /lustre1/home/mass/eskalon/miniconda/bin/activate proteinort
#########################################

### Below you can enter your program job command ###

#!/usr/bin/env bash

set -euo pipefail

STYLO_CDS="${1:-cds_from_genomic.fna}"
THREADS="${2:-16}"

EVALUE="1e-5"
OUTDIR="proteinortho_nucleotide_results"
CLEAN_STYLO="${OUTDIR}/Stylophora_cds_LOC.fna"

QUERY_FILES=(
    "Meyer2011.down_after_CCA.transcripts.fasta"
    "Meyer2011.up_after_CCA.transcripts.fasta"
    "Walker2019.probing.transcripts.fasta"
)

[[ -f "$STYLO_CDS" ]] || {
    echo "Error: Stylophora CDS file not found: $STYLO_CDS" >&2
    exit 1
}

command -v proteinortho >/dev/null 2>&1 || {
    echo "Error: proteinortho is not available in PATH." >&2
    exit 1
}

command -v python3 >/dev/null 2>&1 || {
    echo "Error: python3 is not available in PATH." >&2
    exit 1
}

mkdir -p "$OUTDIR"

echo "Preparing Stylophora CDS headers..."

python3 - "$STYLO_CDS" "$CLEAN_STYLO" <<'PY'
import re
import sys

input_fasta, output_fasta = sys.argv[1:3]

seen = set()
written = 0
missing_gene = 0
duplicate_ids = 0

with open(input_fasta, encoding="utf-8", errors="replace") as source, \
     open(output_fasta, "w", encoding="utf-8") as output:

    for line_number, line in enumerate(source, start=1):

        if not line.startswith(">"):
            output.write(line)
            continue

        header = line[1:].strip()

        # First token in the original header:
        # lcl|NW_019217784.1_cds_XP_022785885.1_1
        first_token = header.split(maxsplit=1)[0]

        # Extract gene from [gene=LOC111326177]
        match = re.search(r"\[gene=([^\]]+)\]", header)

        if match:
            gene = match.group(1).strip()
            new_id = f"{first_token}_{gene}"
        else:
            missing_gene += 1
            new_id = first_token
            print(
                f"Warning: no [gene=...] field on FASTA header line "
                f"{line_number}: {first_token}",
                file=sys.stderr,
            )

        # Proteinortho requires unique sequence identifiers.
        if new_id in seen:
            duplicate_ids += 1
            suffix = 2
            unique_id = f"{new_id}_copy{suffix}"

            while unique_id in seen:
                suffix += 1
                unique_id = f"{new_id}_copy{suffix}"

            print(
                f"Warning: duplicate ID {new_id}; renamed to {unique_id}",
                file=sys.stderr,
            )
            new_id = unique_id

        seen.add(new_id)
        output.write(f">{new_id}\n")
        written += 1

print(f"Stylophora CDS records written: {written}")
print(f"Headers without gene field:     {missing_gene}")
print(f"Duplicate IDs renamed:          {duplicate_ids}")
PY

echo "Cleaned Stylophora CDS: $CLEAN_STYLO"

SUMMARY="${OUTDIR}/proteinortho_nucleotide_summary.tsv"

printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "dataset" \
    "input_transcripts" \
    "input_transcripts_with_match" \
    "input_transcripts_without_match" \
    "unique_stylophora_cds" \
    "unique_stylophora_genes" \
    > "$SUMMARY"

for QUERY in "${QUERY_FILES[@]}"; do

    if [[ ! -f "$QUERY" ]]; then
        echo "Warning: query file not found, skipping: $QUERY" >&2
        continue
    fi

    DATASET="${QUERY%.transcripts.fasta}"
    RUN_DIR="${OUTDIR}/${DATASET}"
    PROJECT="${DATASET}_vs_Stylophora"

    mkdir -p "$RUN_DIR"

    echo
    echo "============================================================"
    echo "Dataset: $DATASET"
    echo "Query:   $QUERY"
    echo "============================================================"

    # Use absolute paths so the command works after changing directory.
    QUERY_ABS="$(readlink -f "$QUERY")"
    STYLO_ABS="$(readlink -f "$CLEAN_STYLO")"

    (
        cd "$RUN_DIR"

        proteinortho \
            -project="$PROJECT" \
            -p=blastn+ \
            -e="$EVALUE" \
            -cpus="$THREADS" \
            "$QUERY_ABS" \
            "$STYLO_ABS"
    )

    PROTEINORTHO_TSV="${RUN_DIR}/${PROJECT}.proteinortho.tsv"

    if [[ ! -s "$PROTEINORTHO_TSV" ]]; then
        echo "Error: expected Proteinortho output not found:" >&2
        echo "  $PROTEINORTHO_TSV" >&2
        exit 1
    fi

    STYLO_CDS_IDS="${RUN_DIR}/${DATASET}.stylophora_cds.txt"
    STYLO_GENES="${RUN_DIR}/${DATASET}.stylophora_genes.txt"
    MATCHED_QUERY="${RUN_DIR}/${DATASET}.matched_input_transcripts.txt"
    UNMATCHED_QUERY="${RUN_DIR}/${DATASET}.unmatched_input_transcripts.txt"
    STATS="${RUN_DIR}/${DATASET}.statistics.tsv"

    python3 - \
        "$PROTEINORTHO_TSV" \
        "$QUERY" \
        "$CLEAN_STYLO" \
        "$STYLO_CDS_IDS" \
        "$STYLO_GENES" \
        "$MATCHED_QUERY" \
        "$UNMATCHED_QUERY" \
        "$STATS" \
        "$SUMMARY" \
        "$DATASET" <<'PY'
import csv
import re
import sys
from pathlib import Path

(
    proteinortho_tsv,
    query_fasta,
    stylophora_fasta,
    cds_output,
    genes_output,
    matched_output,
    unmatched_output,
    statistics_output,
    summary_output,
    dataset,
) = sys.argv[1:]


def read_fasta_ids(path):
    identifiers = set()

    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                identifier = line[1:].strip().split()[0]

                if identifier:
                    identifiers.add(identifier)

    return identifiers


def split_cell(cell):
    cell = cell.strip()

    if not cell or cell in {"*", "-", "NA"}:
        return []

    return [
        value.strip()
        for value in cell.split(",")
        if value.strip() and value.strip() not in {"*", "-", "NA"}
    ]


def find_column(header, fasta_path):
    basename = Path(fasta_path).name

    for index, column in enumerate(header):
        clean = column.strip().lstrip("#")

        if Path(clean).name == basename:
            return index

    for index, column in enumerate(header):
        if basename in column:
            return index

    raise RuntimeError(
        f"Could not identify column for {basename}.\n"
        f"Proteinortho header: {header}"
    )


def write_sorted(path, values):
    with open(path, "w", encoding="utf-8") as output:
        for value in sorted(values):
            output.write(value + "\n")


all_query_ids = read_fasta_ids(query_fasta)

matched_query_ids = set()
stylophora_cds_ids = set()

with open(
    proteinortho_tsv,
    encoding="utf-8",
    errors="replace",
    newline="",
) as handle:

    reader = csv.reader(handle, delimiter="\t")

    try:
        header = next(reader)
    except StopIteration:
        raise RuntimeError(
            f"Proteinortho output is empty: {proteinortho_tsv}"
        )

    query_column = find_column(header, query_fasta)
    stylophora_column = find_column(header, stylophora_fasta)

    for row in reader:

        if len(row) <= max(query_column, stylophora_column):
            continue

        query_ids = split_cell(row[query_column])
        stylophora_ids = split_cell(row[stylophora_column])

        # Only count orthogroups represented in both datasets.
        if query_ids and stylophora_ids:
            matched_query_ids.update(query_ids)
            stylophora_cds_ids.update(stylophora_ids)


matched_query_ids &= all_query_ids
unmatched_query_ids = all_query_ids - matched_query_ids

stylophora_genes = set()
unparseable_ids = set()

for cds_id in stylophora_cds_ids:

    # Expected suffix:
    # ..._LOC111326177
    #
    # Also tolerates duplicate suffixes:
    # ..._LOC111326177_copy2
    match = re.search(r"_(LOC\d+)(?:_copy\d+)?$", cds_id)

    if match:
        stylophora_genes.add(match.group(1))
    else:
        unparseable_ids.add(cds_id)


write_sorted(cds_output, stylophora_cds_ids)
write_sorted(genes_output, stylophora_genes)
write_sorted(matched_output, matched_query_ids)
write_sorted(unmatched_output, unmatched_query_ids)

if unparseable_ids:
    unparseable_path = genes_output.replace(
        ".stylophora_genes.txt",
        ".stylophora_cds_without_LOC.txt",
    )
    write_sorted(unparseable_path, unparseable_ids)

    print(
        f"Warning: {len(unparseable_ids)} Stylophora CDS IDs did not "
        f"contain a terminal LOC identifier.",
        file=sys.stderr,
    )

statistics = {
    "dataset": dataset,
    "input_transcripts": len(all_query_ids),
    "input_transcripts_with_match": len(matched_query_ids),
    "input_transcripts_without_match": len(unmatched_query_ids),
    "unique_stylophora_cds": len(stylophora_cds_ids),
    "unique_stylophora_genes": len(stylophora_genes),
}

with open(statistics_output, "w", encoding="utf-8") as output:
    output.write("metric\tvalue\n")

    for metric, value in statistics.items():
        output.write(f"{metric}\t{value}\n")

with open(summary_output, "a", encoding="utf-8") as output:
    output.write(
        f"{dataset}\t"
        f"{statistics['input_transcripts']}\t"
        f"{statistics['input_transcripts_with_match']}\t"
        f"{statistics['input_transcripts_without_match']}\t"
        f"{statistics['unique_stylophora_cds']}\t"
        f"{statistics['unique_stylophora_genes']}\n"
    )

print(f"Input transcripts:                 {len(all_query_ids)}")
print(f"Input transcripts with matches:    {len(matched_query_ids)}")
print(f"Input transcripts without matches: {len(unmatched_query_ids)}")
print(f"Unique Stylophora CDS matches:      {len(stylophora_cds_ids)}")
print(f"Unique Stylophora genes:            {len(stylophora_genes)}")
PY

done

echo
echo "All analyses finished."
echo "Summary: $SUMMARY"
