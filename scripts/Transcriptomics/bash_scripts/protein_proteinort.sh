set -euo pipefail
shopt -s nullglob

STYLO_FASTA="${1:-GCF_002571385.2_Stylophora_pistillata_v1.1_protein.short.faa}"
MAPPING_FILE="${2:-xp_to_loc.tsv}"
THREADS="${3:-16}"

# Relaxed ProteinOrtho/DIAMOND settings for short or divergent proteins.
EVALUE="1e-5"
IDENTITY="25"
COVERAGE="20"
SIMILARITY="0.8"
DIAMOND_OPTIONS="--more-sensitive"

# Resolve all paths before changing directories inside individual runs.
WORKDIR="$(pwd -P)"
STYLO_FASTA="$(readlink -f "$STYLO_FASTA")"
MAPPING_FILE="$(readlink -f "$MAPPING_FILE")"
OUTDIR="${WORKDIR}/proteinortho_stylophora_results_relaxed"

# Candidate files to process.
CANDIDATE_FILES=( *.proteins.fasta )

if (( ${#CANDIDATE_FILES[@]} == 0 )); then
    echo "Error: no *.proteins.fasta files found in the current directory." >&2
    exit 1
fi

[[ -f "$STYLO_FASTA" ]] || {
    echo "Error: Stylophora protein FASTA not found: $STYLO_FASTA" >&2
    exit 1
}

[[ -f "$MAPPING_FILE" ]] || {
    echo "Error: protein-to-gene mapping file not found: $MAPPING_FILE" >&2
    exit 1
}

command -v proteinortho >/dev/null 2>&1 || {
    echo "Error: proteinortho6 is not installed or is not in PATH." >&2
    exit 1
}

command -v python3 >/dev/null 2>&1 || {
    echo "Error: python3 is not installed or is not in PATH." >&2
    exit 1
}

mkdir -p "$OUTDIR"

SUMMARY="${OUTDIR}/proteinortho_summary.tsv"

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "dataset" \
    "input_proteins" \
    "input_proteins_with_stylophora_match" \
    "input_proteins_without_match" \
    "unique_stylophora_proteins" \
    "unique_stylophora_genes" \
    "unmapped_stylophora_proteins" \
    > "$SUMMARY"

for QUERY_FASTA in "${CANDIDATE_FILES[@]}"; do

    # Do not accidentally process the Stylophora proteome itself if its
    # filename happens to end in .proteins.fasta.
    if [[ "$(readlink -f "$QUERY_FASTA")" == "$(readlink -f "$STYLO_FASTA")" ]]; then
        continue
    fi

    DATASET="${QUERY_FASTA%.proteins.fasta}"

    # Replace characters that could cause inconvenient output paths.
    SAFE_NAME=$(printf '%s' "$DATASET" | sed 's/[^A-Za-z0-9._-]/_/g')

    RUN_DIR="${OUTDIR}/${SAFE_NAME}"
    PROJECT="${SAFE_NAME}_vs_Stylophora"

    mkdir -p "$RUN_DIR"

    # DIAMOND/ProteinOrtho may create databases and cache files beside the
    # FASTA paths supplied on the command line. Use run-local symlinks so all
    # generated files stay inside the designated dataset directory.
    QUERY_ABS="$(readlink -f "$QUERY_FASTA")"
    QUERY_LOCAL="$(basename "$QUERY_FASTA")"
    STYLO_LOCAL="$(basename "$STYLO_FASTA")"

    ln -sfn "$QUERY_ABS" "$RUN_DIR/$QUERY_LOCAL"
    ln -sfn "$STYLO_FASTA" "$RUN_DIR/$STYLO_LOCAL"
    mkdir -p "$RUN_DIR/tmp"

    echo
    echo "============================================================"
    echo "Dataset: $DATASET"
    echo "Query:   $QUERY_FASTA"
    echo "============================================================"

    # Run from the dataset directory using local FASTA symlinks. This keeps
    # ProteinOrtho output, DIAMOND databases, caches, and temporary files here.
    (
        cd "$RUN_DIR"
        export TMPDIR="$PWD/tmp"

        proteinortho \
            -project="$PROJECT" \
            -p=diamond \
            -e="$EVALUE" \
            -identity="$IDENTITY" \
            -cov="$COVERAGE" \
            -sim="$SIMILARITY" \
            -subparaBlast="$DIAMOND_OPTIONS" \
            -cpus="$THREADS" \
            -singles \
            -keep \
            -force \
            "$QUERY_LOCAL" \
            "$STYLO_LOCAL"
    )

    PROTEINORTHO_TSV="${RUN_DIR}/${PROJECT}.proteinortho.tsv"

    if [[ ! -s "$PROTEINORTHO_TSV" ]]; then
        echo "Error: expected output was not produced: $PROTEINORTHO_TSV" >&2
        exit 1
    fi

    STYLO_PROTEINS="${RUN_DIR}/${SAFE_NAME}.stylophora_proteins.txt"
    STYLO_GENES="${RUN_DIR}/${SAFE_NAME}.stylophora_genes.txt"
    UNMAPPED="${RUN_DIR}/${SAFE_NAME}.stylophora_proteins_unmapped.txt"
    MATCHED_QUERY_IDS="${RUN_DIR}/${SAFE_NAME}.matched_input_proteins.txt"
    UNMATCHED_QUERY_IDS="${RUN_DIR}/${SAFE_NAME}.unmatched_input_proteins.txt"
    STATS="${RUN_DIR}/${SAFE_NAME}.statistics.tsv"

    python3 - \
        "$PROTEINORTHO_TSV" \
        "$QUERY_ABS" \
        "$STYLO_FASTA" \
        "$MAPPING_FILE" \
        "$STYLO_PROTEINS" \
        "$STYLO_GENES" \
        "$UNMAPPED" \
        "$MATCHED_QUERY_IDS" \
        "$UNMATCHED_QUERY_IDS" \
        "$STATS" \
        "$SUMMARY" \
        "$DATASET" <<'PY'
import csv
import os
import re
import sys
from pathlib import Path

(
    proteinortho_tsv,
    query_fasta,
    stylo_fasta,
    mapping_file,
    stylo_proteins_out,
    stylo_genes_out,
    unmapped_out,
    matched_query_out,
    unmatched_query_out,
    stats_out,
    summary_out,
    dataset,
) = sys.argv[1:]


def fasta_ids(path):
    """Return unique first-token FASTA identifiers."""
    identifiers = set()

    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                identifier = line[1:].strip().split()[0]
                if identifier:
                    identifiers.add(identifier)

    return identifiers


def split_proteinortho_cell(value):
    """Split comma-separated Proteinortho IDs and discard missing values."""
    value = value.strip()

    if not value or value in {"*", "-", "NA"}:
        return []

    return [
        item.strip()
        for item in value.split(",")
        if item.strip() and item.strip() not in {"*", "-", "NA"}
    ]


def basename_variants(path):
    """Generate filename forms useful for matching Proteinortho headers."""
    name = Path(path).name
    resolved_name = Path(os.path.realpath(path)).name
    return {name, resolved_name, path, os.path.realpath(path)}


def identify_species_column(header, fasta_path):
    """
    Find a Proteinortho species column by matching its header against the
    FASTA filename. Proteinortho normally uses the input filename.
    """
    variants = basename_variants(fasta_path)

    for index, column in enumerate(header):
        clean_column = column.strip().lstrip("#")
        column_name = Path(clean_column).name

        if clean_column in variants or column_name in variants:
            return index

    fasta_basename = Path(fasta_path).name

    for index, column in enumerate(header):
        if fasta_basename in column:
            return index

    raise RuntimeError(
        f"Could not identify the column for {fasta_path}.\n"
        f"Proteinortho columns were: {header}"
    )


def clean_field(value):
    return value.strip().strip('"').strip("'")


def load_mapping(path):
    """
    Read protein-to-gene mapping.

    Expected first two columns:
        protein_accession<TAB>gene_id

    Comma-separated input is also accepted. Header lines are skipped
    automatically when they do not resemble protein accessions.
    """
    exact = {}
    versionless = {}

    with open(path, encoding="utf-8-sig", errors="replace", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)

        # Prefer tab, but permit comma-separated mappings.
        delimiter = "\t" if "\t" in sample else ","
        reader = csv.reader(handle, delimiter=delimiter)

        for row in reader:
            if len(row) < 2:
                continue

            protein = clean_field(row[0])
            gene = clean_field(row[1])

            if not protein or not gene:
                continue

            lower_protein = protein.lower()
            lower_gene = gene.lower()

            if (
                lower_protein in {
                    "protein",
                    "protein_id",
                    "protein id",
                    "accession",
                    "accession.version",
                }
                or lower_gene in {
                    "gene",
                    "gene_id",
                    "gene id",
                    "locus",
                    "locus_tag",
                }
            ):
                continue

            exact[protein] = gene
            versionless.setdefault(re.sub(r"\.\d+$", "", protein), gene)

    return exact, versionless


all_query_ids = fasta_ids(query_fasta)
all_stylo_reference_ids = fasta_ids(stylo_fasta)

matched_query_ids = set()
stylophora_protein_ids = set()

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
        raise RuntimeError(f"Empty Proteinortho table: {proteinortho_tsv}")

    query_column = identify_species_column(header, query_fasta)
    stylo_column = identify_species_column(header, stylo_fasta)

    for row in reader:
        if not row:
            continue

        required_columns = max(query_column, stylo_column)

        if len(row) <= required_columns:
            continue

        query_ids = split_proteinortho_cell(row[query_column])
        stylo_ids = split_proteinortho_cell(row[stylo_column])

        # Count a source protein as matched only when the orthogroup
        # contains at least one Stylophora protein.
        if query_ids and stylo_ids:
            matched_query_ids.update(query_ids)
            stylophora_protein_ids.update(stylo_ids)


# Restrict query statistics to IDs genuinely present in the input FASTA.
matched_query_ids &= all_query_ids
unmatched_query_ids = all_query_ids - matched_query_ids

# This check catches header/parsing problems.
unexpected_stylo_ids = stylophora_protein_ids - all_stylo_reference_ids

if unexpected_stylo_ids:
    print(
        "Warning: some extracted Stylophora IDs were not found in the "
        "Stylophora FASTA:",
        file=sys.stderr,
    )
    for identifier in sorted(unexpected_stylo_ids):
        print(f"  {identifier}", file=sys.stderr)

mapping_exact, mapping_versionless = load_mapping(mapping_file)

stylophora_gene_ids = set()
unmapped_stylophora_ids = set()

for protein in stylophora_protein_ids:
    gene = mapping_exact.get(protein)

    # Permit a versionless fallback:
    # XP_123456.1 can match XP_123456 when necessary.
    if gene is None:
        protein_without_version = re.sub(r"\.\d+$", "", protein)
        gene = mapping_versionless.get(protein_without_version)

    if gene:
        stylophora_gene_ids.add(gene)
    else:
        unmapped_stylophora_ids.add(protein)


def write_list(path, values):
    with open(path, "w", encoding="utf-8") as output:
        for value in sorted(values):
            output.write(value + "\n")


write_list(stylo_proteins_out, stylophora_protein_ids)
write_list(stylo_genes_out, stylophora_gene_ids)
write_list(unmapped_out, unmapped_stylophora_ids)
write_list(matched_query_out, matched_query_ids)
write_list(unmatched_query_out, unmatched_query_ids)

statistics = {
    "dataset": dataset,
    "input_proteins": len(all_query_ids),
    "input_proteins_with_stylophora_match": len(matched_query_ids),
    "input_proteins_without_match": len(unmatched_query_ids),
    "unique_stylophora_proteins": len(stylophora_protein_ids),
    "unique_stylophora_genes": len(stylophora_gene_ids),
    "unmapped_stylophora_proteins": len(unmapped_stylophora_ids),
}

with open(stats_out, "w", encoding="utf-8") as output:
    output.write("metric\tvalue\n")

    for metric, value in statistics.items():
        output.write(f"{metric}\t{value}\n")

with open(summary_out, "a", encoding="utf-8") as output:
    output.write(
        f"{dataset}\t"
        f"{statistics['input_proteins']}\t"
        f"{statistics['input_proteins_with_stylophora_match']}\t"
        f"{statistics['input_proteins_without_match']}\t"
        f"{statistics['unique_stylophora_proteins']}\t"
        f"{statistics['unique_stylophora_genes']}\t"
        f"{statistics['unmapped_stylophora_proteins']}\n"
    )

print(f"Input proteins:                   {len(all_query_ids)}")
print(f"Input proteins with match:        {len(matched_query_ids)}")
print(f"Unique Stylophora proteins:       {len(stylophora_protein_ids)}")
print(f"Unique mapped Stylophora genes:   {len(stylophora_gene_ids)}")
print(f"Unmapped Stylophora proteins:     {len(unmapped_stylophora_ids)}")
PY

done

echo
echo "All analyses finished. Results: $OUTDIR"
