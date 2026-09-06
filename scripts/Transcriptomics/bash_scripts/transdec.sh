#!/bin/bash
#################################################################################################################
#SBATCH --job-name=sbatchTemplate ## Name of your job
#SBATCH --ntasks=8 ## number of cpu's to allocate for a job
#SBATCH --ntasks-per-node=8 ## number of cpu's to allocate per each node
#SBATCH --nodes=1 ## number of nodes to allocate for a job
#SBATCH --mem=36G ## memory to allocate for your job in MB
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

set -euo pipefail

THREADS=8
MIN_PROTEIN_LENGTH=50

FILES=(
#    "Meyer2011.down_after_CCA.transcripts.fasta"
#    "Meyer2011.up_after_CCA.transcripts.fasta"
    "Walker2019.probing.transcripts.clean.fasta"
)

command -v TransDecoder.LongOrfs >/dev/null 2>&1 || {
    echo "Error: TransDecoder.LongOrfs is not available in PATH." >&2
    exit 1
}

command -v TransDecoder.Predict >/dev/null 2>&1 || {
    echo "Error: TransDecoder.Predict is not available in PATH." >&2
    exit 1
}

mkdir -p transdecoder_results

for FASTA in "${FILES[@]}"; do
    if [[ ! -f "$FASTA" ]]; then
        echo "Warning: file not found, skipping: $FASTA" >&2
        continue
    fi

    NAME="${FASTA%.fasta}"
    OUTDIR="transdecoder_results/${NAME}.transdecoder"

    echo "Processing: $FASTA"

    mkdir -p "$OUTDIR"

    TransDecoder.LongOrfs \
        -t "$FASTA" \
        -m "$MIN_PROTEIN_LENGTH" \
        -O "$OUTDIR"

    TransDecoder.Predict \
        -t "$FASTA" \
        --single_best_only \
        -O "$OUTDIR"

    PEP="${OUTDIR}/$(basename "$FASTA").transdecoder.pep"
    CDS="${OUTDIR}/$(basename "$FASTA").transdecoder.cds"

    if [[ -s "$PEP" ]]; then
        cp "$PEP" "transdecoder_results/${NAME}.predicted_proteins.faa"
    else
        echo "Warning: expected protein output not found: $PEP" >&2
    fi

    if [[ -s "$CDS" ]]; then
        cp "$CDS" "transdecoder_results/${NAME}.predicted_CDS.fna"
    fi

    echo "Finished: $FASTA"
done

echo
echo "Results are in: transdecoder_results/"
