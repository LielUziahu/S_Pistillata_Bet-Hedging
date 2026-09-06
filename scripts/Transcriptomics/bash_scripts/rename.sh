#!/bin/bash
#################################################################################################################
#SBATCH --job-name=sbatchTemplate ## Name of your job
#SBATCH --ntasks=1 ## number of cpu's to allocate for a job
#SBATCH --ntasks-per-node=1 ## number of cpu's to allocate per each node
#SBATCH --nodes=1 ## number of nodes to allocate for a job
#SBATCH --mem=12G ## memory to allocate for your job in MB
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
#source /lustre1/home/mass/eskalon/miniconda/bin/activate rnaseq
#########################################

### Below you can enter your program job command ###

#!/usr/bin/env bash

set -euo pipefail

# Read mapping from ids.tsv into an array
declare -A MAP

while IFS=$'\t' read -r key value; do
    MAP["$key"]="$value"
done < ../ids.csv

# Process files
for f in *.nonrrna.fastq.gz; do
    # Extract prefix before the first underscore
    prefix=$(echo "$f" | cut -d"_" -f1)

    # Look up mapping
    newname="${MAP[$prefix]:-}"

    if [[ -z "$newname" ]]; then
        echo "No mapping found for $prefix, skipping..."
        continue
    fi

    # Construct final name
    final="${newname}.nonrrna.fastq.gz"

    echo "Renaming:  $f  →  $final"
    mv "$f" "$final"
done

