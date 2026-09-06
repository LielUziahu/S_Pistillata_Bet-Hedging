#!/bin/bash
#################################################################################################################
#SBATCH --job-name=sbatchTemplate ## Name of your job
#SBATCH --ntasks=64 ## number of cpu's to allocate for a job
#SBATCH --ntasks-per-node=64 ## number of cpu's to allocate per each node
#SBATCH --nodes=1 ## number of nodes to allocate for a job
#SBATCH --mem=128G ## memory to allocate for your job in MB
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
source /lustre1/home/mass/eskalon/miniconda/bin/activate ribodetector
#########################################

INPUT_DIR="./trim"
OUTPUT_DIR="./ribo"

mkdir -p "$OUTPUT_DIR"

for infile in ${INPUT_DIR}/*.gz; do
    # Extract file name without path
    fname=$(basename "$infile")    

    # Remove the final ".gz"
    base=${fname%.trim.fastq.gz}               

    # Add suffix and extension 
    outfile="${OUTPUT_DIR}/${base}.nonrrna.fastq.gz"
    # skip if output already exists
    if [[ -f "$outfile" ]]; then
        echo "Skipping (already done): $outfile"
        continue
    fi

    echo "Processing: $infile"
    echo "Output:    $outfile"

    ribodetector_cpu \
      -t 64 \
      -l 150 \
      -i "$infile" \
      --chunk_size 256 \
      -o "$outfile"

done
