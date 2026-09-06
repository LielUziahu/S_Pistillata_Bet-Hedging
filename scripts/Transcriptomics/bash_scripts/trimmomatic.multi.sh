#!/bin/bash
#################################################################################################################
#SBATCH --job-name=sbatchTemplate ## Name of your job
#SBATCH --ntasks=32 ## number of cpu's to allocate for a job
#SBATCH --ntasks-per-node=32 ## number of cpu's to allocate per each node
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
#########################################
source /lustre1/home/mass/eskalon/miniconda/bin/activate rnaseq

### Below you can enter your program job command ###
for R1 in  ../raw/data/*R2_001.fastq.gz; do
    fname=$(basename "$R1")

    # Replace ".fastq.gz" with ".trim.fastq.gz" and put into ./trim
    out="trim/${fname%.fastq.gz}.trim.fastq.gz"

    echo "Trimming $R1 -> $out"

    trimmomatic SE -threads 32 "$R1" "$out" \
    ILLUMINACLIP:Sequencing_adaptors.fasta:2:30:10 \
    SLIDINGWINDOW:4:5 MAXINFO:50:0.6 MINLEN:30


done

fastqc --noextract -o fastqc_after  -t 32 trim/*.trim.fastq.gz

source /lustre1/home/mass/eskalon/miniconda/bin/activate multiqc
multiqc fastqc_after -o fastqc_after --filename fastqc_after
