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
source /lustre1/home/mass/eskalon/miniconda/bin/activate rnaseq
#########################################

### Below you can enter your program job command ###

fastqc --noextract -o ./fastqc_before -t 64 ../raw/data/*R2*.fastq.gz
fastqc --noextract -o ./ribo -t 64 ./ribo/*

source /lustre1/home/mass/eskalon/miniconda/bin/activate multiqc
multiqc --outdir fastqc_before fastqc_before/ --filename fastqc.before.html
fastqc --noextract -o ./fastqc_after -t 4 trim/*.fastq.gz
multiqc --outdir ribo ribo/ --filename ribo.html
