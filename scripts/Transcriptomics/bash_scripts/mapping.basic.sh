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
#########################################
source /lustre1/home/mass/eskalon/miniconda/bin/activate rnaseq
### Below you can enter your program job command ###

ulimit -n 5000
mkdir -p bam

for R1 in /lustre1/home/mass/eskalon/hiba/pl_align/raw/*nonrrna.fastq.gz; do
R2="${R1%.nonrrna.fastq.gz}"


echo "$R1"
echo "$R2"
STAR --genomeDir index_SP2.gtf \
 --readFilesIn "$R1" \
 --readFilesCommand gunzip -c \
 --runThreadN 16 --outSAMtype BAM SortedByCoordinate \
 --outFileNamePrefix "$R2" \
 --clip3pAdapterSeq AAAAAAAAAA \
 --outReadsUnmapped Fastx \
 --quantMode GeneCounts \
 --limitBAMsortRAM 64424509440  --outBAMsortingBinsN 50

done

mv ./raw/*bam ./bam/
mv ./raw/*out ./bam/
mv ./raw/*tab ./bam/
mv ./raw/*mate1 ./bam

source /lustre1/home/mass/eskalon/miniconda/bin/activate samtools

for i in ./bam/*.bam; do
  samtools index $i
done

source /lustre1/home/mass/eskalon/miniconda/bin/activate multiqc

multiqc bam -o bam --filename star

