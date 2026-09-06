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
#./etc/profile.d/modules.sh # Required line for modules environment to work
#module load openmpi/1.8.4 python/2.7 # Load modules that are required by your program
#source /lustre1/home/mass/eskalon/miniconda/bin/activate ribodetector
#########################################

#!/bin/bash

echo -e "Sample\tTotal_Reads\tNonrRNA_Reads\trRNA_Reads\trRNA_Percent" > rrna_stats.tsv

for raw in *_R2_001.trim.fastq.gz; do
    sample=$(basename "$raw")
    clean="../ribo/${sample/.trim.fastq.gz/.nonrrna.fastq.gz}"

    # count reads (zcat + wc -l / 4)
    total=$(zcat "$raw" | wc -l)
    total=$((total / 4))

    nonrrna=$(zcat "$clean" | wc -l)
    nonrrna=$((nonrrna / 4))

    rrna=$((total - nonrrna))

    # avoid division by zero
    if [ "$total" -gt 0 ]; then
        percent=$(awk -v r="$rrna" -v t="$total" 'BEGIN { printf("%.2f", (r/t)*100) }')
    else
        percent="NA"
    fi

    echo -e "${sample}\t${total}\t${nonrrna}\t${rrna}\t${percent}" >> rrna_stats.tsv
done

