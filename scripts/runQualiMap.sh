#!/bin/bash

#SBATCH --job-name qualimap
#SBATCH -A naiss2025-22-494
#SBATCH --array=1-10
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=16GB
#SBATCH -t 05:30:00
#SBATCH --output=slurm-logs/qualimap/SLURM-%j-%a.out
#SBATCH --error=slurm-logs/qualimap/SLURM-%j-%a.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL


# Start the process
echo "$(date) [INFO]        Starting script execution"

# Load in package
source /cfs/klemming/home/a/andbou/.bashrc
mamba activate qualimap

# SLURM array jobid
JOBID=${SLURM_ARRAY_TASK_ID}
CPU=8
MEM=14G

# Paths and variables
BAM_TXT=doc/qualimap_bam_txt_paths-metagenome.txt
OUTDIR=11-QUALIMAP
BED=07-bin-quality-assessment/CHST-pt_042/bakta/bakta_output/bin.3/bin.3.bed

mkdir -p $OUTDIR

# Read the fastq files this task will work with.
BAM_LIST_IN=$(sed -n "${JOBID}p" $BAM_TXT)
PREFIX=$(basename $BAM_LIST_IN -bams.txt)

echo "$(date) [INFO]        Running qualimap for: $PREFIX"
qualimap --java-mem-size=$MEM multi-bamqc \
    --data $BAM_LIST_IN \
    --outdir $OUTDIR/$PREFIX \
    --feature-file $BED \
    --outformat PDF:HTML \
    --run-bamqc

echo "$(date) [DONE]        Finished qualimap for: $PREFIX"