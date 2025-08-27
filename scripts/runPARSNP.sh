#!/bin/bash

#SBATCH --job-name parsnp
#SBATCH -A naiss2025-22-494
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --mem=36GB
#SBATCH -t 00:30:00
#SBATCH --output=slurm-logs/parsnp/SLURM-%j.out
#SBATCH --error=slurm-logs/parsnp/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Load in the kraken2 mamba environment
source /cfs/klemming/home/a/andbou/.bashrc
mamba activate parsnp

# Start the process
echo "$(date) [INFO]        Starting script execution"



parsnp \
    -r 17-parsnp/ctg4994.fa \
    -d 17-parsnp/mags/*.fa \
    --threads 12 \
    --alignment-program mafft \
    --output-dir 17-parsnp/run1


# Start the process
echo "$(date) [INFO]        Done!"