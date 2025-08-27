#!/bin/bash

#SBATCH --job-name busco
#SBATCH -A naiss2025-22-494
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=24GB
#SBATCH -t 15:30:00
#SBATCH --output=slurm-logs/busco/SLURM-%j.out
#SBATCH --error=slurm-logs/busco/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Start the process
echo "$(date) [INFO]        Starting script execution"



# Folder holding mags.