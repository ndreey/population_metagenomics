#!/bin/bash

#SBATCH --job-name ncbi-downl
#SBATCH -A naiss2024-22-580
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 6
#SBATCH --mem=32GB
#SBATCH -t 06:00:00
#SBATCH --output=slurm-logs/thesis-pipe/kraken2/SLURM-%j-ncbidown.out
#SBATCH --error=slurm-logs/thesis-pipe/kraken2/SLURM-%j-ncbidown.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL


# Load in the metaMDBG mamba environment
source /cfs/klemming/home/a/andbou/.bashrc
mamba activate ncbidown

echo "$(date) [START]     Loading module and starting script execution"

# Slurm variable
CPU=$SLURM_CPUS_ON_NODE

echo "$(date)       Starting download"
# Download datasets
datasets download genome accession \
    --inputfile doc/kraken-custom/bacteria-archaea-fungi-cirsium-human-insecta-2025-04-13.acc \
    --filename data/dehydrated-ncbi-set.zip --dehydrated

mkdir -p data/custom-kraken2-genomes

echo "$(date)       Unzipping dehydrated genome set"
unzip data/dehydrated-ncbi-set.zip -d data/custom-kraken2-genomes

echo "$(date)       Rehydrating it with full genomes"
datasets rehydrate --directory data/custom-kraken2-genomes

echo "$(date)       Done!"
