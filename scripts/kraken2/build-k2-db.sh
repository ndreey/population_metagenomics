#!/bin/bash

#SBATCH --job-name kraken
#SBATCH -A naiss2025-22-494
#SBATCH -p main
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --mem=100GB
#SBATCH -t 18:30:00
#SBATCH --output=slurm-logs/kraken2/SLURM-%j.out
#SBATCH --error=slurm-logs/kraken2/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Load in the kraken2 mamba environment
source /cfs/klemming/home/a/andbou/.bashrc
mamba activate kraken2

# Start the process
echo "$(date) [INFO]        Starting script execution"

# SLURM variables
CPU=24

# Define variables
K2_DB=/cfs/klemming/projects/snic/snic2020-6-222/Projects/Tconura/working/Andre/databases/k2_thistle_db
FASTA_TO_ADD_DIR=/cfs/klemming/projects/snic/snic2020-6-222/Projects/Tconura/working/Andre/databases/mag-host-k2-friendly-fasta

##### BUILD DATABASE ######
# Get taxonomy
echo "$(date) [INFO]        Download taxonomy"
k2 download-taxonomy --db $K2_DB

# Get bacteria library
echo "$(date) [INFO]        Download bacteria library"
k2 download-library --library bacteria --db $K2_DB --threads $CPU

# Add Stammerula and Tephritis
echo "$(date) [INFO]        Add custom fasta $(basename $FASTA_TO_ADD_DIR)"
k2 add-to-library \
    --files $FASTA_TO_ADD_DIR/ctg4994.fa $FASTA_TO_ADD_DIR/filtered-tconura.fasta \
    --threads $CPU \
    --db $K2_DB

echo "$(date) [INFO]        Taxonomy, bacteria and custom fasta added" # 4h
echo "$(date) [INFO]        Building database" # 10h

k2 build \
    --db $K2_DB \
    --threads $CPU

echo "$(date) [COMPLETE]        Done!"
