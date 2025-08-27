#!/bin/bash

#SBATCH --job-name bakta
#SBATCH -A naiss2025-22-494
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=24GB
#SBATCH -t 15:30:00
#SBATCH --output=slurm-logs/bakta/SLURM-%j.out
#SBATCH --error=slurm-logs/bakta/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL



# Start the process
echo "$(date) [INFO]        Starting script execution"

# Load in package
source /cfs/klemming/home/a/andbou/.bashrc
mamba activate bakta

# SLURM array jobid
CPU=8

BAKTA_DB=../databases/bakta-db_v6/db
PREFIX=CHSt1
OUTDIR=04-bakta
GENOME=../MAG-MAYHEM/data/chst1-mag/chst1-mag.fa

echo "$(date) [INFO]        PREFIX: $PREFIX"
echo "$(date) [INFO]        GENOME: $GENOME"

#bakta_db download --output $BAKTA_DB --type full

bakta \
    --db $BAKTA_DB \
    --prefix $PREFIX \
    --output $OUTDIR/$PREFIX \
    --genus Stammerula \
    --species tephritidis \
    --strain CHSt1 \
    --gram - \
    --keep-contig-headers \
    --threads $CPU \
    --force \
    $GENOME


echo "$(date) [FINISH]        Done!"