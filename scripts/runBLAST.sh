#!/bin/bash

#SBATCH --job-name BLAST
#SBATCH -A naiss2025-22-494
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 6
#SBATCH --mem=24GB
#SBATCH -t 01:30:00
#SBATCH --output=slurm-logs/blast/SLURM-%j.out
#SBATCH --error=slurm-logs/blast/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

#ml load PDCOLD/23.12
#ml load blast+/2.15.0-cpeGNU-23.12

module load PDC/24.11
module load blast+/2.16.0-cpeGNU-24.11

# Start the process
echo "$(date) [INFO]        Starting script execution"

# Define variables
OUTDIR=08-BLAST
DB_DIR=$OUTDIR/db

#QUERY=doc/stammerula16s.fasta
QUERY=doc/Erwinia_Stammerula_Wolbachia-16SrRNA.fa
BLAST_OUT=$OUTDIR/blastn-result.tsv

# BLASTn search
blastn \
    -query ${QUERY} \
    -db $DB_DIR/custom-db \
    -outfmt "6 qseqid qlen sseqid slen nident pident length mismatch evalue score sstart send sseq" \
    -out ${BLAST_OUT} \
    -num_threads 6

echo "$(date) [COMPLETE]        Done!"