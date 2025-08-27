#!/bin/bash

#SBATCH --job-name DRAM_main
#SBATCH -A naiss2025-22-494
#SBATCH -p main
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mem=256GB
#SBATCH -t 01:30:00
#SBATCH --output=slurm-logs/dram/SLURM-%j.out
#SBATCH --error=slurm-logs/dram/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Start the process
echo "$(date) [INFO]        Starting script execution"

# Load in package
source /cfs/klemming/home/a/andbou/.bashrc
mamba activate DRAM


MAG=MAG/ctg4994.fa
DB_PATH=/cfs/klemming/projects/snic/snic2020-6-222/Projects/Tconura/working/Andre/databases/DRAM_v1.5.0
OUTDIR=08-DRAM
CPU=8

DRAM.py annotate \
    --input_fasta $MAG \
    --output_dir $OUTDIR \
    --threads $CPU

DRAM.py distill -i $OUTDIR/annotations.tsv -o $OUTDIR/genome_summary --trna_path $OUTDIR/trnas.tsv --rrna_path $OUTDIR/rrnas.tsv

DRAM.py distill -i $OUTDIR/annotations.tsv -o $OUTDIR/genome_summary2 --trna_path $OUTDIR/trnas.tsv --rrna_path $OUTDIR/rrnas.tsv --distillate_gene_names
#DRAM-setup.py prepare_databases --output_dir $DB_PATH --threads 10 --clear_config