#!/bin/bash

#SBATCH --job-name kraken
#SBATCH -A naiss2025-22-494
#SBATCH --array=1-3
#SBATCH -p main
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --mem=125GB
#SBATCH -t 01:00:00
#SBATCH --output=slurm-logs/kraken2/SLURM-%j-%a.out
#SBATCH --error=slurm-logs/kraken2/SLURM-%j-%a.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Load in the kraken2 mamba environment
source /cfs/klemming/home/a/andbou/.bashrc
mamba activate kraken2

# Start the process
echo "$(date) [INFO]        Starting script execution"

# SLURM variables
JOBID=${SLURM_ARRAY_TASK_ID}
CPU=12

# Remove trailing slashes
OUTDIR="${OUTDIR%/}"


# Define variables
K2_DB=/cfs/klemming/projects/snic/snic2020-6-222/Projects/Tconura/working/Andre/databases/k2_thistle_db
READS=doc/lr-raw-paths.txt
OUTDIR=09-KRAKEN2/raw
K2_OUT=$OUTDIR/kraken2_output
K2_MPA=$OUTDIR/kreport2mpa

mkdir -p $K2_OUT $K2_MPA

sed -n "${JOBID}p" "$READS" | while read LR; do

    IDtmp=$(basename $LR .ccsreads.fastq.gz)
    ID=$(echo $IDtmp | cut -f 1,2,4 -d "_" )
    echo "$(date) [INFO]        Classifying $ID"

    k2 classify \
        --db $K2_DB \
        --threads $CPU \
        --confidence 0.5 \
        --output $K2_OUT/$ID.k2.output \
        --report $K2_OUT/$ID.k2.report \
        --use-names \
        $LR
    
    echo "$(date) [INFO]        Converting to mpa format"
    kreport2mpa.py \
        --report-file $K2_OUT/$ID.k2.report \
        --output $K2_MPA/$ID.mpa \
        --percentages

done

echo "$(date) [COMPLETE]        Done!"