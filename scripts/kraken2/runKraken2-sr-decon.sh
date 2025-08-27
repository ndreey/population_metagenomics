#!/bin/bash

#SBATCH --job-name kraken
#SBATCH -A naiss2025-22-494
#SBATCH --array=1-84
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --mem=64GB
#SBATCH -t 02:00:00
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

# Define variables
#K2_DB=/cfs/klemming/projects/snic/snic2020-6-222/Projects/Tconura/working/Andre/databases/k2_thistle_db
K2_DB=/sw/data/Kraken2_data/prebuilt/k2_pluspfp_20250714/
READS=doc/sample-merged-sr-paths-no-SC.txt
OUTDIR=09-KRAKEN2/decon-0_2-pluspfp-full
K2_OUT=$OUTDIR/kraken2_output
K2_MPA=$OUTDIR/kreport2mpa

mkdir -p $K2_OUT $K2_MPA

sed -n "${JOBID}p" "$READS" | while read R1 R2; do

    ID=$(basename $R1 _R1.fq.gz)
    echo "$(date) [INFO]        Classifying $ID"

    k2 classify \
        --db $K2_DB \
        --threads $CPU \
        --confidence 0.2 \
        --output $K2_OUT/$ID.k2.output \
        --report $K2_OUT/$ID.k2.report \
        --paired \
        --memory-mapping \
        --use-names \
        $R1 $R2
done

echo "$(date) [COMPLETE]        Done!"