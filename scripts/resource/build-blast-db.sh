#!/bin/bash

#SBATCH --job-name BLAST
#SBATCH -A naiss2025-22-494
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 6
#SBATCH --mem=24GB
#SBATCH -t 02:30:00
#SBATCH --output=slurm-logs/blast/SLURM-%j.out
#SBATCH --error=slurm-logs/blast/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

#ml load PDCOLD/23.12
#ml load blast+/2.15.0-cpeGNU-23.12

module load PDC/24.11
module spider blast+/2.16.0-cpeGNU-24.11

# Start the process
echo "$(date) [INFO]        Starting script execution"

# Define variables
OUTDIR=08-BLAST
SILVA=../databases/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
DB_DIR=$OUTDIR/db
TMP_DIR=$OUTDIR/tmp-dir

mkdir -p $OUTDIR $TMP_DIR $DB_DIR

echo "$(date) [INFO]        Create unique identifiers of all bins"
# Add unique identifier to all the bins.
for DIR in 06-metaWRAP-refined-bins/*; do 
    
    ID=$(basename $DIR)
    BIN_DIR=06-metaWRAP-refined-bins/$ID/bin_refinement/metawrap*_bins


    for BIN in $BIN_DIR/*.fa; do
        BIN_ID="${ID}_$(basename $BIN .fa)"

        # Add the $ID to each header
        sed "s/^>/>${BIN_ID}_/" $BIN > $TMP_DIR/$BIN_ID.fa
    done

done

echo "$(date) [INFO]        Concatenate all bins and SILVA database"
cat $TMP_DIR/*.fa > $OUTDIR/custom-db.fa
zcat $SILVA >> $OUTDIR/custom-db.fa

echo "$(date) [INFO]        Make BLAST database"
makeblastdb \
    -in $OUTDIR/custom-db.fa \
    -dbtype nucl \
    -parse_seqids \
    -out $DB_DIR/custom-db

echo "$(date) [COMPLETE]        Done!"