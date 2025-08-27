#!/bin/bash

# Start the process
echo "$(date) [INFO]        Starting script execution"

# Load bbmap
# mamba activate qc

WORK=/cfs/klemming/projects/snic/snic2020-6-222/Projects/Tconura/working/Andre/thesis_stammerula
STATS_OUT=$WORK/stats/metagenome-assembly-stats.tsv

# Create summary for number of bins
echo -e "id\tn_scaffolds\tn_contigs\tscaf_bp\tcontig_bp\tgap_pct\tscaf_N50\tscaf_L50\tctg_N50\tctg_L50\tscaf_N90\tscaf_L90\tctg_N90\tctg_L90\tscaf_max\tctg_max\tscaf_n_gt50K\tscaf_pct_gt50K\tgc_avg\tgc_std" > $STATS_OUT

for DIR in $WORK/05-metagenomes/02-megahit/*; do 
    
    ID=$(basename $DIR)
    FASTA=$DIR/*.fa
    
    echo "$(date) [INFO]        Calculating metagenome assembly statistics for $ID"
    STATS=$(stats.sh $FASTA format=6 | grep -v "^#")

    echo -e "$ID\t$STATS" >> $STATS_OUT
done

for METAG in $WORK/05-metagenomes/01-metamdbg/*.gz; do 
    
    ID=$(basename $METAG .contigs.fa.gz)
      
    echo "$(date) [INFO]        Calculating metagenome assembly statistics for $ID"
    STATS=$(stats.sh $METAG format=6 | grep -v "^#")

    echo -e "$ID\t$STATS" >> $STATS_OUT
done

echo "$(date) [COMPLETE]        Stats: $STATS_OUT"