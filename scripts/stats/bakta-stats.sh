#!/bin/bash

echo "$(date) [INFO] Starting BUSCO summary"

# Define variables
WORK=/cfs/klemming/projects/snic/snic2020-6-222/Projects/Tconura/working/Andre/thesis_stammerula
STATS_OUT=$WORK/stats/bakta-summary.tsv

# Header
echo -e "sample\tbin\ttRNAs\trRNAs\t5S.rRNA\t16S.rRNA\t23S.rRNA\tCDSs\tpseudogenes\thypotheticals\toriCs\toriVs\toriTs\tCDS.dens" > $STATS_OUT

# Parse the bakta output files
for DIR in $WORK/07-bin-quality-assessment/*; do
    SAMPLE=$(basename $DIR)
    echo "$(date) [INFO] Summarizing BAKTA for $SAMPLE"

    for BAKTA_DIR in $DIR/bakta/bakta_output/bin.*; do
        BIN=$(basename $BAKTA_DIR)
        SUM_FILE=$BAKTA_DIR/$BIN.txt
        GFF=$BAKTA_DIR/$BIN.gff3

        # Parse the values
        tRNA=$(cat $SUM_FILE | grep tRNAs | cut -f 2 -d ":")
        rRNA=$(cat $SUM_FILE | grep rRNAs | cut -f 2 -d ":")
        CDS=$(cat $SUM_FILE | grep CDSs | cut -f 2 -d ":")
        pseudo=$(cat $SUM_FILE | grep pseudogenes | cut -f 2 -d ":")
        hypo=$(cat $SUM_FILE | grep hypotheticals | cut -f 2 -d ":")
        oriC=$(cat $SUM_FILE | grep oriCs | cut -f 2 -d ":")
        oriV=$(cat $SUM_FILE | grep oriVs | cut -f 2 -d ":")
        oriT=$(cat $SUM_FILE | grep oriTs | cut -f 2 -d ":")
        cdsDens=$(cat $SUM_FILE | grep "coding density" | cut -f 2 -d ":")
        RNA5=$(cat $GFF | grep "Name=5S ribosomal RNA" | wc -l)
        RNA16=$(cat $GFF | grep "Name=16S ribosomal RNA" | wc -l)
        RNA23=$(cat $GFF | grep "Name=23S ribosomal RNA" | wc -l)
      
        
        echo -e "$SAMPLE\t$BIN\t$tRNA\t$rRNA\t$RNA5\t$RNA16\t$RNA23\t$CDS\t$pseudo\t$hypo\t$oriC\t$oriV\t$oriT\t$cdsDens" >> $STATS_OUT
    done
done

echo "$(date)   [COMPLETE] Stats: $STATS_OUT"