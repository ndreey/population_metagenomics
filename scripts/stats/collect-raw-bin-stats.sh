#!/bin/bash

# Start the process
echo "$(date) [INFO]        Starting script execution"

WORK=/cfs/klemming/projects/snic/snic2020-6-222/Projects/Tconura/working/Andre/thesis_stammerula
STATS_OUT=$WORK/stats/number-of-bins.tsv

# Create summary for number of bins
echo -e "id\tmaxbin2\tmetabat2\tconcoct\ttot.raw.bins\trefined" > $STATS_OUT

for DIR in $WORK/06-metaWRAP-refined-bins/*; do 
    
    ID=$(basename $DIR)
    ID_DIR=$WORK/06-metaWRAP-refined-bins/$ID/bin_refinement
    
    echo "$(date) [INFO]        Summarizing number of bins for $ID"
    NUMB_CONCOCT_BINS=$(ls $ID_DIR/concoct_bins/*.fa | wc -l)
    NUMB_MAXBIN2_BINS=$(ls $ID_DIR/maxbin2_bins/*.fa | wc -l)
    NUMB_METABAT2_BINS=$(ls $ID_DIR/metabat2_bins/*.fa | wc -l)
    NUMB_REFINED_BINS=$(ls $ID_DIR/metawrap*_bins/*.fa | wc -l)

    TOTAL_RAW=$(($NUMB_CONCOCT_BINS + $NUMB_MAXBIN2_BINS + $NUMB_METABAT2_BINS))

    echo -e "$ID\t$NUMB_MAXBIN2_BINS\t$NUMB_METABAT2_BINS\t$NUMB_CONCOCT_BINS\t$TOTAL_RAW\t$NUMB_REFINED_BINS" >> $STATS_OUT
done

echo "$(date) [COMPLETE]        Stats: $STATS_OUT"