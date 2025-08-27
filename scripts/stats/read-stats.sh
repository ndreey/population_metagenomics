#!/bin/bash

# Start the process
echo "$(date) [INFO]        Starting script execution"

# Define variables
TMPsr=$(mktemp)
TMPlr=$(mktemp)
STATS_OUT=stats/read-stats.tsv

SR_RAW_STATS=stats/sr-raw-seq-stats.tsv
SR_TRIM_STATS=stats/sr-trim-seq-stats.tsv
SR_DECON_STATS=stats/sr-decon-seq-stats.tsv

LR_RAW_STATS=stats/lr-raw-seq-stats.tsv
LR_DECON_STATS=stats/lr-decon-seq-stats.tsv

# Create tmp
cat $SR_RAW_STATS | cut -f 1 | grep -v "^file" > $TMPsr
cat $LR_DECON_STATS | cut -f 1 | grep -v "^file" > $TMPlr

# Header
echo -e "read\traw\ttrim\tdecon\tfraction.kept" > $STATS_OUT

# Summarize short
while read READ; do

    FILE=$(basename $READ)
    ID=$(echo $FILE | cut -f 1 -d "." | sed "s/_001//")

    RAW=$(cat $SR_RAW_STATS | grep $ID | cut -f 4)

    TRIM=$(cat $SR_TRIM_STATS | grep $ID | cut -f 4)

    DECON=$(cat $SR_DECON_STATS | grep $ID | cut -f 4)

    KEPT=$(echo "scale=4; $DECON / $RAW" | bc -l)


    echo -e "$ID\t$RAW\t$TRIM\t$DECON\t$KEPT" >> $STATS_OUT

done < $TMPsr

# Summarize long
while read READ; do

    ID=$(basename $READ -clean.fq.gz)
    CELL=$(echo $ID | cut -f 3 -d "_")

    RAW=$(cat $LR_RAW_STATS | grep $CELL | cut -f 4)

    TRIM="na"

    DECON=$(cat $LR_DECON_STATS | grep $CELL | cut -f 4)

    KEPT=$(echo "scale=4; $DECON / $RAW" | bc -l)


    echo -e "$ID\t$RAW\t$TRIM\t$DECON\t$KEPT" >> $STATS_OUT

done < $TMPlr

rm $TMPsr $TMPlr

echo "$(date) [COMPLETE]        Done!"