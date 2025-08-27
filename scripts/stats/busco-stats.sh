#!/bin/bash

echo "$(date) [INFO] Starting BUSCO summary"

# Define variables
WORK=/cfs/klemming/projects/snic/snic2020-6-222/Projects/Tconura/working/Andre/thesis_stammerula
STATS_OUT=$WORK/stats/busco-summary.tsv

# Header
echo -e "sample\tbin\tComp.pct\tSC.pct\tComp.Dupe.pct\tFrag.pct\tMiss.pct\tn_markers\tC\tS\tD\tF\tM" > $STATS_OUT

# Parse the BUSCO .json output files
for DIR in $WORK/07-bin-quality-assessment/*; do
    SAMPLE=$(basename $DIR)
    echo "$(date) [INFO] Summarizing BUSCO for $SAMPLE"

    for BUSCO_DIR in $DIR/busco_output/*_busco; do
        BIN=$(basename $BUSCO_DIR _busco)
        JSON_FILE=$BUSCO_DIR/short*busco.json
        
        # Use jq to parse .json
        LINE=$(jq -r '
          [
            .results["Complete percentage"],
            .results["Single copy percentage"],
            .results["Multi copy percentage"],
            .results["Fragmented percentage"],
            .results["Missing percentage"],
            .results["n_markers"],
            .results["Complete BUSCOs"],
            .results["Single copy BUSCOs"],
            .results["Multi copy BUSCOs"],
            .results["Fragmented BUSCOs"],
            .results["Missing BUSCOs"]
          ] | @tsv
        ' $JSON_FILE)
        
        echo -e "$SAMPLE\t$BIN\t$LINE" >> $STATS_OUT
    done
done

echo "$(date)   [COMPLETE] Stats: $STATS_OUT"
