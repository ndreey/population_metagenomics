#!/bin/bash

# Start the process
echo "$(date) [INFO]        Starting script execution"

WORK=/cfs/klemming/projects/snic/snic2020-6-222/Projects/Tconura/working/Andre/thesis_stammerula
STATS_OUT=$WORK/stats/checkm2.tsv

# Header
echo -e "sample\tbin\tCompleteness_General\tContamination\tCompleteness_Specific\tCompleteness_Model_Used\tTranslation_Table_Used\tCoding_Density\tContig_N50\tAverage_Gene_Length\tGenome_Size\tGC_Content\tTotal_Coding_Sequences\tTotal_Contigs\tMax_Contig_Length\tAdditional_Notes" > $STATS_OUT

# Parse the checkm2 output files
for DIR in $WORK/07-bin-quality-assessment/*; do
    SAMPLE=$(basename $DIR)
    CHECKM2=$DIR/checkm2/checkm2_output/quality_report.tsv
    echo "$(date) [INFO] Summarizing CheckM2 for $SAMPLE"

    TMP=$(mktemp)

    cat $CHECKM2 | awk -v s="$SAMPLE" 'NR==1{print "sample\t"$0; next} {print s"\t"$0}' > $TMP
    cat $TMP | grep -v "sample" >> $STATS_OUT

    # Remove tmp
    rm $TMP

done

echo "$(date) [COMPLETE]        Stats: $STATS_OUT"