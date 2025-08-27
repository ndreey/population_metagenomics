#!/bin/bash

DIR=10-ALIGNMENT/dupe-marked-bams
SUM_CSV=stats/qualimap-summary-stats.csv

# Initiate the stats file
echo "sample,n_reads,n_mapped_reads,GC,mean.Cov,std.Cov,bam_file" > $SUM_CSV

for RESULTS in $(ls $DIR/*_stats/genome_results.txt); do
    if [ -f "$RESULTS" ]; then

        # Get bam file
        BAM=$(basename $(cat $RESULTS | grep -A 2 ">>>>>>> Input" | grep "bam file =" | cut -f 2 -d "="))

        # Get sample
        SAMPLE=$(basename $BAM .sorted.marked.bam)

        # Extract number of reads
        n_reads=$(cat $RESULTS | grep -A 4 ">>>>>>> Globals" | grep "number of reads" | cut -f 2 -d "=" | tr -d " " | tr -d ",")

        # Extract number of mapped reads (from Globals inside)
        n_mapped_reads=$(cat $RESULTS | grep -A 3 ">>>>>>> Globals inside" | grep "number of mapped reads" | cut -f 2 -d "=" | tr -d " " | tr -d "," | cut -f 1 -d "(")

        # Extract GC percentage
        GC=$(cat $RESULTS | grep -A 8 ">>>>>>> ACTG content" | grep "GC percentage" | cut -f 2 -d "=" | tr -d " " | tr -d "%")

        # Extract mean coverage
        meanCov=$(cat $RESULTS | grep -A 3 ">>>>>>> Coverage" | grep "mean" | cut -f 2 -d "=" | tr -d " " | tr -d "X")

        # Extract standard deviation of coverage
        stdCov=$(cat $RESULTS | grep -A 3 ">>>>>>> Coverage" | grep "std" | cut -f 2 -d "=" | tr -d " " | tr -d "X")

        # Write to stats file
        echo "$SAMPLE,$n_reads,$n_mapped_reads,$GC,$meanCov,$stdCov,$BAM" >> $SUM_CSV

    else
        echo "File not found for $RESULTS"
    fi

done
