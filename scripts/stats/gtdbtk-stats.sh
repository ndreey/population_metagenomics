#!/bin/bash

# Start the process
echo "$(date) [INFO]        Starting script execution"

WORK=/cfs/klemming/projects/snic/snic2020-6-222/Projects/Tconura/working/Andre/thesis_stammerula
STATS_OUT=$WORK/stats/gtdbtk.tsv

# Header
echo -e "sample\tbin\tclassification\tclosest_genome_reference\tclosest_genome_reference_radius\tclosest_genome_taxonomy\tclosest_genome_ani\tclosest_genome_af\tclosest_placement_reference\tclosest_placement_radius\tclosest_placement_taxonomy\tclosest_placement_ani\tclosest_placement_af\tpplacer_taxonomy\tclassification_method\tnote\tother_related_references(genome_id,species_name,radius,ANI,AF)\tmsa_percent\ttranslation_table\tred_value\twarnings" > $STATS_OUT

# Parse the checkm2 output files
for DIR in $WORK/07-bin-quality-assessment/*; do
    SAMPLE=$(basename $DIR)
    GTDBTK=$DIR/GTDB-Tk/gtdbtk_output/gtdbtk.bac120*.tsv

    echo "$(date) [INFO] Summarizing GTDB-Tk for $SAMPLE"

    TMP=$(mktemp)

    cat $GTDBTK | awk -v s="$SAMPLE" 'NR==1{print "sample\t"$0; next} {print s"\t"$0}' > $TMP
    cat $TMP | grep -v "sample" >> $STATS_OUT

    # Remove tmp
    rm $TMP
done

echo "$(date) [COMPLETE]        Stats: $STATS_OUT"