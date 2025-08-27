#!/bin/bash


STATS_OUT_RAW=stats/k2-stam-prop-raw.tsv
STATS_OUT_DECON=stats/k2-stam-prop-decon.tsv

echo -e "sample\tconura\tstammerula" > $STATS_OUT_RAW
echo -e "sample\tconura\tstammerula" > $STATS_OUT_DECON

for file in $(ls 09-KRAKEN2/decon/kreport2mpa/*.mpa); do
    echo "$file"

    ID=$(basename $file .mpa)

    P_CONURA=$(cat $file | grep "s__Tephritis_conura" | cut -f 2)
    P_STAM=$(cat $file | grep "s__Candidatus_Stammerula_tephritidis" | cut -f 2)

    echo -e "$ID\t$P_CONURA\t$P_STAM" >> $STATS_OUT_RAW

done


for file in $(ls 09-KRAKEN2/raw/kreport2mpa/*.mpa); do
    echo "$file"

    ID=$(basename $file .mpa)

    P_CONURA=$(cat $file | grep "s__Tephritis_conura" | cut -f 2)
    P_STAM=$(cat $file | grep "s__Candidatus_Stammerula_tephritidis" | cut -f 2)

    echo -e "$ID\t$P_CONURA\t$P_STAM" >> $STATS_OUT_DECON

done