#!/bin/bash


#ml load PDC/23.12 
#ml load apptainer/1.3.6-cpeGNU-23.12

#ml load PDC/24.11
#ml load apptainer/1.4.0-cpeGNU-24.11
#ml load nextflow/24.04.2
#ml load bioinfo-tools
#ml load Graphviz/9.0.0


#wget https://raw.githubusercontent.com/EBI-Metagenomics/mettannotator/master/tests/test.csv


nextflow run -r v1.4.0 ebi-metagenomics/mettannotator \
   -profile local,apptainer \
   --input test.csv \
   --outdir test \
   --dbs ../../databases/mettannotator

#nextflow run -r v1.4.0 ebi-metagenomics/mettannotator \
#    -profile <docker/singularity/...> \
#    --input test.csv \
#    --outdir <OUTDIR> \
#    --dbs <PATH/TO/WHERE/DBS/WILL/BE/SAVED>


#nextflow run ebi-metagenomics/mettannotator \
#    -profile apptainer \
#    --input assemblies_sheet.csv \
#    --outdir 00-stam-ctg4994 \
#    --dbs ../../databases/mettannotator
