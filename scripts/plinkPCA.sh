#!/bin/bash 


VCF=vcfs/stam.filt.bial.snps.passed.no-sc.no-outlier.vcf.gz
OUTDIR=pca
PREFIX=stam.filt.bial.snps.passed.no-sc.no-outlier

OUTFILE=$OUTDIR/$PREFIX
mkdir -p $OUTFILE

plink \
    --vcf $VCF \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --make-bed \
    --pca \
    --out $OUTFILE/$PREFIX