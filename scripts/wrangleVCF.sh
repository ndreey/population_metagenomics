#!/bin/bash

#SBATCH --job-name wrangleVCF
#SBATCH -A naiss2025-22-494
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=16GB
#SBATCH -t 00:30:00
#SBATCH --output=slurm-logs/wrangleVCF/SLURM-%j.out
#SBATCH --error=slurm-logs/wrangleVCF/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL


########################################################################
# 1. Combines the GVCFs
# 2. Genotypes them to create a final metagenome VCF
# 3. Subsets the VCF to only have records for MAG
########################################################################



# Parse argument flags
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --vcf)          VCF_F="$2"; shift ;;
        --ref)          REF="$2"; shift ;;
        --outdir)       OUTDIR="$2"; shift ;;
        --prefix)       PREFIX="$2"; shift ;;
        -h|--help)
            echo "Usage: $0 --vcf FILE --REF FILE --outdir DIR --prefix STRING"
            exit 0 ;;
        *) echo "[ERROR] Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done


# Validate required arguments
if [[ -z $VCF_F || -z $OUTDIR || -z $REF || -z $PREFIX ]]; then
    echo "[ERROR] Missing required arguments."
    echo "Usage: $0 --vcf FILE --ref FILE --outdir DIR --prefix STRING"
    exit 1
fi

# Check if required files actually exist
if [ ! -s "$VCF_F" ]; then echo "[ERROR] VCF file not found: $VCF_F"; exit 1; fi
if [ ! -s "$REF" ]; then echo "[ERROR] reference file not found: $REF"; exit 1; fi



# Start the process
echo "$(date) [INFO]        Starting script execution"

# Switch from GATK to VCF processing tools
ml load gatk/4.5.0.0
ml load vcftools/0.1.16
ml load bcftools/1.20




# Remove trailing slash from --outdir, if present
OUTDIR="${OUTDIR%/}"

# Output directories
FINDIR="$OUTDIR/vcf-final"
STATSDIR="$OUTDIR/stats"
WORKDIR="$OUTDIR/work-dir"
mkdir -p $FINDIR $STATSDIR $WORKDIR

# Raw VCFs
RAW_VARIANTS=${FINDIR}/${PREFIX}.raw.variants.vcf.gz
INVARIANTS=${FINDIR}/${PREFIX}.invariants.vcf.gz
RAW_BIALLELIC=${FINDIR}/${PREFIX}.bial.snps.vcf.gz

# Filtering intermediate steps
GATK_FLAGGED=${WORKDIR}/${PREFIX}.GATK-flagged.vcf.gz
GATK_FILT_VCF=${WORKDIR}/${PREFIX}.filt.var.vcf.gz
THINNED=${WORKDIR}/${PREFIX}.GATK.filt.thinned
INDEL_PROX=${WORKDIR}/${PREFIX}.GATK.filt.indel.thinned.vcf.gz
FILT_TMP=${WORKDIR}/${PREFIX}.tmp.vcf.gz

# Final filtered vcfs
FILT_BIAL=${FINDIR}/${PREFIX}.filt.bial.snps.vcf.gz
FILT_MISS=${FINDIR}/${PREFIX}.filt.no-miss.vcf.gz
FIN_VCF=${FINDIR}/${PREFIX}.filt.no-miss.complete.vcf.gz

# Summary files
STATS_OUT=$STATSDIR/${PREFIX}-bcftools.stats
RAW_STATS_OUT=$STATSDIR/${PREFIX}-raw.biallelic.base.stats.tsv
FILT_STATS_OUT=$STATSDIR/${PREFIX}-filtered.biallelic.base.stats.tsv

echo "$(date) [START]       Processing VCF-WRANGLE for: $(basename $VCF_F)"

#########################################################################
#                     Prepare VCF and Reference
#########################################################################
# Check so index of original vcf exists
if [ ! -f "$VCF_F.tbi" ]; then
    echo "$(date) [INFO]        Create index for $(basename $VCF_F)"
    gatk IndexFeatureFile -I $VCF_F
else
    echo "$(date) [SKIP]        Index file for $(basename $VCF_F) exists"
fi

# Check so index of reference exists
if [ ! -f "$REF.fai" ]; then
    echo "$(date) [INFO]        Create samtools index for $(basename $REF)"
    samtools faidx $REF
else
    echo "$(date) [SKIP]        Index file for $(basename $REF) exists"
fi

# Check so dictionary of reference exists
if [ ! -f "$(basename $REF .fa).dict" ]; then
    echo "$(date) [INFO]        Create GATK dictionary for $(basename $REF)"
    gatk CreateSequenceDictionary -R $REF
else
    echo "$(date) [SKIP]        Dictionary file for $(basename $REF) exists"
fi

#########################################################################
#                   Split VCF into variants and invariants
#########################################################################
# Generate VCF file with invariants only
if [ ! -f "$INVARIANTS" ]; then
    echo "$(date) [INFO]        Splitting $(basename $VCF_F) into invariants only: $INVARIANTS"
    bcftools view -C 0 -O z -o "$INVARIANTS" "$VCF_F"
    tabix "$INVARIANTS"
else
    echo "$(date) [SKIP]        Invariants already exist: $INVARIANTS"
fi

# Generate VCF file for all variants (snps, indels, etc)
if [ ! -f "$RAW_VARIANTS" ]; then
    echo "$(date) [INFO]        Splitting $(basename $VCF_F) into variants only: $RAW_VARIANTS"
    bcftools view -c 1 -O z -o "$RAW_VARIANTS" "$VCF_F"
    tabix "$RAW_VARIANTS"
else
    echo "$(date) [SKIP]        Raw variants already exist: $RAW_VARIANTS"
fi

# Generate VCF file with only raw biallelic SNPs
if [ ! -f "$RAW_BIALLELIC" ]; then
    echo "$(date) [INFO]        Splitting $(basename $VCF_F) into only biallelic SNPs: $RAW_BIALLELIC"
    bcftools view -m2 -M2 -v snps -O z -o "$RAW_BIALLELIC" "$RAW_VARIANTS"
    tabix "$RAW_BIALLELIC"
else
    echo "$(date) [SKIP]        Raw biallelic SNPs already exist: $RAW_BIALLELIC"
fi

#########################################################################
#                       Filter VCF
#########################################################################
# Apply GATK filteration
if [ ! -f $GATK_FLAGGED ]; then
    echo "$(date) [EXEC]      Filtering $RAW_BIALLELIC according to GATK best practices"
    echo "$(date) [INFO]      Flagging records below/above the following:"
    echo "$(date) [FILTER]          QD < 2.0"
    echo "$(date) [FILTER]          SOR > 3.0"
    echo "$(date) [FILTER]          FS > 60.0"
    echo "$(date) [FILTER]          MQ < 40.0"
    
    gatk VariantFiltration \
        -V $RAW_VARIANTS \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -O $GATK_FLAGGED
else
    echo "$(date) [SKIP]        Skipping GATK filtering: File already exists - $GATK_FLAGGED"
fi

# Keep only records flagged PASS
if [ ! -f $GATK_FILT_VCF ]; then
    echo -e "\n$(date) [INFO]        Keeping only records that PASSED GATK filtering: $(basename $GATK_FILT_VCF)"
    bcftools view -f PASS -O z -o $GATK_FILT_VCF $GATK_FLAGGED
    tabix $GATK_FILT_VCF
    #rm $GATK_FLAGGED $GATK_FLAGGED.tbi
else
    echo "$(date) [SKIP]        Skipping keeping PASS flags: File already exists - $GATK_FILT_VCF"
fi

# Remove SNPs within 3bp of another SNP
if [ ! -f "$THINNED.recode.vcf" ]; then
    echo "$(date) [FILTER]      Removing SNPs within 3bp of another SNP"
    vcftools --gzvcf $GATK_FILT_VCF --thin 3 --recode --out $THINNED
else
    echo "$(date) [SKIP]        Skipping SNP thinning: File already exists - $THINNED.recode.vcf"
fi

# Remove SNPs within 10bp of INDEL
if [ ! -f $INDEL_PROX ]; then
    echo "$(date) [FILTER]      Removing SNPs within 10bp of INDEL"
    bcftools filter -G 10 -o $INDEL_PROX $THINNED.recode.vcf
else
    echo "$(date) [SKIP]        Skipping INDEL proximity filtering: File already exists - $INDEL_PROX"
fi

# Generate VCF file with only filtered biallelic SNPs
if [ ! -f "$FILT_BIAL" ]; then
    echo "$(date) [INFO]        Splitting $(basename $INDEL_PROX) into only biallelic SNPs: $FILT_BIAL"
    # Get only biallelic snps
    bcftools view -m2 -M2 -v snps -O z -o "$FILT_TMP" "$INDEL_PROX"
    tabix "$FILT_TMP"
    
    # Extract from raw variants to retain the stats info.
    bcftools view -R $FILT_TMP -O z -o $FILT_BIAL $RAW_BIALLELIC
    tabix $FILT_BIAL

else
    echo "$(date) [SKIP]        Raw biallelic SNPs already exist: $RAW_BIALLELIC"
fi


## Discard records with high missigness
#if [ ! -f $FILT_MISS ]; then
#    vcftools \
#        --gzvcf $FILT_BIAL \
#        --max-missing 0.9 \
#        --recode --stdout | bgzip > $FILT_MISS
#
#    tabix $FILT_MISS
#else
#    echo "$(date) [SKIP]        Skipping removing records with high missigness: File already exists - $FILT_MISS"
#fi


##########################################################################
##                  Merge filtered variants with invariants
##########################################################################
## Merge invariants and biallelic SNPs vcf.
#if [ ! -f "$FIN_VCF" ]; then
#    echo "$(date) [INFO]        Merging invariants with filtered variants: $FIN_VCF"
#    echo "$(date) [INFO]        Merging $(basename $INVARIANTS) + $(basename $FILT_MISS)"
#    gatk MergeVcfs \
#    -I $INVARIANTS \
#    -I $FILT_BIAL \
#    -O $FIN_VCF
#
#else
#    echo "$(date) [SKIP]        Filtered merged VCF already exists: $FIN_VCF"
#fi

echo -e "$(date) [FINISH]       VCF-wrangle complete!\n"


#########################################################################
#                  Summarize BCFTOOLS Stats
#########################################################################
# Get quick stats for the files
echo -e "$(date) [INFO]        Summarising stats from each VCF file: $STATS_OUT"
echo -e "VCF\tsamples\trecords\tinvariants\tSNPs\tMNPs\tindels\tothers\tmulti-sites\tmulti-SNPs" > $STATS_OUT

for FILE in $VCF_F $INVARIANTS $RAW_VARIANTS $INDEL_PROX $RAW_BIALLELIC $FILT_BIAL; do

    echo -e "$(date) [INFO]        Summarising $FILE"
    echo -e "$(basename $FILE)\t$(bcftools stats $FILE | grep -v "# SN," | grep -A 9 "# SN" | grep -v "#" | cut -f 4 | paste -sd$'\t')" >> $STATS_OUT

done

# Preview the file
echo -e ">> BCFTOOLS STATS DESCRIPTION << \n"
echo "number of records   .. number of data rows in the VCF"
echo "number of no-ALTs   .. reference-only sites, ALT is either '.' or identical to REF"
echo "number of SNPs      .. number of rows with a SNP"
echo "number of MNPs      .. number of rows with a MNP, such as CC>TT"
echo "number of indels    .. number of rows with an indel"
echo "number of others    .. number of rows with other type, for example a symbolic allele or a complex substitution, such as ACT>TCGA"
echo "number of multiallelic sites     .. number of rows with multiple alternate alleles"
echo "number of multiallelic SNP sites .. number of rows with multiple alternate alleles, all SNPs"
echo "Note that rows containing multiple types will be counted multiple times, in each counter. For example, a row with a SNP and an indel increments both the SNP and the indel counter."
echo ""

# Print out the stats in a clean table
echo -e "\n >> BCFTOOLS STATS <<"
column -t "$STATS_OUT"


#########################################################################
#                  Validate VCF files and Summarize Errors
#########################################################################
# Run GATK ValidateVariants on all output VCFs
echo ""
echo "$(date) [INFO]        Running GATK ValidateVariants on final VCFs with --warn-on-errors"
for FILE in $RAW_BIALLELIC $FILT_BIAL; do
    OUTLOG=$STATSDIR/$(basename "$FILE" .vcf.gz)-validate.log
    if [ ! -f "$OUTLOG" ]; then
        gatk ValidateVariants -V "$FILE" -R $REF --warn-on-errors &> "$OUTLOG"
    else
        echo "$(date) [SKIP]        Validation log already exists: $OUTLOG"
    fi
done

# Prepare summary output
WARN_SUMMARY=$STATSDIR/${PREFIX}-validation-warnings.summary.tsv
echo -e "VCF\twarning_count" > "$WARN_SUMMARY"

# Process each file's validation output
for FILE in $RAW_BIALLELIC $FILT_BIAL; do
    BASENAME=$(basename "$FILE")
    LOGFILE=$STATSDIR/${BASENAME%.vcf.gz}-validate.log
    OUT_TSV=$STATSDIR/${BASENAME%.vcf.gz}-validate.tsv

    # Write header to per-file TSV
    echo -e "file\tposition\ttype\terror" > "$OUT_TSV"

    if [ -f "$LOGFILE" ]; then
        # Extract and write warnings to per-file TSV
        grep "fails strict validation" "$LOGFILE" | \
        sed -E 's/.*fails strict validation of type ([^:]+): (.*) at position ([^ ]+).*/\3\t\1\t\2/' | \
        awk -v f="$BASENAME" '{print f "\t" $0}' >> "$OUT_TSV"
    fi

    # Count how many warnings were written (excluding header)
    COUNT=$(($(wc -l < "$OUT_TSV") - 1))
    echo -e "$BASENAME\t$COUNT" >> "$WARN_SUMMARY"
done

# Print summary to screen
echo -e "\n >> VALIDATION WARNINGS SUMMARY <<\n"
column -t "$WARN_SUMMARY"


#########################################################################
#                                  VCF Stats
#########################################################################


echo "$(date) [INFO]        Collecting site-level QC stats"

# RAW BIALLELIC
if [ ! -f "$RAW_STATS_OUT" ]; then
    echo "$(date) [INFO]        Writing raw biallelic site stats: $RAW_STATS_OUT"

    echo -e "CHROM\tPOS\tQUAL\tDP\tQD\tFS\tSOR\tMQ\tMQRankSum\tReadPosRankSum\tAN\tAC\tAF" > "$RAW_STATS_OUT"

    bcftools query -f '%CHROM\t%POS\t%QUAL\t%DP\t%QD\t%FS\t%SOR\t%MQ\t%MQRankSum\t%ReadPosRankSum\t%AN\t%AC\t%AF\n' \
        "$RAW_BIALLELIC" >> "$RAW_STATS_OUT"
else
    echo "$(date) [SKIP]        Raw biallelic site stats already exist: $RAW_STATS_OUT"
fi

# Filtered biallelic
if [ ! -f "$FILT_STATS_OUT" ]; then
    echo "$(date) [INFO]        Writing raw biallelic site stats: $FILT_STATS_OUT"

    echo -e "CHROM\tPOS\tQUAL\tDP\tQD\tFS\tSOR\tMQ\tMQRankSum\tReadPosRankSum\tAN\tAC\tAF" > "$FILT_STATS_OUT"

    bcftools query -f '%CHROM\t%POS\t%QUAL\t%DP\t%QD\t%FS\t%SOR\t%MQ\t%MQRankSum\t%ReadPosRankSum\t%AN\t%AC\t%AF\n' \
        "$FILT_BIAL" >> "$FILT_STATS_OUT"
else
    echo "$(date) [SKIP]        Raw biallelic site stats already exist: $FILT_STATS_OUT"
fi

# Get vcf stats
for VCF_STAT in $RAW_BIALLELIC $FILT_BIAL; do

    ID=$(basename $VCF_STAT .vcf.gz)
    STAT_OUT=$STATSDIR/${PREFIX}-$ID
    vcftools --gzvcf "$VCF_STAT" --site-mean-depth --out "$STAT_OUT"          # Get the per-site mean depth average across samples (.ldepth.mean)
    vcftools --gzvcf "$VCF_STAT" --missing-site --out "$STAT_OUT"             # Get the proportion of missingness per site (.lmiss)
    vcftools --gzvcf "$VCF_STAT" --depth --out "$STAT_OUT"                    # Get mean depth per sample (.idepth)
    vcftools --gzvcf "$VCF_STAT" --missing-indv --out "$STAT_OUT"             # Get the proportion of missingness per sample (.imiss)
    vcftools --gzvcf "$VCF_STAT" --SNPdensity 10000 --out "$STAT_OUT"         # Get SNP density across the genome (.snpden)
done

echo -e "\n$(date) [FINISH]       Script Complete!\n"