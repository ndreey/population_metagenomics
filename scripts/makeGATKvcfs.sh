#!/bin/bash

#SBATCH --job-name GATK_VCF
#SBATCH -A naiss2025-22-494
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mem=32GB
#SBATCH -t 12:30:00
#SBATCH --output=slurm-logs/GATK_VCF/SLURM-%j.out
#SBATCH --error=slurm-logs/GATK_VCF/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL


########################################################################
# 1. Combines the GVCFs
# 2. Genotypes them to create a final metagenome VCF
# 3. Subsets the VCF to only have records for MAG
########################################################################


# Start the process
echo "$(date) [INFO]        Starting script execution"

# Load in modules
ml load gatk/4.5.0.0
ml load bcftools/1.20


while [[ "$#" -gt 0 ]]; do
    case $1 in
        --gvcf-list)    GVCF_F="$2"; shift ;;
        --ref)          REF="$2"; shift ;;
        --mag-headers)  MAG_HEADERS="$2"; shift ;;
        --outdir)       OUTDIR="$2"; shift ;;
        --prefix)        PREFIX="$2"; shift ;;
        -h|--help)
            echo "Usage: $0 --gvcf-list FILE --ref REF --mag-headers MAG --outdir DIR --prefix NAME"
            exit 0 ;;
        *) echo "[ERROR] Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done


# Validate required arguments
if [[ -z $GVCF_F || -z $REF || -z $MAG_HEADERS || -z $OUTDIR || -z $PREFIX ]]; then
    echo "[ERROR] Missing required arguments."
    echo "Usage: $0 --gvcf-list FILE --ref REF --mag-headers MAG --outdir DIR --PREFIX NAME"
    exit 1
fi

# Check if required files actually exist
if [ ! -s "$GVCF_F" ]; then echo "[ERROR] GVCF list file not found or empty: $GVCF_F"; exit 1; fi
if [ ! -s "$REF" ]; then echo "[ERROR] Reference genome not found: $REF"; exit 1; fi
if [ ! -s "$MAG_HEADERS" ]; then echo "[ERROR] MAG headers file not found: $MAG_HEADERS"; exit 1; fi


# Remove trailing slashes
OUTDIR="${OUTDIR%/}"

# Output directories
FINDIR="$OUTDIR/$PREFIX/vcf-final"
mkdir -p $FINDIR

# GATK outputs
COMBINED_GVCF=$FINDIR/${PREFIX}.metagenome.combined.g.vcf.gz
META_VCF=$FINDIR/${PREFIX}.metagenome.vcf.gz

# Raw VCFs
MAG_VCF=$FINDIR/${PREFIX}.mag.raw.vcf.gz

# Info message
echo "$(date) [START]       Processing VCF-WRANGLE for: $PREFIX"

#########################################################################
#################### COMBINE AND GENOTYPE GVCFs #########################
#########################################################################

# Combine GVCFs if the output does not already exist
if [ ! -f "$COMBINED_GVCF" ]; then
    echo "$(date) [INFO]        Combining GVCFs for GVCFs in $GVCF_F"
    gatk CombineGVCFs \
        -R $REF \
        --variant $GVCF_F \
        -O $COMBINED_GVCF
else
   echo "$(date) [SKIP]     Skipping CombineGVCFs: File already exists - $COMBINED_GVCF"
fi


# Genotype the GVCF to create a VCF file.
if [ ! -f "$META_VCF" ]; then
    echo "$(date) [INFO]      Running GenotypeGVCFs for $PREFIX"
    gatk GenotypeGVCFs \
        -R $REF \
        -V $COMBINED_GVCF \
        -O $META_VCF \
        -all-sites \
        --sample-ploidy 1
    # Index the vcf
    gatk IndexFeatureFile -I $META_VCF
else
    echo "$(date) [SKIP]        Skipping GenotypeGVCFs: File already exists - $META_VCF"
fi


# Extract only the records from the MAG
if [ ! -f $MAG_VCF ]; then 
    echo "$(date) [INFO]        Keeping only records that belong to MAG: $(basename $MAG_HEADERS)"
    REGIONS=$(paste -sd "," $MAG_HEADERS)
    bcftools view -O z -o $MAG_VCF -r "$REGIONS" $META_VCF
    gatk IndexFeatureFile -I $MAG_VCF
else
    echo "$(date) [SKIP]        Skipping MAG extraction: File already exists - $MAG_VCF"
fi

echo "$(date) [FINISH]        Done!"