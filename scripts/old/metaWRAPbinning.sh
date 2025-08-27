#!/bin/bash

#SBATCH --job-name metaWRAP-bin
#SBATCH -A naiss2025-22-494
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mem=56GB
#SBATCH -t 16:30:00
#SBATCH --output=slurm-logs/metaWRAP-bin/del-SLURM-%j.out
#SBATCH --error=slurm-logs/metaWRAP-bin/del-SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Load modules
ml load PDCOLD/23.12
ml load bioinfo-tools
ml load CONCOCT/1.1.0
ml load MaxBin/2.2.7   # Maxbin2
ml load metabat/2.15   # MetaBat2
ml load bwa/0.7.17
ml load samtools/1.20
ml load hmmer/3.4-cpeGNU-23.12
ml load checkm/1.2.2-cpeGNU-23.12

# Arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --id)           ID="$2"; shift ;;
        --ref)          REF="$2"; shift ;;
        --r1)           R1="$2"; shift ;;
        --r2)           R2="$2"; shift ;;
        --outdir)       OUTDIR="$2"; shift ;;
        -h|--help)
            echo "Usage: $0 --id PREFIX --ref METAGENOME.fa --r1 R1 --r2 R2 --outdir DIR"
            exit 0 ;;
        *) echo "[ERROR] Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done


# Validate required arguments
if [[ -z $REF || -z $R1 || -z $R2 || -z $OUTDIR ]]; then
    echo "[ERROR] Missing required arguments."
    echo "Usage: $0 --id PREFIX --ref METAGENOME.fa --r1 R1 --r2 R2 --outdir DIR"
    exit 1
fi

# Check if required files actually exist
if [ ! -s "$REF" ]; then echo "[ERROR] Reference genome not found: $REF"; exit 1; fi
if [ ! -s "$R1" ]; then echo "[ERROR] forward fastq file not found or empty: $R1"; exit 1; fi
if [ ! -s "$R2" ]; then echo "[ERROR] reverse fastq file not found or empty: $R2"; exit 1; fi


# Clean any trailing "/"
OUTDIR="${OUTDIR%/}"

# Define values
CPU=10
MEMORY=56
COMP=50
CONT=10

# Define paths
WORK=$OUTDIR/$ID/tmp
BIN_RES=$OUTDIR/$ID/bins
REFINED_RES=$OUTDIR/$ID

mkdir -p $WORK $BIN_RES $REFINED_RES

echo "$(date)   [INFO] Starting metaWRAP binning for $ID with:"
echo "$(date)       Metagenome: $REF"
echo "$(date)       R1 reads: $R1"
echo "$(date)       R2 reads: $R2"

# Create temporary uncompressed fastq files
R1_TEMP="$WORK/$(basename ${R1} _R1.fq.gz)_1.fastq"
R2_TEMP="$WORK/$(basename ${R2} _R2.fq.gz)_2.fastq"

echo "$(date)   [INFO] metaWRAP friendly reads: $R1_TEMP $R2_TEMP"

echo "$(date)   [INFO] Decompressing reads..."
zcat $R1 > $R1_TEMP
zcat $R2 > $R2_TEMP


echo "$(date)   [INFO] Running metaWRAP binning..."
metawrap binning \
    -a $REF \
    -o $BIN_RES \
    -t $CPU \
    -m $MEMORY \
    -l 2500 \
    --metabat2 \
    --maxbin2 \
    --concoct \
    $R1_TEMP $R2_TEMP

# Clean up temporary files
#rm \$R1_TEMP \$R2_TEMP

echo "$(date)   [FINISH] metaWRAP binning complete for $ID"

echo "$(date)   [INFO] Starting metaWRAP bin_refinement for $ID"

# Set CHECKM data path
CHECKM_DATA_PATH=/cfs/klemming/projects/snic/snic2020-6-222/Projects/Tconura/working/Andre/databases/checkm-db

# Set checkm root
checkm data setRoot $CHECKM_DATA_PATH

# Bin refine
echo "$(date)   [INFO] Running metaWRAP bin_refinement..."
bash ../bin/metaWRAP/bin/metawrap-modules/bin_refinement.sh \
    -o $REFINED_RES/bin_refinement_c${COMP}_x${CONT} \
    -t $CPU \
    -m $MEMORY \
    -A $BIN_RES/concoct_bins \
    -B $BIN_RES/maxbin2_bins \
    -C $BIN_RES/metabat2_bins \
    -c ${COMP} \
    -x ${CONT}

echo "$(date)   [FINISH] metaWRAP bin_refinment complete for $ID"

# Load modules
ml load PDC/24.11
ml load apptainer/1.4.0-cpeGNU-24.11
ml load nextflow/24.04.2

echo "$(date)   [INFO] Starting metaWRAP bin_refinment for $ID"

