#!/bin/bash

#SBATCH --job-name map2meta
#SBATCH -A naiss2025-22-494
#SBATCH --array=1-115
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mem=24GB
#SBATCH -t 08:30:00
#SBATCH --output=slurm-logs/map2meta/fixed-SLURM-%j-%a.out
#SBATCH --error=slurm-logs/map2meta/fixed-SLURM-%j-%a.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL


########################################################################
# 1. Align to metagenome using BWA-MEM
# 2. Marks PCR duplicates with picard
# 3. Call variants to a GVCF with GATK HaplotypeCaller
########################################################################


# Start time and date
echo "$(date) [START]     Starting script execution"

ml load bioinfo-tools
ml load bwa/0.7.17
ml load samtools/1.20
ml load picard/3.3.0
ml load gatk/4.5.0.0

# Get SLURM variables
JOBID=${SLURM_ARRAY_TASK_ID}

# Arguments
META_F=doc/sample-merged-sr-paths.txt
OUTDIR=10-ALIGNMENT
BAM_DIR=$OUTDIR/bams
MARKED_DIR=$OUTDIR/dupe-marked-bams
METRICS_DIR=$OUTDIR/dupe-metrics
GVCF_DIR=11-GVCFs
META_G=metagenome/pt_042.contigs.fa
CPU=10

mkdir -p $BAM_DIR $GVCF_DIR $MARKED_DIR $METRICS_DIR


sed -n "${JOBID}p" "$META_F" | while read R1 R2; do

	# Get unique sample name
	ID=$(basename $R1 _R1.fq.gz)
	POP=$(basename $R1 | cut -f 1 -d "_")
	SAMPLE=$(basename $R1 | cut -f 2,3 -d "_")

	BAM_SORTED=$BAM_DIR/$ID.sorted.bam
	BAM_MARKED=$MARKED_DIR/$ID.sorted.marked.bam
	DUPE_STATS=$METRICS_DIR/$ID.dupe-metrics.txt
	GVCF_OUT=${GVCF_DIR}/${ID}.g.vcf.gz

echo "$(date) [INFO]      Processing: $ID"
echo "$(date) [INFO]      Array number: $JOBID"

#########################################################################
##################### Align Against Metagenome ##########################
#########################################################################

	### Align to Metagenome
	if [ ! -f "$BAM_SORTED" ]; then
	echo "$(date) [EXEC]        Aligning reads to metagenome: $(basename $META_G)"
	    bwa mem \
	        -R "@RG\tID:${ID}\tSM:${ID}\tLB:${SAMPLE}\tPL:ILLUMINA" \
	        -t $CPU $META_G "$R1" "$R2" | \
	        samtools view -h -b -@ $CPU | \
	        samtools sort -@ $CPU --write-index -o $BAM_SORTED -
	
	else
	    echo "$(date) [SKIP]     Skipping Alignment: $BAM_SORTED already exists"
	fi

	# Mark Duplicates
	if [ ! -f "$BAM_MARKED" ]; then
		echo "$(date) [EXEC]        Marking duplicates with picard"
		picard MarkDuplicates \
			I=$BAM_SORTED \
			O=$BAM_MARKED \
			M=$DUPE_STATS \
			REMOVE_DUPLICATES=false \
			ASSUME_SORTED=true \
			CREATE_INDEX=true
	else
		echo "$(date) [SKIP]     Skipping picard: $BAM_MARKED already exists"
	fi

	# Remove temporary files
	echo "$(date) [INFO]          Removing intermediate files"
	rm $BAM_SORTED $BAM_SORTED.csi

		
#########################################################################
######################## Call GVCFs with GATK ###########################
#########################################################################

	# Generate a GVCF for the marked duped bam file
	if [ ! -f "$GVCF_OUT" ]; then
		echo "$(date) [EXEC]        Calling GVCF with GATK HaplotypeCaller for $ID"
		gatk HaplotypeCaller \
			-R $META_G \
			-I $BAM_MARKED\
			-O $GVCF_OUT \
			--emit-ref-confidence GVCF \
			--sample-ploidy 1 \
			--native-pair-hmm-threads $CPU
	else
		echo "$(date) [SKIP]        GVCF already exists for $ID"
	fi

	# Finished
	echo "$(date) [COMPLETE]        map2meta completed for sample $ID"

done