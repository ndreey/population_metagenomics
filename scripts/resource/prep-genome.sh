#!/bin/bash

#SBATCH --job-name index
#SBATCH -A naiss2025-22-494
#SBATCH --array=1-3
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16GB
#SBATCH -t 01:30:00
#SBATCH --output=slurm-logs/map2meta/index-SLURM-%j-%a.out
#SBATCH --error=slurm-logs/map2meta/index-SLURM-%j-%a.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Start time and date
echo "$(date)[START]     Loading module and starting script execution"

JOBID=${SLURM_ARRAY_TASK_ID}
GENOME=$1
PREFIX=$1

if [ $JOBID -eq 1 ]; then
    # Load in modules
    ml load bwa/0.7.17

    echo "$(date)[INFO]      bwa index: $GENOME"

    # Index genome
    bwa index -p ${PREFIX} ${GENOME}

    # End time and date
    echo "$(date)[FINISHED]"
fi


if [ $JOBID -eq 2 ]; then
    # Load in modules
    ml load gatk/4.5.0.0

    echo "$(date)[INFO]      gatk CreateSequenceDictionary: $GENOME"

    # Create dictionary
    gatk CreateSequenceDictionary -R $GENOME 

    # End time and date
    echo "$(date)[FINISHED]"
fi    


if [ $JOBID -eq 3 ]; then
    # Load in modules
    ml load samtools/1.20

    echo "$(date)[INFO]      samtools faidx: $GENOME"

    # Index genome
    samtools faidx ${GENOME}

    # End time and date
    echo "$(date)[End]      Finished indexing genome with samtools"
fi