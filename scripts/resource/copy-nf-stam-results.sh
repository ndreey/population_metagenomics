#!/bin/bash

#SBATCH --job-name copy
#SBATCH -A naiss2025-22-494
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=12GB
#SBATCH -t 02:30:00
#SBATCH --output=slurm-logs/copyResults/SLURM-%j.out
#SBATCH --error=slurm-logs/copyResults/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Start the process
echo "$(date) [INFO]        Starting script execution"

nfSTAM_RES=/cfs/klemming/projects/snic/snic2020-6-222/Projects/Tconura/working/Andre/stammerula2025/results

# Follow the link and copy over the results.
echo "$(date) [INFO]        Copying decontaminated long reads"
rsync -aL --progress $nfSTAM_RES/02-decontamination/clean-reads-lr .
mv clean-reads-lr/ 02-clean-reads-lr

echo "$(date) [INFO]        Copying decontaminated sample merged reads"
rsync -aL --progress $nfSTAM_RES/03-sample-merged-sr .

echo "$(date) [INFO]        Copying decontaminated pop merged reads"
rsync -aL --progress $nfSTAM_RES/04-pop-merged-sr .

echo "$(date) [INFO]        Copying metagenomes"
rsync -aL --progress $nfSTAM_RES/05-metagenomes .

echo "$(date) [INFO]        Copying refined bins"
rsync -aL --progress $nfSTAM_RES/06-metaWRAP-refined-bins .

echo "$(date) [INFO]        Copying bin quality assessment"
rsync -aL --progress $nfSTAM_RES/07-bin-quality-assessment .

echo "$(date) [INFO]        Copying stats"
rsync -aL --progress $nfSTAM_RES/stats .

echo "$(date) [COMPLETE]        Done!"