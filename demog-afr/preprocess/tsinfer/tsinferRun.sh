#!/bin/bash
#SBATCH --job-name=tsinfer_chr
#SBATCH --output=logs/tsinfer_chr_%A_%a.out
#SBATCH --error=logs/tsinfer_chr_%A_%a.err
#SBATCH --array=1-23
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --time=110:00:00


mkdir -p logs


if [ ${SLURM_ARRAY_TASK_ID} -eq 23 ]; then
    CHR="X"
else
    CHR=${SLURM_ARRAY_TASK_ID}
fi


echo "Starting chromosome ${CHR}"


python runInferTskit.py ${CHR}