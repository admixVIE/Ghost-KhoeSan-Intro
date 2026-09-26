#!/bin/bash
#
#SBATCH --job-name=dmred
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=8:00:00
#SBATCH --output=/lisc/home/user/kuhlwilm/logs/wsr_%j.out
#SBATCH --error=/lisc/home/user/kuhlwilm/logs/wsr_%j.err
#SBATCH --array=2025,2026,2044,2047,2051,2053,2058,2069,2073,2078,2080,2082,2083,2085,2086,2102,2195,2261,2314,2366,2367

ID=$SLURM_ARRAY_TASK_ID
R_LIBS=/lisc/data/scratch/admixlab/mk_data/rlib/

## load module
module load R-demografr/0.1.0-foss-2025a-R-4.5.1

## run script
Rscript ~/demog-afr/simulation/model_comp_wss_red.R ${ID}


exit
