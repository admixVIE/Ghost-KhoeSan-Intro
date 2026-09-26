#!/bin/bash
#
#SBATCH --job-name=dmcomp
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=6:00:00
#SBATCH --output=/lisc/home/user/kuhlwilm/logs/mg_%j.out
#SBATCH --error=/lisc/home/user/kuhlwilm/logs/mg_%j.err
#SBATCH --array=1779

ID=$SLURM_ARRAY_TASK_ID
R_LIBS=/lisc/data/scratch/admixlab/mk_data/rlib/

## load module
module load R-demografr/0.1.0-foss-2025a-R-4.5.1

## run script
Rscript ~/demog-afr/simulation/model_comp.R ${ID}


exit
