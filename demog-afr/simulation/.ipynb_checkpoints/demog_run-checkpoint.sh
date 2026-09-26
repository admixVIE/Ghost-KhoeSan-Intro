#!/bin/bash
#
#SBATCH --job-name=dmg
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=14:00:00
#SBATCH --output=/lisc/home/user/kuhlwilm/logs/dm_%j.out
#SBATCH --error=/lisc/home/user/kuhlwilm/logs/dm_%j.err
#SBATCH --array=2544,2881,2886,2996,3017

ID=$SLURM_ARRAY_TASK_ID
R_LIBS=/lisc/data/scratch/admixlab/mk_data/rlib/

## load module
module load R-demografr/0.1.0-foss-2025a-R-4.5.1

## run script
Rscript ~/demog-afr/simulation/demog_sim.R ${ID}


exit
