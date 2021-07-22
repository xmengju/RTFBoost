#!/bin/bash
#SBATCH --job-name=FAM
#SBATCH --output=arrayJob_%A_%a.out
#SBATCH --error=arrayJob_%A_%a.err
#SBATCH --array=1-200
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --mail-user=xiaomeng.ju@stat.ubc.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rrg-matiass


######################
# Begin work section #
######################

# Print this sub-job's task ID
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

# Do some work based on the SLURM_ARRAY_TASK_ID
# For example: 
# ./my_process $SLURM_ARRAY_TASK_ID
# 
# where my_process is you executable

# Run sbatch --export=family='logistic' cedar.shs
# Run sbatch --export=family='logistic' cedar.shs
module load StdEnv/2018.3
module spider r
module load  nixpkgs/16.09  gcc/8.3.0 r/3.6.0

srun Rscript Code/Exp_A/conduct_exp.R  $SLURM_ARRAY_TASK_ID 1