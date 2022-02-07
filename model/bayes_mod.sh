#!/bin/bash

#SBATCH --job-name=bayes_mod
#SBATCH --time=11-0
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=50G

# Load modules
module load r/4.1.0
module load gcc/6.3.0

# Run rscript
Rscript script_model_bayes_trial.R