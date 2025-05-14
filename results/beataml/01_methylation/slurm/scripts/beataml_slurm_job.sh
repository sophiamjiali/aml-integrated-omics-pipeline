#!/bin/bash
#SBATCH --job-name=beataml_dmrs
#SBATCH --output=/scratch/.../beataml_dmrs_generation.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1        
#SBATCH --cpus-per-task=40   
#SBATCH --time=12:00:00

module load NiaEnv/2019b gcc/8.3.0 r/4.1.2

Rscript beataml_slurm_analysis.R
