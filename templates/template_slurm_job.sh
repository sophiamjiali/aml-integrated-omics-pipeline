#!/bin/bash
#SBATCH --job-name={{job_name}}
#SBATCH --output={{output_file}}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1        
#SBATCH --cpus-per-task={{n_cores}}   
#SBATCH --time=12:00:00

module load NiaEnv/2019b gcc/8.3.0 r/4.1.2

Rscript {{r_script}}
