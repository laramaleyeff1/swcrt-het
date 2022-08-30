#!/bin/bash
#SBATCH -c 10
#SBATCH -N 1
#SBATCH -p shared
#SBATCH -t 0-05:00
#SBATCH --array=1-100
#SBATCH --mem=5G
#SBATCH	-o %j.out
#SBATCH	--mail-type=ALL
#SBATCH	--mail-user=laramaleyeff@gmail.com

module load gcc/7.1.0-fasrc01 R/4.0.5-fasrc01
Rscript main.R $1 $2 $3 ${SLURM_ARRAY_TASK_ID}
