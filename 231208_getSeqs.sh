#!/bin/bash 

#SBATCH -p general
#SBATCH -N 1 
#SBATCH --mem=10g
#SBATCH -n 1 
#SBATCH -t 1:00:00
#SBATCH -o 231208_getSeqs.out
#SBATCH -e 231208_getSeqs.err 
#SBATCH --mail-type=end 
#SBATCH --mail-user=askav@email.unc.edu 

module load r/4.3.1
Rscript ~/rna/autoregulation/231208_getSeqs.R
