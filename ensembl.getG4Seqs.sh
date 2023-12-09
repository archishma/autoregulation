#!/bin/bash

#SBATCH -N 1 
#SBATCH -n 1 
#SBATCH -e g4_%j.err
#SBATCH --mem=15g
#SBATCH -t 1:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=askav@email.unc.edu

cd /nas/longleaf/home/askav/rna/autoregulation/g4

for seqfile in "cds.GRCh38.seqs.tsv" "utr5.GRCh38.seqs.tsv" "utr3.GRCh38.seqs.tsv"; do	
	outfile="g4.${seqfile}"
	grep -E "[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}" $seqfile >> $outfile		
done
