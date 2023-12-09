#!/bin/bash

#SBATCH -p general
#SBATCH -N 1 
#SBATCH -n 1 
#SBATCH -e g4_job.err
#SBATCH --mem=2g
#SBATCH -t 5:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=askav@email.unc.edu

REGIONS="/nas/longleaf/home/askav/rna/autoregulation/g4/genome.proteinCoding.subset.gtf"
GENOME="/proj/RNA_lab/Genomic_Indexes/Human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
OUTFILE="/nas/longleaf/home/askav/rna/autoregulation/g4/g4.proteinCodingRegions.tsv"

module load samtools 

while read -r line
do
	chr="`echo ${line} | cut -d " " -f 1| tr -d "chr"`"
	annot="`echo ${line} | cut -d " " -f 3`"
	pos_start="`echo ${line} | cut -d " " -f 4`"
	pos_end="`echo ${line} | cut -d " " -f 5`"
	strand="`echo ${line} | cut -d " " -f 6`"
	
	sequence=""
	# echo $chr $pos_start $pos_end $strand
	if [ $strand == "+" ]; then		
		sequence="`samtools faidx --mark-strand no ${GENOME} ${chr}:${pos_start}-${pos_end} | tr -d "\n"`"
	else # reverse strand 
		sequence="`samtools faidx -i --mark-strand no ${GENOME} ${chr}:${pos_start}-${pos_end} | tr -d "\n"`"
	fi
	
	g4_list=`echo "$sequence" | grep -E -o "[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}"`
	for g4 in ${g4_list}; do
		echo -e "$chr\t$pos_start\t$pos_end\t$strand\t$annot\t$g4" >> $OUTFILE	
	done	
done < $REGIONS
