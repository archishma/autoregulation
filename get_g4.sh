#/bin/bash
REGIONS="/nas/longleaf/home/askav/rna/autoregulation/g4/genome.proteinCoding.subset.gtf"
GENOME="/proj/RNA_lab/Genomic_Indexes/Human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

module load samtools 

while read -r line
do
	chr="`echo ${line} | cut -d " " -f 1| tr -d "chr"`"
	pos_start="`echo ${line} | cut -d " " -f 4`"
	pos_end="`echo ${line} | cut -d " " -f 5`"
	strand="`echo ${line} | cut -d " " -f 6`"
	
	sequence=""
	echo $chr $pos_start $pos_end $strand
	if [ $strand == "+" ]; then		
		sequence="`samtools faidx --mark-strand no ${GENOME} ${chr}:${pos_start}-${pos_end} | tr -d "\n"`"
	else # reverse strand 
		sequence="`samtools faidx -i --mark-strand no ${GENOME} ${chr}:${pos_start}-${pos_end} | tr -d "\n"`"
	fi
	echo $sequence
	echo "$sequence" | grep -E -o "[T]{3,5}.AAT"
	
	break # just to test on the first line 
done < $REGIONS
