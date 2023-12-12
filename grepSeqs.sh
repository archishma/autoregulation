# 231209
# the easier way to find G4s using regex, which can just be run at the command line 
# I used this command to create `*.GRCh38.g4.tsv` files in the g4/ directory
grep -E "[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}" $file.seqs.tsv > $file.g4.tsv

# 231212
# what I used to get all G4 containing transcripts from my previous work
# copy pasted from the history command 
 1009  grep -E "[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}" GRCh38.cds.txseqs.tsv > GRCh38.cds.txseqs.g4.tsv
 1010  grep -E "[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}" GRCh38.utr5.txseqs.tsv > GRCh38.utr5.txseqs.g4.tsv
 1011  grep -E "[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}" GRCh38.utr3.txseqs.tsv > GRCh38.utr3.txseqs.g4.tsv
