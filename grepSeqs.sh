# 231209
# the easier way to find G4s using regex, which can just be run at the command line 
# I used this command to create `*.GRCh38.g4.tsv` files in the g4/ directory
grep -E "[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}" $file.seqs.tsv > $file.g4.tsv
