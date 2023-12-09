library(AnnotationHub)
library(BSgenome)
### if GRCh38 not installed, uncomment the following
# library(BiocManager)
# install("BSgenome.Hsapiens.NCBI.GRCh38")
library(BSgenome.Hsapiens.NCBI.GRCh38)

setwd("~/rna/autoregulation/")
ah <- AnnotationHub()

# ensembl version 110 
# GRCh38.p14 
# release date: July 2023 
ensembl <- ah[["AH113665"]] 

### reference documentation 
### https://rdrr.io/bioc/ensembldb/man/EnsDb-exonsBy.html

organism(ensembl) # output: "Homo sapiens"

cds <- cdsBy(ensembl, by = c("tx", "gene"), columns = NULL,
             filter = AnnotationFilterList(), use.names = FALSE)
# cds <- cdsBy(ensembl, by = c("gene"), columns = NULL, 
#       filter = AnnotationFilterList(), use.names = FALSE)
utr5 <- fiveUTRsByTranscript(ensembl, columns = NULL, 
                             filter = AnnotationFilterList())
utr3 <- threeUTRsByTranscript(ensembl, columns = NULL, 
                              filter = AnnotationFilterList())
proteinCoding <- genes(ensembl, filter = GeneBiotypeFilter("protein_coding"))

length(cds) 
length(utr5)
length(utr3)
length(proteinCoding) # there are about 20,000 protein coding genes 

genome <- getBSgenome(genome = "BSgenome.Hsapiens.NCBI.GRCh38")
metadata(genome)
# $organism
# [1] "Homo sapiens"
# $common_name
# [1] "Human"
# $provider
# [1] "NCBI"
# $genome
# [1] "GRCh38"
# $release_date
# [1] "2013-12-17"
# $source_url
# [1] "ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/"

genome$"1" # get the sequence of the first chromosome 
# expand(cds)

### keeping the data as a df is easier to look at 
### however, we lose the convenience of GRanges objects 
# cds.df <- as.data.frame(cds) # much easier to look at 
# utr5.df <- as.data.frame(utr5)
# utr3.df <- as.data.frame(utr3)

cds <- keepStandardChromosomes(cds, pruning.mode = "coarse")
utr5 <- keepStandardChromosomes(utr5, pruning.mode = "coarse")
utr3 <- keepStandardChromosomes(utr3, pruning.mode = "coarse")

# we now only have the relevant chromosomes 
length(cds) 
length(utr5)
length(utr3)

ex <- cds[[1]] # the first GRanges object in the cds GRangesList 
ex
seq <- getSeq(genome, names = ex) # returns the sequences based on
seq

### 231208 - working with granges objects 
### https://seandavi.github.io/ITR/RangesAndSignal.html
cds <- unlist(cds) # cds is now a GRanges object 
cds[1:10]
as.character(cds[1:10])
as.vector(strand(cds))[1:10]
ranges(cds, use.names = TRUE)[1:10]
names(cds)[1:10]

as.data.frame(getSeq(genome, names = cds[1:10]))
### function to get sequence info from GRanges object with relevant info
### we want: gene name, chr, start, stop, strand

# lapply(grl) is slow: https://support.bioconductor.org/p/77300/

### idea:
# get all seqs and put into a df 
# https://bioinformatics.stackexchange.com/questions/5287/how-to-transform-a-dnastringset-from-the-bioconductor-package-biostrings-to-a-da
cds.head.df <- data.frame(annot=as.character(cds[1:10]), names = names(cds[1:10]), seq=as.character(getSeq(genome, names= cds[1:10])))
cds.head.df
# write.table(cds.head.df, "cds_head.tsv", sep = "\t", row.names = FALSE)

unlist(utr3) # checking output data table dimensions
cds
unlist(utr5) # checking output data table dimensions 
