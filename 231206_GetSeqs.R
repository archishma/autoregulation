library(AnnotationHub)
library(BSgenome)
### if GRCh38 not installed, uncomment the following
# library(BiocManager)
# install("BSgenome.Hsapiens.NCBI.GRCh38")

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

ex <- cds[[1]]
seq <- getSeq(genome, names = ex)
names(seq) <- ex$exon_id
# expand(cds)

cds.df <- as.data.frame(cds) # much easier to look at 
utr5.df <- as.data.frame(utr5)
utr3.df <- as.data.frame(utr3)
