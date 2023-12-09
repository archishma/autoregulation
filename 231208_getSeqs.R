library(AnnotationHub)
library(BSgenome)
library(BSgenome.Hsapiens.NCBI.GRCh38)

setwd("/nas/longleaf/home/askav/rna/autoregulation/g4/")
ah <- AnnotationHub()

# ensembl version 110 
# GRCh38.p14 
# release date: July 2023 
ensembl <- ah[["AH113665"]] 

# get coding sequences, 5' UTR and 3' UTR
cds <- cdsBy(ensembl, by = c("tx", "gene"), columns = NULL,
             filter = AnnotationFilterList(), use.names = FALSE)
utr5 <- fiveUTRsByTranscript(ensembl, columns = NULL, 
                             filter = AnnotationFilterList())
utr3 <- threeUTRsByTranscript(ensembl, columns = NULL, 
                              filter = AnnotationFilterList())

# load genome 
genome <- getBSgenome(genome = "BSgenome.Hsapiens.NCBI.GRCh38")

# keep only standard chromosomes 
cds <- keepStandardChromosomes(cds, pruning.mode = "coarse")
utr5 <- keepStandardChromosomes(utr5, pruning.mode = "coarse")
utr3 <- keepStandardChromosomes(utr3, pruning.mode = "coarse")

# CDS annotations and sequences 
cds <- unlist(cds) # convert to GRanges from GRangesList 
cds.annot <- as.character(cds)
cds.id <- names(cds)
cds.seq <- as.character(getSeq(genome, names = cds))

cds.df <- data.frame(annot = cds.annot, id = cds.id, seq = cds.seq)
write.table(cds.df, 
            file = "GRCh38.cds.seqs.tsv", 
            sep = "\t", 
            row.names = FALSE)

# 5' UTR annotations and sequences 
utr5 <- unlist(utr5)
utr5.annot <- as.character(utr5)
utr5.id <- names(utr5)
utr5.seq <- as.character(getSeq(genome, names = utr5))

utr5.df <- data.frame(annot = utr5.annot, id = utr5.id, seq = utr5.seq)
write.table(utr5.df, 
            file = "GRCh38.utr5.seqs.tsv", 
            sep = "\t", 
            row.names = FALSE)

# 3' UTR annotations and sequences 
utr3 <- unlist(utr3)
utr3.annot <- as.character(utr3)
utr3.id <- names(utr3)
utr3.seq <- as.character(getSeq(genome, names = utr3))

utr3.df <- data.frame(annot = utr3.annot, id = utr3.id, seq = utr3.seq)
write.table(utr3.df, 
            file = "GRCh38.utr3.seqs.tsv", 
            sep = "\t", 
            row.names = FALSE)