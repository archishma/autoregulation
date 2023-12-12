library(AnnotationHub)
library(BSgenome)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(dplyr)

### I modified the 231208 file to see if I can maintain the exon info
### in each transcript 

setwd("/nas/longleaf/home/askav/rna/autoregulation/g4/")
ah <- AnnotationHub()

# ensembl version 110 
# GRCh38.p14 
# release date: July 2023 
ensembl <- ah[["AH113665"]] 

# get EXONS from coding sequences, 5' UTR and 3' UTR
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

# all exons from every transcript in CDS, with corresponding sequences  
cds.df <- data.frame(annot = cds.annot, 
                     id = cds.id, 
                     seq = cds.seq,
                     exon_rank = cds$exon_rank,
                     exon_id = cds$exon_id)
cds.df.tx <- 
  cds.df %>% 
  group_by(id) %>% 
  mutate(tx_seq = paste0(seq, collapse = "")) %>%
  subset(select = c("id", "tx_seq")) %>%
  distinct()

write.table(cds.df.tx,
            file = "GRCh38.cds.txseqs.tsv",
            sep = "\t",
            row.names = FALSE)

# 5' UTR annotations and sequences 
utr5 <- unlist(utr5)
utr5.annot <- as.character(utr5)
utr5.id <- names(utr5)
utr5.seq <- as.character(getSeq(genome, names = utr5))

# all exons from every transcript in 5' UTR, with corresponding sequences  
utr5.df <- data.frame(annot = utr5.annot, 
                     id = utr5.id, 
                     seq = utr5.seq,
                     exon_rank = utr5$exon_rank,
                     exon_id = utr5$exon_id)

utr5.df.tx <- 
  utr5.df %>%
  group_by(id) %>%
  mutate(tx_seq = paste0(seq, collapse = "")) %>%
  subset(select = c("id", "tx_seq")) %>%
  distinct()

write.table(utr5.df.tx,
            file = "GRCh38.utr5.txseqs.tsv",
            sep = "\t",
            row.names = FALSE)

# 3' UTR annotations and sequences 
utr3 <- unlist(utr3)
utr3.annot <- as.character(utr3)
utr3.id <- names(utr3)
utr3.seq <- as.character(getSeq(genome, names = utr3))

# all exons from every transcript in 3' UTR, with corresponding sequences  
utr3.df <- data.frame(annot = utr3.annot, 
                      id = utr3.id, 
                      seq = utr3.seq,
                      exon_rank = utr3$exon_rank,
                      exon_id = utr3$exon_id)

utr3.df.tx <- 
  utr3.df %>%
  group_by(id) %>%
  mutate(tx_seq = paste0(seq, collapse = "")) %>%
  subset(select = c("id", "tx_seq")) %>%
  distinct()

write.table(utr3.df.tx,
            file = "GRCh38.utr3.txseqs.tsv",
            sep = "\t",
            row.names = FALSE) 
