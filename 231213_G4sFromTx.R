### started analysis 231213
# I have all transcripts that contain at least 1 G4 
# I also have start/stop coords for all exons comprising transcripts
library(dplyr)

# load data 
setwd("~/rna/autoregulation/g4/")
cds.exons <- read.table("GRCh38.cds.exons.tsv", 
                        header = TRUE)
cds.tx.g4 <- read.table("GRCh38.cds.txseqs.g4.tsv", 
                        col.names = c("id", "seq"))

utr5.exons <- read.table("GRCh38.utr5.exons.tsv", 
                         header = TRUE)
utr5.tx.g4 <- read.table("GRCh38.utr5.txseqs.g4.tsv", 
                         col.names = c("id", "seq"))

utr3.exons <- read.table("GRCh38.utr3.exons.tsv", 
                         header = TRUE)
utr3.tx.g4 <- read.table("GRCh38.utr3.txseqs.g4.tsv", 
                         col.names = c("id",  "seq"))

# include exons corresponding to G4-containing transcripts 
cds.g4.exons <- cds.exons %>% subset(id %in% cds.tx.g4$id)
utr5.g4.exons <- utr5.exons %>% subset(id %in% utr5.tx.g4$id)
utr3.g4.exons <- utr3.exons %>% subset(id %in% utr3.tx.g4$id)

# get rid of objs we won't need anymore 
rm(cds.exons)
rm(utr5.exons)
rm(utr3.exons)

# trying out grep in R on "ENST00000005340"
example <- cds.tx.g4[cds.tx.g4$id == "ENST00000005340",]
# I know (visually) that 15 exons comprise this transcript 
seq <- example$seq
match <- seq %>% gregexpr(pattern = "[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}")
match
poss <- seq %>% substr(start = 1795, stop = 1795 + 28 - 1)
nchar(poss)
poss

# where is it in the exons 
