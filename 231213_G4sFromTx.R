### started analysis 231213
# I have all transcripts that contain at least 1 G4 
# I also have start/stop coords for all exons comprising transcripts
library(dplyr)
library(tidyr)

# load data 
setwd("~/rna/autoregulation/g4/")
cds.exons <- read.table("GRCh38.cds.exons.tsv", 
                        header = TRUE)
cds.tx.g4 <- read.table("GRCh38.cds.txseqs.g4.tsv", 
                        col.names = c("id", "seq"))
sum(is.na(cds.exons)) # 0
sum(is.na(cds.tx.g4)) # 0

utr5.exons <- read.table("GRCh38.utr5.exons.tsv", 
                         header = TRUE)
utr5.tx.g4 <- read.table("GRCh38.utr5.txseqs.g4.tsv", 
                         col.names = c("id", "seq"))
sum(is.na(utr5.exons)) # 0 
sum(is.na(utr5.tx.g4)) # 0 

utr3.exons <- read.table("GRCh38.utr3.exons.tsv", 
                         header = TRUE)
utr3.tx.g4 <- read.table("GRCh38.utr3.txseqs.g4.tsv", 
                         col.names = c("id",  "seq"))
sum(is.na(utr3.exons)) # 0 
sum(is.na(utr3.tx.g4)) # 0

# include exons corresponding to G4-containing transcripts 
cds.g4.exons <- cds.exons %>% subset(id %in% cds.tx.g4$id)
utr5.g4.exons <- utr5.exons %>% subset(id %in% utr5.tx.g4$id)
utr3.g4.exons <- utr3.exons %>% subset(id %in% utr3.tx.g4$id)

sum(is.na(cds.g4.exons)) # 0 
sum(is.na(utr5.g4.exons)) # 0 
sum(is.na(utr3.g4.exons)) # 0 

# function to split the ranges and account for 1-length sequences 
# TODO: optimize later... im using grepl
split.range <- function(string, num = c(1,2)) {
  if (grepl("-", string)) {
    return(strsplit(string, split = "-")[[1]][num])
  }
  else {
    return(string)
  }
}

### split up chr, start_pos, stop_pos, strand 
# CDS 
cds.tx.exons <- cds.g4.exons %>% 
  separate(col = "annot", into = c("chr", "range", "strand"), sep = ":") 
sum(is.na(cds.tx.exons)) # 0

cds.tx.exons$gen_starts <- cds.tx.exons$range %>% lapply(split.range, num = 1)
cds.tx.exons$gen_stops <- cds.tx.exons$range %>% lapply(split.range, num = 2)
sum(is.na(cds.tx.exons)) # 0. this is not the behavior I get with separate()

# 5' UTR
utr5.tx.exons <- utr5.g4.exons %>%  
  separate(col = "annot", into = c("chr", "range", "strand"), sep = ":")
utr5.tx.exons$gen_starts <- utr5.tx.exons$range %>% lapply(split.range, num = 1)
utr5.tx.exons$gen_stops <- utr5.tx.exons$range %>% lapply(split.range, num = 2)
sum(is.na(utr5.tx.exons)) # 0 

# 3' UTR
utr3.tx.exons <- utr3.g4.exons %>% 
  separate(col = "annot", into = c("chr", "range", "strand"), sep = ":")
utr3.tx.exons$gen_starts <- utr3.tx.exons$range %>% lapply(split.range, num = 1)
utr3.tx.exons$gen_stops <- utr3.tx.exons$range %>% lapply(split.range, num = 2)

# get rid of objs we won't need anymore 
rm(cds.exons, cds.g4.exons)
rm(utr5.exons, utr5.g4.exons)
rm(utr3.exons, utr3.g4.exons)

### make sure that didn't give me NA values 

### TESTING
### trying out grep in R on "ENST00000005340"
example <- cds.tx.g4[cds.tx.g4$id == "ENST00000005340",]
# I know (visually) that 15 exons comprise this transcript 
seq <- example$seq
match <- seq %>% gregexpr(pattern = "[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}")
# https://stackoverflow.com/questions/51875360/accessing-results-of-gregexpr
# this works because there is only one text element: seq
# but I could pass in a list and iterate over match[[1]] ... match[[n]]
starts <- match[[1]]
attributes(starts) <- NULL
lengths <- attr(match[[1]], 'match.length')
starts
lengths

s1 <- 0
s2 <- 0
for (i in length(starts)) {
  s1 <- starts[i]
  s2 <- starts[i] + lengths[i] - 1
  poss <- seq %>% substr(start = starts[i], stop = starts[i] + lengths[i] - 1)
  print(nchar(poss))
  print(poss)
}

# where is it in the transcript. 1-based indexing 
s1
s2
# get exons corresponding to transcript 
exons.5340 <- cds.tx.exons[cds.tx.exons$id == "ENST00000005340",]
exons.5340[,c("gen_starts", "gen_stops")]

# example (naive) method for a single transcript 
# I'm not yet sure of a faster way to do this. 
# make sure chrom start, stop columns are numeric 
exons.5340[,c("gen_starts","gen_stops")] <- sapply(exons.5340[,c("gen_starts","gen_stops")], as.numeric)

# base case 
exons.5340$tx_start <- NA
exons.5340$tx_stop <- NA
exons.5340$tx_start[1] <- 1
exons.5340$tx_stop[1] <- exons.5340$gen_stops[1] - exons.5340$gen_starts[1] + 1
# loop over the exons in the tx 
for (i in seq(2,length(exons.5340$gen_starts))) {
  exons.5340$tx_start[i] <- exons.5340$tx_stop[i-1] + 1
  exons.5340$tx_stop[i] <- exons.5340$tx_start[i] + exons.5340$gen_stops[i] - exons.5340$gen_starts[i]
}
exons.5340$tx_stop[length(exons.5340$tx_stop)] # 2211
nchar(seq) #2211
# my method works 
