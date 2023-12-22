# 231215
library(dplyr)
library(tidyr)
library(data.table)

cds.tx.exons <- read.table("~/rna/autoregulation/231215_G4_outputs/cds.tx.exons.tsv", sep = "\t")
cds.tx.g4 <- read.table("~/rna/autoregulation/231215_G4_outputs/cds.tx.g4.tsv", sep = "\t")
utr5.tx.exons <- read.table("~/rna/autoregulation/231215_G4_outputs/utr5.tx.exons.tsv", sep = "\t")
utr5.tx.g4 <- read.table("~/rna/autoregulation/231215_G4_outputs/utr5.tx.g4.tsv", sep = "\t")
utr3.tx.exons <- read.table("~/rna/autoregulation/231215_G4_outputs/utr3.tx.exons.tsv", sep = "\t")
utr3.tx.g4 <- read.table("~/rna/autoregulation/231215_G4_outputs/utr3.tx.g4.tsv", sep = "\t")

### CDS - calcuate transcript coords
cds.tx.exons$tx_starts <- 0 
cds.tx.exons$tx_stops <- 0 

# group exons by transcript id to loop over transcripts
cds <- cds.tx.exons %>% group_by(id) %>% group_split()

for (i in seq(1,length(cds))) { # accessing elements by index to modify inplace 
  # calculate "transcript coords" based on genomic coords 
  cds[[i]]$tx_starts[1] <- 1
  cds[[i]]$tx_stops[1] <- cds[[i]]$gen_stops[1] - cds[[i]]$gen_starts[1] + 1
  
  if (length(cds[[i]]$gen_starts) > 1) { # if tx is only one exon, don't run the following  
    for (j in seq(2,length(cds[[i]]$gen_starts))) {
      cds[[i]]$tx_starts[j] <- cds[[i]]$tx_stops[j-1] + 1
      cds[[i]]$tx_stops[j] <- 
        cds[[i]]$tx_starts[j] + cds[[i]]$gen_stops[j] - cds[[i]]$gen_starts[j]
    }
  }
}

cds.pos <- rbindlist(cds)

### 5' UTR - calculate transcript coords 
utr5.tx.exons$tx_starts <- 0 
utr5.tx.exons$tx_stops <- 0 

utr5 <- utr5.tx.exons %>% group_by(id) %>% group_split()

for (i in seq(1,length(utr5))) { # accessing elements by index to change in place 
  utr5[[i]]$tx_starts[1] <- 1
  utr5[[i]]$tx_stops[1] <- utr5[[i]]$gen_stops[1] - utr5[[i]]$gen_starts[1] + 1
  
  if (length(utr5[[i]]$gen_starts) > 1) { # if tx is only one exon, don't run the following  
    for (j in seq(2,length(utr5[[i]]$gen_starts))) {
      utr5[[i]]$tx_starts[j] <- utr5[[i]]$tx_stops[j-1] + 1
      utr5[[i]]$tx_stops[j] <- 
        utr5[[i]]$tx_starts[j] + utr5[[i]]$gen_stops[j] - utr5[[i]]$gen_starts[j]
    }
  }
}

utr5.pos <- rbindlist(utr5)

### 3' UTR - calculate transcript coords 
utr3.tx.exons$tx_starts <- 0 
utr3.tx.exons$tx_stops <- 0 

utr3 <- utr3.tx.exons %>% group_by(id) %>% group_split()

for (i in seq(1,length(utr3))) { # accessing elements by index to change in place 
  utr3[[i]]$tx_starts[1] <- 1
  utr3[[i]]$tx_stops[1] <- utr3[[i]]$gen_stops[1] - utr3[[i]]$gen_starts[1] + 1
  
  if (length(utr3[[i]]$gen_starts) > 1) { # if tx is only one exon, don't run the following  
    for (j in seq(2,length(utr3[[i]]$gen_starts))) {
      utr3[[i]]$tx_starts[j] <- utr3[[i]]$tx_stops[j-1] + 1
      utr3[[i]]$tx_stops[j] <- 
        utr3[[i]]$tx_starts[j] + utr3[[i]]$gen_stops[j] - utr3[[i]]$gen_starts[j]
    }
  }
}

utr3.pos <- rbindlist(utr3)

### Now, we get the G4s from each region

# CDS
cds.g4 <- cds.tx.g4$seq %>% 
  gregexpr(pattern = "[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}")
cds.g4.starts <- list() 
cds.g4.stops <- list()

for (i in seq(length(cds.g4))) { 
  # each i corresponds to a transcript 
  starts <- cds.g4[[i]]
  attributes(starts) <- NULL 
  lengths <- attr(cds.g4[[i]], 'match.length')
  stops <- starts + lengths - 1
  
  # input positions next to tx and seq 
  cds.g4.starts[[i]] <- starts
  cds.g4.stops[[i]] <- stops 
}

cds.tx.g4$starts <- cds.g4.starts
cds.tx.g4$stops <- cds.g4.stops

# 5' UTR
utr5.g4 <- utr5.tx.g4$seq %>% 
  gregexpr(pattern = "[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}")
utr5.g4.starts <- list()
utr5.g4.stops <- list()

for (i in seq(length(utr5.g4))) { 
  # each i corresponds to a transcript 
  starts <- utr5.g4[[i]]
  attributes(starts) <- NULL 
  lengths <- attr(utr5.g4[[i]], 'match.length')
  stops <- starts + lengths - 1
  
  
  # input positions next to tx and seq 
  utr5.g4.starts[[i]] <- starts
  utr5.g4.stops[[i]] <- stops 
}

utr5.tx.g4$starts <- utr5.g4.starts
utr5.tx.g4$stops <- utr5.g4.stops

# 3' UTR 
utr3.g4 <- utr3.tx.g4$seq %>% 
  gregexpr(pattern = "[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}[G]{3,5}.{1,7}")
utr3.g4.starts <- list()
utr3.g4.stops <- list()

for (i in seq(length(utr3.g4))) { 
  # each i corresponds to a transcript 
  starts <- utr3.g4[[i]]
  attributes(starts) <- NULL 
  lengths <- attr(utr3.g4[[i]], 'match.length')
  stops <- starts + lengths - 1
  
  # input positions next to tx and seq 
  utr3.g4.starts[[i]] <- starts
  utr3.g4.stops[[i]] <- stops 
}

utr3.tx.g4$starts <- utr3.g4.starts
utr3.tx.g4$stops <- utr3.g4.stops

# reformat the data 
cds.tx.g4 <- unnest(cds.tx.g4, cols = c(starts, stops))
utr5.tx.g4 <- unnest(utr5.tx.g4, cols = c(starts, stops))
utr3.tx.g4 <- unnest(utr3.tx.g4, cols = c(starts, stops))

### STATS 
# how many g4s 
cat("CDS G4s counted: ", dim(cds.tx.g4)[1], "\n") # 4908
cat("UTR5 G4s counted: ", dim(utr5.tx.g4)[1], "\n") # 6248
cat("UTR3 G4s counted: ", dim(utr3.tx.g4)[1], "\n") # 10574

### save the data 
# cds.tx.g4 %>% write.table("~/rna/autoregulation/231215_G4_outputs/cds.tx.g4.coords.tsv", 
#                           sep = "\t")
# utr5.tx.g4 %>% write.table("~/rna/autoregulation/231215_G4_outputs/utr5.tx.g4.coords.tsv", 
#                           sep = "\t")
# utr3.tx.g4 %>% write.table("~/rna/autoregulation/231215_G4_outputs/utr3.tx.g4.coords.tsv", 
#                            sep = "\t")
# cds.pos %>% write.table("~/rna/autoregulation/231215_G4_outputs/cds.tx.exons.coords.tsv",
#                         sep = "\t")
# utr5.pos %>% write.table("~/rna/autoregulation/231215_G4_outputs/utr5.tx.exons.coords.tsv",
#                         sep = "\t")
# utr3.pos %>% write.table("~/rna/autoregulation/231215_G4_outputs/utr3.tx.exons.coords.tsv",
#                          sep = "\t")
