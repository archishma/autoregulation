library(dplyr)
library(tidyr)
library(data.table)

# 231222
# i moved to another file because the other file got clunky 

# load the data I formatted previously
cds.coords <- read.table("~/rna/autoregulation/231215_G4_outputs/cds.tx.exons.coords.tsv")
utr5.coords <- read.table("~/rna/autoregulation/231215_G4_outputs/utr5.tx.exons.coords.tsv")
utr3.coords <- read.table("~/rna/autoregulation/231215_G4_outputs/utr3.tx.exons.coords.tsv")

# G4s 
cds.g4 <- read.table("~/rna/autoregulation/231215_G4_outputs/cds.tx.g4.coords.tsv")
utr5.g4 <- read.table("~/rna/autoregulation/231215_G4_outputs/utr5.tx.g4.coords.tsv")
utr3.g4 <- read.table("~/rna/autoregulation/231215_G4_outputs/utr3.tx.g4.coords.tsv")

# cds.coords[cds.coords$id == "ENST00000005340",]

## making some dummy data for testing 
test.g4 <- data.frame(id = c("a", "a", "a", "b", "b", "b"), 
                     starts = c(2, 19, 23, 1, 32, 1), 
                     stops = c(25, 22, 49, 28, 33, 100)
                     )
# more accurate to what the data actually looks like 
test.g4 <- data.frame(id = c("a", "b"), 
                      starts = c(2, 1), 
                      stops = c(25, 100))
test.coords <- data.frame(id = c("a","a","a", "b", "b", "b"), 
                         gen_starts = c(3001, 4218, 5724, 6001, 7121, 7274), 
                         gen_stops = c(3017, 4223, 5750, 6020, 7173, 7300), 
                         tx_starts = c(1, 18, 24, 1, 21, 74), 
                         tx_stops = c(17, 23, 50, 20, 73, 100), 
                         chr = c(19, 19, 19, 15, 15, 15), 
                         strand = c("+", "+", "+", "-", "-", "-"))

findGenomicCoords <- function(region.g4, region.coords) {
  g4.spans <- vector("list", length = dim(region.g4)[1])
  g4.chr <- vector("character", length = dim(region.g4)[1])
  g4.strand <- vector("character", length = dim(region.g4)[1])
  
  for (i in seq(dim(region.g4)[1])) {
    # get tx id and transcriptomic start/stop coords 
    tx <- region.g4[i,]$id
    start <- region.g4[i,]$starts
    stop <- region.g4[i,]$stops 
    
    # store chromosome and strand info for this G4 
    chr <- region.coords[region.coords$id == tx,]$chr[1]
    strand <- region.coords[region.coords$id == tx,]$strand[1]
    g4.chr[i] <- chr
    g4.strand[i] <- strand
    
    # subset the transcript regions corresponding to current G4 
    table <- region.coords[region.coords$id == tx,]
    table.length <- dim(table)[1]
    
    for (j in seq(table.length)) {
      # start coord is in the current interval
      if (start >= table[j,]$tx_starts && start <= table[j,]$tx_stops) {
        # stop coord is within the same interval
        if (stop <= table[j,]$tx_stops) {
          tx.start <- table[j,]$tx_starts 
          start.offset <- start - tx.start 
          stop.offset <- stop - tx.start 
          
          gen.start <- table[j,]$gen_starts + start.offset
          gen.stop <- table[j,]$gen_starts + stop.offset 
          
          # cat(gen.start, "\n")
          # cat(gen.stop, "\n")
          g4.spans[[i]] <- append(g4.spans[[i]], c(gen.start, gen.stop))
          break 
        }
        # stop coord is in a different interval
        else { 
          tx.start <- table[j,]$tx_starts 
          start.offset <- start - tx.start
          
          gen.start <- table[j,]$gen_starts + start.offset
          
          # cat(gen.start, "\n")
          # cat(table[j,]$gen_stops, "\n")
          g4.spans[[i]] <- append(g4.spans[[i]], c(gen.start, table[j,]$gen_stops))
          for (k in seq(j+1,table.length)) {
            # check if the stop coordinate is in the next interval 
            if (stop >= table[k,]$tx_starts && stop <= table[k,]$tx_stops) {
              tx.start <- table[k,]$tx_starts
              stop.offset <- stop - tx.start 
              gen.stop <- table[k,]$gen_starts + stop.offset 
              # cat(table[k,]$gen_starts, "\n")
              # cat(gen.stop, "\n")
              g4.spans[[i]] <- append(g4.spans[[i]], c(table[k,]$gen_starts,gen.stop))
              break
            }
            else { # the stop coordinate is not in the next interval 
              
              # cat(table[k,]$gen_starts, "\n")
              # cat(table[k,]$gen_stops, "\n")
              g4.spans[[i]] <- append(g4.spans[[i]], c(table[k,]$gen_starts, table[k,]$gen_stops))
            }
          }
          # once the above loop finishes 
          # g4.spans <- append(g4.spans, span)
        }
      }
    }
  }
  out <- list("spans" = g4.spans, "chr" = g4.chr, "strand" = g4.strand)
  return(out)
}

### test the function on dummy data 
example.list <- findGenomicCoords(test.g4, test.coords)
test.g4$genomic_coords <- example.list$spans
test.g4$chr <- example.list$chr 
test.g4$strand <- example.list$strand

### test the function on the first 3 transcripts from CDS 
cds.subset <- cds.g4[1:3,]
cds.subset.example <- findGenomicCoords(cds.subset, cds.coords)
# cds.subset.example$spans
# [[1]]
# [1] 7225897 7225924
# 
# [[2]]
# [1] 45178186 45178218 45165127 45165135
# 
# [[3]]
# [1] 50550089 50550115 50550208 50550216

### 040124
### Run the function on CDS, 5' UTR, and 3' UTR G4 data 

cds.output <- findGenomicCoords(cds.g4, cds.coords)
cds.g4$genomic_coords <- cds.output$spans
cds.g4$chr <- cds.output$chr
cds.g4$strand <- cds.output$strand 

cds.g4 %>% write.table(file = "~/rna/autoregulation/040124_G4GenomicCoords/cds.g4.gencoords.tsv", 
                       sep = "\t")

utr5.output <- findGenomicCoords(utr5.g4, utr5.g4)
utr5.g4$genomic_coords <- utr5.output$spans 
utr5.g4$chr <- utr5.output$chr
utr5.g4$strand <- utr5.output$strand 

utr3.output <- findGenomicCoords(utr3.g4, utr3.g4)
utr3.g4$genomic_coords <- utr3.output$spans 
utr3.g4$chr <- utr3.output$chr
utr3.g4$strand <- utr3.output$strand 
