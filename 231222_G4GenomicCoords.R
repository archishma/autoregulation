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

cds.coords[cds.coords$id == "ENST00000005340",]
## making some dummy data 
cds.g4 <- data.frame(id = c("a", "a", "a"), 
                     starts = c(2, 19, 23), 
                     stops = c(25, 22, 49)
                     )
cds.coords <- data.frame(id = c("a","a","a"), 
                         gen_starts = c(3001, 4218, 5724), 
                         gen_stops = c(3017, 4223, 5750), 
                         tx_starts = c(1, 18, 24), 
                         tx_stops = c(17, 23, 50))


g4.spans <- list()
for (i in seq(dim(cds.g4)[1])) {
  # get tx id and transcriptomic start/stop coords 
  tx <- cds.g4[i,]$id
  start <- cds.g4[i,]$starts
  stop <- cds.g4[i,]$stops 
  span <- list()
  # subset the transcript regions corresponding to current G4 
  table <- cds.coords[cds.coords$id == tx,]
  table.length <- dim(table)[1]
  for (j in seq(table.length)) {
    # if start coord is in a single interval:
    if (start >= table[j,]$tx_starts && start <= table[j,]$tx_stops) {
      if (stop <= table[j,]$tx_stops) {
        tx.start <- table[j,]$tx_starts 
        start.offset <- start - tx.start 
        stop.offset <- stop - tx.start 
        
        gen.start <- table[j,]$gen_starts + start.offset
        gen.stop <- table[j,]$gen_starts + stop.offset 
        
        cat(gen.start, "\n")
        cat(gen.stop, "\n")
        span <- append(span, c(gen.start, gen.stop))
        break 
      }
      else { # stop coord is in a different interval
        tx.start <- table[j,]$tx_starts 
        start.offset <- start - tx.start
        
        gen.start <- table[j,]$gen_starts + start.offset
        
        cat(gen.start, "\n")
        cat(table[j,]$gen_stops, "\n")
        span <- append(span, c(gen.start, table[j,]$gen_stops))
        for (k in seq(j+1,table.length)) {
          # check if the stop coordinate is in the next interval 
          if (stop >= table[k,]$tx_starts && stop <= table[k,]$tx_stops) {
            tx.start <- table[k,]$tx_starts
            stop.offset <- stop - tx.start 
            gen.stop <- table[k,]$gen_starts + stop.offset 
            cat(table[k,]$gen_starts, "\n")
            cat(gen.stop, "\n")
            span <- append(span, c(table[k,]$gen_starts,gen.stop))
            break
          }
          else { # the stop coordinate is not in the next interval 
            
            cat(table[k,]$gen_starts, "\n")
            cat(table[k,]$gen_stops, "\n")
            span <- append(span, c(table[k,]$gen_starts, table[k,]$gen_stops))
          }
        }
        # once the above loop finishes 
        # need to fix this part - not all coords are being added 
        g4.spans <- append(g4.spans, span)
        break
      }
    }
  }
}
