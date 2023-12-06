library(AnnotationHub)

ah <- AnnotationHub()
# human <- subset(ah, species == "Homo sapiens") # learning how AH works 
ensembl <- ah[["AH113665"]]
organism(ensembl) # output: "Homo sapiens"
proteinCoding <- exons(ensembl, filter = GeneBiotypeFilter("protein_coding"))
proteinCoding
listTxbiotypes(ensembl)
listGenebiotypes(ensembl)
