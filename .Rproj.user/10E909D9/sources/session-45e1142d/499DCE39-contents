#Alineamiento 
library (Biostrings)
library (msa)

fasta_seq <- readAAStringSet ("/Users/julie/OneDrive/Escritorio/R/Bioinformática/rep_nov04/data/DivergentGlobins.fasta") 
fasta_seq

globins <- msa (fasta_seq, method = "Muscle") 
globins


#Matriz de distancias con seqinr
library (seqinr)
globins2 <- msaConvert (globins, type = "seqinr::alignment")
globins2

dist_m <- dist.alignment (globins2)

#Árbol filogenético 
install.packages("ape")
library (ape)  
globinstree <- nj (dist_m)
print (globinstree)

#Mandar a results
pdf ("results/globinstree.pdf")
plot (globinstree)
dev.off ()



