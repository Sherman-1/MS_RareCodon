#! /usr/bin/Rscript

# This script takes as input a fasta file and output a csv file containing 
# the codon counts needed for further prediction

# Based on Paul Roginski's script


# In order to pass arguments with bash
args <- commandArgs(trailingOnly = TRUE)

file= args[1]

# Return the codon frequency table
get_codon_count <- function(fasta_path){
    
  library(coRdon)
  #library(coRdon, lib.loc="/home/paul.roginski/R/x86_64-pc-linux-gnu-library/3.5", verbose=TRUE)

  set <- readSet(file = fasta_path)
  codon_Table <- codonTable(set)
  codon_Counts <- codonCounts(codon_Table)
  row.names(codon_Counts) <- codon_Table@ID
    
  return(codon_Counts)
}
 

# For now we always need to get all the features :
output <- get_codon_count(fasta_path = file)

cc_file = sub(x = file, pattern = ".fna", replacement = "_CC.csv")

write.csv(output, cc_file)
