require(Biostrings)
require(dplyr)
require(ComplexHeatmap)
require(circlize)

upstream = read.csv("~/Bureau/Gits/MS_RareCodon/output/upstream_codons.csv", check.names = FALSE)
downstream = read.csv("~/Bureau/Gits/MS_RareCodon/output/downstream_codons.csv", check.names = FALSE)



normalize_codon_frequencies <- function(glob_freq, data) {
  
  require(dplyr)
  bases <- c("A", "T", "G", "C")
  codons <- as.vector(outer(outer(bases, bases, paste0), bases, paste0))  
  codons <- sort(codons)
  
  if (!identical(sort(row.names(glob_freq)), codons)) {
    stop("The row names of df1 should exactly match the 64 codons and be sorted alphabetically.")
  }

  # Validate and set row names for df2
  if (!identical(sort(row.names(data)), codons)) {
    stop("The row names of df2 should exactly match the 64 codons and be sorted alphabetically.")
  }

  merged_df <- merge(data, glob_freq, by = "row.names")
  merged_df <- merged_df[,-ncol(merged_df)]
  output <- merged_df %>%
      mutate_at(vars(-freq, -Row.names), list(~ ./freq))
  
  rownames(output) = output[,"Row.names"]
  output = output[,c(-1,-ncol(output))]
  
  

  
  return(output)
  
}

## GENERAL

codons <- c("AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
            "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
            "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
            "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT")

kyte_doolittle <- c('A' = 1.8, 'R' = -4.5, 'N' = -3.5, 'D' = -3.5, 'C' = 2.5, 
                    'Q' = -3.5, 'E' = -3.5, 'G' = -0.4, 'H' = -3.2, 'I' = 4.5, 
                    'L' = 3.8, 'K' = -3.9, 'M' = 1.9, 'F' = 2.8, 'P' = -1.6, 
                    'S' = -0.8, 'T' = -0.7, 'W' = -0.9, 'Y' = -1.3, 'V' = 4.2)

codon_freq = read.csv("~/Bureau/Gits/MS_RareCodon/input/codon_freq.csv",
                      row.names = 1,
                      col.names = "freq",
)

codon_freq <- codon_freq %>% mutate(Codon = row.names(.))


raw_up = upstream[,-ncol(upstream)]
raw_down = downstream[,-ncol(downstream)]

######

# DATA PREP UPSTREAM 

######

res_up <- data.frame(Column1 = matrix(0, nrow = 64, ncol = 1), 
                     Column2 = matrix(0, nrow = 64, ncol = 1), 
                     Column3 = matrix(0, nrow = 64, ncol = 1))
rownames(res_up) = codons

colnames(res_up) = colnames(raw_up)

######

# DATA PREP DOWNSTREAM 

######

res_down <- data.frame(Column1 = matrix(0, nrow = 64, ncol = 1), 
                       Column2 = matrix(0, nrow = 64, ncol = 1), 
                       Column3 = matrix(0, nrow = 64, ncol = 1))
rownames(res_down) = codons
colnames(res_down) = colnames(raw_down)


### CODONS SORTING HYDROP.
amino_acids <- sapply(row.names(res_up), function(codon) {
  aa <- GENETIC_CODE[[as.character(DNAString(codon))]]
  return(aa)
})

hydrophobicity <- sapply(amino_acids, function(aa) {
  return(kyte_doolittle[aa])
})


##########

## RESULTS UPSTREAM CODONS

#########

for (col in names(raw_up)) {
  
  counts <- table(raw_up[[col]])
  
  frequencies <- counts/sum(counts)
  
  res_up[[col]] <- frequencies[codons]
  
}

# Sort codons by hydrophobicity
res_up <- res_up[order(hydrophobicity), ]

##########

## RESULTS UPSTREAM CODONS

#########

for (col in names(raw_down)) {
  
  counts <- table(raw_down[[col]])
  
  frequencies <- counts/sum(counts)
  
  res_down[[col]] <- frequencies[codons]
  
}

# Sort codons by hydrophobicity
res_down <- res_down[order(hydrophobicity), ]


##### PLOT HEATMAPS

## UPSTREAM
mat <- as.matrix(res_up)

upper_bound <- max(mat)

Heatmap(mat, 
        name = "Codon Frequencies", 
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = colorRamp2(breaks = c(0, upper_bound/2, upper_bound), colors = c("blue", "white","red")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_gp = 
        column_title = "Upstream codons")
        

## DOWNSTREAM

mat <- as.matrix(res_down)

upper_bound <- max(mat)

Heatmap(mat, 
        name = "Codon Frequencies", 
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = colorRamp2(breaks = c(0, upper_bound/2, upper_bound), colors = c("blue", "white","red")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_title = "Downstream codons")



## PULLING RESULTS 

# DOWNSTREAM

vector_down <- as.vector(unlist(raw_down))
frequencies <- table(vector_down) / length(vector_down)
codon_freq <- data.frame("Codon" = names(frequencies), "Frequency" = as.vector(frequencies))
codon_freq <- codon_freq[order(hydrophobicity), ]

mat <- matrix(codon_freq$Frequency, nrow=length(codon_freq$Frequency), ncol=1, 
              dimnames = list(codon_freq$Codon, "Frequency"))

upper_bound <- max(mat)

Heatmap(mat, 
        name = "Codon Frequencies", 
        show_column_names = FALSE,
        show_row_names = FALSE,
        col = colorRamp2(breaks = c(0, upper_bound/2, upper_bound), colors = c("blue", "white","red")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_title = "Pulled downstream codons")


vector_up <- as.vector(unlist(raw_up))
frequencies <- table(vector_up) / length(vector_up)
codon_freq <- data.frame("Codon" = names(frequencies), "Frequency" = as.vector(frequencies))
codon_freq <- codon_freq[order(hydrophobicity), ]

mat <- matrix(codon_freq$Frequency, nrow=length(codon_freq$Frequency), ncol=1, 
              dimnames = list(codon_freq$Codon, "Frequency"))

upper_bound <- max(mat)

Heatmap(mat, 
        name = "Codon Frequencies", 
        show_column_names = FALSE,
        show_row_names = TRUE,
        col = colorRamp2(breaks = c(0, upper_bound/2, upper_bound), colors = c("blue", "white","red")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_title = "Pulled upstream codons")

###################

### NORMALISATION

###################


norm_down <- normalize_codon_frequencies(codon_freq, res_down)
norm_up <- normalize_codon_frequencies(codon_freq, res_up)


## HEATMAP DOWNSTREAM
mat <- as.matrix(norm_down)

upper_bound <- max(mat)

Heatmap(mat, 
        name = "Codon Frequencies", 
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = colorRamp2(breaks = c(0, upper_bound/2, upper_bound), colors = c("blue", "white","red")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_gp = grid::gpar(fontsize = 8),
        column_title = "Downstream codons normalized")

## HEATMAP UPSTREAM
mat <- as.matrix(norm_up)

upper_bound <- max(mat)

Heatmap(mat, 
        name = "Codon Frequencies", 
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = colorRamp2(breaks = c(0, upper_bound/2, upper_bound), colors = c("blue", "white","red")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_gp = grid::gpar(fontsize = 8),
        column_title = "Upstream codons normalized")



