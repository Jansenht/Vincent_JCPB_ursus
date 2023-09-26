# read in normalized expression values for condition
Active37_meta_input <- read_csv("Desktop/Gene Level Round 3/Active37_meta_input.csv")
  dataset <- read.table("Desktop/Gene Level Round 3/Active37_meta_input.csv", sep="\t", header = TRUE, row.names =  NULL, dec = ",")

  # transpose matrix
 trans_a37 = t(Active37_meta_input)
 
 #make 1st row col names
library(janitor)
 a37_complete = row_to_names(trans_a37, 1)


# convert the data set to numeric for WGCNA
dataset_numeric = apply(a37_complete, 2, as.numeric)
rownames(dataset_numeric) = rownames(a37_complete)
colMeans(dataset_numeric)


**************************************************************************** # ROUND 2

# read in normalized expression values for all conditions
normalized_expressions = gene_level_norm_counts_for_WGCNA_2DAYS

# transpose matrix
 trans_ne = t(normalized_expressions)
 
 # make 1st row col names
library(janitor)
 complete_ne = row_to_names(trans_ne, 1)


# convert the data set to numeric for WGCNA
dataset_numeric = apply(complete_ne, 2, as.numeric)
rownames(dataset_numeric) = rownames(complete_ne)
colMeans(dataset_numeric)

