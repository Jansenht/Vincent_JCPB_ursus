

# read in raw counts, skipping first row that just contains info on how featureCounts was run

raw_counts_final = read.delim("/Users/elleryvincent/Desktop/Gene Level Round 3/featureCountsgenelevel_output.txt",  skip = 1) 

# Create matrix containing only gene count columns

counts.mat.final = as.matrix(raw_counts_final[,-c(1,2,3,4,5,6)])

# set gene IDs as row names

row.names(counts.mat.final) = raw_counts_final$Geneid 

# Calculate TPM normalized counts

x = counts.mat.final / raw_counts_final$Length              # Normalize counts by gene length
tpm.final = t( t(x) * 1e6 / colSums(x) )         # Convert to transcripts per million (TPM)