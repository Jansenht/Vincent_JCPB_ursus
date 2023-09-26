library(tidyverse)

# read in all counts
counts <- gene_level_norm_counts_for_WGCNA

# separate by treatment
a34.counts <- counts %>% select(1,contains('A34'))
a37.counts <- counts %>% select(1,contains('A37'))
h34.counts <- counts %>% select(1,contains('H34'))
h37.counts <- counts %>% select(1,contains('H37'))

# define minimum count threshold
min_count <- 1

# calculates the mean of all columns except for the first one (gene ids), filters out rows that don't meet minimum count
a34.counts.filtered <- a34.counts %>% filter(rowMeans(.[,-1]) > min_count) 
write.csv(a34.counts.filtered, "./A34_low_filtered.csv")
a37.counts.filtered <- a37.counts %>% filter(rowMeans(.[,-1]) > min_count) 
write.csv(a37.counts.filtered, "./A37_low_filtered.csv")
h34.counts.filtered <- h34.counts %>% filter(rowMeans(.[,-1]) > min_count) 
write.csv(h34.counts.filtered, "./H34_low_filtered.csv")
h37.counts.filtered <- h37.counts %>% filter(rowMeans(.[,-1]) > min_count) 
write.csv(h37.counts.filtered, "./H37_low_filtered.csv")

# see number of removed genes per treatment
nrow(a34.counts) - nrow(a34.counts.filtered) # 13,689 removed
nrow(a37.counts) - nrow(a37.counts.filtered) # 13,026 removed
nrow(h34.counts) - nrow(h34.counts.filtered) # 13,830 removed
nrow(h37.counts) - nrow(h37.counts.filtered) # 13,716 removed

# see percent removed genes per treatment
nrow(a34.counts.filtered) / nrow(a34.counts) # 47% removed
nrow(a37.counts.filtered) / nrow(a37.counts) # 49% removed
nrow(h34.counts.filtered) / nrow(h34.counts) # 47% removed
nrow(h37.counts.filtered) / nrow(h37.counts) # 47% removed

