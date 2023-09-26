
library(tidyverse)
library(edgeR)

# read in normalized expression values for all conditions
normalized_expressions = read_csv('~/Desktop/PCA_time.csv') %>% 
  as.data.frame() %>% 
  column_to_rownames('...1') # move first column with gene IDs to rownames %>% 

length(unique(rownames(normalized_expressions)))

# Define group/treatment for each sample
group <- str_split_fixed(names(normalized_expressions),'_',2)[,1] # make a list of treatments by splitting each column name at the underscore and keeping the first half (i.e., A34)
group <- factor(group) # convert to factor format

# Build a DGElist data object using count table and group information
y <- DGEList(counts=normalized_expressions,group=group)
y <- calcNormFactors(y) 

# Run plotMDS to run MDS analysis and save to a dataframe. Will actually plot using ggplot2 below (easier and nicer looking)
mds.data <- plotMDS(y, top=1000, gene.selection = "common")

# Extract variance explained for principal component 1 (PC1) and PC2 - this is the percent of variation in the data that the axis "captures" 
pc1_varExplained = mds.data$var.explained[1]
pc2_varExplained = mds.data$var.explained[2]

# Build dataframe of MDS results for plotting with ggplot2. Define columns, then combine into a dataframe
sample <- row.names(mds.data$distance.matrix.squared)
pc1 <- mds.data$eigen.vectors[,1]
pc2 <- mds.data$eigen.vectors[,2]

mds_results <- data.frame(sample,group,pc1,pc2) %>% 
  mutate(temp = str_sub(group,2,3)) %>% # make a column with just the temperature
  mutate(season = str_sub(group,1,1)) %>%  # make a column with just the season
  mutate(time = str_split_fixed(sample,'[_]',3)[,2] %>% as.numeric()) %>% 
  mutate(day = str_split_fixed(sample,'[_]',3)[,3])

# Plot MDS results - included a few variations here

# # Plot with sample names, colored by treatment
# ggplot(mds_results,aes(x=pc1,y=pc2,label=sample,color=group)) +
#   geom_text() +
#   labs(color='orange', 'blue',
#        x=paste('Principal Component 1 (',round(pc1_varExplained*100, digits = 2),'%)',sep=''),
#        y=paste('Principal Component 2 (',round(pc2_varExplained*100, digits = 2),'%)',sep='')) +
#   theme_linedraw(base_size = 16) 
# 
# # Plot points, colored by treatment
# ggplot(mds_results,aes(x=pc1,y=pc2,label=sample,color=group)) +
#   geom_point(size=3) +
#   labs(color='orange', 'blue',
#        x=paste('Principal Component 1 (',round(pc1_varExplained*100, digits = 2),'%)',sep=''),
#        y=paste('Principal Component 2 (',round(pc2_varExplained*100, digits = 2),'%)',sep='')) +
#   theme_linedraw(base_size = 16) 
# 
# # Plot points with shape based on temperature and color based on season
# ggplot(mds_results,aes(x=pc1,y=pc2,label=sample,color=season,shape=temp)) +
#   geom_point(size=3) +
#   labs(shape='Temperature', col = "Season",
#        color='orange', 'blue',
#        x=paste('Principal Component 1 (',round(pc1_varExplained*100, digits = 2),'%)',sep=''),
#        y=paste('Principal Component 2 (',round(pc2_varExplained*100, digits = 2),'%)',sep='')) +
#   scale_color_manual(values = c("orange", "darkblue")) +
#   theme_linedraw(base_size = 16) +
#   theme(text=element_text(family="sans", size=16)) 

# New plot with transparency based on time
ggplot(mds_results,aes(x=pc1,y=pc2,label=sample,color=season,shape=temp,alpha=time)) +
  geom_point(size=3) +
  labs(shape='Temperature', col = "Season",
       alpha='Circadian Time',
       x=paste('Principal Component 1 (',round(pc1_varExplained*100, digits = 2),'%)',sep=''),
       y=paste('Principal Component 2 (',round(pc2_varExplained*100, digits = 2),'%)',sep='')) +
  scale_alpha_continuous(range = c(0.05,1),breaks=c(0,3,6,9,12,15,18,21)) +
  scale_color_manual(values = c("orange", "darkblue")) +
  theme_linedraw(base_size = 16) +
  theme(text=element_text(family="sans", size=16)) 

# Same as above but with day 1, day 2, and day 3 plotted separately
ggplot(mds_results,aes(x=pc1,y=pc2,label=sample,color=season,shape=temp,alpha=time)) +
  geom_point(size=3) +
  labs(shape='Temperature', col = "Season",
       alpha='Circadian Time',
       x=paste('Principal Component 1 (',round(pc1_varExplained*100, digits = 2),'%)',sep=''),
       y=paste('Principal Component 2 (',round(pc2_varExplained*100, digits = 2),'%)',sep='')) +
  scale_alpha_continuous(range = c(0.05,1),breaks=c(0,3,6,9,12,15,18,21)) +
  scale_color_manual(values = c("orange", "darkblue")) +
  facet_wrap(~day) +
  theme_linedraw(base_size = 16) +
  theme(text=element_text(family="sans", size=16)) 





  
