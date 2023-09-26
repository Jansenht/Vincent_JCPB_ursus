library(WGCNA)
library(tidyverse)     # tidyverse will pull in ggplot2, readr, other useful libraries
library(magrittr)      # provides the %>% operator
library(ggplot2)


# allowWGCNAThreads()          # allow multi-threading (optional)
#> Allowing multi-threading with up to 4 threads.


# ****** Added by Blair: making sample info table -------------------------

datTraits <- as.data.frame(row.names(dataset_numeric)) %>%  # make dataframe with sample names
  select(sample = 1) %>% # rename first column as "sample"
  mutate(season = ifelse(str_sub(sample,1,1) == 'A',0,1)) %>%  # New column "season", active = 0, hibernation = 1
  mutate(temp = ifelse(str_sub(sample,2,3) == '34',0,1)) %>%  # New column "temp", 34 = 0, 37 = 1
  column_to_rownames('sample') # move "sample" column to rownames

# Cluster samples by trait
sampleTree2 = hclust(dist(dataset_numeric), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")



# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  dataset_numeric,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
  )


par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

picked_power = 9
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
net <- blockwiseModules(dataset_numeric,                # <= input here

                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",

                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 10000, ### ********** CHANGED BY BLAIR: 4000 -> 10000

                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,

                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",

                          # == Output Options
                          numericLabels = T,
                          verbose = 3)


cor <- temp_cor     # Return cor function to original namespace
######################################################################
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  net$dendrograms[[1]],
  mergedColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )


# ********** Added by Blair: plot dendrograms for other blocks -----------------------

plotDendroAndColors(
  net$dendrograms[[2]],
  mergedColors[net$blockGenes[[2]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

plotDendroAndColors(
  net$dendrograms[[3]],
  mergedColors[net$blockGenes[[3]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )


# export gene modules

module_df <- data.frame(
  gene_id = names(net$colors),
  colors = labels2colors(net$colors)
)

module_df %>% group_by(colors) %>% tally() %>% arrange(-n)

module_df[1:5,]

write_delim(module_df,
            file = "/Users/elleryvincent/Desktop/gene_modules.txt",
            delim = "\t")


# ********* Added by Blair: Correlating modules with traits ---------------

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

# Define numbers of genes and samples
nGenes = ncol(dataset_numeric);
nSamples = nrow(dataset_numeric);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(dataset_numeric, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.2,
               zlim = c(-1,1),
               cex.lab.y = 0.4,
               main = paste("Module-trait relationships"))


# To further explore: -----------------------------------------------------
# Follow section: 3.b Gene relationship to trait and important modules: Gene Significance and Module Membership
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-03-relateModsToExt.pdf







#############################################################################################
# find modules of interest pick out a few modules of interest here
modules_of_interest = c("green", "turquoise", "tan")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes
library(tidyr)
dataset_numeric[1:5,1:10]

subexpr = dataset_numeric
[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = "blue"),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")


# create coexpression network
genes_of_interest = module_df %>%
  subset(colors %in% modules_of_interest)

expr_of_interest =dataset_numeric[genes_of_interest$gene_id, ]
expr_of_interest[1:5,1:5]

# Only recalculate TOM for modules of interest (faster, altho there's some online discussion if this will be slightly off)
TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = picked_power)


# Add gene names to row and columns
row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)

edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )

head(edge_list)


# Export Network file to be read into Cytoscape, VisANT, etc
write_delim(edge_list,
            file = "/Users/elleryvincent/Desktop/edgelist.tsv",
            delim = "\t")

