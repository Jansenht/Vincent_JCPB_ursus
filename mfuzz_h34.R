if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Mfuzz")
library(Mfuzz)

###########################################################H37
# to get files for use below, refer to "cluster2.R'
# line plots gene expression on y, time of day on x (needs to be paired wi/ cluster2.r to work)

h34 = subset(Hibernation34_meta_input,`...1` %in% matchedh34$CycID)
write.csv(h34, "Desktop/for_mfuzz_h34.csv")

# make time points 1 row name in excel

h34 = for_mfuzz_h34

h34 = h34 [ -c(1:1)]


#create temp file, make gene expression set, and filter out NAs
tmp =  tempfile()
write.table(h34,file=tmp, sep='\t', quote = F, col.names=NA)
expr = table2eset(file = tmp)
nona = filter.NA(expr)

# standarize data
data.s = standardise(expr)

# estimate fuzzifier lol
m1 = mestimate(data.s)
m1

# determine cluster number

Dmin(data.s, m=m1, crange=seq(2,22,1), repeats=3, visu=TRUE)
clust=10
c = mfuzz(data.s,c=clust,m=m1)

# visuals

mfuzz.plot(data.s,cl=c,mfrow=c(1,1),time.labels=c(0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69),new.window=FALSE)


# check cluster matrix

h34_cluster_matrix = cor(t(c[[1]]))

write.csv(h34_cluster_matrix, "Desktop/h34_cluster_matrix.csv")

# extract lists of genes in clusters
acore <- acore(data.s,c,min.acore=0)

acore_list <- do.call(rbind, lapply(seq_along(acore), function(i){ data.frame(CLUSTER=i, acore[[i]])}))

# export and line up gene names

write.csv(acore_list, "Desktop/h34_mfuzz_genelist_acore.csv")
