library(dplyr)
library(RVenn)

mmm = readRDS('data/newData/countsInPlants.rds')
foo = readRDS('data/newData/countsInInoculum.rds')

ggvenn(Venn(list(mmm$Combined, foo$Combined)))

mmm2 = inner_join(mmm, foo, by = 'Combined')
mmm2 = mmm2[complete.cases(mmm2), ]
mmm2 = mmm2[!(mmm2$Abundance == 0), ]
mmm2$newAbundance = round((mmm2$Count / mmm2$Abundance) * 1000)

rdf = reshape(data = mmm2[, c(1,3,9)], idvar = 'Isolate.x', timevar = 'Plant', direction = 'wide')
row.names(rdf) = rdf$Isolate.x
rdf = rdf[,-1]
rdf[is.na(rdf)] = 0
colnames(rdf) = sapply(strsplit(colnames(rdf), split = '\\.'), function(x) x[[2]])
saveRDS(rdf, file = 'data/newData/adjustedPlantCounts.rds')
