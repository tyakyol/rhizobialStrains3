library(magrittr)
library(dplyr)
library(vegan)
library(ggplot2)
library(pheatmap)
library(RVenn)

# Count table reshaping =========================
df = read.delim('data/inoculumCounts.csv', sep = ',')
df = df[df$pl_name != 'none', ]
df = df[, c(6, 3, 4)]
rdf = reshape(data = as.data.frame(df), idvar = 'pl_name', timevar = 'ID', direction = 'wide')
rdf[is.na(rdf)] = 0
row.names(rdf) = rdf$pl_name
rdf = rdf[, -1]
colnames(rdf) = colnames(rdf) %>% strsplit(split = 'ID_') %>% sapply(function(x) x[[2]])
dfRA = round(decostand(rdf, method = 'total', MARGIN = 2), 8)
foo = reshape(data = dfRA, idvar = 'ID', ids = row.names(dfRA), time = names(dfRA), timevar = 'isolate', direction = 'long', varying = list(names(dfRA)))
colnames(foo) = c('Batch', 'Abundance', 'Isolate')
foo$Batch = gsub(foo$Batch, pattern = 'Inoculum', replacement = 'Batch')
foo$Combined = paste(foo$Batch, foo$Isolate, sep = '_')
saveRDS(foo, file = 'data/newData/countsInInoculum.rds')
