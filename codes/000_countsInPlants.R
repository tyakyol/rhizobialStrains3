library(magrittr)
library(dplyr)
library(vegan)
library(ggplot2)

# Count table reshaping =========================
df = read.delim('data/counts.csv', sep = ',')
df = df[df$pl_name != 'none', ]
df = df[, c(4, 1, 2)]
ndf = df %>% group_by(pl_name, ID) %>% summarise(count = sum(count))
rdf = reshape(data = as.data.frame(ndf), idvar = 'pl_name', timevar = 'ID', direction = 'wide')
rdf[is.na(rdf)] = 0
row.names(rdf) = rdf$pl_name
rdf = rdf[, -1]
colnames(rdf) = colnames(rdf) %>% strsplit(split = '\\.') %>% sapply(function(x) x[[2]])
rdf = rdf[, colSums(rdf) > 2000]  # 608 reduced to 603
rm(df, ndf)

dfRA = round(decostand(rdf, method = 'total', MARGIN = 2), 8)
foo = reshape(data = dfRA, idvar = 'isolate', ids = row.names(dfRA), time = names(dfRA), timevar = 'ID', direction = 'long', varying = list(names(dfRA)), v.names = 'Abundance')


st = read.delim('data/metadata.csv', sep = ',')
plantsToDiscard = c('Hedin2', 'VF172-3Cv', 'Giza402Gö', 'ILB938/2', 'Melodie', 'Alameda')
st = st[!(st$Line1 %in% plantsToDiscard), ]
row.names(st) = paste0('ID_100', st$Plant)
length(intersect(row.names(st), colnames(rdf)))  # 603
st$Plant = row.names(st)
st__ = st[colnames(rdf), ]
identical(row.names(st__), colnames(rdf))  # T
st = st__[, c('Plant', 'Batch')]


mmm = left_join(foo, st, by = c('ID' = 'Plant'))
mmm$Batch = gsub(mmm$Batch, pattern = ' ', replacement = '_')
mmm$Batch = gsub(mmm$Batch, pattern = 'MMS_', replacement = '')
mmm$Comb = paste(mmm$Batch, mmm$isolate, sep = '_')
colnames(mmm) = c('Plant', 'Count', 'Isolate', 'Batch', 'Combined')
saveRDS(mmm, file = 'data/newData/countsInPlants.rds')
