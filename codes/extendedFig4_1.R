### Cumulative abundance ###

library(magrittr)
library(dplyr)
library(vegan)
library(ggplot2)
library(purrr)
library(scales)
library(seqtime)
library(caret)

# Load the data =================================
counts = readRDS('data/counts.rds')
st = counts$samples
# pgR = counts$samples
isoInfo = read.csv('results/competitiveness/table1.csv', row.names = 1)
grm = read.delim('results/grm.csv', sep = ',')
rm(counts)

# Relative abundance as int in 10e5 =============
isoToKeep = row.names(isoInfo)[isoInfo$class != 'Transients']
dfRA = readRDS('data/newData/adjustedPlantCounts.rds')
dfRA = dfRA[isoToKeep, ]

GG = dfRA[row.names(isoInfo)[isoInfo$class == 'Generalists'], ]
DD = dfRA[row.names(isoInfo)[isoInfo$class == 'Dominants'], ]
SS = dfRA[row.names(isoInfo)[isoInfo$class == 'Specialists'], ]

# Cumulative abundance
raG = colSums(GG)
raD = colSums(DD)
raS = colSums(SS)

# Evenness
eveG = sheldon(t(GG))
eveD = sheldon(t(DD))
eveS = sheldon(t(SS))


# Alpha diversity
alpG = vegan::diversity(t(GG), index = 'shannon')
alpD = vegan::diversity(t(DD), index = 'shannon')
alpS = vegan::diversity(t(SS), index = 'shannon')


# MDS ===========================================
mdsG = as.data.frame(cmdscale(vegdist(t(GG), method = 'bray'), k = 4))
mdsD = as.data.frame(cmdscale(vegdist(t(DD), method = 'bray'), k = 4))
mdsS = as.data.frame(cmdscale(vegdist(t(SS), method = 'bray'), k = 4))

lf = st
lf$Plant = paste0('ID_100', lf$Plant)
mergD = lf
mergD = left_join(mergD, grm[,1:11], by = c('Line1' = 'X'))
mergD = cbind(mergD, raG, raD, raS)
write.csv(mergD, 'results/cumulativeAbundance.csv')
mergD = mergD[complete.cases(mergD), ]
modB = as.data.frame(model.matrix(~mergD$Batch)[,-1])
mergD = cbind(mergD, modB)
mergD = mergD[, c(-1, -3, -4, -2)]

# ranger with top 5 =============================
xg1 = function(seed) {
  set.seed(seed)
  
  sss = mergD[, -1]
  inTrain = createDataPartition(y = mergD$Batch, p = 0.8, list = FALSE)
  training = sss[inTrain,]
  testing = sss[-inTrain,]
  
  fitControl = trainControl(method = 'repeatedcv', number = 6, repeats = 2, 
                            search = 'random')
  modFit = train(Biomass ~ ., method = 'ranger', data = training, 
                 verbose = FALSE,
                 trControl = fitControl)
  
  qqq = predict(modFit, training[,setdiff(colnames(testing), 'Biomass')])
  out1 = list(qqq, training, 
              predict(modFit, testing[,setdiff(colnames(testing), 'Biomass')]),
              testing)
  
  return(out1)
}

set.seed(33); n = sample(1:10000, size = 100, replace = FALSE)

res1_ = c()
res1__ = c()
for(i in 1:100) {
  res1 = xg1(n[i])
  res1_[i] = R2(res1[[1]], res1[[2]][,'Biomass'])
  res1__[i] = R2(res1[[3]], res1[[4]][,'Biomass'])
} 

saveRDS(res1_,  'results/growthPredictionWithDiversity/abundance1.rds')
saveRDS(res1__, 'results/growthPredictionWithDiversity/abundance2.rds')
