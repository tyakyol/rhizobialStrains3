# This code is an example of the RF runs on GWAS results (top 200 markers) for Dominants. We ran the same code for 100 different partitions.
# Output is feature importance.

library(dplyr)
library(ggplot2)
library(caret)

pheno = read.delim('results/dominantsRFShannon.csv', sep = ',')
pheno[,'Dominants'] = log2(pheno[,'Dominants'] + 1e-6)
res = read.delim('results/gwas99429DominantsShannon.csv', sep = ',')
res = res[order(res[,'P.value']),]
res[, 'Label'] = paste(res[, 'Chr'], res[, 'Pos'], sep = '__')
snps = read.delim(unz('data/snps.csv.zip', 'snps.csv'), sep = ',', header = T)
snps[, 'Label'] = paste(snps[, 'CHR'], snps[, 'POS'], sep = '__')
row.names(snps) = snps[, 'Label']
snps2 = snps[res[, 'Label'][1:200],]
snps2 = snps2[,3:214]
snps2 = as.data.frame(t(snps2))
snps2[, 'Taxa'] = row.names(snps2)
snps2 = left_join(snps2, pheno, by = 'Taxa')

# RF with top 200 ===============================
xg1 = function(seed) {
  set.seed(seed)
  sss = snps2[, c(203, 1:200)]
  colnames(sss)[1] = 'abun'
  
  inTrain = createDataPartition(y = sss[,'abun'], p = 0.8, list = FALSE)
  training = sss[inTrain,]
  testing = sss[-inTrain,]
  
  fitControl = trainControl(method = 'repeatedcv', number = 4, repeats = 2, 
                            search = 'random')
  modFit = train(abun ~ ., method = 'ranger', data = training, verbose = FALSE,
                 importance = "permutation",
                 trControl = fitControl)
  
  qqq = predict(modFit, testing)
  out1 = R2(qqq, testing[,1])
  
  return(modFit)
}

res1 = varImp(xg1(99429))[["importance"]]
write.csv(res1, file ='results/dominants99429FeatureImportanceShannon.csv')
