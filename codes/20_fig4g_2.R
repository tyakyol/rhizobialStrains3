# This code doesn't perform any analysis, it only generates the data for GWAS and RF.
# For Generalists.

library(caret)
library(dplyr)

# Numbers for random seeds
set.seed(32); numbers = sample(x = 1:100000, size = 1000, replace = FALSE)

# Generate partitions
partition = function(s) {
  set.seed(s)
  inTrain = createDataPartition(y = 1:212, p = 0.8, list = FALSE)
  return(inTrain)
}

parts = list()
for(n in 1:1000) {
  parts[[n]] = partition(numbers[n])
}

# Diversity data
pheno = read.delim('results/shannonDiversity.csv', sep = ',')
pheno = pheno[, c('Line1', 'alpG')]
pheno = pheno %>% group_by(Line1) %>% summarise(Generalists = mean(alpG))
pheno$Generalists = round(pheno$Generalists, 6)

# To GWAS format
template = read.delim('data/templateForGWAS.csv', sep = ',')
class(pheno) = class(pheno)[3]
row.names(pheno) = pheno$Line1
pheno = pheno[template$Taxa,]
row.names(pheno) = NULL
colnames(pheno)[1] = 'Taxa'
write.csv(pheno, 'results/generalistsForGWAS.csv')

# Write the data tables out
writer = function(x, i) {
  df = x[i, ]
  write.table(df,
              paste0('results/inputForGWAS/generalists/generalists', as.character(numbers[n]), '.csv'),
              quote = TRUE, sep = ',',
              row.names = FALSE, col.names = TRUE)
}

# Only first 100
for(n in 1:100) {
  writer(pheno, parts[[n]])
}
