library(caret)
library(dplyr)

set.seed(32); numbers = sample(x = 1:100000, size = 1000, replace = FALSE)

partition = function(s) {
  set.seed(s)
  inTrain = createDataPartition(y = 1:212, p = 0.8, list = FALSE)
  return(inTrain)
}

parts = list()
for(n in 1:1000) {
  parts[[n]] = partition(numbers[n])
}

pheno = read.delim('results/shannonDiversity.csv', sep = ',')
pheno = pheno[, c('Line1', 'alpS')]
pheno = pheno %>% group_by(Line1) %>% summarise(Specialists = mean(alpS))
pheno$Specialists = round(pheno$Specialists, 6)

template = read.delim('data/templateForGWAS.csv', sep = ',')
class(pheno) = class(pheno)[3]
row.names(pheno) = pheno$Line1
pheno = pheno[template$Taxa,]
row.names(pheno) = NULL
write.csv(pheno, 'results/specialistsForGWAS.csv')

writer = function(x, i) {
  df = x[i, ]
  write.table(df,
              paste0('results/inputForGWAS/specialists/specialists', as.character(numbers[n]), '.csv'),
              quote = TRUE, sep = ',',
              row.names = FALSE, col.names = TRUE)
}

for(n in 1:100) {
  writer(pheno, parts[[n]])
}
