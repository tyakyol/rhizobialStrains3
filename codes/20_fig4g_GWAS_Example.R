# This code is an example of the GWAS runs for Dominants. We ran the same code for 100 different partitions.

source('codes/__gapitFunctions__.R')

p1 = 'results/dominants99429Shannon.csv'
myY = read.delim(p1, header = TRUE, sep = ',')
colnames(myY)[1] = 'Taxa'
myGD = read.delim(unz('data/snps.csv.zip', 'snps.csv'), sep = ',', header = T)
myGM = read.delim('data/ProFaba_chrAll_qual20_depth3_2snp_miss50_maf0.01_imputed_GM.csv', sep = ',', header = TRUE)

dir.create('dominants99429.csv')
setwd('dominants99429.csv')

# Function to do GAPIT GWAS
GP_GBLUP = function() {
  myGAPIT = GAPIT(
    Y = myY,
    GD = myGD,
    GM = myGM,
    SNP.MAF = 0.05,
    PCA.totalv = v3,
    model = 'Blink'
    )
}

# Modify
myGD = myGD[, c(-1, -2)]
myGD = myGD[, myY[, 'Taxa']]
myGD = as.data.frame(t(myGD))
myGD[,'taxa'] = row.names(myGD)
myGD = myGD[, c(ncol(myGD), 1:(ncol(myGD)-1))]

# Run
GP_GBLUP()
