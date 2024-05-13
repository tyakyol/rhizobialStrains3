# https://github.com/cks2903/White_Clover_GenomicPrediction_2020/blob/aeb1147d640cf4326a619ed14806b9132db6ba98/GRM/GRMscript
# https://github.com/cks2903/White_Clover_GenomicPrediction_2020/blob/aeb1147d640cf4326a619ed14806b9132db6ba98/PlotScripts/GRM_PCA.R

library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Loading the data ==============================
geno = read.delim(unz('data/snps.csv.zip', 'snps.csv'), sep = ',', header = T)
geno = geno[, c(-1, -2)]
geno = geno[ , order(colnames(geno))]
genoT = t(geno)

# VanRaden (2009)'s method ======================
M = dim(genoT)[2]
N = dim(genoT)[1]

# Cathrine's function
CALCULATE_P_MULTIPLIED_BY_Q = function(SNP_column) {
  SNP_column = data.matrix(SNP_column)
  countsfirstallele = sum(SNP_column, na.rm = TRUE)
  p_freq = countsfirstallele/(length(na.omit(SNP_column))*2)
  q_freq = 1-p_freq
  value = p_freq*q_freq
  return(value)
}

# Apply function to all columns
# Sum across results
results = apply(genoT, 2, CALCULATE_P_MULTIPLIED_BY_Q)
SUM = sum(results)

#Centered genotype matrix. Every column (SNP) get a mean of 0
Z = scale(genoT, scale = FALSE) 
# Z[is.na(Z)] = 0  # Set NA to zero in scaled matrix. That is setting them equal to avg.

# Calculate GRM matrix
G = (Z%*%t(Z))/(2*SUM)

# PCA ===========================================
d.pca = prcomp(G,
              center = TRUE,
              scale. = TRUE)

percentVar = d.pca$sdev^2 / sum( d.pca$sdev^2)
pV = percentVar
barplot((pV*100)[1:10], names.arg = paste0('PC', 1:10), main = 'Scree plot')
percentVar = round(100 * percentVar)
d.pca.rot = as.data.frame(d.pca$rotation)

ggplot(d.pca.rot, aes(PC1, PC2)) +
  coord_fixed() +
  geom_point(size = 4, shape = 21,
             color = 'dodgerblue3', fill = 'lightblue', stroke = 1.5) +
  xlab(paste0('PC1: ', percentVar[1], '% variance')) +
  ylab(paste0('PC2: ', percentVar[2], '% variance')) + 
  theme_classic() +
  theme(aspect.ratio = 1, axis.ticks = element_blank())
# ggsave(filename = 'results/GRM.pdf', device = 'pdf', width = 5, height = 5)

write.csv(d.pca.rot, file = 'results/grm.csv', quote = FALSE)
