library(qqman)
library(purrr)
library(doBy)
library(dplyr)

# Load the data
do3 = list()
for(i in dir(path = 'results/genomicPredictionForAbundance/specialistVarImp/', full.names = T)) {
  do3[[i]] = read.delim(i, sep = ',')
}

# Merge
ov = reduce(do3, function(x, y) full_join(x, y, by = 'X'))
ov[is.na(ov)] = 0
row.names(ov) = ov$X
ov = ov[,-1]

# Adjust for manhattan plot
lt0 = apply(ov, 1, function(x) sum(x > 0))
lt0 = data.frame(chr = sapply(strsplit(row.names(ov), split = '__'), function(x) x[[1]]),
                 pos = sapply(strsplit(row.names(ov), split = '__'), function(x) x[[2]]),
                 pv = lt0)
lt0$chr = as.numeric(substr(lt0$chr, start = 2, stop = 100))
lt0$pos = as.numeric(gsub(lt0$pos, pattern = '`', replacement = ''))
lt0$snp = row.names(lt0)
write.csv(lt0, 'results/manhattanAbundance/specialists.csv')

# Plot
pdf('results/29_extendedFig6f.pdf', width = 7, height = 5, useDingbats = F)
qqman::manhattan(x = lt0, chr = 'chr', bp = 'pos', p = 'pv', snp = 'snp', logp = F, suggestiveline = 60, genomewideline = F, chrlabs = c('chr1L', 'chr1S', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6'), xlab = '', ylab = 'Number of runs', ylim = c(0, 105))
dev.off()
