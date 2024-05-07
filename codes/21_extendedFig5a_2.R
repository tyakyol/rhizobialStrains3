library(magrittr)
library(dplyr)
library(lme4)
library(vegan)
library(ggplot2)
library(purrr)
library(scales)
library(seqtime)
library(gt)

# Load the data =================================
counts = readRDS('data/counts.rds')
st = counts$samples
pgR = counts$samples
isoInfo = read.csv('results/competitiveness/table1.csv', row.names = 1)
grm = read.delim('results/grm.csv', sep = ',')
rm(counts)

# Relative abundance for generalists ============
isoToKeep = row.names(isoInfo)[isoInfo$class != 'Transients']
dfRA = readRDS('data/newData/adjustedPlantCounts.rds')
dfRA = dfRA[isoToKeep, ]
gene = isoInfo[isoInfo$class == 'Generalists',]
dfRA = dfRA[row.names(gene),]

# MDS ===========================================
mp = cmdscale(vegdist(t(dfRA), method = 'cao'), k = 4)
mp = as.data.frame(mp)
mp$Plant = row.names(mp)

# Modify the growth data ========================
pgR = pgR %>% group_by(Line1) %>% mutate(repCount = sum(!is.na(Biomass)))
pgR = arrange(pgR, Line1, Ref)
pgR = pgR %>% group_by(Line1, Ref) %>% summarize(newBM = mean(Biomass))

# Growth analysis with LMM ======================
lf = left_join(st, grm[, 1:3], by = c('Line1' = 'X'))  # combine with grm
lf$Plant = paste0('ID_100', lf$Plant)
lf = left_join(lf, mp, by = 'Plant')

# Calculate abundance
rela = as.data.frame(colSums(dfRA))

# Merge with biomass
rela$Plant = colnames(dfRA)
colnames(rela)[1] = 'Abundance'
rela = left_join(rela, lf, by = 'Plant')
rela$Comb = paste(rela$Ref, rela$Line1, sep = '_')
pgR$Comb = paste(pgR$Ref, pgR$Line1, sep = '_')

fff = left_join(rela, pgR, by = 'Comb')
fff = fff[complete.cases(fff), ]
fff$Abundance = as.numeric(scale(fff$Abundance))

# Result
m = lmerTest::lmer(newBM ~ Abundance + (1 | Line1.x) + (1 | Batch) + V2 + V3 + V4, data = fff)
mm = summary(m)
br = broom.mixed::tidy(m, conf.int = TRUE)
gt(br)
