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

# Relative abundance as int in 10e5 =============
isoToKeep = row.names(isoInfo)[isoInfo$class != 'Transients']
dfRA = readRDS('data/newData/adjustedPlantCounts.rds')
dfRA = dfRA[isoToKeep, ]
domi = isoInfo[isoInfo$class == 'Dominants',]
dfRA = dfRA[row.names(domi),]

# MDS ===========================================
mp = cmdscale(vegdist(t(dfRA), method = 'cao'), k = 4)
mp = as.data.frame(mp)
mp$Plant = row.names(mp)

# Modify the growth data ========================
# 1. calculate replicate counts
# 2. sort by genotype and batch
# 3. sometimes a batch can include more than 1 replicates of the same genotype, in those cases i took the mean of those values (eg if there are two growth values in batch-1 for Hedin)
pgR = pgR %>% group_by(Line1) %>% mutate(repCount = sum(!is.na(Biomass)))
pgR = arrange(pgR, Line1, Ref)
pgR = pgR %>% group_by(Line1, Ref) %>% summarize(newBM = mean(Biomass))

# Growth analysis with LMM ======================
lf = left_join(st, grm[, 1:3], by = c('Line1' = 'X'))  # combine with grm
lf$Plant = paste0('ID_100', lf$Plant)
lf = left_join(lf, mp, by = 'Plant')

# Calculate evenness
# eve = as.data.frame(sheldon(t(dfRA)))
eve = as.data.frame(colSums(dfRA))

# Merge with biomass
eve$Plant = colnames(dfRA)
colnames(eve)[1] = 'Evenness'
eve = left_join(eve, lf, by = 'Plant')
eve$Comb = paste(eve$Ref, eve$Line1, sep = '_')
pgR$Comb = paste(pgR$Ref, pgR$Line1, sep = '_')
# write.csv(eve, 'results/evenness/dominants.csv')

fff = left_join(eve, pgR, by = 'Comb')
fff = fff[complete.cases(fff), ]
# fff$Evenness = fff$Evenness / max(fff$Evenness)
fff$Evenness = as.numeric(scale(fff$Evenness))

# togwas = read.delim('~/Desktop/806rcbc491_3.csv', sep = ',')
# fff3 = fff %>% group_by(Line1.x) %>% summarise(Evenness = mean(Evenness))
# togwas = left_join(togwas, fff3, by = c('Taxa' = 'Line1.x'))
# togwas = togwas[, -2]
# write.table(togwas, 
#             '~/Desktop/evenness.csv', 
#             quote = TRUE, sep = ',', 
#             row.names = FALSE, col.names = TRUE)

# Test
m = lmerTest::lmer(newBM ~ Evenness + (1 | Line1.x) + (1 | Batch) + V2 + V3 + V4, data = fff)
# m = lmerTest::lmer(newBM ~ (1|newFactor) + (1|Line1.x) + (PC1 || Batch) + Batch + Ref.x, data = fff)
mm = summary(m)
br = broom.mixed::tidy(m, conf.int = TRUE)
# gt(br)
