library(vegan)
library(ggplot2)
library(DOC)

source('codes/__rarefy__.R')

# Load the data =================================
isoInfo = read.csv('results/competitiveness/table1.csv', row.names = 1)

# Prepare the data for DOC analysis =============
isoToKeep = row.names(isoInfo)[isoInfo$class == 'Generalists']
df = readRDS('data/newData/adjustedPlantCounts.rds')
df = df[isoToKeep, ]
df = rarefyFilter(df, min = 1000)$rar
dfRA = vegan::decostand(df, method = 'total', MARGIN = 2)

# DOC
ad4 = DOC(dfRA, R = 100, span = 0.5, cores = 10)

# Plot
pdf('results/09_fig3h.pdf', height = 7, width = 7)
plot(ad4)
dev.off()

# Fns
print(paste0('Fns of generalists = ', mean(ad4$FNS[,1])))
