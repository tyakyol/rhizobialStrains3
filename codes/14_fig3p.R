library(magrittr)
library(dplyr)
library(vegan)
library(ggplot2)
library(purrr)
library(scales)

source('codes/__colors__.R')

# Load the data =================================
counts = readRDS('data/counts.rds')
st = counts$samples
pgR = counts$samples
isoInfo = read.csv('results/competitiveness/table1.csv', row.names = 1)
grm = read.delim('results/grm.csv', sep = ',')
rm(counts)

# Abundance
dfRA = readRDS('data/newData/adjustedPlantCounts.rds')

# Split data by isolates
specialists = row.names(isoInfo)[isoInfo$class == 'Specialists']
dominants__ = row.names(isoInfo)[isoInfo$class == 'Dominants']
generalists = row.names(isoInfo)[isoInfo$class == 'Generalists']
transients = row.names(isoInfo)[isoInfo$class == 'Transients']

# Total abundance of isolates
specialists = colSums(dfRA[specialists,])
dominants__ = colSums(dfRA[dominants__,])
generalists = colSums(dfRA[generalists,])
transients = colSums(dfRA[transients,])

# Making a data frame
fff = data.frame(Specialists = specialists,
                 Generalists = generalists,
                 Dominants__ = dominants__,
                 Transients = transients,
                 DandS = dominants__ + specialists,
                 DandG = dominants__ + generalists,
                 SandG = specialists + generalists)

# Correlation tests
cor.test(fff$DandS, fff$Generalists)
cor.test(fff$DandG, fff$Specialists)
cor.test(fff$SandG, fff$Dominants__)
cor.test(fff$Generalists, fff$Transients)
cor.test(fff$Dominants__, fff$Transients)
cor.test(fff$Specialists, fff$Transients)
cor.test(fff$Specialists+fff$Dominants__+fff$Transients, fff$Generalists)

# Plot ==========================================
ggplot(data = fff, aes(x = DandS, y = Generalists)) +
  geom_point(color = 'black', size = 4, shape = 21, stroke = 1.2) +
  geom_smooth(method = 'lm', se = FALSE, size = 3, color = colors_$yetanotherpurple) +
  xlab('Cumulative abundance of Dominants and\nSpecialists') +
  ylab('Cumulative abundance of\nGeneralists') +
  annotate(geom = 'label',x = 385000, y = 105000, size = 5.5,
  label = 'r = -0.45\np < 1e-3') +
  theme(legend.direction = 'vertical',
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_blank(), panel.grid = element_blank(), 
        axis.ticks = element_blank())
ggsave(filename = 'results/14_fig3p.pdf', device = 'pdf', width = 6, height = 6)
