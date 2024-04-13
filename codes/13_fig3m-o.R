library(magrittr)
library(dplyr)
library(vegan)
library(ggplot2)
library(purrr)
library(scales)

source('codes/__colors__.R')

# Load the data =================================
counts = readRDS('data/counts.rds')
df = counts$counts
st = counts$samples
pgR = counts$samples
isoInfo = read.csv('results/competitiveness/table1.csv', row.names = 1)
grm = read.delim('results/grm.csv', sep = ',')
rm(counts)

# Relative abundance as int in 10e5 =============
isoToKeep = row.names(isoInfo)[isoInfo$class != 'Transients']
dfRA = readRDS('data/newData/adjustedPlantCounts.rds')
dfRA = dfRA[isoToKeep, ]

# Growth analysis with LMM ======================
specialists = row.names(isoInfo)[isoInfo$class == 'Specialists']
dominants__ = row.names(isoInfo)[isoInfo$class == 'Dominants']
generalists = row.names(isoInfo)[isoInfo$class == 'Generalists']

specialists = colSums(dfRA[specialists,])
dominants__ = colSums(dfRA[dominants__,])
generalists = colSums(dfRA[generalists,])

# Any correlations between accumulated abundances?
cor.test(specialists, dominants__)
cor.test(specialists, generalists)
cor.test(dominants__, generalists)

fff = data.frame(Specialists = specialists,
                 Generalists = generalists,
                 Dominants__ = dominants__)

# Nicer plots
ggplot(data = fff, aes(x = Specialists, y = Dominants__)) +
  geom_point(color = 'black', size = 4, shape = 21, stroke = 1.2) +
  geom_smooth(method = 'lm', se = FALSE, size = 3, color = colors_$yetanotherpurple) +
  xlab('Cumulative abundance of Specialists') + 
  ylab('Cumulative abundance of Dominants') +
  # annotate(geom = 'label',x = 15.5, y = 16.75, size = 5.5, 
  #          label = 'r = -0.86\np < 1e-3') +
  annotate(geom = 'label',x = 220000, y = 330000, size = 5.5,
  label = 'r = -0.84\np < 1e-3') +
  theme(legend.direction = 'vertical',
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_blank(), panel.grid = element_blank(), 
        axis.ticks = element_blank())
ggsave(filename = 'results/13_fig3m.pdf', device = 'pdf', width = 6, height = 6)

ggplot(data = fff, aes(x = Dominants__, y = Generalists)) +
  geom_point(color = 'black', size = 4, shape = 21, stroke = 1.2) +
  geom_smooth(method = 'lm', se = FALSE, size = 3, color = colors_$yetanotherpurple) +
  xlab('Cumulative abundance of Dominants') + 
  ylab('Cumulative abundance of Generalists') +
  # annotate(geom = 'label',x = 14.65, y = 4, size = 5.5, 
  #          label = 'r = -0.15\np < 1e-3') +
  annotate(geom = 'label',x = 340000, y = 105000, size = 5.5,
           label = 'r = -0.2\np < 1e-3') +
  theme(legend.direction = 'vertical',
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_blank(), panel.grid = element_blank(), 
        axis.ticks = element_blank())
ggsave(filename = 'results/13_fig3n.pdf', device = 'pdf', width = 6, height = 6)


ggplot(data = fff, aes(x = Specialists, y = Generalists)) +
  geom_point(color = 'black', size = 4, shape = 21, stroke = 1.2) +
  geom_smooth(method = 'lm', se = FALSE, size = 3, color = colors_$yetanotherpurple) +
  xlab('Cumulative abundance of Specialists') + 
  ylab('Cumulative abundance of Generalists') +
  # annotate(geom = 'label',x = 12.3, y = 4, size = 5.5, 
  #          label = 'r = -0.03\nNS') +
  annotate(geom = 'label',x = 215000, y = 105000, size = 5.5, 
           label = 'r = -0.05\nNS') +
  theme(legend.direction = 'vertical',
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_blank(), panel.grid = element_blank(), 
        axis.ticks = element_blank())
ggsave(filename = 'results/13_fig3o.pdf', device = 'pdf', width = 6, height = 6)
