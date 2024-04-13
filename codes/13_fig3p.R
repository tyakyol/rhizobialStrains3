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

# Relative abundance as int in 10e5 =============
dfRA = readRDS('data/newData/adjustedPlantCounts.rds')

# ===============================================
specialists = row.names(isoInfo)[isoInfo$class == 'Specialists']
dominants__ = row.names(isoInfo)[isoInfo$class == 'Dominants']
generalists = row.names(isoInfo)[isoInfo$class == 'Generalists']
transients = row.names(isoInfo)[isoInfo$class == 'Transients']

specialists = colSums(dfRA[specialists,])
dominants__ = colSums(dfRA[dominants__,])
generalists = colSums(dfRA[generalists,])
transients = colSums(dfRA[transients,])

fff = data.frame(Specialists = specialists,
                 Generalists = generalists,
                 Dominants__ = dominants__,
                 Transients = transients,
                 DandS = dominants__ + specialists,
                 DandG = dominants__ + generalists,
                 SandG = specialists + generalists)

cor.test(fff$DandS, fff$Generalists)
cor.test(fff$DandG, fff$Specialists)
cor.test(fff$SandG, fff$Dominants__)
cor.test(fff$Generalists, fff$Transients)
cor.test(fff$Dominants__, fff$Transients)
cor.test(fff$Specialists, fff$Transients)
cor.test(fff$Specialists+fff$Dominants__+fff$Transients, fff$Generalists)

# Nicer plots
ggplot(data = fff, aes(x = DandS, y = Generalists)) +
  geom_point(color = 'black', size = 4, shape = 21, stroke = 1.2) +
  geom_smooth(method = 'lm', se = FALSE, size = 3, color = colors_$yetanotherpurple) +
  xlab('Cumulative abundance of Dominants and\nSpecialists') +
  ylab('Cumulative abundance of\nGeneralists') +
  # annotate(geom = 'label',x = 15.5, y = 16.75, size = 5.5, 
  #          label = 'r = -0.86\np < 1e-3') +
  annotate(geom = 'label',x = 385000, y = 105000, size = 5.5,
  label = 'r = -0.45\np < 1e-3') +
  theme(legend.direction = 'vertical',
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_blank(), panel.grid = element_blank(), 
        axis.ticks = element_blank())
ggsave(filename = 'results/13_fig3p.pdf', device = 'pdf', width = 6, height = 6)

# ggplot(data = fff, aes(x = DandG, y = Specialists)) +
#   geom_point(color = 'black', size = 4, shape = 21, stroke = 1.2) +
#   geom_smooth(method = 'lm', se = FALSE, size = 3, color = colors_$yetanotherpurple) +
#   # xlab('Cumulative abundance of Dominants') + 
#   # ylab('Cumulative abundance of Generalists') +
#   # annotate(geom = 'label',x = 14.65, y = 4, size = 5.5, 
#   #          label = 'r = -0.15\np < 1e-3') +
#   # annotate(geom = 'label',x = 29000, y = 27000, size = 5.5, 
#            # label = 'r = -0.21\np < 1e-3') +
#   theme(legend.direction = 'vertical',
#         axis.text = element_text(size = 14),
#         axis.title = element_text(size = 18),
#         legend.title = element_blank(), panel.grid = element_blank(), 
#         axis.ticks = element_blank())
# # ggsave(filename = 'results/12_fig3n.pdf', device = 'pdf', width = 6, height = 6)
# 
# 
# ggplot(data = fff, aes(x = SandG, y = Dominants__)) +
#   geom_point(color = 'black', size = 4, shape = 21, stroke = 1.2) +
#   geom_smooth(method = 'lm', se = FALSE, size = 3, color = colors_$yetanotherpurple) +
#   # xlab('Cumulative abundance of Specialists') + 
#   # ylab('Cumulative abundance of Generalists') +
#   # annotate(geom = 'label',x = 12.3, y = 4, size = 5.5, 
#   #          label = 'r = -0.03\nNS') +
#   # annotate(geom = 'label',x = 58000, y = 28000, size = 5.5, 
#            # label = 'r = -0.05\nNS') +
#   theme(legend.direction = 'vertical',
#         axis.text = element_text(size = 14),
#         axis.title = element_text(size = 18),
#         legend.title = element_blank(), panel.grid = element_blank(), 
#         axis.ticks = element_blank())
# # ggsave(filename = 'results/12_fig3o.pdf', device = 'pdf', width = 6, height = 6)
