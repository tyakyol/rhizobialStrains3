library(magrittr)
library(dplyr)
library(vegan)
library(ggplot2)
library(purrr)
library(scales)

source('codes/__colors__.R')

# Load the data =================================
counts = readRDS('data/counts.rds')
dfRA = readRDS('data/newData/adjustedPlantCounts.rds')
st = counts$samples
growth = read.delim('results/competitiveness/table1.csv', sep = ',', row.names = 1)
growth = growth[, -1]
rm(counts)

# Niche breadth (Niwa et al 2018) ===============
commonness = function(arr) {
  a = (arr / sum(arr))^2
  a = 1/sum(a)
  return(a)
}

cmm = apply(dfRA, 1, commonness)
growth = growth[names(cmm),]
growth$commonness = cmm

# Plot ==========================================
c1 = c(colors_$laci, colors_$yellow, colors_$turk, colors_$lightpurple)
c2 = c(colors_$coffee, colors_$lightgreen, colors_$anotherred, colors_$white)

growth$newSoil = ifelse(growth$soil %in% c('UG', 'SE', 'NO'), growth$soil, 'Others')
growth$newSoil = factor(growth$newSoil, levels = c('UG', 'SE', 'NO', 'Others'))

ggplot(data = growth, aes(y = commonness, x = class)) +
  geom_boxplot(aes(color = class), fill = NA, outlier.color = NA, size = 3) +
  geom_jitter(aes(fill = newSoil), width = 0.2, size = 5, shape = 21, stroke = 0.3) +
  scale_fill_manual(values = c2) +
  scale_color_manual(values = c1) +
  ylab('Niche breadth') + 
  xlab('') +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 12))
ggsave(filename = 'results/03_fig2d.pdf', device = 'pdf', width = 6, height = 6)
