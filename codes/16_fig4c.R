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
rm(counts)

isoInfo = read.csv('results/competitiveness/table1.csv', row.names = 1)
isoInfo = isoInfo[, 1:2]
growth = readRDS('results/lmmWithMDS.rds')
growth = left_join(growth, isoInfo, by = c('term' = 'isolate'))
row.names(growth) = growth$term
# growth = growth[, -1]

# Relative abundance as int in 10e5 =============
isoToKeep = row.names(growth)[growth$class != 'Transients']
dfRA = readRDS('data/newData/adjustedPlantCounts.rds')
dfRA = dfRA[isoToKeep, ]

# Niche breadth (Niwa et al 2018) ===============
commonness = function(arr) {
  a = (arr / sum(arr))^2
  a = 1/sum(a)
  return(a)
}

cmm = apply(dfRA, 1, commonness)
growth = growth[names(cmm),]
growth$commonness = cmm

# Plots
c1 = c(colors_$laci, colors_$yellow, colors_$turk)

ggplot(data = growth, aes(x = commonness, y = estimate)) +
  geom_point(aes(color = class), size = 4, alpha = 0.8) +
  geom_smooth(method = 'lm', se = FALSE, size = 3, aes(color = class)) +
  scale_color_manual(values = c1) +
  xlab('Niche breadth') + 
  ylab('Effect on plant growth') +
  theme(legend.position = 'none', legend.direction = 'horizontal',
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_blank(), panel.grid = element_blank(), 
        axis.ticks = element_blank())
ggsave(filename = 'results/16_fig4c.pdf', device = 'pdf', width = 6, height = 6)

# Correlation test ==============================
cor.test(growth$estimate, growth$commonness)  # No correlation

# By group
cor.test(growth$estimate[growth$class == 'Dominants'],
         growth$commonness[growth$class == 'Dominants'])  # No correlation
cor.test(growth$estimate[growth$class == 'Generalists'],
         growth$commonness[growth$class == 'Generalists'])  # Correlation
cor.test(growth$estimate[growth$class == 'Specialists'],
         growth$commonness[growth$class == 'Specialists'])  # No correlation
