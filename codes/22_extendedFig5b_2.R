library(ggplot2)
library(agricolae)

full = readRDS('results/growthPredictionWithDiversity/abundance2.rds')
redu = readRDS('results/growthPredictionWithDiversity/reduced2.rds')
rand = readRDS('results/growthPredictionWithDiversity/randomAbundance2.rds')

source('codes/__colors__.R')

x = c(full, redu, rand)
f = rep(c('a', 'b', 'c'), each = 100)

HSD.test(lm(x ~ f), console = TRUE, trt = 'f')

c2 = c(colors_$pink, colors_$gray, colors_$gray)

qplot(x = f, y = x, geom = 'blank') + 
  geom_boxplot(size = 2, aes(color = f), outlier.colour = NA) +
  geom_jitter(width = 0.17, size = 3, alpha = 0.5, aes(color = f)) +
  xlab('') + 
  ylab('Prediction accuracies') +
  scale_fill_manual(values = c2) +
  scale_color_manual(values = c2) +
  scale_x_discrete(labels = c('Genotype\n+\nBatch\n+\nCRA', 'Genotype\n+\nBatch', 'Genotype\n+\nBatch\n+\nRandom CRA')) +
  guides(color = guide_legend(ncol = 3)) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'white', colour = 'black'), 
        legend.position = 'none', 
        legend.title = element_blank(),
        aspect.ratio = 1,
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.ticks = element_blank())
ggsave(filename = 'results/22_extendedFig5b.pdf', device = 'pdf', width = 6, height = 6)
