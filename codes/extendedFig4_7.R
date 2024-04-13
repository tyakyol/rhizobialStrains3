library(ggplot2)
library(agricolae)

shannon = readRDS('results/growthPredictionWithDiversity/shannon2.rds')
fisher_ = readRDS('results/growthPredictionWithDiversity/fisher2.rds')
evennes = readRDS('results/growthPredictionWithDiversity/sheldon2.rds')
betadiv = readRDS('results/growthPredictionWithDiversity/betaDiversity2.rds')
richnes = readRDS('results/growthPredictionWithDiversity/speciesCount2.rds')
simpson = readRDS('results/growthPredictionWithDiversity/simpson2.rds')
reduced = readRDS('results/growthPredictionWithDiversity/reduced2.rds')
cumuabu = readRDS('results/growthPredictionWithDiversity/abundance2.rds')

source('codes/__colors__.R')

x = c(reduced, betadiv, cumuabu, richnes, shannon, simpson, evennes, fisher_)
f = rep(letters[1:8], each = 100)

HSD.test(lm(x ~ f), console = TRUE, trt = 'f')
pairwise.t.test(x = x, g = f, p.adjust.method = 'none')

c2 = c(colors_$gray,colors_$gray,colors_$gray,colors_$gray, colors_$pink, colors_$gray, colors_$pink, colors_$gray)

qplot(x = f, y = x, geom = 'blank') + 
  geom_boxplot(size = 0.7, aes(color = f), outlier.colour = NA) +
  geom_jitter(width = 0.17, size = 2, alpha = 0.25, aes(color = f)) +
  stat_summary(fun.y = mean, geom = "point", shape = '-', size = 36, color = colors_$blue) +
  geom_hline(yintercept = mean(reduced), size = 2, alpha = 0.8, linetype = 'dashed') +
  xlab('') + 
  ylab('Prediction accuracy') +
  scale_fill_manual(values = c2) +
  scale_color_manual(values = c2) +
  scale_x_discrete(labels = c('Null', 'Beta\ndiversity', 'Cumulative\nabundance', 'Richness', 'Shannon', 'Simpson', 'Evenness', 'Fisher')) +
  guides(color = guide_legend(ncol = 3)) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'white', colour = 'black'), 
        legend.position = 'none', 
        legend.title = element_blank(),
        aspect.ratio = 0.4,
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.ticks = element_blank())
ggsave(filename = 'results/extendedFig4.pdf', device = 'pdf', width = 9, height = 5)
