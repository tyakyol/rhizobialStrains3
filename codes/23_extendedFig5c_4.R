library(ggplot2)
library(agricolae)

source('codes/__colors__.R')

# Loading the results
dom = read.delim('results/genomicPredictionForAbundance/dominantScores.txt', header = F)
gen = read.delim('results/genomicPredictionForAbundance/generalistScores.txt', header = F)
spe = read.delim('results/genomicPredictionForAbundance/specialistScores.txt', header = F)

dom2 = read.delim('results/genomicPredictionForAbundance/random/dominantScores.txt', header = F)
gen2 = read.delim('results/genomicPredictionForAbundance/random/generalistScores.txt', header = F)
spe2 = read.delim('results/genomicPredictionForAbundance/random/specialistScores.txt', header = F)

# Comparison
summary(dom$V1)
summary(gen$V1)
summary(spe$V1)

x = c(dom$V1, gen$V1, spe$V1)
f = rep(c('Dominants', 'Generalists', 'Specialists'), each = 100)

HSD.test(lm(x ~ f), console = TRUE, trt = 'f')

# Plot
df = data.frame(Score = c(x, dom2$V1, gen2$V1, spe2$V1), 
                Isolate = c(f, f), 
                Marker = as.character(rep(1:2, each = 300)))

c2 = c(colors_$laci, colors_$yellow, colors_$turk)

ggplot(data = df[1:300,], aes(x = Isolate, y = Score)) + 
  geom_boxplot(size = 2, aes(color = Isolate), outlier.colour = NA) +
  geom_jitter(width = 0.17, size = 3, alpha = 0.5, aes(color = Isolate)) +
  xlab('') + 
  ylab('Prediction scores') +
  scale_fill_manual(values = c2) +
  scale_color_manual(values = c2) +
  geom_boxplot(data = df[301:600,], aes(y = Score, x = Isolate), outlier.colour = NA, 
               size = 2, color = 'gray') + 
  geom_jitter(data = df[301:600,], aes(y = Score, x = Isolate), width = 0.2, color = 'gray',
              size = 3, alpha = 0.5, ) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'white', colour = 'black'), 
        legend.position = 'none', 
        legend.title = element_blank(),
        aspect.ratio = 1,
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.ticks = element_blank())
ggsave(filename = 'results/23_extendedFig5c.pdf', device = 'pdf', width = 6, height = 6)
