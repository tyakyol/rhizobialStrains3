library(magrittr)
library(dplyr)
library(vegan)
library(ggplot2)
library(RVenn)
library(RColorBrewer)

source('codes/__colors__.R')

# Load the data =================================
counts = readRDS('data/counts.rds')
soil = readRDS('data/noduleSizeAndSoil.rds')

# Relative abundance as int in 10e5 =============
dfRA = readRDS('data/newData/adjustedPlantCounts.rds')
st = counts$samples
rm(counts)

# Zero counts of the strains ====================
zeroCounts1 = apply(dfRA, 1, function(x) sum(x == 0))

# Abundance of the strains ======================
ava = apply(dfRA, 1, mean)

# Strain classification by abundance ============
qplot(ava, zeroCounts1, geom = 'blank') + 
  geom_point(size = 3, alpha = 0.7) +
  scale_x_continuous(trans = 'log2', breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000)) +
  xlab('Average abundance in log2') + ylab('Number of zeroes in 603 samples')

abundantAndCommon = names(ava)[ava > quantile(ava, 0.6) & zeroCounts1 < 302]
fewAndCommon = names(ava)[ava <= quantile(ava, 0.6) & zeroCounts1 < 302]
abundantAndSporadic = names(ava)[ava > quantile(ava, 0.6) & zeroCounts1 >= 302]
noise = setdiff(names(ava), c(abundantAndCommon, fewAndCommon, abundantAndSporadic))

# Sanity check (no overlap, sum to 603)
overlap_pairs(Venn(list(
  AC = abundantAndCommon,
  FC = fewAndCommon,
  AS = abundantAndSporadic,
  N = noise
)))
sum(length(abundantAndCommon), length(abundantAndSporadic), length(fewAndCommon), length(noise))

# A strain class vector for plotting
strainClass = rep(c('Dominants', 'Generalists', 'Specialists', 'Transients'), sapply(list(abundantAndCommon, fewAndCommon, abundantAndSporadic, noise), length))

# Adjust names of the parameters
ava1 = ava[c(abundantAndCommon, fewAndCommon, abundantAndSporadic, noise)]
zco1 = zeroCounts1[c(abundantAndCommon, fewAndCommon, abundantAndSporadic, noise)]

# qplot(log2(ava1+1), geom = 'blank') +
#   geom_density(color = NA, fill = colors_$turk, alpha = 0.5) +
#   xlab('') + 
#   ylab('') +
#   geom_vline(xintercept = median(log2(ava1+1))) +
#   geom_vline(xintercept = quantile(log2(ava1+1), 0.6)) +
#   theme(panel.grid = element_blank(),
#         axis.ticks = element_blank(), 
#         axis.text = element_text(size = 14),
#         axis.title = element_text(size = 18), 
#         legend.text = element_text(size = 12))

qplot(ava1, geom = 'blank') +
  geom_density(color = NA, fill = colors_$turk, alpha = 0.5) +
  xlab('Abundance') + 
  ylab('Density') +
  geom_vline(xintercept = median(ava), color = colors_$black, linetype = 'dashed', size = 1.25) +
  geom_vline(xintercept = quantile(ava, 0.6), color = colors_$red, linetype = 'dashed', size = 1.25) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 12))
ggsave(filename = 'results/extendedFig2a__.pdf', device = 'pdf', width = 7.5, height = 5.5)

qplot((603-zco1), geom = 'blank') +
  geom_density(color = NA, fill = colors_$turk, alpha = 0.5) +
  xlab('Occurrence') + 
  ylab('Density') +
  geom_vline(xintercept = 302, color = colors_$red, linetype = 'dashed', size = 1.25) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 12))
ggsave(filename = 'results/extendedFig2b__.pdf', device = 'pdf', width = 7.5, height = 5.5)

# Occupancy-abundance plot colored by arbitrary abundance thresholds
c1 = c(colors_$laci, colors_$yellow, colors_$turk, colors_$lightpurple)
qplot(ava1, (603-zco1), geom = 'blank') + 
  geom_point(aes(fill = strainClass), size = 5, shape = 21, stroke = 0.3) +
  scale_fill_manual(values = c1) +
  scale_x_continuous(trans = 'log2', breaks = c(0.001, 0.1, 10, 1000), labels = c(0.001, 0.1, 10, 1000)) +
  xlab('Abundance') + 
  ylab('Occurrence') +
  guides(fill=guide_legend(ncol = 2)) +
  theme(legend.position = 'none', legend.direction = 'vertical',
        legend.title = element_blank(), panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 12))
ggsave(filename = 'results/02_fig2b.pdf', device = 'pdf', width = 6, height = 6)

# Map to the soil data ==========================
table(soil[abundantAndCommon, 'soilName'])  # Mostly UG
table(soil[fewAndCommon, 'soilName'])  # SEJET and NOR
table(soil[abundantAndSporadic, 'soilName'])  # Mostly UG
table(soil[noise, 'soilName'])  # NOR, SEJET, UG

# Soil vector for plotting
soilForPlot = soil[c(abundantAndCommon, fewAndCommon, abundantAndSporadic, noise),'soilName']
soilForPlot = ifelse(soilForPlot %in% c('UG', 'SEJET', 'Nor'), soilForPlot, 'Others')
soilForPlot[soilForPlot == 'SEJET'] = 'SE'
soilForPlot[soilForPlot == 'Nor'] = 'NO'
soilForPlot = factor(soilForPlot, levels = c('UG', 'SE', 'NO', 'Others'))

# Occupancy-abundance plot colored by soil names
c2 = c(colors_$coffee, colors_$lightgreen, colors_$anotherred, colors_$white)
qplot(ava1, (603-zco1), geom = 'blank') + 
  geom_point(aes(fill = soilForPlot), size = 5, shape = 21, stroke = 0.3) +
  scale_fill_manual(values = c2) +
  scale_x_continuous(trans = 'log2', breaks = c(0.001, 0.1, 10, 1000), labels = c(0.001, 0.1, 10, 1000)) +
  xlab('Abundance') + 
  ylab('Occurrence') +
  guides(fill=guide_legend(ncol = 2)) +
  theme(legend.position = 'none', legend.direction = 'vertical',
        legend.title = element_blank(), panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 12))
ggsave(filename = 'results/02_fig2c.pdf', device = 'pdf', width = 6, height = 6)

# All information together
strainDF = data.frame(
  isolate = names(ava1),
  class = strainClass,
  averageAbundance = ava1,
  occupancyPer = ((603-zco1)/603)*100,
  soil = soilForPlot
)
write.csv(strainDF, file = 'results/competitiveness/table1.csv')

# Abundance-soil comparison
outTab_ = table(soilForPlot, strainClass)
outTab = as.matrix.data.frame(outTab_)
row.names(outTab) = row.names(outTab_)
colnames(outTab) = colnames(outTab_)
outTabPer = t(round(apply(outTab, 1, function(x) x/sum(x))*100))
outTabPer = t(apply(outTabPer, 1, function(x) paste0(x, '%')))
outTabPer = cbind(outTabPer, rowSums(outTab))
colnames(outTabPer) = c(colnames(outTab), 'Total number')
outTabPer = outTabPer[order(-as.integer(outTabPer[,5])), ]
write.csv(outTabPer, file = 'results/competitiveness/table2.csv')
