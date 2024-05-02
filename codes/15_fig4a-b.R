library(magrittr)
library(dplyr)
library(lme4)
library(vegan)
library(ggplot2)
library(purrr)
library(scales)
library(seqtime)

source('codes/__colors__.R')

# Load the data =================================
counts = readRDS('data/counts.rds')
dfRA = readRDS('data/newData/adjustedPlantCounts.rds')
st = counts$samples
pgR = counts$samples
isoInfo = read.csv('results/competitiveness/table1.csv', row.names = 1)
grm = read.delim('results/grm.csv', sep = ',')
rm(counts)

# Abundance without transients ==================
isoToKeep = row.names(isoInfo)[isoInfo$class != 'Transients']
dfRA = dfRA[isoToKeep, ]

# MDS ===========================================
mp = cmdscale(vegdist(t(dfRA), method = 'cao'), k = 4)
mp = as.data.frame(mp)
mp$Plant = row.names(mp)

batch = substr(st$Batch, 5, 100)
  
qplot(mp$V1, mp$V2, geom = 'blank') +
  geom_point(size = 4, aes(color = batch)) +
  xlab('MDS 1') + 
  ylab('MDS 2') +
  scale_color_brewer(palette="Spectral") +
  guides(color = guide_legend(ncol = 4)) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'white', colour = 'black'), 
        legend.position = 'top', 
        legend.title = element_blank(),
        aspect.ratio = 1,
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.ticks = element_blank())
ggsave(filename = 'results/15_extendedFig4a.pdf', device = 'pdf', width = 6, height = 6)

qplot(mp$V3, mp$V4, geom = 'blank') +
  geom_point(size = 4, aes(color = batch)) +
  xlab('MDS 3') + 
  ylab('MDS 4') +
  scale_color_brewer(palette="Spectral") +
  guides(color = guide_legend(ncol = 4)) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'white', colour = 'black'), 
        legend.position = 'top', 
        legend.title = element_blank(),
        aspect.ratio = 1,
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.ticks = element_blank())
ggsave(filename = 'results/15_extendedFig4b.pdf', device = 'pdf', width = 6, height = 6)

# Merge with the growth data ====================
pgR = pgR %>% group_by(Line1) %>% mutate(repCount = sum(!is.na(Biomass)))
pgR = arrange(pgR, Line1, Ref)
pgR = pgR %>% group_by(Line1, Ref) %>% summarize(newBM = mean(Biomass))

dfForLmm = round(log2(dfRA+1), 6)
dfForLmm = scale(t(dfForLmm), center = T, scale = T)
dfForLmm = as.data.frame(t(dfForLmm))

lf = left_join(st, grm[, 1:3], by = c('Line1' = 'X'))
lf$Plant = paste0('ID_100', lf$Plant)
lf = left_join(lf, mp, by = 'Plant')

# Function to split every abundance to data frames with growth
LMM1 = function(r) {
  rh = data.frame(Var1 = t(dfForLmm[r, ]))
  rh$Plant = row.names(rh)
  rh = left_join(rh, lf, by = 'Plant')
  rh = rh[complete.cases(rh), ]
  return(rh)
}

# Apply the function above
inForLMM = list()
for(i in row.names(dfForLmm)) {
  inForLMM[[i]] = LMM1(i)
}

# For the data frames generated with LMM1(), apply growth analysis with the following for loop
br = list()
for(l in 1:length(inForLMM)) {
  df = inForLMM[[l]]
  colnames(df)[1] = 'Abundance'
  m = lmerTest::lmer(Biomass ~ Abundance + (1 | Line1) + (1 | Batch) + V2 + V3 + V4, data = df)
  mm = summary(m)
  br[[l]] = broom.mixed::tidy(m, conf.int = TRUE)
}

# Generate data frame for the abundance related coefficients from LMM and plot
coef_estimates = as.data.frame(matrix(nrow = length(inForLMM), ncol = 5))
coef_estimates$V1 = names(inForLMM)
coef_estimates$V2 = unlist(sapply(br, function(x) x[2, 4]))
coef_estimates$V3 = unlist(sapply(br, function(x) x[2, 8]))
coef_estimates$V4 = unlist(sapply(br, function(x) x[2, 9]))
coef_estimates$V5 = unlist(sapply(br, function(x) x[2, 10]))
colnames(coef_estimates) = c('term', 'estimate', 'pval', 'conf.low', 'conf.high')
coef_estimates$term = strsplit(split = '=', x = coef_estimates$term) %>% sapply(function(x) x[[1]])
coef_estimates$sign = ifelse(coef_estimates$pval < 1e-3, 'Significant', 'NS')
coef_estimates$adjp = p.adjust(coef_estimates$pval, 'BH')
coef_estimates$newSign = ifelse(coef_estimates$adjp < 0.05, 'Significant', 'NS')

for(i in 1:216) {
  if(coef_estimates$newSign[i] == 'NS') {
    coef_estimates$newSign_[i] = 'NS'
  }  
  if(coef_estimates$newSign[i] == 'Significant' & coef_estimates$estimate[i] > 0) {
    coef_estimates$newSign_[i] = 'Positive'
  }
  if(coef_estimates$newSign[i] == 'Significant' & coef_estimates$estimate[i] < 0) {
    coef_estimates$newSign_[i] = 'Negative'
  }
}
coef_estimates$newSign__ = ifelse(coef_estimates$newSign_ == 'Negative', '*', '')
coef_estimates = arrange(coef_estimates, estimate)
coef_estimates$term = factor(coef_estimates$term, levels = unique(coef_estimates$term))
saveRDS(coef_estimates, 'results/lmmWithMDS.rds')

# Plots for the abundance coefficients ==========
# Plot-1
ggplot(coef_estimates, aes(x = term, y = estimate,
                           ymin = conf.low, ymax = conf.high)) +
  geom_linerange(color = colors_$blue, size = 1.7, alpha = 0.4) +
  geom_point(color = colors_$blue, size = 5, alpha = 0.6) +
  geom_text(aes(label = newSign__), nudge_x = -3.5, size = 10) +
  coord_flip() +
  xlab('Rhizobium strains') +
  ylab('Effect on plant growth') +
  geom_hline(yintercept = 0, color = colors_$black, size = 2,
             linetype = 'dashed') +
  guides(color = guide_legend(ncol = 3)) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        legend.position = 'top',
        legend.title = element_blank(),
        aspect.ratio = 1,
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.ticks = element_blank(),
        axis.text.y = element_blank())
ggsave(filename = 'results/15_fig4a.pdf', device = 'pdf', width = 6, height = 6)

# Plot-2
isoInfo_ = isoInfo[isoInfo$class != 'Transients',]
tabOut3 = left_join(isoInfo_, coef_estimates, by = c('isolate' = 'term'))
# tabOut3 %>% group_by(class) %>% summarise(mean(estimate))
# tabOut3 %>% group_by(class) %>% reframe(range(estimate))

c1 = c(colors_$laci, colors_$yellow, '#66666655', colors_$turk)
c2 = c1[-3]

ggplot(tabOut3, aes(x = class, y = estimate)) +
  geom_boxplot(outlier.shape = NA, size = 2, aes(color = class)) +
  geom_jitter(width = 0.17, size = 3, alpha = 0.5, aes(color = class)) +
  xlab('') + 
  ylab('Effect on plant growth') +
  scale_fill_manual(values = c2) +
  scale_color_manual(values = c2) +
  guides(color = guide_legend(ncol = 3)) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'white', colour = 'black'), 
        legend.position = 'none', 
        legend.title = element_blank(),
        aspect.ratio = 1,
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.ticks = element_blank())
ggsave(filename = 'results/15_fig4b.pdf', device = 'pdf', width = 6, height = 6)
