library(magrittr)
library(dplyr)
library(vegan)
library(ggplot2)
library(pheatmap)
library(RVenn)

source('codes/__colors__.R')

# Information from previous results
isolateInfo = read.delim('data/competitivenessFromPreviousAnalysis.csv', sep = ',', row.names = 1, header = TRUE)

# Count table reshaping =========================
df = read.delim('data/inoculumCounts.csv', sep = ',')
df = df[df$pl_name != 'none', ]
df = df[, c(6, 3, 4)]
rdf = reshape(data = as.data.frame(df), idvar = 'pl_name', timevar = 'ID', direction = 'wide')
rdf[is.na(rdf)] = 0
row.names(rdf) = rdf$pl_name
rdf = rdf[, -1]
colnames(rdf) = colnames(rdf) %>% strsplit(split = 'ID_') %>% sapply(function(x) x[[2]])
dfRA = round(decostand(rdf, method = 'total', MARGIN = 2) * 100000)

# Some summary statistics
rawIso = rowMeans(rdf)
relIso = rowMeans(dfRA)
rawIno = colSums(rdf)

# Relative abundance and raw counts correlates perfectly
qplot(rawIso, relIso) + geom_smooth(method = 'lm', se = F)

# Align the isolates
ggvenn(Venn(list(P = row.names(isolateInfo), I = names(relIso)))) # 397 are matching
ggsave('results/inoculum/isolateMatch.pdf')

iii = overlap(Venn(list(row.names(isolateInfo), names(relIso))))
relIso2 = relIso[iii]
relIso2 = sort(relIso2, decreasing = TRUE)
isolateInfo2 = isolateInfo[names(relIso2),]
dfRA2 = dfRA[row.names(isolateInfo2), ]

# Heatmaps ======================================
# Annotation colors
colForAnn = list(
                 class = c('Dominants' = colors_$laci,
                           'Specialists' = colors_$turk,
                           'Generalists' = colors_$yellow,
                           'Transients' = colors_$lightpurple)
)

# Heatmap colors
colForHM = colorRampPalette(c(colors_$black, 'white', 'red2'))(n = 100)

# pheatmap(log2(dfRA2 + 1), border_color = NA, scale = 'none', cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = TRUE, annotation_row = isolateInfo[, 'class', drop = FALSE], annotation_colors = colForAnn, color = colForHM, annotation_names_row = FALSE, 
#          filename = 'results/inoculum/heatmap.pdf', width = 8, height = 7,
#          annotation_legend = TRUE, legend = TRUE
#          )

pheatmap(log2(dfRA2 + 1), border_color = NA, scale = 'none', cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = TRUE, color = colForHM, 
         filename = 'results/inoculum/heatmap2.pdf', width = 7, height = 7,
         annotation_legend = TRUE, legend = TRUE)

# Barplot
qplot(y = relIso2, x = 1:397, geom = 'blank') +
  geom_col(color = colors_$blue, fill = colors_$lightblue) +
  xlab('Isolates') + 
  ylab('Average abundance') +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 12))
ggsave(filename = 'results/inoculum/barplot.pdf', device = 'pdf', width = 6, height = 6)

# Boxplot
qplot(y = relIso2, x = isolateInfo2$class, geom = 'blank') +
  geom_boxplot() +
  xlab('') + 
  ylab('Average abundance') +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 12))
ggsave(filename = 'results/inoculum/boxplot.pdf', device = 'pdf', width = 6, height = 6)

# Density plot
foo = reshape(data = dfRA, idvar = 'ID', ids = row.names(dfRA), time = names(dfRA), timevar = 'isolate', direction = 'long', varying = list(names(dfRA)))
colnames(foo)[2] = 'abundance'

ggplot(data = foo, aes(x = abundance, group = isolate)) +
  geom_density(color = '#32a85288', linewidth = 1) +
  xlab('Abundance') + 
  ylab('Density') +
  theme(legend.position = 'none', legend.direction = 'vertical',
        legend.title = element_blank(), panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 12))
ggsave(filename = 'results/inoculum/density.pdf', device = 'pdf', width = 6, height = 6)


# Correlation with the previous colonisation results
c1 = c(colors_$laci, colors_$yellow, colors_$turk, colors_$lightpurple)

qplot(relIso2, isolateInfo2$averageAbundance, geom = 'blank') +
  geom_point(size = 5, shape = 21, aes(fill = isolateInfo2$class)) +
  scale_fill_manual(values = c1) +
  xlab('Inoculum abundance') + 
  ylab('Previous abundance') +
  guides(fill=guide_legend(ncol = 2)) +
  theme(legend.position = 'none', legend.direction = 'vertical',
        legend.title = element_blank(), panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 12))
ggsave(filename = 'results/inoculum/previousAbundance.pdf', device = 'pdf', width = 6, height = 6)

qplot(relIso2, isolateInfo2$occupancyPer, geom = 'blank') +
  geom_point(size = 5, shape = 21, aes(fill = isolateInfo2$class)) +
  scale_fill_manual(values = c1) +
  xlab('Inoculum abundance') + 
  ylab('Previous occupancy') +
  guides(fill=guide_legend(ncol = 2)) +
  theme(legend.position = 'none', legend.direction = 'vertical',
        legend.title = element_blank(), panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 12))
ggsave(filename = 'results/inoculum/previousOccupancy.pdf', device = 'pdf', width = 6, height = 6)


# Dissimilarity
# The ones more dissimilar than expected by chance.
monteCarlo = function(d) {
  
  sss = combn(1:ncol(d), m = 2, FUN = c, simplify = F)
  out = c()
  
  for(s in 1:length(sss)) {
    first_ = sss[[s]][1]
    second = sss[[s]][2]
    
    real = vegdist(t(data.frame(a1 = t(d)[first_,], a2 = t(d)[second,])), method = 'cao')
    null = c()
    for(i in 1:1000) {
      null[i] = (vegdist(t(data.frame(a1 = t(dfRA)[first_,], a2 = sample(t(dfRA)[second,]))), method = 'cao'))
    }
    out[s] = 1 - ((sum(real < null)) / 1000)
  }
  
  return(out)
}

pv = monteCarlo(dfRA)
pva = p.adjust(pv, method = 'bonferroni')

ddd = as.matrix(vegdist(t(dfRA2), method = 'cao'))
ddd = reshape(data = as.data.frame(ddd), idvar = 'ID', ids = row.names(ddd), time = colnames(ddd), timevar = 'distance', direction = 'long', varying = list(colnames(ddd)))
ddd$sign = ''
ddd['Inoculum_3.Inoculum_8', 4] = '*'
ddd['Inoculum_8.Inoculum_3', 4] = '*'
ddd['Inoculum_4.Inoculum_8', 4] = '*'
ddd['Inoculum_8.Inoculum_4', 4] = '*'
ddd[,4][ddd$Inoculum_1 == 0] = ''

ggplot(data = ddd, aes(x = distance, y = ID, fill = Inoculum_1)) +
  geom_tile(show.legend = F) +
  geom_text(aes(label = sign), size = 10, color = 'white') +
  xlab('') + 
  ylab('') +
  theme(legend.position = 'none', legend.direction = 'vertical',
        legend.title = element_blank(), panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 12))
ggsave(filename = 'results/inoculum/dissimilarity.pdf', device = 'pdf', width = 10, height = 6)
