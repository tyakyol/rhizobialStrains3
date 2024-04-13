library(magrittr)
library(dplyr)
library(vegan)
library(ggplot2)
library(purrr)
library(scales)
library(igraph)
library(pheatmap)
library(psych)

source('codes/__colors__.R')

# Load the data =================================
counts = readRDS('data/counts.rds')
st = counts$samples
pgR = counts$samples
isoInfo = read.csv('results/competitiveness/table1.csv', row.names = 1)
grm = read.delim('results/grm.csv', sep = ',')
rm(counts)

# Relative abundance as int in 10e5 =============
isoToKeep = row.names(isoInfo)[isoInfo$class == 'Dominants']
dfRA = readRDS('data/newData/adjustedPlantCounts.rds')
dfRA = dfRA[isoToKeep, ]
dfRA[dfRA > 0] = 1

# Again
df__ = readRDS('data/newData/adjustedPlantCounts.rds')
df__ = df__[isoToKeep, ]

# Interaction function
funcForTab = function(i, n) {
  
  i1 = as.integer(dfRA[i,])
  i2 = as.integer(dfRA[n,])
  
  t1 = table(i1, i2)
  # a1 = (t1[2,2]) / ((abs(t1[2,1] + t1[1,2]))^2)
  a1 = (t1[2,2]) / ((abs(t1[2,1] + t1[1,2])))
  return(a1)
}

# Compute for all pairs
andFunctionForAll = function(isolate) {
  x = which(row.names(dfRA) == isolate)
  o1 = c()
  for(r in setdiff(seq_len(nrow(dfRA)), x)) {
    o1[row.names(dfRA)[r]] = funcForTab(isolate, r)
  }
  return(o1)
}

outForComp = list()
for(isol in row.names(dfRA)) {
  outForComp[[isol]] = andFunctionForAll(isol)
}

# Convert to data.frame
df = data.frame(num = unlist(outForComp),
                bac = rep(names(outForComp), each = 85))
# df = df[!df$num %in% c(-Inf, Inf),]

net = df
net$to = strsplit(row.names(net), split = '\\.') %>% sapply(function(x) x[[2]])
colnames(net)[2] = 'from'
colnames(net)[1] = 'weight'
net = net[, c(2, 3, 1)]

# The graph object
gr = graph_from_data_frame(d = net, vertices = row.names(dfRA), directed = FALSE)
gr = simplify(gr, edge.attr.comb = 'max')

df2 = igraph::as_data_frame(gr)
th = quantile(df2$weight, 0.9)
rr = range(df2$weight)

qplot(df2$weight, geom = 'blank') + 
  geom_vline(xintercept = 2.72, linetype = 'dashed',
             color = colors_$red, alpha = 0.9, size = 1.3) +
  geom_histogram(bins = 500, fill = colors_$blue, alpha = 0.7) +
  scale_x_continuous(name = '', limits = c(0, 10)) + 
  scale_y_continuous('Frequency') +
  theme(legend.position = 'none',
        panel.background = element_rect(fill = '#66666622'),
        panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 12))
ggsave('results/extendedFig3a.pdf', height = 5, width = 7)

# Only strong and positive interactions 
df3 = df2[df2$weight > 2.72, ]

# The graph object
gr = graph_from_data_frame(d = df3, vertices = row.names(dfRA), directed = FALSE)

# Number of connections
d1 = degree(gr)

# Colors
c1 = isoInfo[names(d1), 'soil']
c1[c1 == 'NO'] = colors_$anotherred
c1[c1 == 'SE'] = colors_$lightgreen
c1[c1 == 'UG'] = colors_$coffee
c1[c1 == 'Others'] = colors_$white

pdf(file = 'results/04_fig3a.pdf', useDingbats = F, width = 5, height = 5)
par(mar = c(0,0,0,0))
plot(gr,
     edge.color = adjustcolor('black', alpha.f = 0.2),
     edge.width = 1.2,
     layout = layout_in_circle(gr, order = order(c1)),
     vertex.color = adjustcolor(c1, alpha.f = 0.6),
     vertex.frame.color = 'black',
     vertex.shape = 'circle',
     vertex.size = (d1/3)+2,
     vertex.label = NA,
     edge.curved = 0.2)
dev.off()

# Correlation
df__ = log2(df__ + 1)
ct = corr.test(t(df__), adjust = 'BH', ci = F)
g_ = graph_from_adjacency_matrix(ct$r, mode = 'upper', weighted = TRUE, diag = F)
g_ = as_long_data_frame(g_)
summary(g_$weight)

pdf(file = 'results/04_fig3c.pdf', useDingbats = F, width = 6, height = 6)
qplot(g_$weight, geom = 'blank') + 
  geom_histogram(bins = 40, fill = colors_$anotherpurple, alpha = 0.7) + 
  scale_x_continuous(limits = c(-1, 1), 'Pearson correlation') + 
  scale_y_continuous(limits = c(0, 1000), 'Frequency') +
  theme(legend.position = 'none',
        panel.background = element_rect(fill = '#66666622'),
        panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 12))
dev.off()
