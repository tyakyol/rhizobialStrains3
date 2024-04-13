library(magrittr)
library(dplyr)
library(vegan)
library(ggplot2)
library(purrr)
library(scales)
library(RColorBrewer)
library(pheatmap)

source('codes/__colors__.R')

# Load the data =================================
counts = readRDS('data/counts.rds')
dfRA = readRDS('data/newData/adjustedPlantCounts.rds')
st = counts$samples
pgR = counts$samples
rm(counts)
isoInfo = read.csv('results/competitiveness/table1.csv', row.names = 1)

# Average values
raav = as.data.frame(t(dfRA))
raav = raav %>% group_by(st$Line1) %>% summarise_all(mean)
colnames(raav)[1] = 'Genotype'
class(raav) = class(raav)[3]
row.names(raav) = raav$Genotype
raav = raav[,-1]

# Heatmaps ======================================
# Annotation colors
colForAnn = list(soil = c('4Others' = colors_$white,
                           '2NO' = colors_$anotherred,
                           '3SE' = colors_$lightgreen,
                           '1UG' = colors_$coffee),
                 class = c('Dominants' = colors_$laci,
                           'Specialists' = colors_$turk,
                           'Generalists' = colors_$yellow,
                           'Transients' = colors_$lightpurple)
                 )

isoInfo$soil[isoInfo$soil == 'UG'] = '1UG'
isoInfo$soil[isoInfo$soil == 'NO'] = '2NO'
isoInfo$soil[isoInfo$soil == 'SE'] = '3SE'
isoInfo$soil[isoInfo$soil == 'Others'] = '4Others'

# Heatmap colors
colForHM = colorRampPalette(c(colors_$gray, 'white', 'red2'))(n = 100)

# Heatmap without dendrograms and annotations combined
isoInfo2 = isoInfo
isoInfo2 = arrange(isoInfo2, class, soil)
raav2 = raav[,row.names(isoInfo2)]

pheatmap(t(log2(raav2+1)), border_color = NA, scale = 'none', cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE, annotation_row = isoInfo2[, c('soil', 'class')], gaps_row = cumsum(table(isoInfo$class))[1:3], annotation_colors = colForAnn, legend_breaks = c(0, 4, 8, 12, 16), color = colForHM, annotation_names_row = FALSE, filename = 'results/01_fig2a.pdf', width = 8, height = 7, annotation_legend = TRUE, legend = TRUE)
