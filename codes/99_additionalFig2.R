library(dplyr)
library(RVenn)
library(ggplot2)

# Dominants =====================================
dom1 = read.delim('results/manhattanDiversity/dominants.csv', sep = ',', row.names = 1)
dom2 = read.delim('results/differentMembership/manhattanDiversity/dominants.csv', sep = ',', row.names = 1)

# Venn diagram
ggvenn(Venn(list('Grouping 1' = dom1$snp, 'Grouping 2' = dom2$snp)))
ggsave(filename = 'results/99_additionalFig2a.pdf', device = 'pdf', width = 6, height = 6)

# Correlation
mixDom = inner_join(dom1, dom2, by = 'snp')
qplot(mixDom$pv.x, mixDom$pv.y) + xlab('Number of runs (grouping 1)') + ylab('Number of runs (grouping 2)') + theme(text = element_text(size = 20))
ggsave(filename = 'results/99_additionalFig2b.pdf', device = 'pdf', width = 6, height = 6)
cor.test(mixDom$pv.x, mixDom$pv.y)

# Generalists ===================================
gen1 = read.delim('results/manhattanDiversity/generalists.csv', sep = ',', row.names = 1)
gen2 = read.delim('results/differentMembership/manhattanDiversity/generalists.csv', sep = ',', row.names = 1)

# Venn diagram
ggvenn(Venn(list('Grouping 1' = gen1$snp, 'Grouping 2' = gen2$snp)))
ggsave(filename = 'results/99_additionalFig2c.pdf', device = 'pdf', width = 6, height = 6)

# Correlation
mixGen = inner_join(gen1, gen2, by = 'snp')
qplot(mixGen$pv.x, mixGen$pv.y) + xlab('Number of runs (grouping 1)') + ylab('Number of runs (grouping 2)') + theme(text = element_text(size = 20))
ggsave(filename = 'results/99_additionalFig2d.pdf', device = 'pdf', width = 6, height = 6)
cor.test(mixGen$pv.x, mixGen$pv.y)

# Specialists ===================================
spe1 = read.delim('results/manhattanDiversity/specialists.csv', sep = ',', row.names = 1)
spe2 = read.delim('results/differentMembership/manhattanDiversity/specialists.csv', sep = ',', row.names = 1)

# Venn diagram
ggvenn(Venn(list('Grouping 1' = spe1$snp, 'Grouping 2' = spe2$snp)))
ggsave(filename = 'results/99_additionalFig2e.pdf', device = 'pdf', width = 6, height = 6)

# Correlation
mixSpe = inner_join(spe1, spe2, by = 'snp')
qplot(mixSpe$pv.x, mixSpe$pv.y) + xlab('Number of runs (grouping 1)') + ylab('Number of runs (grouping 2)') + theme(text = element_text(size = 20))
ggsave(filename = 'results/99_additionalFig2f.pdf', device = 'pdf', width = 6, height = 6)
cor.test(mixSpe$pv.x, mixSpe$pv.y)
