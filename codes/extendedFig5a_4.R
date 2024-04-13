library(ggplot2)

source('codes/__colors__.R')

c2 = c(colors_$laci, colors_$yellow, colors_$turk)

df = data.frame(
  term = c('Dominants', 'Generalists', 'Specialists'),
  estimate = c(-0.0294, 0.0290, -0.0267), 
  ymin =     c(-0.0729, -0.0169, -0.0704), 
  ymax =     c(0.0141, 0.0749, 0.0170),
  sign = c('', '', '')
)

ggplot(df, aes(x = term, y = estimate,
               ymin = ymin, ymax = ymax)) +
  geom_linerange(size = 3.5, aes(color = term)) + 
  geom_point(size = 8.5, aes(color = term)) + 
  # geom_text(aes(label = sign), nudge_x = -0.055, nudge_y = 0.068, size = 13) +
  coord_flip() + 
  ylab('Effect on plant growth') + 
  xlab('') +
  scale_color_manual(values = c2) +
  guides(color = guide_legend(ncol = 3)) +
  theme(text = element_text(size = 14),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        legend.position = 'none', 
        legend.title = element_blank(),
        aspect.ratio = 1,
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.ticks = element_blank())
ggsave(filename = 'results/extendedFig5a.pdf', device = 'pdf', width = 6, height = 6)
