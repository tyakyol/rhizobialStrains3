library(ggplot2)

source('codes/__colors__.R')

c2 = c(colors_$laci, colors_$yellow, colors_$turk)

df = data.frame(
  term = c('Dominants', 'Generalists', 'Specialists'),
  estimate = c(0.0823, -0.0834, 0.0703), 
  ymin =     c(0.0166, -0.137, 0.00509), 
  ymax =     c(0.148, -0.03, 0.136),
  sign = c('*', '**', '*')
)

ggplot(df, aes(x = term, y = estimate,
               ymin = ymin, ymax = ymax)) +
  geom_linerange(size = 3.5, aes(color = term)) + 
  geom_point(size = 8.5, aes(color = term)) + 
  geom_text(aes(label = sign), nudge_x = -0.055, nudge_y = 0.068, size = 13) +
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
ggsave(filename = 'results/18_fig4e.pdf', device = 'pdf', width = 6, height = 6)
