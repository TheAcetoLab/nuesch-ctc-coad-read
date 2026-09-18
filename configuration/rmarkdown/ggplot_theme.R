require(ggplot2)
require(cowplot)
use_font_size <- 5
# Size 1pt in illustrator: https://stackoverflow.com/questions/47519624/how-is-the-line-width-size-defined-in-ggplot2
one_pt <- 1/2.141959
linewidth_default <- one_pt/2
theme_project <- theme_half_open(
  font_size = use_font_size,
  rel_small = use_font_size/use_font_size,
  rel_tiny = use_font_size/use_font_size,
  rel_large = use_font_size/use_font_size) +
  theme (
    plot.title = element_text(size = use_font_size, face = 'plain'),
    plot.caption = element_text(size = use_font_size),

    # axis.line = element_line(linewidth = rel(0.5)),
    # Axis 1pt in illustrator: https://stackoverflow.com/questions/47519624/how-is-the-line-width-size-defined-in-ggplot2
    axis.line = element_line(linewidth = linewidth_default, color = 'black'),
    axis.ticks = element_line(linewidth = linewidth_default, color = 'black'),

    text = element_text(size = use_font_size),
    legend.text = element_text(size = use_font_size),
    legend.title = element_text(size = use_font_size),

    legend.key.width = unit(0.1, "cm"),

    strip.background = element_blank()
  )

theme_set(theme_project)

theme_mimimal_ygrid <- theme(
  strip.text.x = element_text(size=use_font_size, hjust = 0),
  strip.background = element_blank(),
  panel.grid.major.y = element_line(linewidth = (0.2), colour="grey")
)

theme_facet <-
  theme(
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = linewidth_default/2, color = 'black'),
    panel.border = element_rect(fill =NULL,
                                color = "black",
                                # double the size : https://github.com/tidyverse/ggplot2/issues/5382
                                linewidth = linewidth_default)
  )

geom_text_size <- use_font_size * 0.352777778 # https://stackoverflow.com/questions/17311917/ggplot2-the-unit-of-size
