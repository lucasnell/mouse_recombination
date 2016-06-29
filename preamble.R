library('ggplot2')
library('reshape2')
library('dplyr')
library('grid')
library('boot')

options(stringsAsFactors = FALSE)


plotTheme <- function(base_size = 10, base_family = 'Helvetica') {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      strip.text = element_text(face = 'bold'),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray50"),
      axis.ticks = element_line(color = "gray50"),
      axis.ticks.length = unit(2, 'points'),
      legend.position = 'none'
    )
}


# Loading `subunit` data frame
load(file = 'subunit.RData')

