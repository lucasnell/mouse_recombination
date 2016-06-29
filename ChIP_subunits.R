## @knitr packages
library('ggplot2')
library('dplyr')
library('tidyr')
library('readr')
library('readxl')


## @knitr ggplot_theme
plot_theme <- function(base_size = 10, base_family = 'Helvetica') {
    theme_minimal(base_size = base_size, base_family = base_family) %+replace%
        theme(
            strip.text = element_text(face = 'bold'),
            panel.grid = element_blank(),
            panel.border = element_rect(fill = NA, color = "gray50"),
            axis.ticks = element_line(color = "gray50"),
            axis.ticks.length = unit(2, 'points')
        )
}

## @knitr load_RData
load(file = 'ChIP_subunits.RData')



## @knitr input_peaks
one_peak <- function(sample_name, subunit_data_frame) {
    
    if (! sample_name %in% c('B6', '9R', '13R')){
        stop('Please input a valid sample_name: "B6", "9R", or "13R".')
    }
    
    # (Paths to *.bsites files should work on any Mac connected to MW_18TB)
    peak_df <- paste0('/Volumes/MW_18TB/Lucas_Nell/lan/musChIP/bed/', 
                      sample_name, '_Y_Dq_Pv01.bsites') %>%
        read_tsv(., skip = 57, comment = '=', 
                 col_names = c('chr', 'start', 'end', 'tags', 'fold', 'p')) %>%
        mutate(sample = sample_name) %>%
        # Filter for peaks within MSY long arm:
        filter(start >= min(subunit_data_frame$start), 
               end <= max(subunit_data_frame$end)) %>%
        select(sample, start, end, tags, fold, p)
    
    return(peak_df)
}

peak_df <- lapply(c('9R', '13R', 'B6'), one_peak, subunit_data_frame = subunit_df) %>% 
    bind_rows




## @knitr comparison_plot
subunit_df %>%
    ggplot(aes(x = start/1e6)) + 
    plot_theme() +
    geom_histogram(aes(y = ..density..), binwidth = 1.000000, 
                   fill = 'gray80', color = NA) + 
    geom_rug(data = peak_df, aes(color = sample), size = 0.75, sides = 'b') +
    ylab('Density of subunits') +
    xlab('Location on chromosome (Mb)') +
    scale_color_brewer(type = "qual", palette = 'Dark2')





## @knitr over_counts
over_counts <- function(centers, flank_size = 1e6) {
    overlaps <- sapply(centers, function(x) {
        nrow(subunit_df[subunit_df$start <= (x + flank_size) & 
                            subunit_df$end >= (x - flank_size),])
        }, USE.NAMES = FALSE)
    return(overlaps)
}






## @knitr compute_overs
peak_df <- peak_df %>% 
    mutate(center = (start + end) / 2,
           overs = over_counts(center),
           sample = factor(sample, levels = c('B6', '9R', '13R')))

summary_df <- peak_df %>% 
    group_by(sample) %>% 
    summarize(n = length(overs),
              sd = sd(overs), 
              mean = mean(overs))
summary_df



## @knitr input_boot_results
# save(pos_boot_matrix, file = 'pos_boot_matrix.RData', compress = TRUE)
load('pos_boot_matrix.RData')




