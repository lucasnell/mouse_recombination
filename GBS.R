
library('ggplot2')
library('tidyr')
library('dplyr')
library('readr')

# library('grid')
# library('boot')

options(stringsAsFactors = FALSE)


plotTheme <- function(base_size = 10, base_family = 'Helvetica') {
    theme_minimal(base_size = base_size, base_family = base_family) %+replace%
        theme(
            strip.text = element_text(face = 'bold'),
            panel.grid = element_blank(),
            panel.border = element_rect(fill = NA, color = "gray50"),
            axis.ticks = element_line(color = "gray50"),
            axis.ticks.length = unit(2, 'points')
        )
}


readBED <- function(suffix){
    df <- read_tsv(paste0('/Volumes/MW_18TB/Lucas_Nell/lan/musGBS/bam_mapped_comb/', 
                          'run_132571', suffix, '.bed'),
                   col_names = c('chrom', 'start', 'end', 'reads')) %>%
        mutate(start = start + 1) %>% # Turning into 1-based indexing
        gather('pos_type', 'pos', start, end) %>%
        arrange(chrom, pos)
    return(df)
}

# Forward strand, MAPQ >= 20, "properly paired"
r1q <- readBED('_1_qPP') %>% 
    mutate(strand = 1)
# Reverse strand, MAPQ >= 20, "properly paired"
r2q <- readBED('_2_qPP') %>% 
    mutate(strand = 2)

rqPP <- bind_rows(r1q, r2q) %>% 
    mutate(strand = factor(strand, levels = c(1, 2), labels = c('+', '−')))

# Forward strand, MAPQ >= 20
r1q <- readBED('_1_q') %>% 
    mutate(strand = 1)
# Reverse strand, MAPQ >= 20
r2q <- readBED('_2_q') %>% 
    mutate(strand = 2)

rq <- bind_rows(r1q, r2q) %>% 
    mutate(strand = factor(strand, levels = c(1, 2), labels = c('+', '−')))

rm(r1q, r2q)



# # Both reads, MAPQ >= 20
# rq <- readBED('_q')

# Filtering for reads < `threshold`
threshold <- 96 * 1
rqPP_filt <- rqPP %>%
    mutate(reads = sapply(reads, function(x){ifelse(x > threshold, x, 0)}))
rq_filt <- rq %>%
    mutate(reads = sapply(reads, function(x){ifelse(x > threshold, x, 0)}))


rqPP %>%
    filter(chrom == 'chrX', pos >= (166410000 - 1e6)) %>%
    ggplot(aes(x = pos/1e6, y = reads, color = strand)) + 
    geom_vline(xintercept = 166410000/1e6, linetype = 3, color = 'gray60') +
    geom_line() + 
    scale_color_brewer(type = 'qual', palette = 'Set1') + 
    ylab('Reads') +
    xlab('Position on chromosome (Mb)') +
    plotTheme() + 
    theme(legend.direction = 'horizontal', legend.title = element_blank(),
          legend.position = c(0.5, 0.90))

rq %>%
    filter(chrom == 'chrX', pos >= (166410000 - 1e6)) %>%
    ggplot(aes(x = pos/1e6, y = reads, color = strand)) + 
    geom_vline(xintercept = 166410000/1e6, linetype = 3, color = 'gray60') +
    geom_line() + 
    scale_color_brewer(type = 'qual', palette = 'Set1') + 
    ylab('Reads') +
    xlab('Position on chromosome (Mb)') +
    plotTheme() + 
    theme(legend.direction = 'horizontal', legend.title = element_blank(),
          legend.position = c(0.5, 0.90))




rq %>%
    filter(chrom == 'chrX', pos >= 166.434e6, pos <= 166.436e6) %>%
    ggplot(aes(x = pos, y = reads, color = strand)) + 
    geom_line() + 
    scale_color_brewer(type = 'qual', palette = 'Set1') + 
    facet_grid(strand ~ .) + 
    theme(legend.position = 'none')






rqPP %>%
    filter(chrom == 'chrX', reads > 96 * 1) %>%
    ggplot(aes(x = reads, y = ..count..)) + 
    geom_density()





