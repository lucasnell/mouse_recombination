
library('ggplot2')
library('tidyr')
library('dplyr')
library('readr')
library('scales')
library('grid')
# library('boot')



plot_theme <- function(base_size = 10, base_family = 'Helvetica') {
    theme_minimal(base_size = base_size, base_family = base_family) %+replace%
        theme(
            strip.text = element_text(face = 'plain', size = 8),
            panel.grid = element_blank(),
            panel.border = element_rect(fill = NA, color = "gray50"),
            axis.ticks = element_line(color = "gray50"),
            axis.ticks.length = unit(2, 'points'),
            legend.position = 'none'
        )
}


read_bed <- function(filename, do_gather = TRUE){
    # filename <- paste0('/Volumes/MW_18TB/Lucas_Nell/lan/musGBS/dep/', name, '.bed.gz')
    sample_name <- (filename %>% 
                        gsub(pattern = '.gz|.bed', replacement = '', x = .) %>%
                        strsplit('/'))[[1]] %>% tail(1)
    df <- read_tsv(filename,
                   col_names = c('chrom', 'start', 'end', 'reads')) %>%
        mutate(start = start + 1, # Turning into 1-based indexing
               sample = sample_name)
    if (do_gather) {
        df <- df %>%
            gather('pos_type', 'pos', start, end) %>%
            arrange(chrom, pos)
    }
    return(df)
}





combine_beds <- function(dir, do_gather = TRUE){
    files <- list.files(dir, '.bed')
    one_bed <- function(x) {read_bed(paste(dir, x, sep = '/'), do_gather)}
    out_df <- lapply(files, one_bed) %>% 
        bind_rows %>%
        mutate(sample = factor(sample))
    return(out_df)
}




dir_1 <- '/Users/lucasnell/Desktop/dep/run_132571'
dir_2 <- '/Users/lucasnell/Desktop/dep/run_2'

run_1 <- combine_beds(dir_1)
run_2 <- combine_beds(dir_2)


run_1 %>% 
    filter(chrom == 'chrX') %>%
    ggplot(aes(x = pos/1e6, y = reads)) +
    plot_theme() +
    ggtitle('Chromosome X coverage for 10 random samples, run # 1') +
    geom_line(aes(color = sample)) +
    scale_x_continuous('Location on chromosome (Mb)', breaks = seq(0, 150, 25)) +
    scale_y_continuous('Reads', breaks = seq(0, 20, 10), limits = c(0, 30)) +
    facet_grid(sample ~ ., scales = 'free_y')



run_2 %>% 
    filter(chrom == 'chrX') %>%
    ggplot(aes(x = pos/1e6, y = reads)) +
    plot_theme() +
    ggtitle('Chromosome X coverage for 10 random samples, run # 2') +
    geom_line(aes(color = sample)) +
    scale_x_continuous('Location on chromosome (Mb)', breaks = seq(0, 150, 25)) +
    scale_y_continuous('Reads', breaks = seq(0, 20, 10), limits = c(0, 30)) +
    facet_grid(sample ~ .)


# rm(run_1, run_2)

# Here, 'wreads' --> not `gather`ed and only locations with > 0 reads
wreads_1 <- combine_beds(dir_1, do_gather = FALSE) %>% 
    filter(reads > 0)
wreads_2 <- combine_beds(dir_2, do_gather = FALSE) %>% 
    filter(reads > 0)

# Adding column of bp gap between locations
wreads_1 <- wreads_1 %>% 
    group_by(sample, chrom) %>% 
    mutate(gap = start - lag(end) - 1) %>% 
    ungroup
wreads_2 <- wreads_2 %>% 
    group_by(sample, chrom) %>% 
    mutate(gap = start - lag(end) - 1) %>% 
    ungroup


gap_plot_1 <- wreads_1 %>%
    filter(gap > 0) %>%
    ggplot(aes(x = gap, color = sample)) +
    plot_theme() +
    theme(axis.text.y = element_blank()) +
    ggtitle('Gaps between GBS peaks, run # 1') +
    geom_freqpoly(aes(y = ..density..), bins = 50) +
    scale_y_continuous("Density", breaks = seq(0, 0.5, 0.1)) +
    scale_x_log10('Gap size (bp)', breaks = 10^(seq(1, 7, 2)), 
                  labels = trans_format("log10", math_format(10^.x)))

gap_plot_2 <- wreads_2 %>%
    filter(gap > 0) %>%
    ggplot(aes(x = gap, color = sample)) +
    plot_theme() +
    theme(axis.text.y = element_blank()) +
    ggtitle('Gaps between GBS peaks, run # 2') +
    geom_freqpoly(aes(y = ..density..), bins = 50) +
    scale_y_continuous("Density", breaks = seq(0, 0.5, 0.1)) +
    scale_x_log10('Gap size (bp)', breaks = 10^(seq(1, 7, 2)), 
                  labels = trans_format("log10", math_format(10^.x)))

gap_plot <- rbind(ggplotGrob(gap_plot_1), 
                    ggplotGrob(gap_plot_2), 
                    size = "last")

grid.newpage()
grid.draw(gap_plot)



summ_1 <- wreads_1 %>%
    filter(gap > 0) %>%
    group_by(sample, chrom) %>%
    summarize(mean = mean(gap), 
              median = median(gap),
              sd = sd(gap)) %>%
    ungroup %>%
    arrange(chrom, sample) %>%
    select(chrom, sample, mean, median, sd)

summ_2 <- wreads_2 %>%
    filter(gap > 0) %>%
    group_by(sample, chrom) %>%
    summarize(mean = mean(gap), 
              median = median(gap),
              sd = sd(gap)) %>%
    ungroup %>%
    arrange(chrom, sample) %>%
    select(chrom, sample, mean, median, sd)


# write_tsv(summ_1, '~/Desktop/run1_gaps.txt')
# write_tsv(summ_2, '~/Desktop/run2_gaps.txt')







# , breaks = 10^(seq(1, 7, 2)), 
# labels = trans_format("log10", math_format(10^.x))
paste0('10^', seq(1, 7, 2))



# # Forward strand, MAPQ >= 20, "properly paired"
# r1q <- readBED('_1_qPP') %>% 
#     mutate(strand = 1)
# # Reverse strand, MAPQ >= 20, "properly paired"
# r2q <- readBED('_2_qPP') %>% 
#     mutate(strand = 2)
# 
# rqPP <- bind_rows(r1q, r2q) %>% 
#     mutate(strand = factor(strand, levels = c(1, 2), labels = c('+', '−')))
# 
# # Forward strand, MAPQ >= 20
# r1q <- readBED('_1_q') %>% 
#     mutate(strand = 1)
# # Reverse strand, MAPQ >= 20
# r2q <- readBED('_2_q') %>% 
#     mutate(strand = 2)
# 
# rq <- bind_rows(r1q, r2q) %>% 
#     mutate(strand = factor(strand, levels = c(1, 2), labels = c('+', '−')))
# 
# rm(r1q, r2q)
# 
# # # Both reads, MAPQ >= 20
# # rq <- readBED('_q')
# 
# # Filtering for reads < `threshold`
# threshold <- 96 * 1
# rqPP_filt <- rqPP %>%
#     mutate(reads = sapply(reads, function(x){ifelse(x > threshold, x, 0)}))
# rq_filt <- rq %>%
#     mutate(reads = sapply(reads, function(x){ifelse(x > threshold, x, 0)}))
# 
# 
# rqPP %>%
#     filter(chrom == 'chrX', pos >= (166410000 - 1e6)) %>%
#     ggplot(aes(x = pos/1e6, y = reads, color = strand)) + 
#     geom_vline(xintercept = 166410000/1e6, linetype = 3, color = 'gray60') +
#     geom_line() + 
#     scale_color_brewer(type = 'qual', palette = 'Set1') + 
#     ylab('Reads') +
#     xlab('Position on chromosome (Mb)') +
#     plotTheme() + 
#     theme(legend.direction = 'horizontal', legend.title = element_blank(),
#           legend.position = c(0.5, 0.90))
# 
# rq %>%
#     filter(chrom == 'chrX', pos >= (166410000 - 1e6)) %>%
#     ggplot(aes(x = pos/1e6, y = reads, color = strand)) + 
#     geom_vline(xintercept = 166410000/1e6, linetype = 3, color = 'gray60') +
#     geom_line() + 
#     scale_color_brewer(type = 'qual', palette = 'Set1') + 
#     ylab('Reads') +
#     xlab('Position on chromosome (Mb)') +
#     plotTheme() + 
#     theme(legend.direction = 'horizontal', legend.title = element_blank(),
#           legend.position = c(0.5, 0.90))
# 
# 
# 
# 
# rq %>%
#     filter(chrom == 'chrX', pos >= 166.434e6, pos <= 166.436e6) %>%
#     ggplot(aes(x = pos, y = reads, color = strand)) + 
#     geom_line() + 
#     scale_color_brewer(type = 'qual', palette = 'Set1') + 
#     facet_grid(strand ~ .) + 
#     theme(legend.position = 'none')
# 
# 
# 
# 
# 
# 
# rqPP %>%
#     filter(chrom == 'chrX', reads > 96 * 1) %>%
#     ggplot(aes(x = reads, y = ..count..)) + 
#     geom_density()
# 
# 



