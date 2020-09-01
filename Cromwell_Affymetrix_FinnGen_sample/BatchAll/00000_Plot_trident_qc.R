
# 
setwd("/Users/aoxliu/Documents/Project2_Finngen_mCA/Analysis_FinnGen_mCA/FinnGenBatchAll/Result")
df_calls <- read.table("FinnGenBatchAll.filtered.calls.tsv",header=T)

# df_calls = Var_toKeep
df_calls$SV <- factor(df_calls$chrom, levels = c('0-2 Mbp', '2-10 Mbp', '10-50 Mbp', '50-250 Mbp'))
df_calls$SV[df_calls$length <= 250e6] <- '50-250 Mbp'
df_calls$SV[df_calls$length <= 50e6] <- '10-50 Mbp'
df_calls$SV[df_calls$length <= 10e6] <- '2-10 Mbp'
df_calls$SV[df_calls$length <= 2e6] <- '0-2 Mbp'
df_calls$bdev[is.na(df_calls$bdev)] <- -.05
​
library(optparse)
library(ggplot2)
options(bitmapType = 'cairo')
​
pdf("/Users/aoxliu/Documents/Project2_Finngen_mCA/Analysis_FinnGen_mCA/FinnGenBatchAll/Result/init_plot.noCFfilter.pdf", width = 6, height = 6)
​
idx <- !( df_calls$chrom %in% c('chrX', 'chrY', 'MT') )
p <- ggplot(df_calls[idx,], aes(x=bdev, y=rel_cov, color=type, shape=type)) +
  geom_hline(yintercept = c(1.0, 2.0, 3.0), color = 'gray', size = .5, linetype = 'dashed') +
  geom_segment(aes(x = 0.0, y = 2.0, xend = 1.0/6.0, yend = 3.0), color = 'gray', size = .5, linetype = 'dashed') +
  geom_segment(aes(x = 0.0, y = 2.0, xend = 1.0/6.0, yend = 1.5), color = 'gray', size = .5, linetype = 'dashed') +
  geom_point(size = 1, alpha = 1/2) +
  scale_color_manual('', values = c('CN-LOH' = 'orange', 'CNP Deletion' = 'lightblue', 'CNP Duplication' = 'violetred', 'Loss' = 'blue', 'Gain' = 'red', 'Undetermined' = 'gray50')) +
  scale_shape_manual('', values = 0:5) +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(color=NA, fill=NA), legend.position = 'bottom', legend.box = 'horizontal') +
  facet_wrap(~SV)
print(p + scale_x_continuous('BAF deviation (Bdev)', breaks = c(-0.05, 0.0, 0.05, 0.1, 0.15, 0.2)) + scale_y_continuous(expression(paste('Relative coverage (rescaled ', 2^LRR, ')')), breaks = c(0, 1, 2, 3, 4)) + coord_cartesian(xlim = c(-0.05, 0.2), ylim = c(0, 4)))
print(p + scale_x_continuous('BAF deviation (Bdev)', breaks = c(0.0, 0.01, 0.02, 0.03, 0.04, 0.05)) + scale_y_continuous(expression(paste('Relative coverage (rescaled ', 2^LRR, ')')), breaks = c(1.8, 1.9, 2.0, 2.1, 2.2)) + coord_cartesian(xlim = c(0.0, 0.05), ylim = c(1.8, 2.2)))
dev.off()
  

  
## what's the summary(df$cf) for the variants with length2-50Mbp with relcov>2.75?
summary(df_calls[df_calls$length > 2e6 & df_calls$relcov >2.75,"cf"])




## can you figure out what is causing females to cluster as two different Y nonPAR LRR groups?




## there are many XXY individuals classified as unknown sex and there are individuals that are clearly mLOY that are classified as unknown sex ...
 this is how Affymetrix Axiom classifies them ... , the latter could be fixed if you want to include expanded mLOY in your analyses











