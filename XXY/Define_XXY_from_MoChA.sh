

#######################################################
#          Define XXY from MoChA output               #
#######################################################

# gcloud compute ssh "cromwell-dsge" --project "finngen-refinery-dsgelab" --zone europe-west1-b 


# You can identify XXY simply by thresholding on the mLRR-X and mLRR-Y. 
# As you can see below there is clear separation of the XXY (orange) from the XY and XX clusters.


cd /home/aoxliu/LOX
gcsfuse --implicit-dirs  dsge-cromwell  dsge-cromwell


setwd("/home/aoxliu/LOX/")

# install.packages("data.table")
library(data.table)
library(dplyr)
library(ggplot2)
'%!in%' <- function(x,y)!('%in%'(x,y))


# read in data -----------------
df_stats <- data.frame(read.table("dsge-cromwell/mocha/43ba66a5-1bed-4179-94ad-f9fe9379c889/call-mocha_stats_tsv/FinnGenR9Affymetrix.stats.tsv", sep="\t", header=T)) 
dim(df_stats)  # 316,920     20
df_stats <- df_stats %>% mutate(sample_id=as.character(sample_id), computed_gender=as.character(computed_gender))
df_stats %>% group_by(computed_gender) %>% count()
#  computed_gender      n
#1 F               178537
#2 M               136221
#3 U                 2162



# adjust by autosomes, based on https://github.com/freeseek/mocha/blob/master/summary_plot.R  -----------------
df_stats$x_nonpar_adj_lrr_median <- df_stats$x_nonpar_lrr_median - df_stats$lrr_median
df_stats$y_nonpar_adj_lrr_median <- df_stats$y_nonpar_lrr_median - df_stats$lrr_median
df_stats$mt_adj_lrr_median <- df_stats$mt_lrr_median - df_stats$lrr_median


# XXY with and without adjust -----------------
df_stats %>% filter(y_nonpar_lrr_median>-0.5 & x_nonpar_lrr_median>-0.125) %>% nrow()   # 211
df_stats %>% filter(y_nonpar_adj_lrr_median>-0.5 & x_nonpar_adj_lrr_median>-0.125) %>% nrow()   # 211, x_nonpar_lrr_median are small
df_stats %>% filter(y_nonpar_adj_lrr_median>-0.5 & y_nonpar_adj_lrr_median<0.5 & x_nonpar_adj_lrr_median>-0.125 & x_nonpar_adj_lrr_median<0.125) %>% nrow()   # 209, x_nonpar_lrr_median are small


# defined XXY -----------------
df_stats <- df_stats %>% mutate(Type=ifelse(y_nonpar_adj_lrr_median>-0.5 & y_nonpar_adj_lrr_median<0.5 & x_nonpar_adj_lrr_median>-0.125 & x_nonpar_adj_lrr_median<0.125, "K", computed_gender))
df_stats %>% group_by(Type) %>% count()
#  Type       n
#1 F     178,537
#2 K        209
#3 M     136,217
#4 U       1,957


# write out a list for individuals with XXY
XXY_lst <- df_stats %>% filter(Type=="K") %>% select(sample_id)
write.table(XXY_lst, "FinnGenR9.klinefelter.lines", append=F, quote=F, sep="\t", row.names=F, col.names=F)
system("gsutil cp  FinnGenR9.klinefelter.lines  gs://dsge-aoxing/mocha/FinnGenXXY/")


keep_lst <- df_stats %>% filter(y_nonpar_adj_lrr_median>-0.5 & x_nonpar_adj_lrr_median>-0.125) %>% select(sample_id)
write.table(keep_lst, "FinnGenR9.klinefelter.keep.lines", append=F, quote=F, sep="\t", row.names=F, col.names=F)
system("gsutil cp  FinnGenR9.klinefelter.keep.lines  gs://dsge-aoxing/mocha/FinnGenXXY/")

remove_lst <- df_stats %>% filter(y_nonpar_adj_lrr_median>-0.5 & x_nonpar_adj_lrr_median>-0.125) %>% filter(Type!="K") %>% select(sample_id)
write.table(remove_lst, "FinnGenR9.klinefelter.remove.lines", append=F, quote=F, sep="\t", row.names=F, col.names=F)
system("gsutil cp  FinnGenR9.klinefelter.remove.lines  gs://dsge-aoxing/mocha/FinnGenXXY/")



# plot -----------------
p <- ggplot(df_stats, aes(x=x_nonpar_adj_lrr_median, y=y_nonpar_adj_lrr_median, color=Type)) +
  geom_point(data=df_stats[df_stats$call_rate < 0.97 | df_stats$baf_auto > 0.03,], color='black', shape=1, size=0.5, alpha=1/2) +
  geom_point(shape=20, size=0.5, alpha=1/2) +
  scale_x_continuous('X nonPAR median LRR (autosome corrected)') +
  scale_y_continuous('Y nonPAR median LRR (autosome corrected)') +
  scale_color_manual('', values=c('M'='blue', 'F'='orchid', 'K'='orange', 'U'='gray'), labels = c('M'='Male', 'F'='Female', 'K'='Klinefelter', 'U'='Undetermined')) +
  theme_bw(base_size = 16) +
  theme(legend.position = 'bottom', legend.box = 'horizontal')
print(p)
ggsave("DefineXXY_by_LRR.png", p, width=8, height=8)
system("gsutil cp  DefineXXY_by_LRR.png   gs://dsge-aoxing/mocha/FinnGenXXY")


# check prevalence ------------
round(100*209/136221,2)    # 0.15
round(100*2208/5000000,2)  # 0.04



###########################################################################################
#    histogram for 209 Klinefelter cases with respect to the x_nonpar_n_hets variable     # 
###########################################################################################

# what you can do quickly is to make a histogram for those 211 with respect to the x_nonpar_n_hets variable

df_stats_K_cut <- df_stats %>% filter(Type=="K") %>% mutate(grp=ifelse(x_nonpar_n_hets<100,"<100", 
                                        ifelse(x_nonpar_n_hets>=100 & x_nonpar_n_hets<500, "100-500", 
                                        ifelse(x_nonpar_n_hets>=500 & x_nonpar_n_hets<1000, "500-1000", 
                                        ifelse(x_nonpar_n_hets>=1000 & x_nonpar_n_hets<1500, "1000-1500",
                                        ifelse(x_nonpar_n_hets>=1500 & x_nonpar_n_hets<2000, "1500-2000",
                                        ifelse(x_nonpar_n_hets>=2000 & x_nonpar_n_hets<2500, "2000-2500",
                                        ifelse(x_nonpar_n_hets>=2500 & x_nonpar_n_hets<3000, "2500-3000",
                                        ifelse(x_nonpar_n_hets>=3000 & x_nonpar_n_hets<3500, "3000-3500", ">3500")))))))))

df_stats_K_cut_n <- df_stats_K_cut %>% group_by(grp) %>% count()
df_stats_K_cut_n <- df_stats_K_cut_n %>% mutate(grp=factor(grp,levels=c("<100","100-500", "500-1000", "1000-1500", "1500-2000", "2000-2500", "2500-3000", "3000-3500",">3500")))

p_Color <- ggplot(data=df_stats_K_cut_n, aes(y=n, x=grp)) +
  # geom_bar(aes(color=grp_text, fill=grp_text), stat="identity", width=0.5) + 
  geom_bar(stat="identity", width=0.5) + 
  labs(y="N of individuals", x="N of heterozygous sites in the X nonPAR region", title="FinnGen Klinefelter cases (N=209)") + 
  # scale_color_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
  # scale_fill_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
  theme(axis.text.x=element_text(hjust=0, vjust=0, size=6, angle=45, face="bold", color="black"), axis.text.y=element_text(hjust=0, vjust=0, size=6, angle=0, face="bold", color="black"), axis.title=element_text(size=6, face="bold"), plot.title=element_text(size=10, face="bold")) +
  theme_classic() + 
  theme(legend.position="none")
ggsave(paste0("FinnGen209KlinefelterCases_NumberHeterozygousSitesInXnonPAR.png"), p_Color, width=7, height=5)
system("gsutil cp  FinnGen209KlinefelterCases_NumberHeterozygousSitesInXnonPAR.png   gs://dsge-aoxing/mocha/FinnGenXXY")



