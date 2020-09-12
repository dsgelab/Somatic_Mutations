# 

library(data.table)
library(optparse)
library(grid)
library(reshape2)
library(scales)
library(ggrepel)
library(plotly)
library(meta)
library(tidyr)
library(ggsci)
library(binom)
library(ggplot2)
# options(bitmapType = 'cairo')

'%!in%' <- function(x,y)!('%in%'(x,y))





#------------------------------------------------------------------------------------------
#-------- Loading everything in -----------------------------------------------------------
#------------------------------------------------------------------------------------------
          
mCAs <- as.data.frame(fread("FinnGenBatchAll.calls.tsv"))
nrow(mCAs)        # 183,579  


mCAs_QC <- as.data.frame(fread("FinnGenBatchAll.filtered_default_GainCF.calls.tsv"))
nrow(mCAs_QC)     # 38,321   


mCAs_stats <- as.data.frame(fread("FinnGenBatchAll.stats.tsv"))
nrow(mCAs_stats)  # 201,458   


pheno <- read.table("finngen_R6_pheno_cov.txt", header=T, stringsAsFactors=F)
nrow(pheno)       # 260,405  
length(unique(pheno$FINNGENID))  # 260,405
hemeCancer <- read.table("finngen_R6_pheno_hemeCancer_id.lst", header=F)
nrow(hemeCancer)  # 3,526
pheno <- pheno[pheno[,1] %!in% hemeCancer[,1], ]
nrow(pheno)       # 257,011
pheno <- pheno[pheno[,"DEATH"]==0, ]
nrow(pheno)       # 235,640


# Only keep those with phenotypes
mCAs_phe <- merge(mCAs, pheno, by=c(1))
nrow(mCAs_phe)     # 153,310
length(unique(mCAs_phe[,1]))  # 97,291


mCAs_QC_phe <- merge(mCAs_QC, pheno, by=c(1))
nrow(mCAs_QC_phe)  # 29,991
length(unique(mCAs_QC_phe[,1]))  # 24,880


mCAs_stats_phe <- merge(mCAs_stats, pheno, by=c(1))
nrow(mCAs_stats_phe)   # 175,988
length(unique(mCAs_stats_phe[,1]))  # 175,988





#------------------------------------------------------------------------------------------
#-------Sample QC and plots ---------------------------------------------------------------
#------------------------------------------------------------------------------------------

# Sample QC
BAF_Auto_toRemove <- mCAs_stats_phe[which(mCAs_stats_phe$baf_auto>0.05),1] 
length(BAF_Auto_toRemove)   # 222 



# Plot of Baf_plot without any QC (Supplementary Figure 1-A)
# pdf("FinnGenBatchAll.baf_plot.pdf", width=4, height=4)
# tiff("FinnGenBatchAll.baf_plot.tiff", width = 6, height = 6, unit="in", res=300)
tiff("FinnGenBatchAll_nodeath_nohemeCancer.baf_plot.tiff", width = 6, height = 6, unit="in", res=300)
if ('lrr_auto' %in% colnames(mCAs_stats_phe)) { col_x <- 'lrr_auto'; lbl_x <- 'GC-adjusted LRR auto-correlation'} else if ('COV_AUTO' %in% colnames(mCAs_stats_phe)) { col_x <- 'cov_auto'; lbl_x <- 'GC-adjusted coverage auto-correlation' }
p <- ggplot(mCAs_stats_phe, aes_string(x=col_x, y='baf_auto', color='computed_gender')) +
  geom_point(size=.5, alpha=1/2) +
  scale_x_continuous(lbl_x) +
  scale_y_continuous('BAF auto-correlation') +
  scale_color_manual('', values=c('M'='blue', 'F'='orchid', 'U'='gray'), labels = c('M'='Male', 'F'='Female', 'U'='Undetermined')) +
  theme_bw(base_size = 14) +
  theme(legend.position= 'bottom', legend.box='horizontal')
print(p)
dev.off()



# Plot of genderMismatches without any QC (Supplementary Figure 1-B)
mergedStatsPhenos <- mCAs_stats_phe[mCAs_stats_phe$computed_gender!="U", ]
dim(mergedStatsPhenos)   # 175,143     43

mergedStatsPhenos[,"MochaGender_PhenoGender_Mismatch"] <- ifelse((mergedStatsPhenos$computed_gender=="F" & mergedStatsPhenos$SEX=="male") | (mergedStatsPhenos$computed_gender=="M" & mergedStatsPhenos$SEX=="female"), "Sex mismatch", "Sex match")
mismatch_Sex_ids <- mergedStatsPhenos[which(mergedStatsPhenos$MochaGender_PhenoGender_Mismatch=="Sex mismatch"), 1]   # 76 samples

# pdf("FinnGenBatchAll.genderMismatches.pdf", width=9, height=4)
# tiff("FinnGenBatchAll.genderMismatches.tiff", width=9, height=4, unit="in", res=300)
tiff("FinnGenBatchAll_nodeath_nohemeCancer.genderMismatches.tiff", width=9, height=4, unit="in", res=300)
p <- ggplot(mergedStatsPhenos[mergedStatsPhenos$x_nonpar_lrr_median> -1.5 & !is.na(mergedStatsPhenos$MochaGender_PhenoGender_Mismatch),], aes(x=x_nonpar_lrr_median, y=y_nonpar_lrr_median, color=`computed_gender`)) + geom_point() + theme_bw(base_size = 14)+xlab("x_nonpar_lrr_median")+ylab("y_nonpar_lrr_median")  
pp <- p + facet_wrap( ~ MochaGender_PhenoGender_Mismatch, nrow = 1) + theme(legend.position = "right")
print(pp)
dev.off()




  
#------------------------------------------------------------------------------------------
#-------Variant QC and plots---------------------------------------------------------------
#------------------------------------------------------------------------------------------

# Basic QC for variants, including removal of likely germline variants (lod_baf_phase <20 for autosomal variants,lod_baf_phase<5 for ChrX, or annotated as germline CNP)
Var_toKeep <- mCAs_phe[which(!(mCAs_phe$sample_id %in% BAF_Auto_toRemove) & (mCAs_phe$lod_baf_phase > 20 & mCAs_phe$chrom!="chrX")|(mCAs_phe$lod_baf_phase>5 & mCAs_phe$chrom=="chrX") & mCAs_phe$type!="CNP_Gain" & mCAs_phe$type!="CNP_Loss"),]
# Var_toKeep <- Var_toKeep[-which(Var_toKeep$CHROM=="chr6" & Var_toKeep$BEG_GRCh38>27518932 & Var_toKeep$END_GRCh38<33480487),]  # 0 
nrow(Var_toKeep)   # 48,572



df_calls <- Var_toKeep
df_calls$SV <- factor(df_calls$chrom, levels = c('0-2 Mbp', '2-10 Mbp', '10-50 Mbp', '50-250 Mbp'))
df_calls$SV[df_calls$length <= 250e6] <- '50-250 Mbp'
df_calls$SV[df_calls$length <= 50e6] <- '10-50 Mbp'
df_calls$SV[df_calls$length <= 10e6] <- '2-10 Mbp'
df_calls$SV[df_calls$length <= 2e6] <- '0-2 Mbp'
df_calls$bdev[is.na(df_calls$bdev)] <- -.05



# Trident plot with basic QC (Supplementary Figure 2-A)
# pdf("FinnGenBatchAll.init_plot.trident.noCFfilter.pdf", width = 6, height = 6)
# tiff("FinnGenBatchAll.init_plot.trident.noCFfilter.tiff", width = 6, height = 6, unit="in", res=300)
tiff("FinnGenBatchAll_nodeath_nohemeCancer.init_plot.trident.noCFfilter.tiff", width = 6, height = 6, unit="in", res=300)
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
# print(p + scale_x_continuous('BAF deviation (Bdev)', breaks = c(0.0, 0.01, 0.02, 0.03, 0.04, 0.05)) + scale_y_continuous(expression(paste('Relative coverage (rescaled ', 2^LRR, ')')), breaks = c(1.8, 1.9, 2.0, 2.1, 2.2)) + coord_cartesian(xlim = c(0.0, 0.05), ylim = c(1.8, 2.2)))
dev.off()



# Trident plot of with basic QC (Supplementary Figure 2-B)
df_calls$chrom <- as.factor(gsub('^chr', '', gsub('^chrM', 'MT', df_calls$chrom)))
ord <- order(as.numeric(gsub('MT', '26', gsub('Y', '24', gsub('X', '23', levels(df_calls$chrom))))))
df_calls$chrom <- factor(df_calls$chrom, levels(df_calls$chrom)[ord])

# pdf("FinnGenBatchAll.init_plot.CountsbyChrom.noCFfilter.pdf", width = 6, height = 6)
# tiff("FinnGenBatchAll.init_plot.CountsbyChrom.noCFfilter.tiff", width = 6, height = 6, unit="in", res=300)
tiff("FinnGenBatchAll_nodeath_nohemeCancer.init_plot.CountsbyChrom.noCFfilter.tiff", width = 6, height = 6, unit="in", res=300)
idx <- !( df_calls$chrom %in% c('X') ) 
p <- ggplot(df_calls[idx,], aes(x=chrom, fill=type)) +
  geom_bar(stat = 'count', color = 'black') +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual('', values = c('CN-LOH' = 'orange', 'Loss' = 'blue', 'Gain' = 'red', 'Undetermined' = 'gray50')) +
  theme_bw(base_size = 14) +
  theme(legend.position = 'bottom', legend.box = 'horizontal')
print(p)
dev.off()




#------------------------------------------------------------------------------------------
#-------More Variant QC and plots----------------------------------------------------------
#------------------------------------------------------------------------------------------

# More QC for variants, including mocha default QC + length 2-50 Mbp with relcov>2.75 + basic QC
Var_toKeep_Filtered <- mCAs_QC_phe[which(!(mCAs_QC_phe$sample_id %in% BAF_Auto_toRemove) & (mCAs_QC_phe$lod_baf_phase > 20 & mCAs_QC_phe$chrom!="chrX")|(mCAs_QC_phe$lod_baf_phase>5 & mCAs_QC_phe$chrom=="chrX")),]
nrow(Var_toKeep_Filtered)   # 26,229
write.table(Var_toKeep_Filtered, "FinnGenBatchAll_nodeath_nohemeCancer.MoreQC.withoutCFfilter.txt", col.names = T, row.names = F, quote = F, sep = "\t")


df_calls = Var_toKeep_Filtered
df_calls$SV <- factor(df_calls$chrom, levels = c('0-2 Mbp', '2-10 Mbp', '10-50 Mbp', '50-250 Mbp'))
df_calls$SV[df_calls$length <= 250e6] <- '50-250 Mbp'
df_calls$SV[df_calls$length <= 50e6] <- '10-50 Mbp'
df_calls$SV[df_calls$length <= 10e6] <- '2-10 Mbp'
df_calls$SV[df_calls$length <= 2e6] <- '0-2 Mbp'
df_calls$bdev[is.na(df_calls$bdev)] <- -.05



# Trident plot with more QC (Supplementary Figure 2-C)
# pdf("FinnGenBatchAll.init_plot.trident.filter.pdf", width = 6, height = 6)
# tiff("FinnGenBatchAll.init_plot.trident.filter.tiff", width = 6, height = 6, unit="in", res=300)
tiff("FinnGenBatchAll_nodeath_nohemeCancer.init_plot.trident.filter.tiff", width = 6, height = 6, unit="in", res=300)
idx <- !( df_calls$chrom %in% c('chrX', 'Y', 'MT') )
p <- ggplot(df_calls[idx,], aes(x=bdev, y=rel_cov, color=type, shape=type)) +
  geom_hline(yintercept = c(1.0, 2.0, 3.0), color = 'gray', size = .5, linetype = 'dashed') +
  geom_segment(aes(x = 0.0, y = 2.0, xend = 1.0/6.0, yend = 3.0), color = 'gray', size = .5, linetype = 'dashed') +
  geom_segment(aes(x = 0.0, y = 2.0, xend = 1.0/6.0, yend = 1.5), color = 'gray', size = .5, linetype = 'dashed') +
  geom_point(size = 1, alpha = 1/2) +
  scale_color_manual('', values = c('CN-LOH' = 'orange', 'Loss' = 'blue', 'Gain' = 'red', 'Undetermined' = 'gray50')) +
  scale_shape_manual('', values = 0:5) +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(color=NA, fill=NA), legend.position = 'bottom', legend.box = 'horizontal') +
  facet_wrap(~SV)
print(p + scale_x_continuous('BAF deviation (Bdev)', breaks = c(-0.05, 0.0, 0.05, 0.1, 0.15, 0.2)) + scale_y_continuous(expression(paste('Relative coverage (rescaled ', 2^LRR, ')')), breaks = c(0, 1, 2, 3, 4)) + coord_cartesian(xlim = c(-0.05, 0.2), ylim = c(0, 4)))
# print(p + scale_x_continuous('BAF deviation (Bdev)', breaks = c(0.0, 0.01, 0.02, 0.03, 0.04, 0.05)) + scale_y_continuous(expression(paste('Relative coverage (rescaled ', 2^LRR, ')')), breaks = c(1.8, 1.9, 2.0, 2.1, 2.2)) + coord_cartesian(xlim = c(0.0, 0.05), ylim = c(1.8, 2.2)))
dev.off()



# Trident plot with more QC (Supplementary Figure 2-D)
df_calls$chrom <- as.factor(gsub('^chr', '', gsub('^chrM', 'MT', df_calls$chrom)))
ord <- order(as.numeric(gsub('MT', '26', gsub('Y', '24', gsub('X', '23', levels(df_calls$chrom))))))
df_calls$chrom <- factor(df_calls$chrom, levels(df_calls$chrom)[ord])

# pdf("FinnGenBatchAll.init_plot.CountsbyChrom.filter.pdf", width = 6, height = 6)
# tiff("FinnGenBatchAll.init_plot.CountsbyChrom.filter.tiff", width = 6, height = 6, unit="in", res=300)
tiff("FinnGenBatchAll_nodeath_nohemeCancer.init_plot.CountsbyChrom.filter.tiff", width = 6, height = 6, unit="in", res=300)
idx <- !( df_calls$chrom %in% c('X') ) 
p <- ggplot(df_calls[idx,], aes(x=chrom, fill=type)) +
  geom_bar(stat = 'count', color = 'black') +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual('', values = c('CN-LOH' = 'orange', 'Loss' = 'blue', 'Gain' = 'red', 'Undetermined' = 'gray50')) +
  theme_bw(base_size = 14) +
  theme(legend.position = 'bottom', legend.box = 'horizontal')
print(p)
dev.off()





#------------------------------------------------------------------------------------------
#-------Now turning QCed table into 0/1 counts for analysis--------------------------------
#------------------------------------------------------------------------------------------

print(as.data.frame(table(mCAs_stats_phe[,c("computed_gender","SEX")])))
#  computed_gender    SEX   Freq
#1               F female 108752
#2               M female     34  # remove
#3               U female      0
#4               F   male     49  # remove
#5               M   male  80026
#6               U   male   1029


allSamples <- mCAs_stats_phe[mCAs_stats_phe$sample_id %!in% mismatch_Sex_ids, c("sample_id","computed_gender","SEX")]  # remove those with mismatched sex 
nrow(allSamples)         # 189,890
print(as.data.frame(table(allSamples[,c("computed_gender","SEX")])))
#  computed_gender    SEX   Freq
#1               F female 108752
#2               M female      0
#3               U female      0
#4               F   male      0
#5               M   male  80026
#6               U   male   1029


CNVs <- merge(allSamples, Var_toKeep_Filtered[,c("sample_id","chrom","cf")], by=c(1), all.x=TRUE)
nrow(CNVs)               # 180,071
length(unique(CNVs$sample_id))   # 175,912



# Define 4 expousures without considering CF
CNVs$has_MosaicCNV <- ifelse(!is.na(CNVs$chrom), 1, 0)

CNVs$ChrX <- ifelse(CNVs$has_MosaicCNV==0 & CNVs$SEX=="female", 0,
             ifelse(CNVs$chrom=="chrX" & CNVs$SEX=="female", 1, NA))

CNVs$ChrY <- ifelse(CNVs$has_MosaicCNV == 0 & CNVs$SEX=="male", 0,
             ifelse(CNVs$chrom == "chrX" & CNVs$SEX=="male", 1, NA))

CNVs$Autosomes <- ifelse(CNVs$has_MosaicCNV == 0, 0,
                  ifelse(CNVs$chrom %in% paste("chr", 1:22, sep=""), 1, NA))

sum(CNVs$has_MosaicCNV==1,na.rm=T)==sum(sum(CNVs$ChrX==1,na.rm=T), sum(CNVs$ChrY==1,na.rm=T), sum(CNVs$Autosomes==1,na.rm=T))   # 31,264



# Define 4 expousures with CF>0.1
CNVs$cf <- as.numeric(as.character(CNVs$cf))

CNVs$Large_CNV_Clonev2 <- ifelse(CNVs$has_MosaicCNV==0, 0,
                          ifelse(CNVs$cf>=0.1, 1, 
                          ifelse(CNVs$has_MosaicCNV==1 & CNVs$cf<0.1 | is.na(CNVs$CF), NA, 0)))

CNVs$Large_ChrX <- ifelse(CNVs$has_MosaicCNV == 0 & CNVs$SEX=="female", 0,
                   ifelse(CNVs$Large_CNV_Clone == 1 & CNVs$ChrX == 1, 1, 
                   ifelse(CNVs$has_MosaicCNV == 1 & CNVs$Large_CNV_Clone != 1, NA, 0)))
CNVs$Large_ChrX <- ifelse(is.na(CNVs$ChrX), NA, CNVs$Large_ChrX)


CNVs$Large_ChrY <- ifelse(CNVs$has_MosaicCNV == 0 & CNVs$SEX=="male", 0,
                   ifelse(CNVs$Large_CNV_Clone == 1 & CNVs$ChrY == 1, 1, 
                   ifelse(CNVs$has_MosaicCNV == 1 & CNVs$Large_CNV_Clone != 1, NA, 0)))
CNVs$Large_ChrY <- ifelse(is.na(CNVs$ChrY), NA, CNVs$Large_ChrY)


CNVs$Large_Autosomes <- ifelse(CNVs$has_MosaicCNV == 0, 0,
                        ifelse(CNVs$Large_CNV_Clone == 1 & CNVs$Autosomes == 1, 1, 
                        ifelse(CNVs$has_MosaicCNV == 1 & CNVs$Large_CNV_Clone != 1, NA, 0)))

sum(CNVs$Large_CNV_Clonev2==1,na.rm=T)==sum(sum(CNVs$Large_ChrX==1,na.rm=T), sum(CNVs$Large_ChrY==1,na.rm=T), sum(CNVs$Large_Autosomes==1,na.rm=T))   # 15,960



# QC and format exposure
variation <- c("has_MosaicCNV", "Large_CNV_Clonev2", "ChrX", "ChrY", "Autosomes", "Large_ChrX", "Large_ChrY", "Large_Autosomes")
uniqueSamples <- unique(CNVs[-which(CNVs[,1] %in% c(BAF_Auto_toRemove)),1])
length(uniqueSamples)  # 175,690

max_DF <- data.frame(ourSid = uniqueSamples)
for (i in 1:length(variation)){
	tmp <- aggregate(CNVs[,variation[i]]~sample_id, data=CNVs, max)
	colnames(tmp)[2] <- variation[i]
	max_DF <- merge(max_DF, tmp, by=c(1), all.x=TRUE)
}

max_DF <- merge(max_DF, mCAs_stats_phe[, c("sample_id","SEX","BL_AGE")], by=c(1))
dim(max_DF)    # 175,690     11

# write.table(max_DF, "mCA_01_Counts.afterQC.withoutCF.txt", col.names=T, row.names=F, quote=F, sep="\t")
write.table(max_DF, "FinnGenBatchAll_nodeath_nohemeCancer.mCA_01_Counts.afterQC.withoutCF.txt", col.names=T, row.names=F, quote=F, sep="\t")






#------------------------------------------------------------------------------------------
#--------- Making plots with age-----------------------------------------------------------
#------------------------------------------------------------------------------------------

max_DF$Age <- max_DF$BL_AGE
max_DF$Sex <- ifelse(max_DF$SEX=="female", "Female", ifelse(max_DF$SEX=="male", "Male", NA))
max_DF$age_bins <- ifelse(max_DF$Age<40, "<40",
                   ifelse(max_DF$Age>=40 & max_DF$Age<45, "[40,45)",
                   ifelse(max_DF$Age>=45 & max_DF$Age<50, "[45,50)",
                   ifelse(max_DF$Age>=50 & max_DF$Age<55, "[50,55)",
                   ifelse(max_DF$Age>=55 & max_DF$Age<60, "[55,60)",
                   ifelse(max_DF$Age>=60 & max_DF$Age<65, "[60,65)",
                   ifelse(max_DF$Age>=65 & max_DF$Age<70, "[65,70)",
                   ifelse(max_DF$Age>=70 & max_DF$Age<75, "[70,75)",
                   ifelse(max_DF$Age>=75 & max_DF$Age<80, "[75,80)",
                   ifelse(max_DF$Age>=80 & max_DF$Age<85, "[80,85)",
                   ifelse(max_DF$Age>=85 , ">85", NA)))))))))))
max_DF$age_bins <- ordered(max_DF$age_bins , levels=c("<40", "[40,45)", "[45,50)", "[50,55)", "[55,60)", "[60,65)", "[65,70)", "[70,75)", "[75,80)", "[80,85)", ">85"))
max_DF$Sex <- ordered(max_DF$Sex, levels=c("Male", "Female"))



# Age  
Mosaic_Names <- c('has_MosaicCNV', 'Large_CNV_Clonev2', "Autosomes", "Large_Autosomes") 
FemaleMaleMosaics <- c('ChrX','Large_ChrX','ChrY','Large_ChrY')
max_DF <- max_DF[which(!is.na(max_DF$Sex) & !is.na(max_DF$age_bins)),]
nrow(max_DF)   # 175,690
 
# mean_cl_binom_tol <- function(x) {out <- binom.bayes(sum(x>0, na.rm=TRUE),sum(!is.na(x)),tol=0.65/10)[c("mean","lower","upper")]; names(out) <- c("y","ymin","ymax"); return(out)}
mean_cl_binom_tol <- function(x) {out <- binom.bayes(sum(x>0, na.rm=TRUE),sum(!is.na(x)))[c("mean","lower","upper")]; names(out) <- c("y","ymin","ymax"); return(out)}

for (i in Mosaic_Names){
	# pdf(paste0("FinnGenBatchAll.mCAAgeassoc.", i, ".pdf"), width = 5, height= 4)
	# tiff(paste0("FinnGenBatchAll.mCAAgeassoc", i, ".tiff"), width = 5, height= 4, units="in", res=300)
	tiff(paste0("FinnGenBatchAll_nodeath_nohemeCancer.mCAAgeassoc", i, ".tiff"), width = 5, height= 4, units="in", res=300)
	p4 <- ggplot(max_DF, aes(x=as.numeric(age_bins), y=max_DF[,i], linetype=Sex)) +
	      stat_summary(fun=mean, geom="line", size=1.7) + 
	      stat_summary(fun.data=mean_cl_binom_tol, geom="ribbon", alpha=0.25, color="transparent", fill="gray") +
	      xlab("Age (years)") + 
	      ylab("Proportion") + theme_minimal() + guides(colour = FALSE)+scale_x_continuous(breaks=1:length(levels(max_DF$age_bins)),labels=levels(max_DF$age_bins)) +
	      theme(axis.text.x=element_text(size=14, angle=45, hjust=1), strip.text.y=element_text(size=14), axis.text.y=element_text(size=14, hjust=1, color="black"), axis.ticks=element_line(colour="black"), 
	            axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),
	            legend.text=element_text(size=14),legend.position="right",legend.key = element_blank(),legend.title = element_blank(),legend.background = element_blank(), strip.text = element_text(size=14)) 
	print(p4)
	dev.off()
}


for (i in FemaleMaleMosaics){
	# pdf(paste0("FinnGenBatchAll.mCAAgeassoc", i, ".pdf"), width = 4, height= 4)
	# tiff(paste0("FinnGenBatchAll.mCAAgeassoc", i, ".tiff"), width = 4, height= 4, units="in", res=300)
	tiff(paste0("FinnGenBatchAll_nodeath_nohemeCancer.mCAAgeassoc", i, ".tiff"), width = 4, height= 4, units="in", res=300)
	p4 <- ggplot(max_DF, aes(x=as.numeric(age_bins), y=max_DF[,i])) +
	      stat_summary(fun=mean,geom="line", size=1.7) +
          stat_summary(fun.data=mean_cl_binom_tol, geom="ribbon",alpha=0.25,color="transparent",fill="gray") +
          xlab("Age (years)") + 
          ylab("Proportion") + theme_minimal() + guides(colour = FALSE)+scale_x_continuous(breaks=1:length(levels(max_DF$age_bins)),labels=levels(max_DF$age_bins)) +
          theme(axis.text.x = element_text(size=14, angle = 45, hjust = 1),
          strip.text.y = element_text(size = 14),
          axis.text.y = element_text(size=14,hjust=1, color = "black"),
          axis.ticks =  element_line(colour = "black"), 
          axis.title.y= element_text(size=14),  axis.title.x= element_text(size=14),
          legend.text=element_text(size=14),legend.position="right",legend.key = element_blank(),legend.title = element_blank(),legend.background = element_blank(), strip.text = element_text(size=14)) 
	print(p4)
	dev.off()
}


## Autosomes   Large_ChrX  Large_Autosomes
# In binom.bayes(sum(x > 0, na.rm = TRUE), sum(!is.na(x))) :
#  1 confidence interval failed to converge (marked by '*').
#  Try changing 'tol' to a different value.






#------------------------------------------------------------------------------------------
#---------Summary/Count of mCA by age_bin-----------------------------------------------------------
#------------------------------------------------------------------------------------------

mca_age <- matrix(NA, ncol=15, nrow=11)
colnames(mca_age) <- c("Age", "male_N",   "male_mCA_AnyChr",   "male_mCA_LargeAnyChr",   "male_mCA_Autosomes",   "male_mCA_LargeAutosomes", "male_mCA_ChrY",   "male_mCA_LargeChrY",
                              "female_N", "female_mCA_AnyChr", "female_mCA_LargeAnyChr", "female_mCA_Autosomes", "female_mCA_LargeAutosomes", "female_mCA_ChrX", "female_mCA_LargeChrX")
bins <- c("<40", "[40,45)", "[45,50)", "[50,55)", "[55,60)", "[60,65)", "[65,70)", "[70,75)", "[75,80)", "[80,85)", ">85")

	
for (i in 1:length(bins)){
	mca_age[i,"Age"] <- bins[i]
	d <- max_DF[max_DF$age_bins==bins[i],]
	
	mca_age[i,"male_N"] <- sum(d$SEX=="male", na.rm=T)
	mca_age[i,"female_N"] <- sum(d$SEX=="female", na.rm=T)
	
	mca_age[i,"male_mCA_AnyChr"] <- sum(d$SEX=="male" & d$has_MosaicCNV=="1", na.rm=T)		
	mca_age[i,"female_mCA_AnyChr"] <- sum(d$SEX=="female" & d$has_MosaicCNV=="1", na.rm=T)			
	
	mca_age[i,"male_mCA_LargeAnyChr"] <- sum(d$SEX=="male" & d$Large_CNV_Clonev2=="1", na.rm=T)						
	mca_age[i,"female_mCA_LargeAnyChr"] <- sum(d$SEX=="female" & d$Large_CNV_Clonev2=="1", na.rm=T)		
			
	mca_age[i,"male_mCA_Autosomes"] <- sum(d$SEX=="male" & d$Autosomes=="1", na.rm=T)				
	mca_age[i,"female_mCA_Autosomes"] <- sum(d$SEX=="female" & d$Autosomes=="1", na.rm=T)			
	
	mca_age[i,"male_mCA_LargeAutosomes"] <- sum(d$SEX=="male" & d$Large_Autosomes=="1", na.rm=T)	
	mca_age[i,"female_mCA_LargeAutosomes"] <- sum(d$SEX=="female" & d$Large_Autosomes=="1", na.rm=T)	

	mca_age[i,"male_mCA_ChrY"] <- sum(d$SEX=="male" & d$ChrY=="1", na.rm=T)				
	mca_age[i,"female_mCA_ChrX"] <- sum(d$SEX=="female" & d$ChrX=="1", na.rm=T)	

	mca_age[i,"male_mCA_LargeChrY"] <- sum(d$SEX=="male" & d$Large_ChrY=="1", na.rm=T)					
	mca_age[i,"female_mCA_LargeChrX"] <- sum(d$SEX=="female" & d$Large_ChrX=="1", na.rm=T)		
}

write.table(mca_age, "FinnGenBatchAll_nodeath_nohemeCancer.Age_Summary_Count.txt", col.names=T, row.names=F, quote=F, sep="\t")



