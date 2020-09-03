
# Copy files to working directory
# cd /Users/aoxliu/Documents/Project2_Finngen_mCA/Analysis_FinnGen_mCA/FinnGenBatchAll/Result
# wf_id="3bb7e444-d05c-49d7-9e48-d5e012240aa9"
# gsutil cp  gs://dsge-cromwell/mocha/${wf_id}/call-mocha_calls_tsv/*tsv .
# gsutil cp  gs://dsge-cromwell/mocha/${wf_id}/call-mocha_stats_tsv/FinnGenBatchAll.stats.tsv  .
# gsutil cat gs://finngen-production-library-red/finngen_R6/phenotype_2.0/data/finngen_R6_v2_minimum.txt|awk '{print $1, $2,$3,$4,$5,$6,$7,$8,$9,$12,$13,$15}' > finngen_R6_v2_min.txt


setwd("/Users/aoxliu/Documents/Project2_Finngen_mCA/Analysis_FinnGen_mCA/FinnGenBatchAll/Result")

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

pheno <- read.table("finngen_R6_v2_min.txt", header=T, stringsAsFactors=F)
dim(pheno)   # 271343     12              
                  

mCAs <- as.data.frame(fread("FinnGenBatchAll.calls.tsv"))
dim(mCAs)     # 183579    22

# mCAs <- as.data.frame(fread("FinnGenBatchAll.filtered_default.calls.tsv"))
mCAs_QC <- as.data.frame(fread("FinnGenBatchAll.filtered_default_GainCF.calls.tsv"))
dim(mCAs_QC)  # 38321    22


mCAs_stats <- as.data.frame(fread("FinnGenBatchAll.stats.tsv"))
dim(mCAs_stats)  # 201458     20


# merged <- merge(mCAs, pheno, by=c(1), all.x=TRUE)
# dim(merged)   # 183,579    33
merged <- merge(mCAs, pheno, by=c(1))
dim(merged)     # 175,050    33




#------------------------------------------------------------------------------------------
#-------Sample QC and plots ---------------------------------------------------------------
#------------------------------------------------------------------------------------------

# Sample QC
BAF_Auto_toRemove <- mCAs_stats[which(mCAs_stats$baf_auto > 0.05),1] 
length(BAF_Auto_toRemove)   # 664 



# Plot of Baf_plot without any QC (Supplementary Figure 1-A)
pdf("FinnGenBatchAll.baf_plot.pdf", width=4, height=4)
if ('lrr_auto' %in% colnames(mCAs_stats)) { col_x <- 'lrr_auto'; lbl_x <- 'GC-adjusted LRR auto-correlation'} else if ('COV_AUTO' %in% colnames(mCAs_stats)) { col_x <- 'cov_auto'; lbl_x <- 'GC-adjusted coverage auto-correlation' }
p <- ggplot(mCAs_stats, aes_string(x=col_x, y='baf_auto', color='computed_gender')) +
  geom_point(size=.5, alpha=1/2) +
  scale_x_continuous(lbl_x) +
  scale_y_continuous('BAF auto-correlation') +
  scale_color_manual('', values=c('M'='blue', 'F'='orchid', 'U'='gray'), labels = c('M'='Male', 'F'='Female', 'U'='Undetermined')) +
  theme_bw(base_size = 14) +
  theme(legend.position= bottom', legend.box='horizontal')
print(p)
dev.off()



# Plot of genderMismatches without any QC (Supplementary Figure 1-B)
pheno[,"gender1"] <- toupper(substr(pheno$SEX,1,1))
mergedStatsPhenos <- merge(mCAs_stats[mCAs_stats$computed_gender!="U", ], pheno, by.x=c(1), by.y=c(1))
dim(mergedStatsPhenos)   # 194854     31
mergedStatsPhenos[,"MochaGender_PhenoGender_Mismatch"] <- ifelse((mergedStatsPhenos$computed_gender == "F" & mergedStatsPhenos$gender1=="M") | (mergedStatsPhenos$computed_gender == "M" & mergedStatsPhenos$gender1 == "F"), "Sex mismatch", "Sex match")
mismatch_Sex_ids <- mergedStatsPhenos[which(mergedStatsPhenos$MochaGender_PhenoGender_Mismatch == "Sex mismatch"),1]   # 84 samples

pdf("FinnGenBatchAll.genderMismatches.pdf", width = 9, height = 4)
p=ggplot(mergedStatsPhenos[mergedStatsPhenos$x_nonpar_lrr_median> -1.5 & !is.na(mergedStatsPhenos$MochaGender_PhenoGender_Mismatch),], aes(x=x_nonpar_lrr_median, y=y_nonpar_lrr_median, color=`computed_gender`)) + geom_point() + theme_bw(base_size = 14)+xlab("x_nonpar_lrr_median")+ylab("y_nonpar_lrr_median")  
pp <- p + facet_wrap( ~ MochaGender_PhenoGender_Mismatch, nrow = 1) + theme(legend.position = "right")
print(pp)
dev.off()



  
#------------------------------------------------------------------------------------------
#-------Variant QC and plots---------------------------------------------------------------
#------------------------------------------------------------------------------------------

# Basic QC for variants, including removal of likely germline variants (lod_baf_phase <20 for autosomal variants,lod_baf_phase<5 for ChrX, or annotated as germline CNP)
Var_toKeep <- merged[which(!(merged$sample_id %in% BAF_Auto_toRemove) & (merged$lod_baf_phase > 20 & merged$chrom!="chrX")|(merged$lod_baf_phase>5 & merged$chrom=="chrX") & merged$type != "CNP_Gain" & merged$type != "CNP_Loss"),]
# Var_toKeep <- Var_toKeep[-which(Var_toKeep$CHROM == "chr6" & Var_toKeep$BEG_GRCh38 >28510120 & Var_toKeep$END_GRCh38 <33480577),]  # 0 ########## any QC for MHC could be added here
nrow(Var_toKeep)   # 33,662

df_calls <- Var_toKeep
df_calls$SV <- factor(df_calls$chrom, levels = c('0-2 Mbp', '2-10 Mbp', '10-50 Mbp', '50-250 Mbp'))
df_calls$SV[df_calls$length <= 250e6] <- '50-250 Mbp'
df_calls$SV[df_calls$length <= 50e6] <- '10-50 Mbp'
df_calls$SV[df_calls$length <= 10e6] <- '2-10 Mbp'
df_calls$SV[df_calls$length <= 2e6] <- '0-2 Mbp'
df_calls$bdev[is.na(df_calls$bdev)] <- -.05



# Trident plot with basic QC (Supplementary Figure 2-A)
pdf("FinnGenBatchAll.init_plot.trident.noCFfilter.pdf", width = 6, height = 6)
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



# Trident plot of with basic QC (Supplementary Figure 2-B)
df_calls$chrom <- as.factor(gsub('^chr', '', gsub('^chrM', 'MT', df_calls$chrom)))
ord <- order(as.numeric(gsub('MT', '26', gsub('Y', '24', gsub('X', '23', levels(df_calls$chrom))))))
df_calls$chrom <- factor(df_calls$chrom, levels(df_calls$chrom)[ord])

pdf("FinnGenBatchAll.init_plot.CountsbyChrom.noCFfilter.pdf", width = 6, height = 6)
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
Var_toKeep_Filtered <- mCAs_QC[which(!(mCAs_QC$sample_id %in% BAF_Auto_toRemove) & (mCAs_QC$lod_baf_phase > 20 & mCAs_QC$chrom!="chrX")|(mCAs_QC$lod_baf_phase>5 & mCAs_QC$chrom=="chrX")),]
nrow(Var_toKeep_Filtered)   # 33,340
write.table(Var_toKeep_Filtered, "MoreQC.withoutCFfilter.txt", col.names = T, row.names = F, quote = F, sep = "\t")


df_calls = Var_toKeep_Filtered
df_calls$SV <- factor(df_calls$chrom, levels = c('0-2 Mbp', '2-10 Mbp', '10-50 Mbp', '50-250 Mbp'))
df_calls$SV[df_calls$length <= 250e6] <- '50-250 Mbp'
df_calls$SV[df_calls$length <= 50e6] <- '10-50 Mbp'
df_calls$SV[df_calls$length <= 10e6] <- '2-10 Mbp'
df_calls$SV[df_calls$length <= 2e6] <- '0-2 Mbp'
df_calls$bdev[is.na(df_calls$bdev)] <- -.05



# Trident plot with more QC (Supplementary Figure 2-C)
pdf("FinnGenBatchAll.init_plot.trident.filter.pdf", width = 6, height = 6)
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
print(p + scale_x_continuous('BAF deviation (Bdev)', breaks = c(0.0, 0.01, 0.02, 0.03, 0.04, 0.05)) + scale_y_continuous(expression(paste('Relative coverage (rescaled ', 2^LRR, ')')), breaks = c(1.8, 1.9, 2.0, 2.1, 2.2)) + coord_cartesian(xlim = c(0.0, 0.05), ylim = c(1.8, 2.2)))
dev.off()



# Trident plot with more QC (Supplementary Figure 2-D)
df_calls$chrom <- as.factor(gsub('^chr', '', gsub('^chrM', 'MT', df_calls$chrom)))
ord <- order(as.numeric(gsub('MT', '26', gsub('Y', '24', gsub('X', '23', levels(df_calls$chrom))))))
df_calls$chrom <- factor(df_calls$chrom, levels(df_calls$chrom)[ord])

pdf("FinnGenBatchAll.init_plot.CountsbyChrom.filter.pdf", width = 6, height = 6)
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

allSamples <- mCAs_stats[, c("sample_id","computed_gender")]
nrow(allSamples)         # 200,018

CNVs <- merge(allSamples, Var_toKeep_Filtered, by=c(1), all.x=TRUE)
nrow(CNVs)  # 207533
length(unique(CNVs$sample_id)) # 201458



# Define 4 expousures without considering CF
CNVs$has_MosaicCNV <- ifelse(!is.na(CNVs$chrom), 1, 0)

CNVs$ChrX <- ifelse(CNVs$has_MosaicCNV==0 & CNVs$computed_gender.x=="F", 0,
             ifelse(CNVs$chrom=="chrX" & CNVs$computed_gender.x=="F", 1, NA))

CNVs$ChrY <- ifelse(CNVs$has_MosaicCNV == 0 & CNVs$computed_gender.x %in% c("M", "U"), 0,
             ifelse(CNVs$chrom == "chrX" & CNVs$computed_gender.x%in% c("M", "U"), 1, NA))

CNVs$Autosomes <- ifelse(CNVs$has_MosaicCNV == 0, 0,
                  ifelse(CNVs$chrom %in% paste("chr", 1:22, sep=""), 1, NA))

sum(CNVs$has_MosaicCNV==1,na.rm=T)==sum(sum(CNVs$ChrX==1,na.rm=T), sum(CNVs$ChrY==1,na.rm=T), sum(CNVs$Autosomes==1,na.rm=T))  # 33340



# Define 4 expousures with CF>0.1
CNVs$cf <- as.numeric(as.character(CNVs$cf))

CNVs$Large_CNV_Clonev2 <- ifelse(CNVs$has_MosaicCNV==0, 0,
                          ifelse(CNVs$cf>=0.1, 1, 
                          ifelse(CNVs$has_MosaicCNV==1 & CNVs$cf<0.1 | is.na(CNVs$CF), NA, 0)))

CNVs$Large_ChrX <- ifelse(CNVs$has_MosaicCNV == 0 & CNVs$computed_gender.x == "F", 0,
                   ifelse(CNVs$Large_CNV_Clone == 1 & CNVs$ChrX == 1, 1, 
                   ifelse(CNVs$has_MosaicCNV == 1 & CNVs$Large_CNV_Clone != 1, NA, 0)))

CNVs$Large_ChrY <- ifelse(CNVs$has_MosaicCNV == 0 & CNVs$computed_gender.x %in% c("M", "U"), 0,
                   ifelse(CNVs$Large_CNV_Clone == 1 & CNVs$ChrY == 1, 1, 
                   ifelse(CNVs$has_MosaicCNV == 1 & CNVs$Large_CNV_Clone != 1, NA, 0)))

CNVs$Large_Autosomes <- ifelse(CNVs$has_MosaicCNV == 0, 0,
                        ifelse(CNVs$Large_CNV_Clone == 1 & CNVs$Autosomes == 1, 1, 
                        ifelse(CNVs$has_MosaicCNV == 1 & CNVs$Large_CNV_Clone != 1, NA, 0)))

sum(CNVs$Large_CNV_Clonev2==1,na.rm=T)==sum(sum(CNVs$Large_ChrX==1,na.rm=T), sum(CNVs$Large_ChrY==1,na.rm=T), sum(CNVs$Large_Autosomes==1,na.rm=T))  # 17,399



# QC and format exposure
colnames(CNVs)[1] <- "x"
max_DF <- data.frame(ourSid = uniqueSamples)
variation <- c("has_MosaicCNV", "Large_CNV_Clonev2", "ChrX", "ChrY", "Autosomes", "Large_ChrX", "Large_ChrY", "Large_Autosomes")
for (i in 1:length(variation)){
	indx <- which(colnames(CNVs)==variation[i])
	tmp <- aggregate(CNVs[,indx]~ x, data=CNVs, max)
	colnames(tmp)[2] <- variation[i]
	max_DF <- merge(max_DF, tmp, by=c(1), all.x=TRUE)
}

uniqueSamples <- unique(CNVs[-which(CNVs[,"x"] %in% c(BAF_Auto_toRemove)),"x"])
allSamples <- data.frame(sample_id = uniqueSamples)
allSamples <- merge(allSamples, mCAs_stats[, c("sample_id","computed_gender")], by=c(1), all.x=TRUE)

max_DF <- merge(max_DF, allSamples, by=c(1), all.x=TRUE)
length(max_DF)  # 200,794

write.table(max_DF, "mCA_01_Counts.afterQC.withoutCF.txt", col.names = T, row.names=F, quote=F, sep="\t")





#------------------------------------------------------------------------------------------
#--------- Making plots with age-----------------------------------------------------------
#------------------------------------------------------------------------------------------

mermdf <- merge(max_DF, pheno, by=c(1))
nrow(mermdf)   # 195,650
mermdf$Age <- mermdf$BL_AGE
mermdf$Sex <- ifelse(mermdf$SEX=="female", "Female", ifelse(mermdf$SEX=="male", "Male", NA))
mermdf$age_bins <- ifelse(mermdf$Age>=40 & mermdf$Age<45, "[40,45)",
                   ifelse(mermdf$Age>=45 & mermdf$Age<50, "[45,50)",
                   ifelse(mermdf$Age>=50 & mermdf$Age<55, "[50,55)",
                   ifelse(mermdf$Age>=55 & mermdf$Age<60, "[55,60)",
                   ifelse(mermdf$Age>=60 & mermdf$Age<65, "[60,65)",
                   ifelse(mermdf$Age>=65 & mermdf$Age<70, "[65,70)",
                   ifelse(mermdf$Age>=70 & mermdf$Age<75, "[70,75)",
                   ifelse(mermdf$Age>=75 & mermdf$Age<80, "[75,80)",
                   ifelse(mermdf$Age>=80 & mermdf$Age<85, "[80,85)",
                   ifelse(mermdf$Age>=85 , ">85", NA))))))))))
mermdf$age_bins <- ordered(mermdf$age_bins , levels = c("[45,50)","[50,55)", "[55,60)", "[60,65)","[65,70)","[70,75)","[75,80)","[80,85)",">85"))
mermdf$Sex <- ordered(mermdf$Sex, levels = c("Male", "Female"))



# Age  
Mosaic_Names <- c('has_MosaicCNV', 'Large_CNV_Clonev2', "Autosomes", "Large_Autosomes") 
FemaleMaleMosaics <- c('ChrX','Large_ChrX','ChrY','Large_ChrY')
mermdf <- mermdf[which(!is.na(mermdf$Sex) & !is.na(mermdf$age_bins)),]
nrow(mermdf)   # 131,677 since samples with age <40 were considred as NA and have beenfurther removed
 
# mean_cl_binom_tol <- function(x) {out <- binom.bayes(sum(x>0, na.rm=TRUE),sum(!is.na(x)),tol=0.65/10)[c("mean","lower","upper")]; names(out) <- c("y","ymin","ymax"); return(out)}
mean_cl_binom_tol <- function(x) {out <- binom.bayes(sum(x>0, na.rm=TRUE),sum(!is.na(x)))[c("mean","lower","upper")]; names(out) <- c("y","ymin","ymax"); return(out)}

for (i in Mosaic_Names){
	# pdf(paste0("FinnGenBatchAll.mCAAgeassoc.", i, ".pdf"), width = 5, height= 4)
	png(paste0("FinnGenBatchAll.mCAAgeassoc", i, ".png"), width = 4, height= 4, units="in", res=300)
	p4 <- ggplot(mermdf, aes(x=as.numeric(age_bins), y=mermdf[,i], linetype=Sex)) +
	      stat_summary(fun=mean, geom="line", size=1.7) + 
	      stat_summary(fun.data=mean_cl_binom_tol, geom="ribbon", alpha=0.25, color="transparent", fill="gray") +
	      xlab("Age (years)") + 
	      ylab("Proportion") + theme_minimal() + guides(colour = FALSE)+scale_x_continuous(breaks=1:length(levels(mermdf$age_bins)),labels=levels(mermdf$age_bins)) +
	      theme(axis.text.x=element_text(size=14, angle=45, hjust=1), strip.text.y=element_text(size=14), axis.text.y=element_text(size=14, hjust=1, color="black"), axis.ticks=element_line(colour="black"), 
	            axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),
	            legend.text=element_text(size=14),legend.position="right",legend.key = element_blank(),legend.title = element_blank(),legend.background = element_blank(), strip.text = element_text(size=14)) 
	print(p4)
	dev.off()
}


for (i in FemaleMaleMosaics){
	# pdf(paste0("FinnGenBatchAll.mCAAgeassoc", i, ".pdf"), width = 4, height= 4)
	png(paste0("FinnGenBatchAll.mCAAgeassoc", i, ".png"), width = 4, height= 4, units="in", res=300)
	p4 <- ggplot(mermdf, aes(x=as.numeric(age_bins), y=mermdf[,i])) +
	      stat_summary(fun=mean,geom="line", size=1.7) +
          stat_summary(fun.data=mean_cl_binom_tol, geom="ribbon",alpha=0.25,color="transparent",fill="gray") +
          xlab("Age (years)") + 
          ylab("Proportion") + theme_minimal() + guides(colour = FALSE)+scale_x_continuous(breaks=1:length(levels(mermdf$age_bins)),labels=levels(mermdf$age_bins)) +
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
  



