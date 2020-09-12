# Association between COVID-19 and mCA events
setwd("/home/aoxliu/mCA/input/dsge-aoxing/mocha/covid")

library(data.table)

'%!in%' <- function(x,y)!('%in%'(x,y))




#------------------------------------------------------------------------------------
#---Loading data---------------------------------------------------------------------
#------------------------------------------------------------------------------------

max_DF <- read.table("FinnGenBatchAll_nodeath_nohemeCancer.mCA_01_Counts.afterQC.withoutCF.txt", header=T)  # basic QC has been done (remove gender missmatch, first and second degree relatives, died_before_COVID(All death events are censored at 31.12.2018), and Hematologic_Cancer before genotyping)
max_DF <- max_DF[, !(names(max_DF) %in% c("SEX","BL_AGE"))]
dim(max_DF)   # 175,690     11


covariates <- read.table("finngen_R6_pheno_cov.txt", header=T, stringsAsFactors=F)
dim(covariates)       # 260405     24
length(unique(covariates$FINNGENID))  # 260,405


covid <- read.table("finngen_R6_infectious_disease_register_corona.txt", sep="\t", header=T, stringsAsFactors=F)
dim(covid)


pheno <- merge(max_DF, covariates, by=c(1))
dim(pheno)   # 175690     32
pheno$BL_AGE_2 <- as.numeric(as.character(pheno$BL_AGE))^2
pheno$SMOKE <- ifelse(pheno$SMOKE3=="never",0,1)
pheno$covid_P <- ifelse(pheno$ourSid %in% covid$FINNGENID,1,0)
pheno$covid_H <- ifelse(pheno$covid_P==1 & (pheno$ourSid %!in% covid[covid$Hospital_treatment==1,"FINNGENID"]),NA,pheno$covid_P)
table(pheno$covid_P)
#      0      1 
# 175408    282 
table(pheno$covid_H)
#      0      1 
# 175408     69 
sum(is.na(pheno$covid_P))  # 0
sum(is.na(pheno$covid_H))  # 213
colnames(pheno)



#------------------------------------------------------------------------------------
#---Loading data---------------------------------------------------------------------
#------------------------------------------------------------------------------------

summaryDF <- data.frame()
outc <- c("covid_P", "covid_H")
CNV_x <- c("has_MosaicCNV", "Large_CNV_Clonev2", "Autosomes", "Large_Autosomes")
popn <- "ALL"

for (i in 1:length(outc)){
	indx <- which(colnames(pheno)==outc[i])
	
	for(j in 1:length(CNV_x)){
		CNV_cur <- which(colnames(pheno)==CNV_x[j])	
		temp <- pheno
		temp$mCA <- temp[,CNV_cur]
		print(paste0(i,":",j))

		# df <- as.data.frame(t(as.data.frame(summary(glm(temp[,indx] ~ mCA +SEX +SMOKE +BL_AGE +BL_AGE_2 +J10_ASTHMA +J10_COPD +I9_ISCHHEART +T2D +C3_CANCER +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10, family="binomial", data=temp))$coeff[2,])))		
		mylogit <- glm(temp[,indx] ~ mCA +SEX +SMOKE +BL_AGE +BL_AGE_2 +J10_ASTHMA +J10_COPD +I9_ISCHHEART +T2D +C3_CANCER +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10, family="binomial", data=temp)
		df <- as.data.frame(t(as.data.frame(summary(mylogit)$coeff[2,])))
		df$OR <- exp(df[1,1])
		ci <- exp(confint(mylogit,"mCA"))
		df$OR_25 <- ci[1]
		df$OR_975 <- ci[2]				
		df$N_Cases <- length(which(temp[,indx]==1 & !is.na(temp[,CNV_cur]) & !is.na(temp$SEX) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$J10_ASTHMA) & !is.na(temp$J10_COPD) & !is.na(temp$I9_ISCHHEART) & !is.na(temp$T2D) & !is.na(temp$C3_CANCER) & !is.na(temp$PC1)))
		df$N_Controls <- length(which(temp[,indx]==0 & !is.na(temp[,CNV_cur]) & !is.na(temp$SEX) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$J10_ASTHMA) & !is.na(temp$J10_COPD) & !is.na(temp$I9_ISCHHEART) & !is.na(temp$T2D) & !is.na(temp$C3_CANCER) & !is.na(temp$PC1)))
		df$N_Cases_withVar <- length(which(temp[,indx]==1 & temp[,CNV_cur]==1 & !is.na(temp$SEX) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$J10_ASTHMA) & !is.na(temp$J10_COPD) & !is.na(temp$I9_ISCHHEART) & !is.na(temp$T2D) & !is.na(temp$C3_CANCER) & !is.na(temp$PC1)))
		df$N_Controls_withVar <- length(which(temp[,indx]==0 & temp[,CNV_cur]==1 & !is.na(temp$SEX) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$J10_ASTHMA) & !is.na(temp$J10_COPD) & !is.na(temp$I9_ISCHHEART) & !is.na(temp$T2D) & !is.na(temp$C3_CANCER) & !is.na(temp$PC1)))
		df$y <- outc[i]
		df$x <- CNV_x[j]
		df$population <- popn
		df$Adjustment <- "Fully_Adjusted"
		summaryDF <- rbind(summaryDF, df)

		# df <- as.data.frame(t(as.data.frame(summary(glm(temp[,indx] ~ mCA +SEX +BL_AGE +BL_AGE_2 +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10, family="binomial", data=temp))$coeff[2,])))			
		mylogit <- glm(temp[,indx] ~ mCA +SEX +BL_AGE +BL_AGE_2 +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10, family="binomial", data=temp)
		df <- as.data.frame(t(as.data.frame(summary(mylogit)$coeff[2,])))
		df$OR <- exp(df[1,1])
		ci <- exp(confint(mylogit,"mCA"))
		df$OR_25 <- ci[1]
		df$OR_975 <- ci[2]				
		df$N_Cases <- length(which(temp[,indx]==1 & !is.na(temp[,CNV_cur]) & !is.na(temp$SEX) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
		df$N_Controls <- length(which(temp[,indx]==0 & !is.na(temp[,CNV_cur]) & !is.na(temp$SEX) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
		df$N_Cases_withVar <- length(which( temp[,indx]==1 & temp[,CNV_cur]==1 & !is.na(temp$SEX) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
		df$N_Controls_withVar <- length(which( temp[,indx]==0 & temp[,CNV_cur]==1 & !is.na(temp$SEX) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
		df$y <- outc[i]
		df$x <- CNV_x[j]
		df$population <- popn
		df$Adjustment <- "Sparsely_Adjusted"
		summaryDF <- rbind(summaryDF, df)

		# df <- as.data.frame(t(as.data.frame(summary(glm(temp[,indx] ~ mCA +SEX +SMOKE +BL_AGE +BL_AGE_2 +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10, family="binomial", data=temp))$coeff[2,])))
		mylogit <- glm(temp[,indx] ~ mCA +SEX +SMOKE +BL_AGE +BL_AGE_2 +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10, family="binomial", data=temp)
		df <- as.data.frame(t(as.data.frame(summary(mylogit)$coeff[2,])))
		df$OR <- exp(df[1,1])
		ci <- exp(confint(mylogit,"mCA"))
		df$OR_25 <- ci[1]
		df$OR_975 <- ci[2]				
		df$N_Cases <- length(which(temp[,indx]==1 & !is.na(temp[,CNV_cur]) & !is.na(temp$SEX) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
		df$N_Controls <- length(which(temp[,indx]==0 & !is.na(temp[,CNV_cur]) & !is.na(temp$SEX) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
		df$N_Cases_withVar <- length(which( temp[,indx]==1 & temp[,CNV_cur]==1 & !is.na(temp$SEX) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
		df$N_Controls_withVar <- length(which( temp[,indx]==0 & temp[,CNV_cur]==1 & !is.na(temp$SEX) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
		df$y <- outc[i]
		df$x <- CNV_x[j]
		df$population <- popn
		df$Adjustment <- "Sparsely_Adjusted_SMOKE"
		summaryDF <- rbind(summaryDF, df)

		# df <- as.data.frame(t(as.data.frame(summary(glm(temp[,indx] ~ mCA, family="binomial", data=temp))$coeff[2,])))
		mylogit <- glm(temp[,indx] ~ mCA, family="binomial", data=temp)
		df <- as.data.frame(t(as.data.frame(summary(mylogit)$coeff[2,])))
		df$OR <- exp(df[1,1])
		ci <- exp(confint(mylogit,"mCA"))
		df$OR_25 <- ci[1]
		df$OR_975 <- ci[2]						
		df$N_Cases <- length(which( temp[,indx]==1 & !is.na(temp[,CNV_cur])))
		df$N_Controls <- length(which(temp[,indx]==0 & !is.na(temp[,CNV_cur])))
		df$N_Cases_withVar <- length(which( temp[,indx]==1 & temp[,CNV_cur]==1))
		df$N_Controls_withVar <- length(which( temp[,indx]==0 & temp[,CNV_cur]==1))
		df$y <- outc[i]
		df$x <- CNV_x[j]
		df$population <- popn
		df$Adjustment <- "Unadjusted"
		summaryDF <- rbind(summaryDF, df)
	}
}



for (popn in c("male","female")){
	if (popn=="male"){
		CNV_x <- c("has_MosaicCNV", "Large_CNV_Clonev2", "Autosomes", "Large_Autosomes", "ChrY", "Large_ChrY")
	} else {
		CNV_x <- c("has_MosaicCNV", "Large_CNV_Clonev2", "Autosomes", "Large_Autosomes", "ChrX", "Large_ChrX")
	}
	
	for (i in 1:length(outc)){
		indx <- which(colnames(pheno)==outc[i])
		
		for(j in 1:length(CNV_x)){
			CNV_cur <- which(colnames(pheno)==CNV_x[j])
			temp <- pheno[pheno$SEX==popn,]
			temp$mCA <- temp[,CNV_cur]
			print(paste0(i,":",j))
			
			# df <- as.data.frame(t(as.data.frame(summary(glm(temp[,indx] ~ mCA +SMOKE +BL_AGE +BL_AGE_2 +J10_ASTHMA +J10_COPD +I9_ISCHHEART +T2D +C3_CANCER +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10, family="binomial", data=temp))$coeff[2,])))
			mylogit <- glm(temp[,indx] ~ mCA +SMOKE +BL_AGE +BL_AGE_2 +J10_ASTHMA +J10_COPD +I9_ISCHHEART +T2D +C3_CANCER +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10, family="binomial", data=temp)
			df <- as.data.frame(t(as.data.frame(summary(mylogit)$coeff[2,])))
			df$OR <- exp(df[1,1])
			ci <- exp(confint(mylogit,"mCA"))
			df$OR_25 <- ci[1]
			df$OR_975 <- ci[2]				
			df$N_Cases <- length(which(temp[,indx]==1 & !is.na(temp[,CNV_cur]) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$J10_ASTHMA) & !is.na(temp$J10_COPD) & !is.na(temp$I9_ISCHHEART) & !is.na(temp$T2D) & !is.na(temp$C3_CANCER) & !is.na(temp$PC1)))
			df$N_Controls <- length(which(temp[,indx]==0 & !is.na(temp[,CNV_cur]) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$J10_ASTHMA) & !is.na(temp$J10_COPD) & !is.na(temp$I9_ISCHHEART) & !is.na(temp$T2D) & !is.na(temp$C3_CANCER) & !is.na(temp$PC1)))
			df$N_Cases_withVar <- length(which(temp[,indx]==1 & temp[,CNV_cur]==1 & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$J10_ASTHMA) & !is.na(temp$J10_COPD) & !is.na(temp$I9_ISCHHEART) & !is.na(temp$T2D) & !is.na(temp$C3_CANCER) & !is.na(temp$PC1)))
			df$N_Controls_withVar <- length(which(temp[,indx]==0 & temp[,CNV_cur]==1 & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$J10_ASTHMA) & !is.na(temp$J10_COPD) & !is.na(temp$I9_ISCHHEART) & !is.na(temp$T2D) & !is.na(temp$C3_CANCER) & !is.na(temp$PC1)))
			df$y <- outc[i]
			df$x <- CNV_x[j]
			df$population <- popn
			df$Adjustment <- "Fully_Adjusted"
			summaryDF <- rbind(summaryDF, df)
	
			# df <- as.data.frame(t(as.data.frame(summary(glm(temp[,indx] ~ mCA +BL_AGE +BL_AGE_2 +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10, family="binomial", data=temp))$coeff[2,])))	
			mylogit <- glm(temp[,indx] ~ mCA +BL_AGE +BL_AGE_2 +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10, family="binomial", data=temp)
			df <- as.data.frame(t(as.data.frame(summary(mylogit)$coeff[2,])))
			df$OR <- exp(df[1,1])
			ci <- exp(confint(mylogit,"mCA"))
			df$OR_25 <- ci[1]
			df$OR_975 <- ci[2]								
			df$N_Cases <- length(which(temp[,indx]==1 & !is.na(temp[,CNV_cur]) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
			df$N_Controls <- length(which(temp[,indx]==0 & !is.na(temp[,CNV_cur]) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
			df$N_Cases_withVar <- length(which( temp[,indx]==1 & temp[,CNV_cur]==1 & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
			df$N_Controls_withVar <- length(which( temp[,indx]==0 & temp[,CNV_cur]==1 & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
			df$y <- outc[i]
			df$x <- CNV_x[j]
			df$population <- popn
			df$Adjustment <- "Sparsely_Adjusted"
			summaryDF <- rbind(summaryDF, df)
	
			df <- as.data.frame(t(as.data.frame(summary(glm(temp[,indx] ~ mCA +SMOKE +BL_AGE +BL_AGE_2 +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10, family="binomial", data=temp))$coeff[2,])))	
			mylogit <- glm(temp[,indx] ~ mCA +SMOKE +BL_AGE +BL_AGE_2 +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10, family="binomial", data=temp)
			df <- as.data.frame(t(as.data.frame(summary(mylogit)$coeff[2,])))
			df$OR <- exp(df[1,1])
			ci <- exp(confint(mylogit,"mCA"))
			df$OR_25 <- ci[1]
			df$OR_975 <- ci[2]										
			df$N_Cases <- length(which(temp[,indx]==1 & !is.na(temp[,CNV_cur]) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
			df$N_Controls <- length(which(temp[,indx]==0 & !is.na(temp[,CNV_cur]) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
			df$N_Cases_withVar <- length(which( temp[,indx]==1 & temp[,CNV_cur]==1 & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
			df$N_Controls_withVar <- length(which( temp[,indx]==0 & temp[,CNV_cur]==1 & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
			df$y <- outc[i]
			df$x <- CNV_x[j]
			df$population <- popn
			df$Adjustment <- "Sparsely_Adjusted_SMOKE"
			rownames(df) <- c("mCA")
			summaryDF <- rbind(summaryDF, df)
	
			df <- as.data.frame(t(as.data.frame(summary(glm(temp[,indx] ~ mCA, family="binomial", data=temp))$coeff[2,])))
			mylogit <- glm(temp[,indx] ~ mCA, family="binomial", data=temp)
			df <- as.data.frame(t(as.data.frame(summary(mylogit)$coeff[2,])))
			df$OR <- exp(df[1,1])
			ci <- exp(confint(mylogit,"mCA"))
			df$OR_25 <- ci[1]
			df$OR_975 <- ci[2]										
			df$N_Cases <- length(which( temp[,indx]==1 & !is.na(temp[,CNV_cur])))
			df$N_Controls <- length(which(temp[,indx]==0 & !is.na(temp[,CNV_cur])))
			df$N_Cases_withVar <- length(which( temp[,indx]==1 & temp[,CNV_cur]==1))
			df$N_Controls_withVar <- length(which( temp[,indx]==0 & temp[,CNV_cur]==1))
			df$y <- outc[i]
			df$x <- CNV_x[j]
			df$population <- popn
			df$Adjustment <- "Unadjusted"
			rownames(df) <- c("mCA")
			summaryDF <- rbind(summaryDF, df)
		}
	}
}

# write.table(summaryDF, "COVID_Mosaic_Associations_Smoke.txt", col.names=T, row.names=F, quote=F, sep="\t")
# write.table(summaryDF[summaryDF[,4]<0.05,], "COVID_Mosaic_Associations_Smoke_Pval005.txt", col.names=T, row.names=F, quote=F, sep="\t")




#------------------------------------------------------------------------------------
#---Format the results---------------------------------------------------------------------
#------------------------------------------------------------------------------------


summaryDF[,"Pheno_Grouping"] <- ifelse(summaryDF[,"y"]=="covid_H",1, 
                                ifelse(summaryDF[,"y"]=="covid_P",2, NA))
                                
summaryDF[,"Mosaic_Variant_Grouping"] <- ifelse(summaryDF[,"x"] %in% c("has_MosaicCNV","Large_CNV_Clonev2"), 1,
                                         ifelse(summaryDF[,"x"] %in% c("Autosomes","Large_Autosomes"), 2,
                                         ifelse(summaryDF[,"x"] %in% c("ChrX","Large_ChrX"), 3,
                                         ifelse(summaryDF[,"x"] %in% c("ChrY","Large_ChrY"), 4, NA))))

summaryDF[,"Clone_size_Filter"] <- ifelse(summaryDF[,"x"] %in% c("has_MosaicCNV","Autosomes","ChrX","ChrY"), 1,
                                   ifelse(summaryDF[,"x"] %in% c("Large_CNV_Clonev2","Large_Autosomes","Large_ChrX","Large_ChrY"), 2, NA))

summaryDF[,"Popn"] <- ifelse(summaryDF[,"population"]=="ALL", 1,
	                  ifelse(summaryDF[,"population"]=="male", 2,
                      ifelse(summaryDF[,"population"]=="female", 3, NA)))
                      
summaryDF[,"Model"] <- ifelse(summaryDF[,"Adjustment"]=="Unadjusted", 1,
	                  ifelse(summaryDF[,"Adjustment"]=="Sparsely_Adjusted", 2,
                      ifelse(summaryDF[,"Adjustment"]=="Fully_Adjusted", 3, 
                       ifelse(summaryDF[,"Adjustment"]=="Sparsely_Adjusted_SMOKE", 4, NA))))
                    
summaryDF[,"OR"] <- round(summaryDF[,"OR"],3)                                        
summaryDF[,"OR_CI"] <- paste0(round(summaryDF[,"OR_25"],3),"-",round(summaryDF[,"OR_975"],3))
summaryDF[,"Beta"] <- summaryDF[,"Estimate"] 
summaryDF[,"SE"] <- summaryDF[,"Std. Error"] 
summaryDF[,"P"] <- summaryDF[,"Pr(>|z|)"] 
summaryDF[,"N_Cases"] <- summaryDF[,"N_Cases"]
summaryDF[,"N_ctls"] <- summaryDF[,"N_Controls"]
summaryDF[,"N_Case_with_variant"] <- summaryDF[,"N_Cases_withVar"]
summaryDF[,"N_Ctls_with_variant"] <- summaryDF[,"N_Controls_withVar"]


res <- summaryDF[,c("Pheno_Grouping", "Mosaic_Variant_Grouping", "Clone_size_Filter", "Popn", "Model", "OR", "OR_CI", "Beta", "SE", "P", "N_Cases", "N_ctls", "N_Case_with_variant","N_Ctls_with_variant")]
res <- res[order(res$Pheno_Grouping, res$Mosaic_Variant_Grouping,res$Clone_size_Filter,res$Popn,res$Model),] 

write.csv(res, "COVID_Mosaic_Associations_FORMAT.csv", row.names=F, quote=F)
write.csv(res[res[,"P"]<0.05,], "COVID_Mosaic_Associations_FORMAT_Pval005.csv", row.names=F, quote=F)

