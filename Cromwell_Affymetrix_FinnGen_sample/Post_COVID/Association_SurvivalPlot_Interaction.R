# Association between Infection and mCA events (Using "aoxing-mocha" VM)

# mount the data from Google bucket to VM instance
# cd  /home/aoxliu/mCA/input
# gcsfuse --implicit-dirs  from-fg-datateam  from-fg-datateam 


setwd("/home/aoxliu/mCA/input/dsge-aoxing/mocha/covid")

library(data.table)
# install.packages("survival")
library("survival")
# install.packages("cmprsk")
library(cmprsk)
# install.packages("survminer")
library(survminer)


'%!in%' <- function(x,y)!('%in%'(x,y))




#------------------------------------------------------------------------------------
#---Loading data---------------------------------------------------------------------
#------------------------------------------------------------------------------------

pheno <- read.table("INFECTION_DATA.txt", header=T)

infect <- read.table("Infection_endpoints.tsv", header=T, stringsAsFactors=F)
head(infect)
colnames(infect)   # "Organ_System", "Category", "FinnGen_Endpoint"

Organ <- unique(infect$Organ_System)
Cate <- unique(infect$Category)
grp <- unique(c(Organ,Cate,"EVERY"))



#------------------------------------------------------------------------------------
#---Making CoxPH plots---------------------------------------------------------------
#------------------------------------------------------------------------------------


pheno$mCACancer <- ifelse(pheno[,"Large_Autosomes"]==1 & pheno[,"C3_CANCER"]==1, "+mCA +Cancer", 
				   ifelse(pheno[,"Large_Autosomes"]==0 & pheno[,"C3_CANCER"]==1, "-mCA +Cancer",
				   ifelse(pheno[,"Large_Autosomes"]==1 & pheno[,"C3_CANCER"]==0, "+mCA -Cancer",
				   ifelse(pheno[,"Large_Autosomes"]==0 & pheno[,"C3_CANCER"]==0, "-mCA -Cancer", NA))))
				   
pheno$mCACancer <- ordered(pheno$mCACancer, levels=c("+mCA +Cancer", "-mCA +Cancer","+mCA -Cancer", "-mCA -Cancer"))
pheno <- pheno[order(pheno$mCACancer),]


# Sepsis, Pneumonia, and Digestive_System
tp <- pheno[pheno$BL_YEAR>=1989,]
CNV_x <- c("Large_Autosomes")
grp_N <- (1:length(grp))[grp %in% c("Sepsis","Pneumonia","Digestive_System")]

for (i in grp_N){
	tmp <- tp[is.na(tp[,grp[i]])==F & is.na(tp[,paste0(grp[i],"_ftime")])==F & is.na(tp[,paste0(grp[i],"_incid")])==F & is.na(tp[,paste0(grp[i],"_preval")])==0, ]	
	tmp[,"mCA"] <- tmp[,"Large_Autosomes"]	
	tmp[,"ftime"] <- tmp[,paste0(grp[i],"_ftime")]	
	tmp[,"incid"] <- tmp[,paste0(grp[i],"_incid")]	
	nrow(tmp)  # 140735
	temp <- tmp[!is.na(tmp$ftime) & !is.na(tmp$incid) & !is.na(tmp$mCACancer),]
	nrow(temp) # 124462
	p <- ggsurvplot(
		fit = survfit(Surv(ftime, incid) ~ mCACancer, data=temp), 
		xlab = "Years", 
		ylab = "Cumulative Proportion",
		break.time.by=2.5,
		conf.int = FALSE, font.x=14, font.y=14, font.tickslab=14,
		risk.table = TRUE, risk.table.y.text.col = TRUE,
		risk.table.col = "strata",  # Change risk table color by groups
		linetype = "strata",        # Change line type by groups
		ggtheme = theme_bw(),       # Change ggplot2 theme
		palette = c('red', 'purple', 'royalblue', 'springgreen4'), fun = "event", font.legend=14, font.xtickslab=14, fontsize=5,
		ylim = c(0, 0.15), xlim=c(-0.5,6.5), legend="none" ,legend.labs =c("+mCA +Cancer",  "-mCA +Cancer","+mCA -Cancer",  "-mCA -Cancer")
	)
	p$table <- ggpubr::ggpar(p$table, font.x = 14, font.y=14, font.xtickslab=14, font.ytickslab=14)
	pdf(paste("FinnGen_",grp[i],".Everyone.ggsurvplot.pdf",sep=""), width = 7, height= 7)
	p
	dev.off()
}



#------------------------------------------------------------------------------------
#---Survivial analysis using Cox proportional hazards model--------------------------
#------------------------------------------------------------------------------------
summaryDF <- data.frame()

for (popn in c("ALL","male","female")){
	if (popn=="male"){
		tp <- pheno[pheno$BL_YEAR>=1989 & pheno$SEX==popn,]
		CNV_x <- c("has_MosaicCNV", "Large_CNV_Clonev2", "Autosomes", "Large_Autosomes", "ChrY", "Large_ChrY")
		# grp_N <- (1:length(grp))[grp %!in% c("Female_Gynecological_Infections","Pregnancy_childbirth_and_puerperial_infections")]
		grp_N <- (1:length(grp))[grp %in% Organ]
	} else if (popn=="female"){
		tp <- pheno[pheno$BL_YEAR>=1989 & pheno$SEX==popn,]
		CNV_x <- c("has_MosaicCNV", "Large_CNV_Clonev2", "Autosomes", "Large_Autosomes", "ChrX", "Large_ChrX")
		# grp_N <- (1:length(grp))[grp %!in% c("Male_Genitourinary_Infections")]
		grp_N <- (1:length(grp))[grp %in% Organ]
	} else if (popn=="ALL"){
		tp <- pheno[pheno$BL_YEAR>=1989,]
		CNV_x <- c("has_MosaicCNV", "Large_CNV_Clonev2", "Autosomes", "Large_Autosomes")
		# grp_N <- (1:length(grp))[grp %!in% c("Female_Gynecological_Infections","Male_Genitourinary_Infections","Pregnancy_childbirth_and_puerperial_infections")]
		grp_N <- (1:length(grp))[grp %in% Organ]
	}

	for (i in grp_N){
		tmp <- tp[is.na(tp[,grp[i]])==F & is.na(tp[,paste0(grp[i],"_ftime")])==F & is.na(tp[,paste0(grp[i],"_incid")])==F & is.na(tp[,paste0(grp[i],"_preval")])==0, ]	

		for (j in 1:length(CNV_x)){	
			tmp[,"mCA"] <- tmp[,CNV_x[j]]	
			tmp[,"ftime"] <- tmp[,paste0(grp[i],"_ftime")]	
			tmp[,"incid"] <- tmp[,paste0(grp[i],"_incid")]			
			tmp[,"PriorCancer"] <- ifelse(tmp$C3_CANCER_AGE < tmp[,paste0(grp[i],"_AGE")],1,0)
			temp <- tmp
			for (k in c("SMOKE_YES")){
				print(paste0(popn,": ",i,": ",j,": ",k))
				if(k=="SMOKE_YES" & popn=="ALL"){
					mycoxph <- coxph(Surv(ftime, incid) ~ mCA*PriorCancer +SEX +BL_AGE +BL_AGE_2 +SMOKE +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10, data=temp)	
					df <- as.data.frame(t(as.data.frame(summary(mycoxph)$coeff["mCA",])))				
					df$N_Cases <- length(which(temp[,"incid"]==1 & !is.na(temp[,"mCA"]) & !is.na(temp$SEX) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
					df$N_Controls <- length(which(temp[,"incid"]==0 & !is.na(temp[,"mCA"]) & !is.na(temp$SEX) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
					df$N_Cases_withVar <- length(which(temp[,"incid"]==1 & temp[,"mCA"]==1 & !is.na(temp$SEX) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
					df$N_Controls_withVar <- length(which(temp[,"incid"]==0 & temp[,"mCA"]==1 & !is.na(temp$SEX) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
				} else if (k=="SMOKE_YES" & popn!="ALL"){
					if (length(which(temp[,"incid"]==1 & temp[,"mCA"]==1 & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))>0){
						mycoxph <- coxph(Surv(ftime, incid) ~ mCA*PriorCancer +BL_AGE +BL_AGE_2 +SMOKE +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10, data=temp)	
						df <- as.data.frame(t(as.data.frame(summary(mycoxph)$coeff["mCA",])))
					} else {
						df[,1:5] <- NA
					}
					df$N_Cases <- length(which(temp[,"incid"]==1 & !is.na(temp[,"mCA"]) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
					df$N_Controls <- length(which(temp[,"incid"]==0 & !is.na(temp[,"mCA"]) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
					df$N_Cases_withVar <- length(which(temp[,"incid"]==1 & temp[,"mCA"]==1 & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
					df$N_Controls_withVar <- length(which(temp[,"incid"]==0 & temp[,"mCA"]==1 & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
				} 

				df[,"mod"] <- k
				df$HR_25 <- exp(confint(mycoxph,"mCA"))[1]
				df$HR_975 <- exp(confint(mycoxph,"mCA"))[2]
				df[,"y"] <- grp[i]
				df[,"x"] <- CNV_x[j]
				df[,"Population"] <- popn
				summaryDF <- rbind(summaryDF, df)
			}
		}
	}
}



#------------------------------------------------------------------------------------
#---Format the results---------------------------------------------------------------------
#------------------------------------------------------------------------------------

summaryDF[,"Pheno_Grouping"] <- ifelse(summaryDF[,"y"] %in% Organ, 2, 
								ifelse(summaryDF[,"y"]=="ALL",3, 
                                ifelse(summaryDF[,"y"]=="EVERY",4,
                            	ifelse(summaryDF[,"y"] %in% Cate, 1, NA))))
                            	                            	
summaryDF[,"Phenotype"] <- ifelse(summaryDF[,"y"]=="ALL", "All_AB",
                           ifelse(summaryDF[,"y"]=="EVERY", "All_AB_EXCEL", summaryDF[,"y"]))

summaryDF[,"Mosaic_Variant_Grouping"] <- ifelse(summaryDF[,"x"] %in% c("has_MosaicCNV","Large_CNV_Clonev2"), 1,
                                         ifelse(summaryDF[,"x"] %in% c("Autosomes","Large_Autosomes"), 2,
                                         ifelse(summaryDF[,"x"] %in% c("ChrX","Large_ChrX"), 3,
                                         ifelse(summaryDF[,"x"] %in% c("ChrY","Large_ChrY"), 4, NA))))

summaryDF[,"Clone_size_Filter"] <- ifelse(summaryDF[,"x"] %in% c("has_MosaicCNV","Autosomes","ChrX","ChrY"), 1,
                                   ifelse(summaryDF[,"x"] %in% c("Large_CNV_Clonev2","Large_Autosomes","Large_ChrX","Large_ChrY"), 2, NA))

summaryDF[,"Sex_Popn"] <- ifelse(summaryDF[,"Population"]=="ALL", 1,
	                      ifelse(summaryDF[,"Population"]=="male", 2,
                          ifelse(summaryDF[,"Population"]=="female", 3, NA)))                         

summaryDF[,"Model"] <- summaryDF[,"mod"]
summaryDF[,"HR"] <- ifelse(summaryDF[,"N_Cases_withVar"]==0, NA, round(summaryDF[,"exp(coef)"],3))
summaryDF[,"HR_CI"] <- ifelse(summaryDF[,"N_Cases_withVar"]==0, NA, paste0(round(summaryDF[,"HR_25"],3),"-",round(summaryDF[,"HR_975"],3)))
summaryDF[,"Beta"] <- ifelse(summaryDF[,"N_Cases_withVar"]==0, NA, summaryDF[,"coef"])
summaryDF[,"SE"] <- ifelse(summaryDF[,"N_Cases_withVar"]==0, NA, summaryDF[,"se(coef)"])
summaryDF[,"P"] <- ifelse(summaryDF[,"N_Cases_withVar"]==0, NA, summaryDF[,"Pr(>|z|)"])
summaryDF[,"N_Cases"] <- summaryDF[,"N_Cases"]
summaryDF[,"N_ctls"] <- summaryDF[,"N_Controls"]
summaryDF[,"N_Case_with_variant"] <- summaryDF[,"N_Cases_withVar"]
summaryDF[,"N_Ctls_with_variant"] <- summaryDF[,"N_Controls_withVar"]


res <- summaryDF[,c("Pheno_Grouping", "Phenotype", "Mosaic_Variant_Grouping", "Clone_size_Filter", "Sex_Popn", "Model", "HR", "HR_CI", "Beta", "SE", "P", "N_Cases", "N_ctls", "N_Case_with_variant","N_Ctls_with_variant")]
res <- res[order(res$Pheno_Grouping, res$Phenotype, res$Mosaic_Variant_Grouping,res$Clone_size_Filter,res$Sex_Popn,res$Model),] 
res_sig <- res[res[,"P"]<0.05 & is.na(res[,"P"])==F,]
write.csv(res, "INFECTION_Organ_mcaInteractCancer_Mosaic_Associations_FORMAT.csv", row.names=F, quote=F)







