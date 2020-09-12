# Association between Infection and mCA events (Using "aoxing-mocha" VM)

# mount the data from Google bucket to VM instance
# cd  /home/aoxliu/mCA/input
# gcsfuse --implicit-dirs  from-fg-datateam  from-fg-datateam 


setwd("/home/aoxliu/mCA/input/dsge-aoxing/mocha/covid")

library(data.table)
# install.packages("survival")
library("survival")

'%!in%' <- function(x,y)!('%in%'(x,y))




#------------------------------------------------------------------------------------
#---Loading data---------------------------------------------------------------------
#------------------------------------------------------------------------------------

max_DF <- read.table("FinnGenBatchAll_nodeath_nohemeCancer.mCA_01_Counts.afterQC.withoutCF.txt", header=T)  # basic QC has been done (remove gender missmatch, first and second degree relatives, died_before_COVID(All death events are censored at 31.12.2018), and Hematologic_Cancer before genotyping)
max_DF <- max_DF[, !(names(max_DF) %in% c("SEX","BL_AGE"))]
dim(max_DF)   # 175,690     11


infect <- read.table("Infection_endpoints.tsv", header=T, stringsAsFactors=F)
head(infect)
colnames(infect)   # "Organ_System", "Category", "FinnGen_Endpoint"


covariates <- read.table("finngen_R6_pheno_cov.txt", header=T, stringsAsFactors=F)
dim(covariates)       # 260405     24
head(covariates)


phe_1 <- read.table("INFECT_1.txt", sep=" ", header=T)
phe_1 <- phe_1[, !(names(phe_1) %in% c("BL_AGE", "BL_YEAR", "SEX", "DEATH", "C3_CANCER"))]
dim(phe_1)     # 270534     18
phe_2 <- read.table("INFECT_2.txt", sep=" ", header=T)
dim(phe_2)     # 270534     40
phe_3 <- read.table("INFECT_3.txt", sep=" ", header=T)
dim(phe_3)     # 270534     40
phe_4 <- read.table("INFECT_4.txt", sep=" ", header=T)
dim(phe_4)     # 270534     46
phe <- cbind(phe_1, phe_2, phe_3, phe_4)
dim(phe)       # 270534    144


pheno <- merge(max_DF, covariates, by=c(1))
dim(pheno)   # 175690     32
pheno <- merge(pheno, phe, by=c(1))
dim(pheno)   # 175690    175
pheno$BL_AGE_2 <- as.numeric(as.character(pheno$BL_AGE))^2
pheno$SMOKE <- ifelse(pheno$SMOKE3=="never",0,1)




#------------------------------------------------------------------------------------
#---Create phenotype for Organ, Category, ALL (all AB), and EVERY (all AB plus nonAB)
#------------------------------------------------------------------------------------

Organ <- unique(infect$Organ_System)
Cate <- unique(infect$Category)
grp <- unique(c(Organ,Cate,"EVERY"))


for (i in 1:length(grp)){
	if (grp[i] %in% Organ){
		ep <- infect[infect$Organ_System==grp[i],"FinnGen_Endpoint"]
	} else if (grp[i] %in% Cate){
		ep <- infect[infect$Category==grp[i],"FinnGen_Endpoint"]
	} else if (grp[i]=="EVERY"){
		ep <- infect[,"FinnGen_Endpoint"]	
	} else {
		print("wrong")
	}
	ep <- colnames(pheno)[colnames(pheno) %in% ep]	
	print(i)
	
	if (length(ep)>1){
		pheno[,grp[i]] <- apply(pheno[,ep], 1, max, na.rm=T)
		pheno[,grp[i]] <- ifelse(pheno[,grp[i]]=="-Inf", NA, pheno[,grp[i]])
		pheno[,paste0(grp[i],"_AGE")] <- apply(pheno[,paste0(ep, "_AGE")], 1, min, na.rm=T)
		pheno[,paste0(grp[i],"_AGE")] <- ifelse(pheno[,paste0(grp[i],"_AGE")]=="-Inf", NA, pheno[,paste0(grp[i],"_AGE")])
	} else {
		pheno[,grp[i]] <- pheno[,ep]
		pheno[,paste0(grp[i],"_AGE")] <- pheno[,paste0(ep, "_AGE")]
	}
	pheno[,paste0(grp[i],"_ftime")] <- ifelse(pheno[,paste0(grp[i],"_AGE")]>pheno$BL_AGE, pheno[,paste0(grp[i],"_AGE")]-pheno$BL_AGE, NA)
	pheno[,paste0(grp[i],"_incid")] <- ifelse((pheno[,paste0(grp[i],"_AGE")]>pheno$BL_AGE & pheno[,grp[i]]==1) ,1,0)
	pheno[,paste0(grp[i],"_preval")] <- ifelse(pheno[,paste0(grp[i],"_AGE")]<=pheno$BL_AGE,1,0)
}


summary(pheno[,grp])
summary(pheno[,paste0(grp,"_AGE")])
summary(pheno[,paste0(grp,"_ftime")])
summary(pheno[,paste0(grp,"_incid")])
summary(pheno[,paste0(grp,"_preval")])
summary(pheno[,c("EVERY","EVERY_AGE","EVERY_ftime","EVERY_incid","EVERY_preval")])
summary(pheno[,c("ALL","ALL_AGE","ALL_ftime","ALL_incid","ALL_preval")])

# pheno <- pheno[,c("ourSid",grp,paste0(grp,"_AGE"),paste0(grp,"_ftime"),paste0(grp,"_incid"),paste0(grp,"_preval"),"C3_CANCER","C3_CANCER_AGE","BL_AGE","BL_YEAR","BL_AGE_2","SEX","SMOKE",paste0("PC",1:10),"has_MosaicCNV","Large_CNV_Clonev2","ChrX","ChrY","Autosomes","Large_ChrX","Large_ChrY","Large_Autosomes")]
write.table(pheno, "INFECTION_DATA.txt", col.names=T, row.names=F, quote=F, sep="\t")




#------------------------------------------------------------------------------------
#---Survivial analysis using Cox proportional hazards model--------------------------
#------------------------------------------------------------------------------------

# pheno <- read.table("INFECTION_DATA.txt", header=T)

summaryDF <- data.frame()
for (popn in c("ALL","male","female")){
	if (popn=="male"){
		tp <- pheno[pheno$BL_YEAR>=1989 & pheno$SEX==popn,]
		CNV_x <- c("has_MosaicCNV", "Large_CNV_Clonev2", "Autosomes", "Large_Autosomes", "ChrY", "Large_ChrY")
		grp_N <- (1:length(grp))[grp %!in% c("Female_Gynecological_Infections","Pregnancy_childbirth_and_puerperial_infections")]
	} else if (popn=="female"){
		tp <- pheno[pheno$BL_YEAR>=1989 & pheno$SEX==popn,]
		CNV_x <- c("has_MosaicCNV", "Large_CNV_Clonev2", "Autosomes", "Large_Autosomes", "ChrX", "Large_ChrX")
		grp_N <- (1:length(grp))[grp %!in% c("Male_Genitourinary_Infections")]
	} else if (popn=="ALL"){
		tp <- pheno[pheno$BL_YEAR>=1989,]
		CNV_x <- c("has_MosaicCNV", "Large_CNV_Clonev2", "Autosomes", "Large_Autosomes")
		grp_N <- (1:length(grp))[grp %!in% c("Female_Gynecological_Infections","Male_Genitourinary_Infections","Pregnancy_childbirth_and_puerperial_infections")]
	}

	for (i in grp_N){
		tmp <- tp[is.na(tp[,grp[i]])==F & is.na(tp[,paste0(grp[i],"_ftime")])==F & is.na(tp[,paste0(grp[i],"_incid")])==F & is.na(tp[,paste0(grp[i],"_preval")])==0, ]	
	
		for (j in 1:length(CNV_x)){
			tmp[,"mCA"] <- tmp[,CNV_x[j]]	
			tmp[,"ftime"] <- tmp[,paste0(grp[i],"_ftime")]	
			tmp[,"incid"] <- tmp[,paste0(grp[i],"_incid")]	
		
			for (l in c("Everyone","WithPriorCancer","WithoutPriorCancer")){
				if(l=="Everyone"){
					temp <- tmp
				} else if (l=="WithPriorCancer"){
					temp <- tmp[tmp$C3_CANCER_AGE < tmp[,paste0(grp[i],"_AGE")], ]
				} else if (l=="WithoutPriorCancer"){
					temp <- tmp[tmp$C3_CANCER_AGE >= tmp[,paste0(grp[i],"_AGE")], ]			
				} 
			
				for (k in c("SMOKE_YES","SMOKE_NO")){
					print(paste0(popn,": ",i,": ",j,": ",l,": ",k))
					if(k=="SMOKE_YES" & popn=="ALL"){
						mycoxph <- coxph(Surv(ftime, incid) ~ mCA +SEX +BL_AGE +BL_AGE_2 +SMOKE +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10, data=temp)	
						df <- as.data.frame(t(as.data.frame(summary(mycoxph)$coeff["mCA",])))				
						df$N_Cases <- length(which(temp[,"incid"]==1 & !is.na(temp[,"mCA"]) & !is.na(temp$SEX) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
						df$N_Controls <- length(which(temp[,"incid"]==0 & !is.na(temp[,"mCA"]) & !is.na(temp$SEX) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
						df$N_Cases_withVar <- length(which(temp[,"incid"]==1 & temp[,"mCA"]==1 & !is.na(temp$SEX) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
						df$N_Controls_withVar <- length(which(temp[,"incid"]==0 & temp[,"mCA"]==1 & !is.na(temp$SEX) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
					} else if (k=="SMOKE_NO" & popn=="ALL"){
						mycoxph <- coxph(Surv(ftime, incid) ~ mCA +SEX +BL_AGE +BL_AGE_2 +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10, data=temp)
						df <- as.data.frame(t(as.data.frame(summary(mycoxph)$coeff["mCA",])))
						df$N_Cases <- length(which(temp[,"incid"]==1 & !is.na(temp[,"mCA"]) & !is.na(temp$SEX) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
						df$N_Controls <- length(which(temp[,"incid"]==0 & !is.na(temp[,"mCA"]) & !is.na(temp$SEX) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
						df$N_Cases_withVar <- length(which(temp[,"incid"]==1 & temp[,"mCA"]==1 & !is.na(temp$SEX) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
						df$N_Controls_withVar <- length(which(temp[,"incid"]==0 & temp[,"mCA"]==1 & !is.na(temp$SEX) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
					} else if (k=="SMOKE_YES" & popn!="ALL"){
						if (length(which(temp[,"incid"]==1 & temp[,"mCA"]==1 & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))>0){
							mycoxph <- coxph(Surv(ftime, incid) ~ mCA +BL_AGE +BL_AGE_2 +SMOKE +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10, data=temp)	
							df <- as.data.frame(t(as.data.frame(summary(mycoxph)$coeff["mCA",])))
						} else {
							df[,1:5] <- NA
						}
						df$N_Cases <- length(which(temp[,"incid"]==1 & !is.na(temp[,"mCA"]) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
						df$N_Controls <- length(which(temp[,"incid"]==0 & !is.na(temp[,"mCA"]) & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
						df$N_Cases_withVar <- length(which(temp[,"incid"]==1 & temp[,"mCA"]==1 & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
						df$N_Controls_withVar <- length(which(temp[,"incid"]==0 & temp[,"mCA"]==1 & !is.na(temp$SMOKE) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
					} else if (k=="SMOKE_NO" & popn!="ALL"){
						if (length(which(temp[,"incid"]==1 & temp[,"mCA"]==1 & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))>0) {
							mycoxph <- coxph(Surv(ftime, incid) ~ mCA +BL_AGE +BL_AGE_2 +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10, data=temp)
							df <- as.data.frame(t(as.data.frame(summary(mycoxph)$coeff["mCA",])))
						} else {
							df[,1:5] <- NA
						}						
						df$N_Cases <- length(which(temp[,"incid"]==1 & !is.na(temp[,"mCA"]) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
						df$N_Controls <- length(which(temp[,"incid"]==0 & !is.na(temp[,"mCA"]) & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
						df$N_Cases_withVar <- length(which(temp[,"incid"]==1 & temp[,"mCA"]==1 & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
						df$N_Controls_withVar <- length(which(temp[,"incid"]==0 & temp[,"mCA"]==1 & !is.na(temp$BL_AGE) & !is.na(temp$PC1)))
					}
				
					df[,"mod"] <- k
					df$HR_25 <- exp(confint(mycoxph,"mCA"))[1]
					df$HR_975 <- exp(confint(mycoxph,"mCA"))[2]
					df[,"y"] <- grp[i]
					df[,"x"] <- CNV_x[j]
					df[,"Population"] <- popn
					df[,"Can_Popn"] <- l
					summaryDF <- rbind(summaryDF, df)
				}
			}
		}
	}
}




#------------------------------------------------------------------------------------
#---Format the results---------------------------------------------------------------------
#------------------------------------------------------------------------------------

summaryDF[,"Pheno_Grouping"] <- ifelse(summaryDF[,"y"]=="ALL",3, 
                                ifelse(summaryDF[,"y"]=="EVERY",4,
                                ifelse(summaryDF[,"y"] %in% Cate, 1, 
                            	ifelse(summaryDF[,"y"] %in% Organ, 2, NA))))
                            	
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
                          
summaryDF[,"Cancer_Popn"] <- ifelse(summaryDF[,"Can_Popn"]=="Everyone", 1,
	                         ifelse(summaryDF[,"Can_Popn"]=="WithPriorCancer", 2,
                             ifelse(summaryDF[,"Can_Popn"]=="WithoutPriorCancer", 3, NA)))


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


res <- summaryDF[,c("Pheno_Grouping", "Phenotype", "Mosaic_Variant_Grouping", "Clone_size_Filter", "Sex_Popn", "Cancer_Popn", "Model", "HR", "HR_CI", "Beta", "SE", "P", "N_Cases", "N_ctls", "N_Case_with_variant","N_Ctls_with_variant")]
res <- res[order(res$Pheno_Grouping, res$Phenotype, res$Mosaic_Variant_Grouping,res$Clone_size_Filter,res$Sex_Popn,res$Cancer_Popn,res$Model),] 
res_sig <- res[res[,"P"]<0.05 & is.na(res[,"P"])==F,]


write.csv(res, "INFECTION_Mosaic_Associations_FORMAT.csv", row.names=F, quote=F)
write.csv(res_sig[order(res_sig$P),], "INFECTION_Mosaic_Associations_FORMAT_Pval005.csv", row.names=F, quote=F)



