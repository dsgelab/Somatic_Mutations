# Association between CHRNA5 and mCA events (Using "aoxing-mocha" VM)

# mount the data from Google bucket to VM instance
# cd  /home/aoxliu/mCA/input
# gcsfuse --implicit-dirs  from-fg-datateam  from-fg-datateam 


setwd("/home/aoxliu/mCA/input/dsge-aoxing/mocha/SPDL1/")

library(data.table)



#------------------------------------------------------------------------------------
#---Loading data---------------------------------------------------------------------
#------------------------------------------------------------------------------------

# mCA events
mca <- read.table("/home/aoxliu/mCA/input/dsge-aoxing/mocha/output/mCA_01_Counts.afterQC.withoutCF_1.txt", header=T)
mca <- mca[, !(names(mca) %in% c("SEX","BL_AGE"))]
colnames(mca) <- c("ourSid","AnyChr","Large_AnyChr","ChrX","ChrY","Autosomes","Large_ChrX","Large_ChrY","Large_Autosomes")
dim(mca)    # 189,557     11
length(unique(mca$ourSid))   # 189,557


# SPDL1
# system("awk 'NF==2' finngen_R6_b01_b51.rs72740964.vcf > finngen_R6_b01_b51.rs72740964_format.vcf")
# snp <- read.table("finngen_R6_b01_b51.rs72740964_format.vcf", header=F)
snp <- read.table("/home/aoxliu/mCA/output/SPDL1/finngen_R6_b01_b51.rs72740964_format.vcf", header=F)
colnames(snp) <- c("id","geno")
dim(snp)    # 201,323      2
length(unique(snp$id))  # 197,242
snp <- unique(snp)
dim(snp)    # 197,242      2


# Co-variates
covariates <- read.table("/home/aoxliu/mCA/input/dsge-aoxing/mocha/output/finngen_R6_pheno_cov.txt", header=T, stringsAsFactors=F)
dim(covariates)   # 260,405     37


# Combine
pheno <- merge(mca, snp, by=c(1))
dim(pheno)  # 189,556     10

pheno <- merge(pheno, covariates, by=c(1))
dim(pheno)  # 189,556     46

length(unique(pheno$ourSid))   # 189,556
colnames(pheno)

 

#------------------------------------------------------------------------------------
#---Association analysis-------------------------------------------------------------
#------------------------------------------------------------------------------------

pheno$geno <- as.character(pheno$geno)
pheno$snp <- ifelse(pheno$geno=="0|0",0, ifelse(pheno$geno=="1|1",2,1))
pheno$BL_AGE_2 <- as.numeric(as.character(pheno$BL_AGE))^2


# All and without PCA--------------------------------
summaryDF <- data.frame()
temp <- pheno
temp[ ,"SEX"] <- ifelse(temp[ ,"SEX"]=="male",1,ifelse(temp[ ,"SEX"]=="female",2,NA))

for (i in c("Autosomes","Large_Autosomes","AnyChr","Large_AnyChr")){
	for (k in c("ALL")){
		# for (j in 1:4){
		for (j in c(1,3)){
			print(paste(i,k,j,sep=":"))
			
			if (j==1){
				df <- summary(glm(temp[ ,i]~snp+BL_AGE+BL_AGE_2+SEX+SMOKE2, family="binomial", data=temp))$coeff[,]
			} else if (j==2) {
				df <- summary(glm(temp[ ,i]~snp+BL_AGE+BL_AGE_2+SEX+SMOKE2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, family="binomial", data=temp))$coeff[,]			
			} else if (j==3){
				df <- summary(glm(temp[ ,i]~snp+BL_AGE+BL_AGE_2+SEX, family="binomial", data=temp))$coeff[,]			
			} else {
				df <- summary(glm(temp[ ,i]~snp+BL_AGE+BL_AGE_2+SEX+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, family="binomial", data=temp))$coeff[,]			
			}
			
			mod <- matrix(NA, nrow=nrow(df), ncol=8)
			colnames(mod) <- c("Model", "Population", "y", "Variable", "N_Cases", "N_Controls", "N_Cases_withVar", "N_Controls_withVar")
			
			if (j==1){			
				mod[,"Model"] <- paste0(i, "=rs72740964+AGE+AGE_2+SEX+SMOKE2")
				mod[,"Variable"] <- c("intercept","snp","BL_AGE","BL_AGE_2","SEX","SMOKE2")	
			} else if (j==2){
				mod[,"Model"] <- paste0(i, "=rs72740964+AGE+AGE_2+SEX+SMOKE2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")
				mod[,"Variable"] <- c("intercept","snp","BL_AGE","BL_AGE_2","SEX","SMOKE2","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")	
			} else if (j==3){
				mod[,"Model"] <- paste0(i, "=rs72740964+AGE+AGE_2+SEX")
				mod[,"Variable"] <- c("intercept","snp","BL_AGE","BL_AGE_2","SEX")	
			} else {
				mod[,"Model"] <- paste0(i, "=rs72740964+AGE+AGE_2+SEX+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")
				mod[,"Variable"] <- c("intercept","snp","BL_AGE","BL_AGE_2","SEX","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")	
			}		
							
			mod[,"Population"] <- k
			mod[,"y"] <- i	
			
			if (j<=2){
				mod[,"N_Cases"] <- length(which(temp[,i]==1 & !is.na(temp$SMOKE2)))
				mod[,"N_Controls"] <- length(which(temp[,i]==0 & !is.na(temp$SMOKE2)))
				mod[,"N_Cases_withVar"] <- length(which(temp[,i]==1 & temp[,"snp"]==1 & !is.na(temp$SMOKE2)))
				mod[,"N_Controls_withVar"] <- length(which(temp[,i]==0 & temp[,"snp"]==1 & !is.na(temp$SMOKE2)))	
			} else {
				mod[,"N_Cases"] <- length(which(temp[,i]==1))
				mod[,"N_Controls"] <- length(which(temp[,i]==0))
				mod[,"N_Cases_withVar"] <- length(which(temp[,i]==1 & temp[,"snp"]==1))
				mod[,"N_Controls_withVar"] <- length(which(temp[,i]==0 & temp[,"snp"]==1))			
			}
			
			res <- cbind(df,mod)
			summaryDF <- rbind(summaryDF, res)
		}
	}
}




# Male(ChrY) and without PCA--------------------------------
for (i in c("ChrY","Large_ChrY","Autosomes","Large_Autosomes","AnyChr","Large_AnyChr")){
	for (k in c("male")){
		# for (j in 1:4){
		for (j in c(1,3)){
			print(paste(i,k,j,sep=":"))
			temp <- pheno[pheno$SEX==k,]
			
			if (j==1){
				df <- summary(glm(temp[ ,i]~snp+BL_AGE+BL_AGE_2+SMOKE2, family="binomial", data=temp))$coeff[,]
			} else if (j==2) {
				df <- summary(glm(temp[ ,i]~snp+BL_AGE+BL_AGE_2+SMOKE2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, family="binomial", data=temp))$coeff[,]			
			} else if (j==3){
				df <- summary(glm(temp[ ,i]~snp+BL_AGE+BL_AGE_2, family="binomial", data=temp))$coeff[,]			
			} else {
				df <- summary(glm(temp[ ,i]~snp+BL_AGE+BL_AGE_2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, family="binomial", data=temp))$coeff[,]			
			}
			
			mod <- matrix(NA, nrow=nrow(df), ncol=8)
			colnames(mod) <- c("Model", "Population", "y", "Variable", "N_Cases", "N_Controls", "N_Cases_withVar", "N_Controls_withVar")
			
			if (j==1){			
				mod[,"Model"] <- paste0(i, "=rs72740964+AGE+AGE_2+SMOKE2")
				mod[,"Variable"] <- c("intercept","snp","BL_AGE","BL_AGE_2","SMOKE2")	
			} else if (j==2){
				mod[,"Model"] <- paste0(i, "=rs72740964+AGE+AGE_2+SMOKE2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")
				mod[,"Variable"] <- c("intercept","snp","BL_AGE","BL_AGE_2","SMOKE2","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")	
			} else if (j==3){
				mod[,"Model"] <- paste0(i, "=rs72740964+AGE+AGE_2")
				mod[,"Variable"] <- c("intercept","snp","BL_AGE","BL_AGE_2")	
			} else {
				mod[,"Model"] <- paste0(i, "=rs72740964+AGE+AGE_2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")
				mod[,"Variable"] <- c("intercept","snp","BL_AGE","BL_AGE_2","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")	
			}		

			mod[,"Population"] <- k
			mod[,"y"] <- i
			
			if (j<=2){
				mod[,"N_Cases"] <- length(which(temp[,i]==1 & !is.na(temp$SMOKE2)))
				mod[,"N_Controls"] <- length(which(temp[,i]==0 & !is.na(temp$SMOKE2)))
				mod[,"N_Cases_withVar"] <- length(which(temp[,i]==1 & temp[,"snp"]==1 & !is.na(temp$SMOKE2)))
				mod[,"N_Controls_withVar"] <- length(which(temp[,i]==0 & temp[,"snp"]==1 & !is.na(temp$SMOKE2)))	
			} else {
				mod[,"N_Cases"] <- length(which(temp[,i]==1))
				mod[,"N_Controls"] <- length(which(temp[,i]==0))
				mod[,"N_Cases_withVar"] <- length(which(temp[,i]==1 & temp[,"snp"]==1))
				mod[,"N_Controls_withVar"] <- length(which(temp[,i]==0 & temp[,"snp"]==1))			
			}

			res <- cbind(df,mod)
			summaryDF <- rbind(summaryDF, res)
		}
	}
}




# Female(ChrX) and without PCA--------------------------------
for (i in c("ChrX","Large_ChrX","Autosomes","Large_Autosomes","AnyChr","Large_AnyChr")){
	for (k in c("female")){
		# for (j in 1:4){
		for (j in c(1,3)){
			print(paste(i,k,j,sep=":"))
			temp <- pheno[pheno$SEX==k,]
			
			if (j==1){
				df <- summary(glm(temp[ ,i]~snp+BL_AGE+BL_AGE_2+SMOKE2, family="binomial", data=temp))$coeff[,]
			} else if (j==2) {
				df <- summary(glm(temp[ ,i]~snp+BL_AGE+BL_AGE_2+SMOKE2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, family="binomial", data=temp))$coeff[,]			
			} else if (j==3){
				df <- summary(glm(temp[ ,i]~snp+BL_AGE+BL_AGE_2, family="binomial", data=temp))$coeff[,]			
			} else {
				df <- summary(glm(temp[ ,i]~snp+BL_AGE+BL_AGE_2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, family="binomial", data=temp))$coeff[,]			
			}
			
			mod <- matrix(NA, nrow=nrow(df), ncol=8)
			colnames(mod) <- c("Model", "Population", "y", "Variable", "N_Cases", "N_Controls", "N_Cases_withVar", "N_Controls_withVar")
			
			if (j==1){			
				mod[,"Model"] <- paste0(i, "=rs72740964+AGE+AGE_2+SMOKE2")
				mod[,"Variable"] <- c("intercept","snp","BL_AGE","BL_AGE_2","SMOKE2")	
			} else if (j==2){
				mod[,"Model"] <- paste0(i, "=rs72740964+AGE+AGE_2+SMOKE2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")
				mod[,"Variable"] <- c("intercept","snp","BL_AGE","BL_AGE_2","SMOKE2","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")	
			} else if (j==3){
				mod[,"Model"] <- paste0(i, "=rs72740964+AGE+AGE_2")
				mod[,"Variable"] <- c("intercept","snp","BL_AGE","BL_AGE_2")	
			} else {
				mod[,"Model"] <- paste0(i, "=rs72740964+AGE+AGE_2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")
				mod[,"Variable"] <- c("intercept","snp","BL_AGE","BL_AGE_2","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")	
			}		

			mod[,"Population"] <- k
			mod[,"y"] <- i
			
			if (j<=2){
				mod[,"N_Cases"] <- length(which(temp[,i]==1 & !is.na(temp$SMOKE2)))
				mod[,"N_Controls"] <- length(which(temp[,i]==0 & !is.na(temp$SMOKE2)))
				mod[,"N_Cases_withVar"] <- length(which(temp[,i]==1 & temp[,"snp"]==1 & !is.na(temp$SMOKE2)))
				mod[,"N_Controls_withVar"] <- length(which(temp[,i]==0 & temp[,"snp"]==1 & !is.na(temp$SMOKE2)))	
			} else {
				mod[,"N_Cases"] <- length(which(temp[,i]==1))
				mod[,"N_Controls"] <- length(which(temp[,i]==0))
				mod[,"N_Cases_withVar"] <- length(which(temp[,i]==1 & temp[,"snp"]==1))
				mod[,"N_Controls_withVar"] <- length(which(temp[,i]==0 & temp[,"snp"]==1))			
			}

			res <- cbind(df,mod)
			summaryDF <- rbind(summaryDF, res)
		}
	}
}



summaryDF <- summaryDF[,c("Model","Population","y","Variable","Estimate","Std. Error","Pr(>|z|)","N_Cases","N_Controls","N_Cases_withVar","N_Controls_withVar")]
colnames(summaryDF) <- c("Model","Population","y","Variable","Estimate","SE(Estimate)","P-val","N_Cases","N_Controls","N_Cases_withVar","N_Controls_withVar")
write.csv(summaryDF, "summaryDF_rs72740964_noPC.csv", row.names=F)
write.csv(summaryDF[summaryDF$Variable=="snp", ], "summaryDF_rs72740964_snp_noPC.csv", row.names=F)


