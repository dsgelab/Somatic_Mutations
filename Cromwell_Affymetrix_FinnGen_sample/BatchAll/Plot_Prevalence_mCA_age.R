### Plot to investigate correlation between mCA events and age (local terminal)

# gsutil cp gs://dsge-aoxing/mocha/covid/${prefix}.mca_age_sex.txt /Users/aoxliu/Documents/Project2_Finngen_mCA/Analysis_FinnGen_mCA/FinnGenBatchAll/Result


prefix <- "FinnGenBatchAll" 
setwd(paste0("/Users/aoxliu/Documents/Project2_Finngen_mCA/Analysis_FinnGen_mCA/",prefix,"/Result"))

library(ggplot2)


# Function
Right <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

std <- function(x) sd(x)/sqrt(length(x))




# Read in mca, age, and sex
m_a_s <- read.table(paste0(prefix,".mca_age_sex.txt"), header=T)

head(m_a_s)
dim(m_a_s)


# Add age bin
m_a_s$age_g <- ifelse(m_a_s$age<40,"<40", 
                ifelse(m_a_s$age>=40 & m_a_s$age <45,"40-45", ifelse(m_a_s$age>=45 & m_a_s$age <50,"45-50", ifelse(m_a_s$age>=50 & m_a_s$age <55,"50-55",
                ifelse(m_a_s$age>=55 & m_a_s$age <60,"55-60", ifelse(m_a_s$age>=60 & m_a_s$age <65,"60-65", ifelse(m_a_s$age>=65 & m_a_s$age <70,"65-70", 
                ifelse(m_a_s$age>=70 & m_a_s$age <75,"70-75", ifelse(m_a_s$age>=75 & m_a_s$age <80,"75-80", 
                ifelse(m_a_s$age>85,">85", NA))))))))))


# Types of mca                 
type_nm <- matrix(c("All mCAs", "Autosomal mCAs", "ChrX mCAs", "ChrY mCAs",
                    "any_1", "auto_2", "X_3", "Y_4",
                    "anyCF_5", "autoCF_6", "XCF_7", "YCF_8"), nrow=3, byrow=T)



# Plot
for (n in 1:ncol(type_nm)){ # 4 mCA type
	
	for (t in 1:2){  # either qc for cell fraction or not		
		m <- aggregate(m_a_s[, type_nm[(t+1), n]], list(m_a_s$age_g,m_a_s$sex), mean)
		colnames(m) <- c("age_g", "sex", "prevalence")
	
		s <- aggregate(m_a_s[, type_nm[(t+1), n]], list(m_a_s$age_g,m_a_s$sex), std)
		colnames(s) <- c("age_g", "sex", "prevalence_se")
		
		dat <- cbind(m,s)[,c(1,2,3,6)]
		dat$age_x <- as.numeric(Right(dat$age_g,2))
		dat <- dat[order(dat$sex,dat$age_x),]
		
		# Add both points and line
		p <- ggplot(dat, aes(x = age_x, y = prevalence,colour = sex)) +geom_point()
		p + stat_smooth(method = "loess", formula = y ~ x,se=F, size = 1)+theme_bw()+ scale_x_continuous(breaks=dat$age_x, label=dat$age_g) +labs(y="Prevalence")+labs(x="Age") +
			labs(title="Prevalence of mCA on sex chromosome with cell fraction >=10%", face="bold")
			# labs(title="Prevalence of mCA on sex chromosome", face="bold")
			# labs(title="Prevalence of autosomal mCA", face="bold")
			# labs(title="Prevalence of autosomal mCA with cell fraction >=10%", face="bold")	
	}
}




