### what is causing females to cluster as two different Y nonPAR LRR groups?


#----------------------------------------
## Any correlation between y_nonpar_lrr_median and variables in "stats_tsv" and "affy_tsv" for females


prefix="FinnGenBatchAll"
wf_id="3bb7e444-d05c-49d7-9e48-d5e012240aa9"
sample_tsv="gs://from-fg-datateam/check/BatchAll_sort_dup.sample.tsv"
stats_tsv="gs://dsge-cromwell/mocha/${wf_id}/call-mocha_stats_tsv/FinnGenBatchAll.stats.tsv"
affy_tsv="gs://dsge-cromwell/mocha/${wf_id}/call-affy_tsv/FinnGenBatchAll.affy.tsv"
prefix="FinnGenBatchAll"

cd /Users/aoxliu/Documents/Project2_Finngen_mCA/Analysis_FinnGen_mCA/FinnGenBatchAll/Result/nonPAR_LRR


# merge multiple files 
gsutil cat ${affy_tsv}|wc -l     # 201,459
gsutil cat ${stats_tsv}|wc -l    # 201,459
gsutil cat ${sample_tsv}|wc -l   # 201,459
paste -d '\t' <(gsutil cat ${affy_tsv}|cut -f 1|awk 'NR>1') <(gsutil cat ${stats_tsv}|cut -f 1|awk 'NR>1')|sed 's/-1//g'|awk '$1!=$2'
paste -d '\t' <(gsutil cat ${affy_tsv}|cut -f 1|awk 'NR>1') <(gsutil cat ${sample_tsv}|cut -f 1|awk 'NR>1')|sed 's/-1//g'|awk '$1!=$2'

paste -d '\t' <(gsutil cat ${sample_tsv}|cut -f -2) <(gsutil cat ${stats_tsv}|cut -f 2-) <(gsutil cat ${affy_tsv}|cut -f 4-) > nonPAR_LRR.txt
awk 'NR==1{print NF}' nonPAR_LRR.txt   # 40
wc -l nonPAR_LRR.txt             # 20,1459



# Plot by R
echo "
setwd("/Users/aoxliu/Documents/Project2_Finngen_mCA/Analysis_FinnGen_mCA/FinnGenBatchAll/Result/nonPAR_LRR")
nonPAR <- read.table("nonPAR_LRR.txt", sep="\t", header=T)

var <- colnames(nonPAR)[c(4:15,17:21,22:38)]
for (k in 1:length(var)){
	print(paste0(k,": ",var[k]))
	tiff(paste0("Batch_b01_09_",var[k],".tiff"))
	par(mfrow=c(3,3))

	for (i in 1:9){
	# for (i in 1:30){  # FinnGen SNP array V1
 	# for (i in 1:30){  # FinnGen SNP array V2 
	# for (i in 1:51){  # FinnGen SNP array V1 & V2
		nonPAR_i <- nonPAR[as.numeric(substr(nonPAR$batch_id,2,3))==i & nonPAR$computed_gender=="F", ]
		nrow(nonPAR_i)
		plot(as.numeric(as.character(nonPAR_i[, var[k]])), as.numeric(as.character(nonPAR_i[, "y_nonpar_lrr_median"])), xlab=var[k], ylab="y_nonpar_lrr_median",col="red", main=paste0("Batch b",i), ylim=c(-3, 0.5))
	}
	dev.off()
} 

# t-test to check whether lower and higher groups of y_nonpar_lrr_median have difference in lrr_auto
nonPAR_6 <- nonPAR[as.numeric(substr(nonPAR$batch_id,2,3))==6 & nonPAR$computed_gender=="F", ]
t.test(nonPAR_6[nonPAR_6$y_nonpar_lrr_median> -1.8, "lrr_auto"], nonPAR_6[nonPAR_6$y_nonpar_lrr_median< -2.2, "lrr_auto"], alternative = "two.sided", var.equal = FALSE)

" > nonPAR_LRR_Plot.R

R nonPAR_LRR_Plot.R






#----------------------------------------
## chrX plots from two females from two groups using mocha_plot.R

# local
# select one individual from each cluster
less /Users/aoxliu/Documents/Project2_Finngen_mCA/Analysis_FinnGen_mCA/FinnGenBatchAll/Result/nonPAR_LRR/nonPAR_LRR.txt|head -n 1 |tr '\t' '\n'|cat -n|grep y_nonpar_lrr_median # 16
less /Users/aoxliu/Documents/Project2_Finngen_mCA/Analysis_FinnGen_mCA/FinnGenBatchAll/Result/nonPAR_LRR/nonPAR_LRR.txt|awk '$2=="b01" && $16> -1.8'|cut -f1
less /Users/aoxliu/Documents/Project2_Finngen_mCA/Analysis_FinnGen_mCA/FinnGenBatchAll/Result/nonPAR_LRR/nonPAR_LRR.txt|awk '$2=="b01" && $16< -2.2'|cut -f1 


# GCP aoxing-mocha
prefix="FinnGenBatchAll"
wf_id="3bb7e444-d05c-49d7-9e48-d5e012240aa9"
stats_tsv="/home/aoxliu/mCA/input/dsge-cromwell/mocha/${wf_id}/call-mocha_stats_tsv/FinnGenBatchAll.stats.tsv"
vcf="/home/aoxliu/mCA/input/dsge-cromwell/mocha/${wf_id}/call-vcf_mocha/shard-0/FinnGenBatchAll.b01.mocha.bcf"
cyto="/home/aoxliu/mCA/input/dsge-aoxing/mocha/GRCh38/mocha.GRCh38/cytoBand.hg38.txt.gz"



# mount the data from Google bucket to VM instance
cd  /home/aoxliu/mCA/input

gcsfuse --implicit-dirs  from-fg-datateam  from-fg-datateam 
gcsfuse --implicit-dirs  dsge-aoxing  dsge-aoxing


cd /home/aoxliu/mCA/output/BatchFirst30/plot

./mocha_plot.R \
  --mocha \
  --stats ${stats_tsv} \
  --vcf ${vcf}\
  --png X_b01_C2_FG22XDY8K2.png \
  --samples FG22XDY8K2 \
  --regions chrX:1-115077367 \
  --cytoband ${cyto}

gsutil cp X_b01_C2_FG22XDY8K2.png gs://dsge-cromwell/mocha/3bb7e444-d05c-49d7-9e48-d5e012240aa9/call-vcf_mocha/shard-0/



./mocha_plot.R \
  --mocha \
  --stats ${stats_tsv} \
  --vcf ${vcf}\
  --png X_b01_C1_FG236F7GCC.png \
  --samples FG236F7GCC \
  --regions chrX:1-115077367 \
  --cytoband ${cyto}

gsutil cp X_b01_C1_FG236F7GCC.png gs://dsge-cromwell/mocha/3bb7e444-d05c-49d7-9e48-d5e012240aa9/call-vcf_mocha/shard-0/


# unmount the data from Google bucket to VM instance
fusermount -u  /home/aoxliu/mCA/input/from-fg-datateam 
fusermount -u  /home/aoxliu/mCA/input/dsge-aoxing


