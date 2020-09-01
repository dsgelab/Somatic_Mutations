### COVID-19 somatic mutation project for https://yale.app.box.com/s/crwbkzqhxzyhjxzrnzzon6ftyruwt2cj
# 8 mCA covariates * 3 models * 2 COVID-19 phenotypes (Severe COVID+, COVID-19 test reporting)



#-----------------------------------------------------------
## Step 0: set up the working environments

# login VM instance
gcloud beta compute --project "finngen-refinery-dsgelab" ssh --zone "europe-west1-b" "aoxing-mocha"     # pathword: blackfriday

# mount the data from Google bucket to VM instance
cd  /home/aoxliu/mCA/input

gcsfuse --implicit-dirs  from-fg-datateam  from-fg-datateam 
gcsfuse --implicit-dirs  dsge-aoxing  dsge-aoxing

# fusermount -u  /home/aoxliu/mCA/input/from-fg-datateam 
# fusermount -u  /home/aoxliu/mCA/input/dsge-aoxing


cd /home/aoxliu/mCA/output/BatchFirst30     # store all outputs of BatchFirst30 here

wf_id="3bb7e444-d05c-49d7-9e48-d5e012240aa9"
prefix="FinnGenBatchAll" 
sample_tsv="gs://from-fg-datateam/check/BatchAll_sort_dup.sample.tsv"





#----------------------------------------------------------------------------------------------------------------------
## Step 1: run MoChA and filtering to high-quality calls

# Filtering the calls 
gsutil cat gs://dsge-cromwell/mocha/${wf_id}/call-mocha_calls_tsv/${prefix}.calls.tsv|wc -l   # 183,580


awk -F "\t" 'NR==FNR && FNR==1 {for (i=1; i<=NF; i++) f[$i] = i} \
  NR==FNR && FNR>1 {sample_id=$(f["sample_id"]); call_rate=$(f["call_rate"]); baf_auto=$(f["baf_auto"])} \
  NR==FNR && FNR>1 && (call_rate<.97 || baf_auto>.03) {xcl[sample_id]++} \
  NR>FNR && FNR==1 {for (i=1; i<=NF; i++) g[$i] = i; print} \
  NR>FNR && FNR>1 {sample_id=$(g["sample_id"]); len=$(g["length"]); p_arm=$(g["p_arm"]); q_arm=$(g["q_arm"]); \
  bdev=$(g["bdev"]); rel_cov=$(g["rel_cov"]); lod_baf_phase=$(g["lod_baf_phase"]); type=$(g["type"]); \
  if (lod_baf_phase=="nan") lod_baf_phase=0} \
  NR>FNR && FNR>1 && !(sample_id in xcl) && rel_cov>0.5 && type!~"^CNP" && \
  (len>5e6 + 5e6 * (p_arm!="N" && q_arm!="N") || len>5e5 && (bdev<1/10 && rel_cov<2.5) && lod_baf_phase>10 || rel_cov<2.1 && lod_baf_phase>10 )' \
   <(gsutil cat gs://dsge-cromwell/mocha/${wf_id}/call-mocha_stats_tsv/${prefix}.stats.tsv) \
   <(gsutil cat gs://dsge-cromwell/mocha/${wf_id}/call-mocha_calls_tsv/${prefix}.calls.tsv) > ${prefix}.filtered.calls.tsv


wc -l ${prefix}.filtered.calls.tsv   # 38,650
 
gsutil cp ${prefix}.filtered.calls.tsv  gs://dsge-cromwell/mocha/${wf_id}/call-mocha_calls_tsv/
gsutil cat gs://dsge-cromwell/mocha/${wf_id}/call-mocha_calls_tsv/${prefix}.filtered.calls.tsv|wc -l    # 38650




# call >1Mbp on chrX in a male (look at the second column for the gender in .calls.tsv) then it is most likely a mLOY event


# Explore the chrY
less ${prefix}.filtered.calls.tsv|grep chrX|awk '$2=="M"'|wc -l
# 16298

less ${prefix}.filtered.calls.tsv|grep chrX|awk '$2=="M" && $6>=1000000'|wc -l
# 16276

less ${prefix}.filtered.calls.tsv|grep chrX|awk '$2=="M" && $6<1000000'|wc -l
# 16276

less ${prefix}.filtered.calls.tsv|grep chrX|awk '$2=="M" && $6<1000000'|wc -l
# 22





#----------------------------------------------------------------------------------------------------------------------
## Step 2: Separate out mosaic variant calls into 8 categories:

gsutil cat ${sample_tsv}|awk 'NR>1{print $1}'|wc -l   # 201,458

# 1) Any mosaic variant (0/1)
cat <(comm -12 <(awk 'NR>1{print $1}' ${prefix}.filtered.calls.tsv|sort|uniq) <(gsutil cat ${sample_tsv}|awk 'NR>1{print $1}'|sort)|awk '{print $1,1}') \
    <(comm -13 <(awk 'NR>1{print $1}' ${prefix}.filtered.calls.tsv|sort|uniq) <(gsutil cat ${sample_tsv}|awk 'NR>1{print $1}'|sort)|awk '{print $1,0}') > ${prefix}.mca_any_1.txt


# 2) Any autosomal mosaic variant (0/1)
cat <(comm -12 <(grep -v "chrX\|chrY" ${prefix}.filtered.calls.tsv|awk 'NR>1{print $1}'|sort|uniq) <(gsutil cat ${sample_tsv}|awk 'NR>1{print $1}'|sort)|awk '{print $1,1}') \
    <(comm -13 <(grep -v "chrX\|chrY" ${prefix}.filtered.calls.tsv|awk 'NR>1{print $1}'|sort|uniq) <(gsutil cat ${sample_tsv}|awk 'NR>1{print $1}'|sort)|awk '{print $1,0}') > ${prefix}.mca_auto_2.txt


# 3) Any ChrX mosaic variant (0/1)
cat <(comm -12 <(awk '$3=="chrX"{print $1}' ${prefix}.filtered.calls.tsv|sort|uniq) <(gsutil cat ${sample_tsv}|awk 'NR>1{print $1}'|sort)|awk '{print $1,1}') \
    <(comm -13 <(awk '$3=="chrX"{print $1}' ${prefix}.filtered.calls.tsv|sort|uniq) <(gsutil cat ${sample_tsv}|awk 'NR>1{print $1}'|sort)|awk '{print $1,0}') > ${prefix}.mca_X_3.txt


# 4) Any ChrY mosaic variant (0/1)
cat <(comm -12 <(awk '$3=="chrY"{print $1}' ${prefix}.filtered.calls.tsv|sort|uniq) <(gsutil cat ${sample_tsv}|awk 'NR>1{print $1}'|sort)|awk '{print $1,1}') \
    <(comm -13 <(awk '$3=="chrY"{print $1}' ${prefix}.filtered.calls.tsv|sort|uniq) <(gsutil cat ${sample_tsv}|awk 'NR>1{print $1}'|sort)|awk '{print $1,0}') > ${prefix}.mca_Y_4.txt


# 5) Any mosaic variant (0/1) with CELL_FRAC > 10%
cat <(comm -12 <(awk 'NR>1 && $NF>0.1{print $1}' ${prefix}.filtered.calls.tsv|sort|uniq) <(gsutil cat ${sample_tsv}|awk 'NR>1{print $1}'|sort)|awk '{print $1,1}') \
    <(comm -13 <(awk 'NR>1 && $NF>0.1{print $1}' ${prefix}.filtered.calls.tsv|sort|uniq) <(gsutil cat ${sample_tsv}|awk 'NR>1{print $1}'|sort)|awk '{print $1,0}') > ${prefix}.mca_anyCF_5.txt


# 6) Any autosomal mosaic variant (0/1) with CELL_FRAC > 10%
cat <(comm -12 <(grep -v "chrX\|chrY" ${prefix}.filtered.calls.tsv|awk 'NR>1 && $NF>0.1{print $1}'|sort|uniq) <(gsutil cat ${sample_tsv}|awk 'NR>1{print $1}'|sort)|awk '{print $1,1}') \
    <(comm -13 <(grep -v "chrX\|chrY" ${prefix}.filtered.calls.tsv|awk 'NR>1 && $NF>0.1{print $1}'|sort|uniq) <(gsutil cat ${sample_tsv}|awk 'NR>1{print $1}'|sort)|awk '{print $1,0}') > ${prefix}.mca_autoCF_6.txt


# 7) Any ChrX mosaic variant (0/1) with CELL_FRAC > 10%
cat <(comm -12 <(awk '$3=="chrX" && $NF>0.1{print $1}' ${prefix}.filtered.calls.tsv|sort|uniq) <(gsutil cat ${sample_tsv}|awk 'NR>1{print $1}'|sort)|awk '{print $1,1}') \
    <(comm -13 <(awk '$3=="chrX" && $NF>0.1{print $1}' ${prefix}.filtered.calls.tsv|sort|uniq) <(gsutil cat ${sample_tsv}|awk 'NR>1{print $1}'|sort)|awk '{print $1,0}') > ${prefix}.mca_XCF_7.txt


# 8) Any ChrY mosaic variant (0/1)  with CELL_FRAC > 10%
cat <(comm -12 <(awk '$3=="chrY" && $NF>0.1{print $1}' ${prefix}.filtered.calls.tsv|sort|uniq) <(gsutil cat ${sample_tsv}|awk 'NR>1{print $1}'|sort)|awk '{print $1,1}') \
    <(comm -13 <(awk '$3=="chrY" && $NF>0.1{print $1}' ${prefix}.filtered.calls.tsv|sort|uniq) <(gsutil cat ${sample_tsv}|awk 'NR>1{print $1}'|sort)|awk '{print $1,0}') > ${prefix}.mca_YCF_8.txt



# Check for all 8 categories
for i in any_1  auto_2  X_3  Y_4  anyCF_5  autoCF_6  XCF_7  YCF_8
do
echo ${i}
wc -l ${prefix}.mca_${i}.txt   # 122,252
awk '{print $2}' ${prefix}.mca_${i}.txt|sort|uniq -c 
done
# any_1:      30988
# auto_2:     7177
# X_3:        25509
# Y_4:        0
# anyCF_5:    13759
# autoCF_6:   3465
# XCF_7:      10675 
# YCF_8:      0


paste -d ' ' <(sort -k1,1 ${prefix}.mca_any_1.txt) <(sort -k1,1 ${prefix}.mca_auto_2.txt) <(sort -k1,1 ${prefix}.mca_X_3.txt) <(sort -k1,1 ${prefix}.mca_Y_4.txt) \
             <(sort -k1,1 ${prefix}.mca_anyCF_5.txt) <(sort -k1,1 ${prefix}.mca_autoCF_6.txt) <(sort -k1,1 ${prefix}.mca_XCF_7.txt) <(sort -k1,1 ${prefix}.mca_YCF_8.txt) > tmp_${prefix}.mca_1_8.txt

awk '$1==$3 && $1=$5 && $1==$7 && $1=$9 && $1==$11 && $1=$13 && $1==$15'  tmp_${prefix}.mca_1_8.txt|wc -l   # 122,252
cat <(echo "id any_1  auto_2  X_3  Y_4  anyCF_5  autoCF_6  XCF_7  YCF_8") \
    <(awk '{print $1, $2, $4, $6, $8, $10, $12, $14, $16}' tmp_${prefix}.mca_1_8.txt) > ${prefix}.mca_1_8.txt 
rm tmp_${prefix}.mca_1_8.txt





#----------------------------------------------------------------------------------------------------------------------
# Step 2S: Add covariates


# extract minimum phenotypes (BL_YEAR|Year of DNA sample collection, BL_AGE|Age at DNA sample collection) and other covariate !!!!
gsutil cat gs://finngen-production-library-red/finngen_R6/phenotype_2.0/data/finngen_R6_v2_minimum.txt|awk '{print $1, $2,$3,$4,$5,$6,$7,$8,$9,$12,$13,$15}' > /home/aoxliu/mCA/input/dsge-aoxing/mocha/covid/finngen_R6_v2_min.txt

### /finngen/library-red/finngen_R6/phenotype_2.0/data/finngen_R6_cov_pheno_1.0.txt.gz  



# only in mCA calling but not in FinnGen R6 (could be due to the removal of related inliers)
comm -13  <(gsutil cat gs://dsge-aoxing/mocha/covid/finngen_R6_v2_min.txt|awk '{print $1}'|sort -k1,1) \
          <(awk '{print $1}' ${prefix}.mca_1_8.txt|sort -k1,1) |grep -v id|grep -v '-'| wc -l     # 5,414


# both in mCA calling and FinnGen R6
join -1 1 -2 1 <(sort -k1,1 ${prefix}.mca_1_8.txt) \
               <(gsutil cat gs://dsge-aoxing/mocha/covid/finngen_R6_v2_min.txt|sort -k1,1) > ${prefix}.mca_1_8.phe_min.txt
               
wc -l ${prefix}.mca_1_8.phe_min.txt   # 195,908


cat <(echo "any_1  auto_2  X_3  Y_4  anyCF_5  autoCF_6  XCF_7  YCF_8  age  sex")  \ 
    <(awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$11,$12}' ${prefix}.mca_1_8.phe_min.txt) > ${prefix}.mca_age_sex.txt

wc -l ${prefix}.mca_age_sex.txt   # 195,909 

gsutil cp ${prefix}.mca*.txt gs://dsge-aoxing/mocha/covid/
# gsutil cp ${prefix}.mca_age_sex.txt gs://dsge-aoxing/mocha/covid/











#----------------------------------------------------------------------------------------------------------------------
## Step 3: Create 2 COVID-19 phenotypes

# covid_file="gs://finngen-production-library-red/finngen_R5/coron/2020-06-11/data/finngen_R5_v7_infectious_disease_register_corona_2020-06-11.txt"
covid_file="gs://finngen-production-library-red/finngen_R6/corona/2020-07-31/data/finngen_R6_v8_infectious_disease_register_corona_2020-07-31.txt"

gsutil cat ${covid_file}|wc -l    # 358 (305 for finngen_R5)






# In the UK Biobank, Severe COVID+ is someone who has tested positive for COVID at least once and has been hospitalized.


# Basic quality control, remove related individuals from the dataset (1st or 2nd degree relatives)



# and remove individuals with any prevalent hematologic cancer history at time of blood draw for genotyping (ICD C81-96, D45-47), 



# and remove individuals who died prior to 3/1/2020.


join -1 1 -2 1 <(gsutil cat ${covid_file}|awk '{print $1}'|sort) <(gsutil cat gs://dsge-aoxing/mocha/covid/${prefix}.mca_any_1.txt|sort -k1,1)  



#-----------------------------------------------------------
Step 4: Run the following logistic regression analyses for each of the 8 variant groupings and 2 COVID-19 phenotypes above: 
A.Unadjusted: a.COVID_Pheno ~ has_MosaicVar 
B.Sparsely adjusted: a.COVID_Pheno ~ has_MosaicVar + age +age2+sex + ever_smoking_status+ PC1-10 of genomic ancestry
C.Fully adjusted: a.COVID_Pheno ~ has_MosaicVar + age +age2+ sex + ever_smoking_status + scale(Townsend) + PC1-10 of genomic ancestry + scale(BMI)+ Asthma + COPD + CAD + HTN + T2D+ anyCancer



#-----------------------------------------------------------
## Batch 72 and 73 are still running






