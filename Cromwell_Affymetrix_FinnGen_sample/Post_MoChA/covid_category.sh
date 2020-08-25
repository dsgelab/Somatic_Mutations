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

wf_id="54024b98-f66d-4eab-b6db-fecd717b5896"
prefix="BatchFirst30" 
sample_tsv="gs://from-fg-datateam/check/BatchFirst30_sort_dup.sample.tsv"



#-----------------------------------------------------------
## Step 1: run MoChA and filtering to high-quality calls

# Filtering the calls 
awk -F "\t" 'NR==FNR && FNR==1 {for (i=1; i<=NF; i++) f[$i] = i} \
  NR==FNR && FNR>1 {sample_id=$(f["sample_id"]); call_rate=$(f["call_rate"]); baf_auto=$(f["baf_auto"])} \
  NR==FNR && FNR>1 && (call_rate<.97 || baf_auto>.03) {xcl[sample_id]++} \
  NR>FNR && FNR==1 {for (i=1; i<=NF; i++) g[$i] = i; print} \
  NR>FNR && FNR>1 {sample_id=$(g["sample_id"]); len=$(g["length"]); p_arm=$(g["p_arm"]); q_arm=$(g["q_arm"]); \
  bdev=$(g["bdev"]); rel_cov=$(g["rel_cov"]); lod_baf_phase=$(g["lod_baf_phase"]); type=$(g["type"]); \
  if (lod_baf_phase=="nan") lod_baf_phase=0} \
  NR>FNR && FNR>1 && !(sample_id in xcl) && rel_cov>0.5 && type!~"^CNP" && \
  (len>5e6 + 5e6 * (p_arm!="N" && q_arm!="N") || len>5e5 && (bdev<1/10 && rel_cov<2.5) && lod_baf_phase>10 || rel_cov<2.1 && lod_baf_phase>10 )' \
   <(gsutil cat gs://dsge-cromwell-runeberg/mocha/${wf_id}/call-mocha_stats_tsv/${prefix}.stats.tsv) \
   <(gsutil cat gs://dsge-cromwell-runeberg/mocha/${wf_id}/call-mocha_calls_tsv/${prefix}.calls.tsv) > ${prefix}.filtered.calls.tsv




#-----------------------------------------------------------
## Step 2: Separate out mosaic variant calls into 8 categories:

gsutil cat ${sample_tsv}|awk 'NR>1{print $1}'|wc -l   # 122,252

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

# any_1: 38393 1
# auto_2: 4302 1
# X_3: 35725 1
# anyCF_5: 33407 1
# autoCF_6: 2261 1
# XCF_7: 31870 1


gsutil cp BatchFirst30.mca*.txt gs://dsge-aoxing/mocha/covid


#-----------------------------------------------------------
## Step 3: Create 2 COVID-19 phenotypes

covid_file="gs://finngen-production-library-red/finngen_R5/coroa/2020-06-11/data/finngen_R5_v7_infectious_disease_register_corona_2020-06-11.txt"

# 1) Severe COVID+ vs. All (all with available COVID-19 test reporting), names of endpoints in FinnGen: "U22_COVID19_CONFIRMED" & "U22_COVID19_SUSPECTED"
# not in the FinnGen R6 endpoints data
gsutil cat gs://finngen-production-library-red/finngen_R6/phenotype_2.0/documentation/finngen_R6_v2_endpoint_definitions.txt|grep COVID
gsutil cat gs://finngen-production-library-red/finngen_R6/phenotype_2.0/data/finngen_R6_v2_endpoint.gz|zcat|head -n 1|tr '\t' '\n'|grep COVID
gsutil cat gs://finngen-production-library-red/finngen_R6/phenotype_2.0/data/finngen_R6_cov_pheno_1.0.txt.gz|zcat|head -n 1|grep COVID


# R5, dataset for coroa
gsutil cat gs://finngen-production-library-red/finngen_R5/coroa/2020-06-11/data/finngen_R5_v7_infectious_disease_register_corona_2020-06-11.txt|wc -l





Note: In the UK Biobank, Severe COVID+ is someone who has tested positive for COVID at least once and has been hospitalized.
Note: as part of basic quality control, remove related individuals from the dataset (1st or 2nd degree relatives)
 
# and remove individuals with any prevalent hematologic cancer history at time of blood draw for genotyping (ICD C81-96, D45-47), 
# and remove individuals who died prior to 3/1/2020.


join -1 1 -2 1 <(gsutil cat ${covid_file}|awk '{print $1}'|sort) <(gsutil cat gs://dsge-aoxing/mocha/covid/${prefix}.mca_any_1.txt|sort -k1,1)  



#-----------------------------------------------------------
Step 4: Run the following logistic regression analyses for each of the 8 variant groupings and 2 COVID-19 phenotypes above: 
A.Unadjusted: a.COVID_Pheno ~ has_MosaicVar 
B.Sparsely adjusted: a.COVID_Pheno ~ has_MosaicVar + age +age2+sex + ever_smoking_status+ PC1-10 of genomic ancestry
C.Fully adjusted: a.COVID_Pheno ~ has_MosaicVar + age +age2+ sex + ever_smoking_status + scale(Townsend) + PC1-10 of genomic ancestry + scale(BMI)+ Asthma + COPD + CAD + HTN + T2D+ anyCancer


