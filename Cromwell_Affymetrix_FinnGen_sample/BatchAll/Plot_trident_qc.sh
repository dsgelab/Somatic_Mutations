

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


