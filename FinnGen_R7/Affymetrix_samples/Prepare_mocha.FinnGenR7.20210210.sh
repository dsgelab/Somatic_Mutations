
## Prepare inputs for WDL pipelines

gcloud compute ssh "cromwell-dsge" --project "finngen-refinery-dsgelab" --zone europe-west1-b -- -fN -L localhost:8000:localhost:80

gcloud compute ssh "cromwell-dsge" --project "finngen-refinery-dsgelab" --zone europe-west1-b 
cd  /home/aoxliu/prepare_mocha_FinnGenR7



#----------------------------------------------------------------------------------
## Prepare sample_tsv_file
cd  /home/aoxliu/prepare_mocha_FinnGenR7
echo "sample_id batch_id cel"|tr ' ' '\t' > temp_FinnGenR7_Affymetrix.sample.tsv   # header

for i in b{01..63}
do
echo $i
gsutil cat gs://from-fg-datateam/R7_cnv_intensity_data/AxiomGT1_${i}/AxiomGT1_${i}*_V*.calls.mapped_selected.txt|grep -v ^#|head -n1|tr '\t' '\n'|tail -n+2|awk -v b=${i} '{print $1,b,$1}'|tr ' ' '\t' >> temp_FinnGenR7_Affymetrix.sample.tsv 
done

awk -F"\t" -v OFS="\t" '{x[$1]++; if (x[$1]>1) $1=$1"-"(x[$1]-1); print}'  temp_FinnGenR7_Affymetrix.sample.tsv  >  FinnGenR7_Affymetrix.sample.tsv  # add suffix(-1) to sample_id to make it unique
awk '{print $1}' FinnGenR7_Affymetrix.sample.tsv|cut -d'-' -f1|sort|uniq -c|awk '$1!=1'|wc -l  # 136

gsutil cp  FinnGenR7_Affymetrix.sample.tsv  gs://dsge-aoxing/mocha/FinnGenR7/

gsutil cat gs://dsge-aoxing/mocha/FinnGenR7/FinnGenR7_Affymetrix.sample.tsv|wc -l # 249,162 for FinnGen R7
gsutil cat gs://from-fg-datateam/check/BatchAll_sort_dup.sample.tsv|wc -l   # 201,459 for FinnGen R6



#----------------------------------------------------------------------------------
## Prepare batch_tsv_file
echo "batch_id csv snp calls summary report"|tr ' ' '\t' > FinnGenR7_Affymetrix.batch.tsv
paste -d '\t' <(echo b{01..63}|tr ' ' '\n') \
              <(cat <(printf 'gs://dsge-aoxing/mocha/input/Axiom_FinnGen1.na36.r1.a1.annot.csv\n%.0s' {1..30}) <(printf 'gs://dsge-aoxing/mocha/input/Axiom_FinnGen2.na36.r2.a2.annot.csv\n%.0s' {31..52}) <(printf 'gs://dsge-aoxing/mocha/input/Axiom_FinnGen2.na36.r2.a2.annot.csv\n%.0s' {53..63})) \
              <(gsutil ls gs://from-fg-datateam/R7_cnv_intensity_data/AxiomGT1_b*/AxiomGT1*.snp-posteriors.txt) \
              <(gsutil ls gs://from-fg-datateam/R7_cnv_intensity_data/AxiomGT1_b*/AxiomGT1*.calls.mapped_selected.txt) \
              <(gsutil ls gs://from-fg-datateam/R7_cnv_intensity_data/AxiomGT1_b*/AxiomGT1*.summary.mapped_selected.txt) \
              <(gsutil ls gs://from-fg-datateam/R7_cnv_intensity_data/AxiomGT1_b*/AxiomGT1*.report_mapped_selected.txt) >> FinnGenR7_Affymetrix.batch.tsv

cat FinnGenR7_Affymetrix.batch.tsv|awk '{print NF}'|sort|uniq -c   # check whether all rows have 6 row
gsutil cp  FinnGenR7_Affymetrix.batch.tsv  gs://dsge-aoxing/mocha/FinnGenR7/


# !!! with OTV for batch 01..30, without OTV for batch 31..51
gsutil cat gs://from-fg-datateam/R7_cnv_intensity_data/AxiomGT1_b11/AxiomGT1_b11_V2P2.snp-posteriors.txt|grep -v ^#|head -n 2|cut -c -100
# id	BB	AB	AA	CV	OTV
# AFFX-SP-000001	-2.940359,0.050627,1397.733934,1397.733934,9.386564,0.076372,-0.004950	0.004702,0.050

gsutil cat gs://from-fg-datateam/R7_cnv_intensity_data/AxiomGT1_b41/AxiomGT1_b41_V2.snp-posteriors.txt|grep -v ^#|head -n 2|cut -c -100
# id	BB	AB	AA	CV
# AFFX-SP-000001	-2.780723,0.043241,1163.999380,1163.999380,9.152170,0.084921,-0.030299	-0.016655,0.04


gsutil cat gs://from-fg-datateam/R7_cnv_intensity_data/AxiomGT1_b55/AxiomGT1_b55_update1_V3.snp-posteriors.txt|grep -v ^#|head -n 2|cut -c -100
# id	BB	AB	AA	CV
# AFFX-SP-000001	-2.844141,0.042691,1303.000324,1303.000324,9.192695,0.079138,-0.034940	-0.042299,0.04




#----------------------------------------------------------------------------------
## check whether batch with V3 used the same manifest file with V2
gsutil cat  gs://from-fg-datateam/R7_cnv_intensity_data/AxiomGT1_b43/AxiomGT1_b43_V2.snp-posteriors.txt|grep -v ^#|awk '{print $1}' > b43.snplst
gsutil cat  gs://from-fg-datateam/R7_cnv_intensity_data/AxiomGT1_b63/AxiomGT1_b63_V3.snp-posteriors.txt|grep -v ^#|awk '{print $1}' > b63.snplst
paste -d ' '  b43.snplst   b63.snplst|awk '$1!=$2'|wc -l



#----------------------------------------------------------------------------------
## Check whether every required file are complete 
# check report
for i in b{01..63}; do echo $i; gsutil cat gs://from-fg-datateam/R7_cnv_intensity_data/AxiomGT1_${i}/AxiomGT1_${i}*_V*.report_mapped_selected.txt|grep -v ^#|head -n 5 ; done
# b31..b61, b63


# check snp-posteriors
for i in b{01..63}; do echo $i; gsutil cat gs://from-fg-datateam/R7_cnv_intensity_data/AxiomGT1_${i}/AxiomGT1_${i}*_V*.snp-posteriors.txt|grep -v ^#|head -n 5|cut -c -100 ; done


# check calls
for i in b{01..63}; do echo $i; gsutil cat gs://from-fg-datateam/R7_cnv_intensity_data/AxiomGT1_${i}/AxiomGT1_${i}*_V*.calls.mapped_selected.txt|grep -v ^#|head -n 5|cut -c -100 ; done


# check summary
for i in b{01..63}; do echo $i; gsutil cat gs://from-fg-datateam/R7_cnv_intensity_data/AxiomGT1_${i}/AxiomGT1_${i}*_V*.summary.mapped_selected.txt|grep -v ^#|head -n 5|cut -c -100 ; done


