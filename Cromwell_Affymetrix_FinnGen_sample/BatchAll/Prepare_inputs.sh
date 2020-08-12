## Prepare inputs for WDL pipelines


#-----------------------------------------
# Prepare sample_tsv_file

cd  /home/aoxliu/mocha/BatchAll

echo "sample_id batch_id cel"|tr ' ' '\t' > BatchAll_sort.sample.tsv   # header

for i in b{01..51}
do
echo $i
gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_${i}/AxiomGT1_${i}_V2*.calls.mapped_selected.txt|grep -v ^#|head -n1|tr '\t' '\n'|tail -n+2|awk -v b=${i} '{print $1,b,$1}'|tr ' ' '\t' >> BatchAll_sort.sample.tsv
done

awk -F"\t" -v OFS="\t" '{x[$1]++; if (x[$1]>1) $1=$1"-"(x[$1]-1); print}' BatchAll_sort.sample.tsv > BatchAll_sort_dup.sample.tsv  # add suffix(-1) to sample_id to make it unique
gsutil cp BatchAll_sort_dup.sample.tsv  gs://from-fg-datateam/check/


#-----------------------------------------
# Prepare batch_tsv_file

echo "batch_id csv snp calls confidences summary report"|tr ' ' '\t' > BatchAll_sort.batch.tsv
paste -d '\t' <(echo b{01..51}|tr ' ' '\n') <(printf 'gs://dsge-aoxing/mocha/input/Axiom_FinnGen1.na36.r1.a1.annot.csv\n%.0s' {1..51}) \
              <(gsutil ls gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b*/AxiomGT1*.snp-posteriors.txt) \
              <(gsutil ls gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b*/AxiomGT1*.calls.mapped_selected.txt) \
              <(gsutil ls gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b*/AxiomGT1*.confidences.mapped_selected.txt) \
              <(gsutil ls gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b*/AxiomGT1*.summary.mapped_selected.txt) \
              <(gsutil ls gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b*/AxiomGT1*.report_mapped_selected.txt|sed 's/'.txt'/'_sorted.txt'/g') >> BatchAll_sort.batch.tsv
gsutil cp BatchAll_sort_dup.sample.tsv  gs://from-fg-datateam/check/


# !!! with OTV for batch 01..30, without OTV for batch 31..51
gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b11/AxiomGT1_b11_V2P2.snp-posteriors.txt|grep -v ^#|head -n 2|cut -c -100
# id	BB	AB	AA	CV	OTV
# AFFX-SP-000001	-2.940359,0.050627,1397.733934,1397.733934,9.386564,0.076372,-0.004950	0.004702,0.050

gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b41/AxiomGT1_b41_V2.snp-posteriors.txt|grep -v ^#|head -n 2|cut -c -100
# id	BB	AB	AA	CV
# AFFX-SP-000001	-2.780723,0.043241,1163.999380,1163.999380,9.152170,0.084921,-0.030299	-0.016655,0.04



#-----------------------------------------
# Prepare sorted .report

for i in {01..51}
do
if [ $i -lt 31 ] ; then v="V2P2"; else v="V2"; fi
echo "b${i}_${v}"
gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b${i}/AxiomGT1_b${i}_${v}.calls.mapped_selected.txt|grep -v ^# | head -n1 | tr '\t' '\n' | tail -n+2 | \
       awk -F"\t" 'NR==FNR {x[NR]=$0} NR>FNR {if ($0~"^#" || $1=="cel_files") print; else y[$1]=$0} END {for (i in x) print y[x[i]]}' - <(gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b${i}/AxiomGT1_b${i}_${v}.report_mapped_selected.txt) > AxiomGT1_b${i}_${v}.report_mapped_selected_sorted.txt
gsutil cp AxiomGT1_b${i}_${v}.report_mapped_selected_sorted.txt  gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b${i}/
diff <(gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b${i}/AxiomGT1_b${i}_V2*.calls.mapped_selected.txt|grep -v ^#|head -n1|tr '\t' '\n'|tail -n+2) \
     <(gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b${i}/AxiomGT1_b${i}_V2*.report_mapped_selected_sorted.txt|grep -v ^#| cut -f1 | tail -n+2)|wc -l
rm AxiomGT1_b${i}_${v}.report_mapped_selected_sorted.txt
done



#-----------------------------------------
# Check whether all .txt files have the same orders

for i in b{01..51}
do
echo $i
# calls and summary
diff <(gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_${i}/AxiomGT1_${i}_V2*.calls.mapped_selected.txt|grep -v ^#|head -n1|tr '\t' '\n'|tail -n+2) \
     <(gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_${i}/AxiomGT1_${i}_V2*.summary.mapped_selected.txt|grep -v ^#|head -n1|tr '\t' '\n'|tail -n+2)|wc -l
# calls and confidences
diff <(gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_${i}/AxiomGT1_${i}_V2*.calls.mapped_selected.txt|grep -v ^#|head -n1|tr '\t' '\n'|tail -n+2) \
     <(gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_${i}/AxiomGT1_${i}_V2*.confidences.mapped_selected.txt|grep -v ^#|head -n1|tr '\t' '\n'|tail -n+2)|wc -l
# calls and report
diff <(gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_${i}/AxiomGT1_${i}_V2*.calls.mapped_selected.txt|grep -v ^#|head -n1|tr '\t' '\n'|tail -n+2) \
     <(gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_${i}/AxiomGT1_${i}_V2*.report_mapped_selected_sorted.txt|grep -v ^#| cut -f1 | tail -n+2)|wc -l
done



#-----------------------------------------
# Additional codes used to check the header (but not used for the pipeline)

# for i in b{01..51}
# do
# join <(gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_${i}/AxiomGT1_${i}_V2*.report_mapped_selected.txt|grep -v '#'|grep -v cel_files|awk '{print $1}'|sort -k1,1) \
#      <(gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_${i}/AxiomGT1_${i}_V2*.calls.mapped_selected.txt|grep '#'|grep cel-[1-9]|sed 's/=/ /g'|awk -v b=${i} '{print $3,b,$2}'|grep -v ^NA|sort -k1,1)|tr ' ' '\t' >> BatchAll_sort.sample.tsv
# n_1=$(gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_${i}/AxiomGT1_${i}_V2*.report_mapped_selected.txt|grep -v '#'|grep -v cel_files|wc -l)
# n_2=$(awk -v b=${i} '$2==b' BatchAll_sort.sample.tsv|wc -l)
# echo $i ${n_1} echo ${n_2}
# if [ ${n_1} -eq ${n_2} ] ; then echo "All samples can be found"; else echo "Some samples cannot be found"; fi
# done


