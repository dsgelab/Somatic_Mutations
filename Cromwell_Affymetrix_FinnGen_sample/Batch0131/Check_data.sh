### This script is to check the difference between b01 and b31, since b31 (FinnGen chip v2) keeps reporting errors while b01 (FinnGen chip v1) is able to get the expected output


#----------------------------------------------------------------
## Whether all probes in .txt files (.snp-posteriors.txt and .calls.txt ) are in .csv, there are some probes in chrM/chrX/chrY had -1 as suffix in .txt files   
## b01
# extract probes in .snp-posteriors.txt
gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b01/AxiomGT1_b01_V2*.snp-posteriors.txt|grep -v ^#|awk '{print $1}' > b01_snp_posterior

# b01_snp_calls

# b01_snp_summary

# probes unique in .snp-posteriors.txt compared to .csv
comm -13 <(gsutil cat gs://dsge-aoxing/mocha/input/Axiom_FinnGen1.na36.r1.a1.annot.csv|grep -v ^#|cut -d, -f1|sed 's/"//g'|sort) <(sort b01_snp_posterior)|grep -v probeset_id|head

# probes unique in .calls.txt compared to .csv
comm -13 <(gsutil cat gs://dsge-aoxing/mocha/input/Axiom_FinnGen1.na36.r1.a1.annot.csv|grep -v ^#|cut -d, -f1|sed 's/"//g'|sort) <(sort b01_snp_calls)|grep -v probeset_id|head   

# probes unique in .calls.txt compared to .snp-posteriors.txt
comm -13 <(sort b01_snp_posterior) <(sort b01_snp_calls)|grep -v probeset_id|head


## b31
# extract probes in .snp-posteriors.txt
gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b31/AxiomGT1_b31_V2.snp-posteriors.txt|grep -v ^#|awk '{print $1}' > b31_snp_posterior

# extract probes in .calls.txt
# b31_snp_calls

# extract probes in .summary.txt
# b31_snp_summary

# probes unique in .snp-posteriors.txt compared to .csv
comm -13 <(gsutil cat gs://dsge-aoxing/mocha/input/Axiom_FinnGen2.na36.r2.a2.annot.csv|grep -v ^#|cut -d, -f1|sed 's/"//g'|sort) <(sort b31_snp_posterior)|grep -v probeset_id|head

# probes unique in .calls.txt compared to .csv
comm -13 <(gsutil cat gs://dsge-aoxing/mocha/input/Axiom_FinnGen2.na36.r2.a2.annot.csv|grep -v ^#|cut -d, -f1|sed 's/"//g'|sort) <(sort b31_snp_calls)|grep -v probeset_id|head

# probes unique in .calls.txt compared to .snp-posteriors.txt
comm -13 <(sort b31_snp_posterior) <(sort b31_snp_calls)|grep -v probeset_id|head




#----------------------------------------------------------------
## Any probe with more than two alleles?
# count of genotype calls for the first 50 FinnGen samples, in order to better understand how the genotypes were coded in .calls.txt
gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b31/AxiomGT1_b31_V2.calls.mapped_selected.txt| grep -v "^#\|^probeset_id" |cut -f2-50 | tr '\t' '\n' | sort | uniq -c  

# Check why the probes which have genotypes coded as 7, 8 ,9, and 12
gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b31/AxiomGT1_b31_V2.calls.mapped_selected.txt| grep -v "^#\|^probeset_id" |cut -f2-50 |grep "7\|8\|9\|12"| tr '\t' '\n' | sort | uniq -c

# further check: taking probes with genotype as 7 as an example
paste -d '\t' <(gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b31/AxiomGT1_b31_V2.calls.mapped_selected.txt| grep -v "^#\|^probeset_id" |cut -f1) \ 
              <(gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b31/AxiomGT1_b31_V2.calls.mapped_selected.txt| grep -v "^#\|^probeset_id" |cut -f2-50|sed -e 's/[\t]//g') |awk '$2 ~ /7/ { print }' 

# AX-178458608	7777700070000070777790070070070779977077777770077
# AX-183851963	0000000000000000000000000000000000000007000000000
# AX-183852366	0000000000000070000000000000000000000000000000000
# AX-183854222	0000000000000000007000070007000000000000000070000

# count of alleles in .summary data
# cut -d'-' -f3 b01_snp_summary|sort|uniq -c   # need to double check


#------------------------------------------------------------------
## Test the codes to remove probes from .summary.txt which are not in .calls.txt 
gsutil cp /Users/aoxliu/Downloads/gtc2vcf_1.10.2-dev.zip   gs://dsge-aoxing/mocha/software/    # http://software.broadinstitute.org/software/gtc2vcf/gtc2vcf_1.10.2-dev.zip
gsutil ls -l gs://dsge-aoxing/mocha/software/gtc2vcf_1.10.2-dev.zip


# login VM instance
gcloud beta compute --project "finngen-refinery-dsgelab" ssh --zone "europe-west1-b" "aoxing-mocha"

gcsfuse --implicit-dirs  from-fg-datateam  from-fg-datateam 

cd /home/aoxliu/mocha/mocha_build/software
gsutil cp gs://dsge-aoxing/mocha/software/gtc2vcf_1.10.2-dev.zip .


# Set input files
annot_file="/home/aoxliu/mCA/input/from-fg-datateam/AxiomReference/Axiom_FinnGen2/V3/Axiom_FinnGen2.na36.r2.a2.annot.csv" 
ref="/home/aoxliu/mCA/software/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
ps_file="/home/aoxliu/mCA/input/from-fg-datateam/AxiomReference/Axiom_FinnGen2/V2/Axiom_FinnGen2.r1.step2.ps"

i_dir="/home/aoxliu/mCA/input/from-fg-datateam/cnv_intensity_data/AxiomGT1_b31"
call="AxiomGT1_b31_V2.calls.mapped_selected.txt"
summary="AxiomGT1_b31_V2.summary.mapped_selected.txt"   
# posteriors="AxiomGT1_b31_V2.snp-posteriors.txt"
# report="AxiomGT1_b31_V2.report_mapped_selected_sorted.txt"

o_dir="/home/aoxliu/mCA/output/b31"
o_prefix="fg_31"  


bcftools/bcftools +$PWD/affy2vcf.so \
--tags GT,BAF,LRR \
--csv ${annot_file}  \
--fasta-ref ${ref} \
--probeset-ids  ${ps_file} \
--calls ${i_dir}/${call} \
--summary ${i_dir}/${summary} \
--output-type b \
--output ${o_dir}/${o_prefix}.bcf


# double checking it did not lose pieces 
bcftools view /home/aoxliu/mCA/output/b31/fg_31.bcf|grep -v "^#\|^probeset_id"|wc -l
# 643360
less /home/aoxliu/mCA/input/from-fg-datateam/AxiomReference/Axiom_FinnGen2/V2/Axiom_FinnGen2.r1.step2.ps|grep -v "^#\|^probeset_id"|wc -l
# 643408
join -1 1 -2 1 <( bcftools view /home/aoxliu/mCA/output/b31/fg_31.bcf|grep -v "^#\|^probeset_id"|cut -f3|sort) <(less /home/aoxliu/mCA/input/from-fg-datateam/AxiomReference/Axiom_FinnGen2/V2/Axiom_FinnGen2.r1.step2.ps|grep -v "^#\|^probeset_id"|cut -f1|sort)|wc -l
# 643360


