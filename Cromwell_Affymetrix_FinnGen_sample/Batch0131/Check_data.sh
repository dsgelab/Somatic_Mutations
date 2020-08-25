

##----------------------------------------------------------------
gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b01/AxiomGT1_b01_V2*.snp-posteriors.txt|grep -v ^#|awk '{print $1}' > b01_snp_posterior
comm -13 <(gsutil cat gs://dsge-aoxing/mocha/input/Axiom_FinnGen1.na36.r1.a1.annot.csv|grep -v ^#|cut -d, -f1|sed 's/"//g'|sort) <(sort b01_snp_calls)|grep -v probeset_id|head
comm -13 <(gsutil cat gs://dsge-aoxing/mocha/input/Axiom_FinnGen1.na36.r1.a1.annot.csv|grep -v ^#|cut -d, -f1|sed 's/"//g'|sort) <(sort b01_snp_posterior)|grep -v probeset_id|head

comm -13 <(sort b01_snp_posterior) <(sort b01_snp_calls)|grep -v probeset_id|head




gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b31/AxiomGT1_b31_V2.snp-posteriors.txt|grep -v ^#|awk '{print $1}' > b31_snp_posterior
comm -13 <(gsutil cat gs://dsge-aoxing/mocha/input/Axiom_FinnGen2.na36.r2.a2.annot.csv|grep -v ^#|cut -d, -f1|sed 's/"//g'|sort) <(sort b31_snp_calls)|grep -v probeset_id|head
comm -13 <(gsutil cat gs://dsge-aoxing/mocha/input/Axiom_FinnGen2.na36.r2.a2.annot.csv|grep -v ^#|cut -d, -f1|sed 's/"//g'|sort) <(sort b31_snp_posterior)|grep -v probeset_id|head

comm -13 <(sort b31_snp_posterior) <(sort b31_snp_calls)|grep -v probeset_id|head


gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b31/AxiomGT1_b31_V2.calls.mapped_selected.txt| grep -v "^#\|^probeset_id" |cut -f2-50 |grep "7\|8\|9\|12"| tr '\t' '\n' | sort | uniq -c

paste -d '\t' <(gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b31/AxiomGT1_b31_V2.calls.mapped_selected.txt| grep -v "^#\|^probeset_id" |cut -f1) \ 
              <(gsutil cat gs://from-fg-datateam/cnv_intensity_data/AxiomGT1_b31/AxiomGT1_b31_V2.calls.mapped_selected.txt| grep -v "^#\|^probeset_id" |cut -f2-50|sed -e 's/[\t]//g') |awk '$2 ~ /7/ { print }' 

# AX-178458608	7777700070000070777790070070070779977077777770077
# AX-183851963	0000000000000000000000000000000000000007000000000
# AX-183852366	0000000000000070000000000000000000000000000000000
# AX-183854222	0000000000000000007000070007000000000000000070000

