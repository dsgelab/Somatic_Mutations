## This script is to check data (Phase 1: general inclusion list of R6 pass samples and single map file) for FinnGen mCA detection prepared by Timo Sipilä

# Login VM instance
gcloud beta compute --project "finngen-refinery-dsgelab" ssh --zone "europe-west1-b" "aoxing-mocha"  



# Mount Cloud Storage buckets as file systems in VM instance
cd /home/aoxliu/mCA/input
gcsfuse --implicit-dirs  from-fg-datateam  from-fg-datateam 
mkdir check



# N of batches
ls -r -lh /home/aoxliu/mCA/input/from-fg-datateam/cnv_intensity_data/|grep AxiomGT1|wc -l   # 49



# Extract ID and add BatchNo_LineNoWithinBatch as label
cd /home/aoxliu/mCA/input/from-fg-datateam/cnv_intensity_data

rm ../check/id.lst

for i in {01..51}
do
echo $i
less AxiomGT1_b${i}/AxiomGT1_b${i}_V2*.report_mapped_selected.txt|grep -v '#'|grep -v cel_files|awk -v b=${i} '{print $1,b"_"NR}' >> ../check/id.lst
done

wc -l ../check/id.lst                                             # 197,364 samples in total
awk '{print $1}' ../check/id.lst|sort|uniq -c|wc -l               # 197,228 unique sample id in total
awk '{print $1}' ../check/id.lst|sort|uniq -c|awk '$1!=1'|wc -l   # 136 id with duplicated samples



# N of samples for each batch
cat <(echo "batch N") <(awk '{print $2}' ../check/id.lst|cut -c 1-2|sort|uniq -c|awk '{print $2,$1}') > ../check/batch_N.lst



# N of samples only available in previous delivery (validation data for checking Timo's data extraction pipeline)
comm -23  <(less ../validation_b01_intensityd_somatic_cnv/AxiomGT1_b01.report_mapped_selected.txt|grep -v '#'|grep -v cel_files|awk '{print $1}') \
          <(awk '{print $1}' ../check/id.lst|sort)|wc -l          # 14, Timo said the validation data are included in the new delivery and should be removed
        
        
        
# Whether all outputs from Affymetrix genotype calling (inputs for mCA detection) are available for all batches
rm ../check/batch.lst

for t in summary snp-posteriors report confidences calls
do
echo ${t}
for i in {01..51}      # assuming batch numbers are from 01 to 51, but actually there are no b32 and b34 currently due to these two batches have too less samples after removing samples with restricted consent
do
ls -alh AxiomGT1_b${i}/AxiomGT1_b${i}_V2*.${t}*.txt >> ../check/batch.lst
done
done

       



