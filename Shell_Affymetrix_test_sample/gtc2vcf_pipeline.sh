## AIM: data preparation for mCAs detection using MoChA (https://github.com/freeseek/mocha) 
## INPUT: Affymetrix Finngen_v1 genotype calls and intensities provided by Affymetrix
## OUTPUT: VCF with BAF and LRR
## TOOLS: bcftools +affy2vcf (https://github.com/freeseek/gtc2vcf), samtools, bwa




##########################
##       Preparation    ##
##########################
gcloud init                             # initialize the gcloud environment
source /homes/aliu/.bashrc              # important!
gcloud beta compute --project "mocha-pipeline" ssh --zone "europe-west1-b" "instance-2"     # pathword: blackfriday
gsutil cp -r gs://affymetrix-test-data/PRO100104_FinnGen_SAX_070518/ .                      # 278 re-genotyped hapmap samples using Affmetrix Finngen_v1 array


##########################
##      Installation    ##
##########################
#Install basic tools and HTSlib librairies 
sudo apt install wget autoconf zlib1g-dev bwa gzip unzip samtools msitools cabextract mono-devel libgdiplus libbz2-dev libssl-dev liblzma-dev libgsl0-dev plink1.9


# Preparation steps
mkdir -p $HOME/bin $HOME/res && cd /tmp


# Download latest version of HTSlib and BCFtools (if not downloaded already)
git clone --branch=develop git://github.com/samtools/htslib.git
git clone --branch=develop git://github.com/samtools/bcftools.git


# Download plugins code
/bin/rm -f bcftools/plugins/{gtc2vcf.{c,h},affy2vcf.c}
wget -P bcftools/plugins https://raw.githubusercontent.com/freeseek/gtc2vcf/master/{gtc2vcf.{c,h},affy2vcf.c}


# Compile latest version of HTSlib and BCFtools (require gcc version 5 or newer; gcc version could be checked by "gcc -v")
cd htslib && autoheader && (autoconf || autoconf) && ./configure --disable-bz2 --disable-gcs --disable-lzma && make && cd ..
cd bcftools && make && cd ..
/bin/cp bcftools/{bcftools,plugins/{gtc,affy}2vcf.so} $HOME/bin/


# Make sure the directory with the plugins is available to bcftools
export PATH="$HOME/bin:$PATH"
export BCFTOOLS_PLUGINS="$HOME/bin"


# Install the GRCh38 human genome reference (following the suggestion from Heng Li)
wget -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | \
  gzip -d > $HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
samtools faidx $HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna  # index the genome reference in the FASTA format 
bwa index $HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna       # index the genome reference



#   ##########################
#   ##    Remap manifest    ##
#   ##########################
#   # This step could be skipped since most of alternate contigs variants did not map anywhere without alternate contigs
#   # Generate an alignment file for the source sequences from a CSV manifest file (-M option to mark shorter split hits as secondary)
#   annot_file="$HOME/PRO100104_FinnGen_SAX_070518/Validation_new/AxiomReference/Axiom_FinnGen1.na36.r1.a1.annot.csv" 
#   ref="$HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
#   bam_alignment_file="$HOME/PRO100104_FinnGen_SAX_070518/Validation_new/AxiomReference/Axiom_FinnGen1.na36.r1.a1.GRCh38.annot.sam"
#   bcftools +affy2vcf \
#     -c $annot_file \
#     --fasta-flank | \
#     bwa mem -M $ref - | \
#     samtools view -bS \
#     -o $bam_alignment_file
#   
#   
#   # Compute the coordinates according to the reference used to align the manifest source sequences
#   annot_file="$HOME/PRO100104_FinnGen_SAX_070518/Validation_new/AxiomReference/Axiom_FinnGen1.na36.r1.a1.annot.csv" 
#   bam_alignment_file="$HOME/PRO100104_FinnGen_SAX_070518/Validation_new/AxiomReference/Axiom_FinnGen1.na36.r1.a1.GRCh38.annot.sam"
#   csv_realigned_file="$HOME/PRO100104_FinnGen_SAX_070518/Validation_new/AxiomReference/Axiom_FinnGen1.na36.r1.a1.GRCh38.annot.csv"
#   bcftools +affy2vcf \
#     -c $annot_file \
#     -s $bam_alignment_file \
#     -o $csv_realigned_file


  
##########################
##    convert to VCF    ##
##########################
# Convert Affymetrix genotype calls (genotype calling from raw CEL files was done by Affymetrix) and intensities to VCF  
annot_file="$HOME/PRO100104_FinnGen_SAX_070518/Validation_new/AxiomReference/Axiom_FinnGen1.na36.r1.a1.annot.csv"       
ref="$HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
out_prefix="test_GRCh38"
path_to_output_folder="$HOME/PRO100104_FinnGen_SAX_070518/Validation_new/cluster_final"

bcftools +affy2vcf \
  --no-version -Ou \
  --csv $annot_file \
  --fasta-ref $ref \
  --calls $path_to_output_folder/AxiomGT1.calls.txt \
  --confidences $path_to_output_folder/AxiomGT1.confidences.txt \
  --summary $path_to_output_folder/AxiomGT1.summary.txt \
  --snp-posteriors $path_to_output_folder/AxiomGT1.snp-posteriors.txt \
  --report $path_to_output_folder/AxiomGT1.report.txt \
  --sex $out_prefix.sex | \
  bcftools sort -Ou -T ./bcftools-sort.XXXXXX | \
  bcftools norm --no-version -Ob -o $out_prefix.bcf -c x -f $ref && \
  bcftools index -f $out_prefix.bcf

  

#########################
##    Plot variants    ##
#########################

# Install basic tools
sudo apt install r-cran-ggplot2 r-cran-data.table r-cran-gridextra


# Download R scripts
/bin/rm -f $HOME/bin/gtc2vcf_plot.R
wget -P $HOME/bin https://raw.githubusercontent.com/freeseek/gtc2vcf/master/gtc2vcf_plot.R
chmod a+x $HOME/bin/gtc2vcf_plot.R


# Plot single variant
out_prefix="test"
chrom="chr1"
pos=3251716
gtc2vcf_plot.R \
  --affymetrix \
  --vcf $out_prefix.bcf \
  --chrom ${chrom} \
  --pos ${pos} \
  --png ${chrom}_${pos}_22.png


# Copy the plot to local path
gcloud compute scp instance-2:/home/aoxliu/PRO100104_FinnGen_SAX_070518/Validation_new/cluster_final/chr1_3251716.png  /Users/aoxliu/Documents/Finngen_mCA --zone=europe-west1-b



###################################
##  genotype consistency check   ## Results showed that the vcf conversion by bcftools +affy2vcf and Affymetrix were both OK
###################################

# Comments from Giulio Genovese: the software that Affymetrix provides to convert to VCF is buggy and gets the coordinates for the indels wrong.
# if you don't have the IDs you could use something like KING to find duplicates when merging with 1000 Genomes data or something like that


# Check 1: compare with vcf provided by Affmetrix using PLINK
out_prefix="test_GRCh38"
bcftools annotate --no-version -x ID -I +'%CHROM:%POS:%REF:%ALT' $out_prefix.bcf \
  |plink1.9 --vcf /dev/stdin --biallelic-only --keep-allele-order --const-fid --make-bed --out $out_prefix.geno   # 278 people with 759,177 loci 
plink1.9 --allow-extra-chr --bfile $out_prefix.geno --bmerge $out_prefix.geno --merge-mode 6   # check data quality by merging with the data itself


affy_vcf="/home/aoxliu/PRO100104_FinnGen_SAX_070518/Validation_new/cluster_final/AxiomGT1.calls.vcf"
plink1.9 --vcf $affy_vcf --biallelic-only --keep-allele-order --const-fid --allow-extra-chr --make-bed --out affy_vcf.geno  # 278 people with 655,492 loci 
plink1.9 --allow-extra-chr --bfile affy_vcf.geno  --bmerge affy_vcf.geno  --merge-mode 6


for i in $out_prefix.geno affy_vcf.geno; do mv ${i}.bim ${i}_original.bim; awk '{print $1,$1":"$4,$3,$4,$5,$6}'  ${i}_original.bim > ${i}.bim; head ${i}.bim; done  # format loci name
plink1.9 --allow-extra-chr --bfile $out_prefix.geno --bmerge affy_vcf.geno --merge-mode 6   # snps with multiple alternative alleles
plink1.9 --allow-extra-chr --bfile $out_prefix.geno --exclude plink.missnp  -make-bed --out source1_tmp
plink1.9 --allow-extra-chr --bfile affy_vcf.geno --exclude plink.missnp  -make-bed --out source2_tmp
plink1.9 --allow-extra-chr --bfile source1_tmp --bmerge source2_tmp --merge-mode 6 && rm source1_tmp.* && rm source2_tmp.*
# 181721372 overlapping calls, 180767174 nonmissing in both filesets. 179206137 concordant, for a concordance rate of 0.991364.


#-----------------------------------
# Check 2: compare with hapmap vcf 
# Link of CEL_ID and hapmap ID provided by Timo
cat <(awk -F '\t' 'match($4, /^[HG|NA]/){print $1,$4}' /home/aoxliu/PRO100104_FinnGen_SAX_070518/Validation_new/cluster_final/additional_sample_data.txt) \
    <(awk -F '\t' 'match($11, /^[HG|NA]/){print $1,$11}' /home/aoxliu/PRO100104_FinnGen_SAX_070518/Validation_new/cluster_final/additional_sample_data.txt)|sort|uniq > CEL_hapmap_Affy.id_lst
wc -l CEL_hapmap_Affy.id_lst && awk '{print $1}' CEL_hapmap_Affy.id_lst|sort|uniq|wc -l    # 285 lines for both


# create links of CEL_ID and hapmap ID using KING
# Install KING
wget http://people.virginia.edu/~wc9c/KING/KINGcode.tar.gz
tar -xzvf KINGcode.tar.gz
c++ -lm -lz -O2 -fopenmp -o king *.cpp


# Prepare 1000 genome data of chromosome 1 for KING (require plink binary format)
bcftools annotate --no-version -x ID -I +'%CHROM:%POS:%REF:%ALT' $HOME/res/kgp/ALL.chr1_GRCh38.genotypes.20170504.bcf \
  |plink1.9 --vcf /dev/stdin --biallelic-only --keep-allele-order --const-fid --make-bed --out ref_chr1.geno    


plink1.9 --allow-extra-chr --bfile ${out_prefix}.geno --chr 1 -make-bed --out ${out_prefix}_chr1.geno  &&  awk '{print $2}' ${out_prefix}_chr1.geno.bim > chr1.snp_lst    
mv ref_chr1.geno.bim  ref_chr1.geno_original.bim && awk '{print $1,$1":"$4,$3,$4,$5,$6}' ref_chr1.geno_original.bim > ref_chr1.geno.bim    # format loci name
plink1.9 --allow-extra-chr --bfile ref_chr1.geno  --extract chr1.snp_lst -make-bed --out ref_chr1_extract.geno     # only keep loci available in Finngen array


# Use KING to find duplicates when merging with 1000 Genomes data (only using chr 1)
/tmp/king -b ${out_prefix}_chr1.geno.bed,ref_chr1_extract.geno.bed --duplicate && wc -l king.con   # 204 pairs
/tmp/king -b ${out_prefix}_chr1.geno.bed,ref_chr1_extract.geno.bed --kinship && awk '$9>0.354{print $2,$3,$9}' king.kin|wc -l  # 206 pairs, kinship >0.354 is MZ (equivalent to --duplicate). However, two pairs ("a550934-4363868-032220-932_E06.CEL-NA19204","NA19331-NA19334") detected here are missing from king.con
awk 'match($4, /^[HG|NA]/){print $2,$4}' king.con > CEL_hapmap_KING.id_lst  &&  wc -l CEL_hapmap_KING.id_lst    # 187 lines


# check links provided by Timo and generated from KING and make a final link list
join -1 1 -2 1 <(sort -dk1,1 CEL_hapmap_KING.id_lst) <(sort -dk1,1 CEL_hapmap_Affy.id_lst) > CEL_hapmap.id_lst
wc -l CEL_hapmap.id_lst  && awk '$2==$3' CEL_hapmap.id_lst|wc -l   # 187 lines for both
awk '{print $1}' CEL_hapmap.id_lst|sort|uniq|wc -l  && awk '{print $2}' CEL_hapmap.id_lst|sort|uniq|wc -l  # 187 lines and 178 lines
join -1 2 -2 1 <(sort -dk2,2 CEL_hapmap_KING.id_lst) <(awk '{print $2}' ref_chr1_extract.geno.fam|sort -dk1,1)|wc -l  # 187 lines 
join -1 2 -2 1 <(sort -dk2,2 CEL_hapmap_Affy.id_lst) <(awk '{print $2}' ref_chr1_extract.geno.fam|sort -dk1,1)|wc -l  # 191 lines (3 pairs was due to test sample has no sequence results, and 1 pair due to detected by --kinship but not by --duplicate). Generally, KING is good when there is no ID links available.



# Use PLINK to check the concordance rates
# For test samples, change CEL_ID to hapmap_ID for
awk '{print 0,$1}' CEL_hapmap.id_lst > CEL.id_lst
plink1.9 --allow-extra-chr --bfile ${out_prefix}.geno --keep CEL.id_lst -make-bed --out ${out_prefix}_keep.geno
mv ${out_prefix}_keep.geno.fam ${out_prefix}_keep.geno.original.fam  &&  join -1 1 -2 1 <(awk '{print $2,NR}' ${out_prefix}_keep.geno.original.fam|sort -k1,1) <(sort -k1,1 CEL_hapmap.id_lst)|sort -nk2,2|awk '{print 0,$3,0,0,0,-9}' > ${out_prefix}_keep.geno.fam
awk '{print 0,$2}' ${out_prefix}_keep.geno.fam > hapmap.id_lst

# For 1000 genome hapmap, extract Finngen loci for Finngen test samples
awk '{print $2}' ${out_prefix}_keep.geno.bim > test_all.snp_lst
for i in {1..22} X Y
do
echo $i
bcftools annotate --no-version -x ID -I +'%CHROM:%POS:%REF:%ALT' $HOME/res/kgp/ALL.chr${i}_GRCh38.genotypes.20170504.bcf \
  |plink1.9 --vcf /dev/stdin --biallelic-only --keep-allele-order --const-fid --keep hapmap.id_lst --make-bed --out ref_chr${i}.geno
mv ref_chr${i}.geno.bim  ref_chr${i}.original.geno.bim  &&  awk '{print $1,$1":"$4,$3,$4,$5,$6}' ref_chr${i}.original.geno.bim > ref_chr${i}.geno.bim    # format loci name
plink1.9 --allow-extra-chr --bfile ref_chr${i}.geno --extract test_all.snp_lst -make-bed --out ref_chr${i}_extract.geno     # only keep loci available in Finngen array  
done

for i in {1..22} X Y
do
echo ref_chr${i}_extract.geno
#echo ref_chr${i}_extract.geno >> merge.lst
echo ref_chr${i}_extract_flip.geno >> merge_flip.lst
done

plink1.9 --allow-extra-chr --merge-list merge.lst --make-bed --out ref_chrall_extract.geno

for i in {1..22} X Y
do
plink1.9 --allow-extra-chr --bfile ref_chr${i}_extract.geno --exclude ref_chrall_extract.geno-merge.missnp  -make-bed --out ref_chr${i}_extract_flip.geno
done

plink1.9 --allow-extra-chr --merge-list merge_flip.lst --make-bed --out ref_chrall_extract.geno


# bmerge of test samples and 1000 genome hapmap
plink1.9 --allow-extra-chr --bfile  ${out_prefix}_keep.geno  --exclude plink.missnp  -make-bed --out ${out_prefix}_keep_flip.geno 

plink1.9 --allow-extra-chr --bfile  ${out_prefix}_keep_flip.geno  --exclude plink.missnp  -make-bed --out ${out_prefix}_keep_flip.geno 
plink1.9 --allow-extra-chr --bfile  ref_chrall_extract.geno  --exclude plink.missnp  -make-bed --out ref_chrall_extract.geno


plink1.9 --allow-extra-chr --bfile ${out_prefix}_keep_flip.geno --bmerge ref_chrall_extract.geno --merge-mode 6

# Performing 1-pass diff (mode 6), writing results to plink.diff .
# 93058044 overlapping calls, 92592149 nonmissing in both filesets.
# 91473528 concordant, for a concordance rate of 0.987919.


for i in {1..22} X Y
do
echo $i
plink1.9 --allow-extra-chr --bfile ${out_prefix}_keep.geno --bmerge ref_chr${i}_extract.geno --merge-mode 6 
done



