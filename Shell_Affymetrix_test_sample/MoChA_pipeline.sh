## AIM: call mosaic chromosomal alterations (mCA) from phased VCF and LRR/BAF for hundreds of FinnGen test samples, in order to build mCA dection pipeline for FinnGen
## INPUT: phased VCF and LRR/BAF prepared by gtc2vcf_pipeline.sh
## OUTPUT: both somatic and germline chromosomal alterations (CA)
## TOOLS: bcftools +mocha (https://github.com/freeseek/mocha), samtools, bedtools, eagle



# Setup GCP working environment
gcloud init                             # initialize the gcloud environment
source /homes/aliu/.bashrc              # important!
gcloud beta compute --project "mocha-pipeline" ssh --zone "europe-west1-b" "instance-2"

gsutil cp -r gs://affymetrix-test-data/PRO100104_FinnGen_SAX_070518/ .     
cd /home/aoxliu/test_data/call2vcf




##########################
##      Installation    ##
##########################

# Install basic tools and HTSlib librairies 
sudo apt install wget autoconf zlib1g-dev gzip unzip samtools bedtools libbz2-dev libssl-dev liblzma-dev libgsl0-dev


# Preparation steps
mkdir -p $HOME/bin $HOME/res $HOME/res/kgp && cd /tmp


# Download latest version of HTSlib and BCFtools (if not downloaded already)
# git clone --branch=develop git://github.com/samtools/htslib.git
# git clone --branch=develop git://github.com/samtools/bcftools.git


# Add patches and code for plugin
#/bin/rm -f bcftools/{{Makefile,main}.patch,vcfmocha.c,{beta_binom,genome_rules}.{c,h}} bcftools/plugins/{trio-phase,mochatools,importFMT,extendFMT}.c
/bin/rm -f bcftools/,vcfmocha.c,{beta_binom,genome_rules}.{c,h}} bcftools/plugins/{trio-phase,mochatools,importFMT,extendFMT}.c
wget -P bcftools https://raw.githubusercontent.com/freeseek/mocha/master/{{Makefile,main}{,patch},vcfmocha.c,{beta_binom,genome_rules}.{c,h}}
#wget -P bcftools https://raw.githubusercontent.com/freeseek/mocha/master/{{Makefile,main},vcfmocha.c,{beta_binom,genome_rules}.{c,h}}
wget -P bcftools/plugins https://raw.githubusercontent.com/freeseek/mocha/master/{trio-phase,mochatools,importFMT,extendFMT}.c
cd bcftools && patch < Makefile.patch && patch < main.patch && cd ..


# Compile latest version of HTSlib (optionally disable bz2, gcs, and lzma) and BCFtools (make sure you are using gcc version 5 or newer)
cd htslib && autoheader && (autoconf || autoconf) && ./configure --disable-bz2 --disable-gcs --disable-lzma && make && cd ..
cd bcftools && make && cd ..
/bin/cp bcftools/{bcftools,plugins/{fill-tags,fixploidy,split,trio-phase,mochatools,importFMT,extendFMT}.so} $HOME/bin/


# /bin/cp bcftools/{plugins/{fill-tags,fixploidy,split,trio-phase,mochatools,importFMT,extendFMT}.so} $HOME/bin/
bcftools plugin -l      # check all available plugin for bcftools


# Notice that you will need some functionalities missing from the base version of bcftools to run the pipeline
# Make sure the directory with the plugins is available to bcftools
export PATH="$HOME/bin:$PATH"
export BCFTOOLS_PLUGINS="$HOME/bin"


# Install latest version of Eagle
wget -O $HOME/bin/eagle https://data.broadinstitute.org/alkesgroup/Eagle/downloads/dev/eagle_v2.4.1
chmod a+x $HOME/bin/eagle




##########################
##   Download GRCh38    ##
##########################

# Download Genetic map
wget -P $HOME/res https://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz


# 1000 Genomes project phase 3 (fixing contig names, removing duplicate variants, removing incomplete variants)
cd $HOME/res/kgp
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr{{1..22},X,Y}_GRCh38.genotypes.20170504.vcf.gz{,.tbi}
ref="$HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
for chr in {1..22} X Y; do
  (bcftools view --no-version -h ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz | \
    grep -v "^##contig=<ID=[GNh]" | sed 's/^##contig=<ID=MT/##contig=<ID=chrM/;s/^##contig=<ID=\([0-9XY]\)/##contig=<ID=chr\1/'; \
  bcftools annotate --no-version -x INFO/END ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz | \   # newly added for debug
  bcftools view --no-version -H -c 2 | \
  grep -v "[0-9]|\.\|\.|[0-9]" | sed 's/^/chr/') | \
  bcftools norm --no-version -Ou -m -any | \
  bcftools norm --no-version -Ob -o ALL.chr${chr}_GRCh38.genotypes.20170504.bcf -d none -f $ref && \
  bcftools index -f ALL.chr${chr}_GRCh38.genotypes.20170504.bcf
done


# List of common germline duplications and deletions
wget -P $HOME/res ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/supporting/GRCh38_positions/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz{,.tbi}
bcftools query -i 'AC>1 && END-POS+1>10000 && TYPE!="INDEL" && (SVTYPE=="CNV" || SVTYPE=="DEL" || SVTYPE=="DUP")' \
  -f "chr%CHROM\t%POS0\t%END\t%SVTYPE\n" $HOME/res/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz > $HOME/res/cnp.grch38.bed


# Minimal divergence intervals from segmental duplications (make sure your bedtools version is not affected by the groupby bug)
wget -O- http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz | gzip -d |
  awk '!($2=="chrX" && $8=="chrY" || $2=="chrY" && $8=="chrX") {print $2"\t"$3"\t"$4"\t"$30}' > $HOME/res/genomicSuperDups.bed
  
awk '{print $1,$2; print $1,$3}' $HOME/res/genomicSuperDups.bed | \
  sort -k1,1 -k2,2n | uniq | \
  awk 'chrom==$1 {print chrom"\t"pos"\t"$2} {chrom=$1; pos=$2}' | \
  bedtools intersect -a $HOME/res/genomicSuperDups.bed -b - | \
  bedtools sort | \
  bedtools groupby -c 4 -o min | \
  awk 'BEGIN {i=0; s[0]="+"; s[1]="-"} {if ($4!=x) i=(i+1)%2; x=$4; print $0"\t0\t"s[i]}' | \
  bedtools merge -s -c 4 -o distinct | \
  grep -v "GL\|KI" | bgzip > $HOME/res/dup.grch38.bed.gz && \
  tabix -f -p bed $HOME/res/dup.grch38.bed.gz


# Download cytoband file
wget -O $HOME/res/cytoBand.hg38.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz


# Setup variables
ref="$HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
map="$HOME/res/genetic_map_hg38_withX.txt.gz"
kgp_pfx="$HOME/res/kgp/ALL.chr"
kgp_sfx="_GRCh38.genotypes.20170504"
rule="GRCh38"
cnp="$HOME/res/cnp.grch38.bed"
dup="$HOME/res/dup.grch38.bed.gz"
cyto="$HOME/res/cytoBand.hg38.txt.gz"




##########################
##  Data preparation    ##
##########################

# Preparation steps
vcf="$HOME/test_data/call2vcf/test_GRCh38.bcf"  # input VCF file with phased GT, LRR, and BAF
thr="1"        # number of threads to use
pfx="test"     # output prefix
out_prefix="test_GRCh38"
sex="$HOME/test_data/call2vcf/$out_prefix.sex" 
#xcl="..."     # VCF file with additional list of variants to exclude (optional)
#ped="..."     # pedigree file to use if parent child duos are present
dir="$HOME/test_data/call2vcf"  # directory where output files will be generated
mkdir -p $dir


# Create a minimal binary VCF (with ALLELE_A, ALLELE_B, GT, BAF, and LRR)
bcftools annotate --no-version -Ob -o $dir/$pfx.unphased.test.bcf $vcf \
  -x ID,QUAL,^INFO/ALLELE_A,^INFO/ALLELE_B,^FMT/GT,^FMT/BAF,^FMT/LRR && \
  bcftools index -f $dir/$pfx.unphased.test.bcf


# Perform basic quality control (the generated list of variants will be excluded from modeling by both eagle and mocha)
# This command will create a list of variants falling within segmental duplications with low divergence (<2%), high levels of missingness (>2%), variants with excess heterozygosity (p<1e-6), and variants that correlate with sex in an unexpected way (p<1e-6). 
n=$(bcftools query -l $dir/$pfx.unphased.bcf|wc -l); \
ns=$((n*98/100)); \
echo '##INFO=<ID=JK,Number=1,Type=Float,Description="Jukes Cantor">' | \
  bcftools annotate --no-version -Ou -a $dup -c CHROM,FROM,TO,-,JK -h /dev/stdin $dir/$pfx.unphased.bcf | \
  bcftools +fill-tags --no-version -Ou -- -t NS,ExcHet | \
  bcftools +mochatools --no-version -Ou -- -x $sex -G | \
  bcftools annotate --no-version -Ob -o $dir/$pfx.xcl.bcf \
    -i 'FILTER!="." && FILTER!="PASS" || JK<.02 || NS<'$ns' || ExcHet<1e-6 || AC_Sex_Test>6' \
    -x FILTER,^INFO/JK,^INFO/NS,^INFO/ExcHet,^INFO/AC_Sex_Test && \
  bcftools index -f $dir/$pfx.xcl.bcf
  



##########################
##    Phasing pipeline   #
##########################

# Phase VCF file by chromosome with Eagle
for chr in {1..22} X; do
  eagle \
    --geneticMapFile $map \
    --outPrefix $dir/$pfx.chr$chr \
    --numThreads $thr \
    --vcfRef $kgp_pfx${chr}$kgp_sfx.bcf \
    --vcfTarget $dir/$pfx.unphased.test.bcf \
    --vcfOutFormat b \
    --noImpMissing \
    --outputUnphased \
    --vcfExclude $dir/$pfx.xcl.bcf \
    --chrom $chr \
    --pbwtIters 3 && \
  bcftools index -f $dir/$pfx.chr$chr.bcf
done


# Extract chromosomes that do not require phasing  # Y and M
bcftools view --no-version -Ob -o $dir/$pfx.other.bcf $dir/$pfx.unphased.test.bcf \
  -t ^$(seq -s, 1 22),X,$(seq -f chr%.0f -s, 1 22),chrX && \
bcftools index -f $dir/$pfx.other.bcf


# Concatenate eagle output into a single VCF file and add GC/CpG content information
bcftools concat --no-version -Ou --threads $thr $dir/$pfx.{chr{{1..22},X},other}.bcf | \
bcftools +mochatools --no-version -Ob -o $dir/$pfx.bcf -- -f $ref && \
bcftools index -f $dir/$pfx.bcf


# Remove unphased VCF and single chromosome files (optional)
/bin/rm $dir/$pfx.{unphased,chr{{1..22},X},other}.bcf{,.csi} $dir/$pfx.chr{{1..22},X}.dose.bcf




#######################################
##  Chromosomal alterations pipeline  #
#######################################

# Preparation steps
pfx="test"    # output prefix
thr="1"       # number of extra threads to use
lst="..."     # file with list of samples to analyze for asymmetries (e.g. samples with 1p CNN-LOH)


# Call mosaic chromosomal alterations with MoChA
# For array data, MoChA's memory requirements will depend on the number of samples (N) and the number of variants (M) in the largest contig and will amount to 9NM bytes. For example, if you are running 4,000 samples and chromosome 1 has ~80K variants, you will need approximately 2-3GB to run MoChA.
bcftools mocha \
  --rules $rule \
  --no-version \
  --output-type b \
  --output $dir/$pfx.mocha.bcf \
  --threads $thr \
  --variants ^$dir/$pfx.xcl.bcf \
  --mosaic-calls $dir/$pfx.mocha.tsv \
  --genome-stats $dir/$pfx.stats.tsv \
  --ucsc-bed $dir/$pfx.ucsc.bed \
  --cnp $cnp \
  --no-BAF-flip \
  $dir/$pfx.bcf && \
bcftools index -f $dir/$pfx.mocha.bcf


# Filter the calls from MoChA, removing samples with BAF_CONC (from stats.tsv, since BAF_CONC from mocha.tsv is for that mCA instead of that genotyped sample) greater than 0.51, 
# removing calls smaller than 100kbp, removing calls with less than a LOD score of 10 for the model based on BAF and genotype phase, 
# removing calls flagged as germline copy number polymorphisms (CNPs), and removing calls with an estimated cell fraction larger than 50%
awk 'NR==FNR && FNR>1 && $6>.51 {x[$1]++}
  NR>FNR && (FNR==1 || !($1 in x) && $6>1e5 && $17>10 && $21!~"CNP" && $22<.5) {print}' \
  $dir/$pfx.stats.tsv $dir/$pfx.mocha.tsv > $dir/$pfx.mocha.filter.tsv




##########################
##     Post analysis     #
##########################

# Install basic tools (Debian/Ubuntu specific if you have admin privileges):
sudo apt install r-cran-ggplot2 r-cran-data.table


# Download R scripts
/bin/rm -f $HOME/bin/{plot_summary,mocha_plot}.R
wget -P $HOME/bin https://raw.githubusercontent.com/freeseek/mocha/master/{plot_summary,mocha_plot}.R
chmod a+x $HOME/bin/{plot_summary,mocha_plot}.R


# N of CAs (chromosomal alterations: somatic + germline) and filtered mCAs (mosaic chromosomal alterations: somatic only)
less $dir/$pfx.mocha.tsv|awk -F '\t' 'NR>1{print $21}'|sort|uniq -c           # all detected CAs, including CNN-LOH, Deletion, Duplication, CNP Deletion, CNP Duplication, and Undetermined
less  $dir/$pfx.mocha.filter.tsv|awk -F '\t' 'NR>1{print $21}'|sort|uniq -c   # flitered mCAs, including CNN-LOH, Deletion, and Duplication


# Generate summary plot
for i in fliter all
do
plot_summary.R \
  --pdf $dir/$pfx.${i}.pdf \
  --stats $dir/$pfx.stats.tsv \
  --calls $dir/$pfx.mocha.${i}.tsv
done


# Plot for each mosaic chromosomal alterations
for i in {1..125}  
do
SAMPLE=$(awk -v i=$i 'NR==i{print $1}' $dir/$pfx.mocha.filter.tsv)
TYPE=$(awk -v i=$i 'NR==i{print $21}' $dir/$pfx.mocha.filter.tsv)
CHROM=$(awk -v i=$i 'NR==i{print $3}' $dir/$pfx.mocha.filter.tsv)
BEG_GRCh38=$(awk -v i=$i 'NR==i{print $4}' $dir/$pfx.mocha.filter.tsv)
END_GRCh38=$(awk -v i=$i 'NR==i{print $5}' $dir/$pfx.mocha.filter.tsv)
echo -e "i:${i}\nSAMPLE:${SAMPLE}\nTYPE:${TYPE}\nCHROM:${CHROM}\nBEG_GRCh38:${BEG_GRCh38}\nEND_GRCh38:${END_GRCh38}"
mocha_plot.R \
  --mocha \
  --png $dir/mCA_plot/${TYPE}_${CHROM}_${BEG_GRCh38}_${END_GRCh38}.png \
  --vcf $dir/$pfx.mocha.bcf \
  --samples ${SAMPLE} \
  --regions ${CHROM}:${BEG_GRCh38}-${END_GRCh38}\
  --cytoband $HOME/res/cytoBand.hg38.txt.gz
done


# Copy plots to local path
gcloud compute scp instance-2:/home/aoxliu/test_data/call2vcf/test*pdf   /Users/aoxliu/Documents/Finngen_mCA/results --zone=europe-west1-b
gcloud compute scp instance-2:/home/aoxliu/test_data/call2vcf/mCA_plot/*.png  /Users/aoxliu/Documents/Finngen_mCA/results --zone=europe-west1-b


