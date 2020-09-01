## Generate pile-up plot for each chromosome using Plot_pileup.R wrriten by Giulio
# https://github.com/freeseek/mocha/blob/master/pileup_plot.R

dir="/Users/aoxliu/Documents/Project2_Finngen_mCA/Analysis_FinnGen_mCA/FinnGenBatchAll/Result/"
pfx="FinnGenBatchAll"
wf_id="3bb7e444-d05c-49d7-9e48-d5e012240aa9"
cyto="cytoBand.hg38.txt.gz"

cd $dir


gsutil cp gs://dsge-aoxing/mocha/GRCh38/mocha.GRCh38/cytoBand.hg38.txt.gz  .
gsutil cp gs://dsge-cromwell/mocha/${wf_id}/call-mocha_stats_tsv/${pfx}.stats.tsv .
gsutil cp gs://dsge-cromwell/mocha/${wf_id}/call-mocha_calls_tsv/${pfx}.calls.tsv .
gsutil cp gs://dsge-cromwell/mocha/${wf_id}/call-mocha_calls_tsv/${prefix}.filtered.calls.tsv .


../Script/pileup_plot.R \
  --pdf ${pfx}.pileup.new.pdf \
  --stats ${pfx}.stats.tsv \
  --calls ${pfx}.filtered.calls.tsv \
  --cytoband ${cyto}

