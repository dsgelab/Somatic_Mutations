cd  /Users/aoxliu/Documents/Project2_Finngen_mCA/Analysis_FinnGen_mCA/FinnGenBatchAll/Result/mhc_qc


## further QC on BTA6
less ../FinnGenBatchAll.filtered.calls.tsv|grep chr6|awk '$4>20000000 && $5<40000000{print "gs://dsge-cromwell/mocha/3bb7e444-d05c-49d7-9e48-d5e012240aa9/call-mocha_plot/shard-*/pngs/",$1,".",$3,"_",$4,"_",$5,".png"}'|sed 's/ //g' >  potential_mhc.lst 
less ../FinnGenBatchAll.filtered.calls.tsv|grep chr6|awk '$4>20000000 && $5<40000000{print "gs://dsge-cromwell/mocha/3bb7e444-d05c-49d7-9e48-d5e012240aa9/call-mocha_plot/shard-*/attempt-2/pngs/",$1,".",$3,"_",$4,"_",$5,".png"}'|sed 's/ //g' >>  potential_mhc.lst 

wc -l potential_mhc.lst   # 208




for i in {1..208}
do
echo ${i}
png=$(awk -v i=${i} 'NR==i' potential_mhc.lst) 
gsutil cp ${png} .

done



