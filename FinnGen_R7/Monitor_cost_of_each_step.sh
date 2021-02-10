# For runs in Google cloud, there is a way to monitor the costs of each step

# it requires passing an option for a monitoring script
# first of all, these four commands should get you the scripts needed
cd  ~/Documents/WDL_pipelines/MoChA/MonitorCost/    # or atlas "cd /homes/aliu/mCA/MonitorCost"
wget https://raw.githubusercontent.com/broadinstitute/gatk-sv/master/scripts/cromwell/cromwell_monitoring_script2.sh
wget https://raw.githubusercontent.com/broadinstitute/gatk-sv/master/scripts/cromwell/get_cromwell_resource_usage2.sh
wget https://raw.githubusercontent.com/broadinstitute/gatk-sv/master/scripts/cromwell/analyze_monitoring_logs2.py
sed -i 's/^  data = pd.read_table(log_file, usecols=lambda x: x not in (\x27IORead\x27, \x27IOWrite\x27))$/  data = pd.read_table(log_file, sep=\x27\\s+\x27, usecols=lambda x: x not in (\x27IORead\x27, \x27IOWrite\x27), comment=\x27H\x27)/' analyze_monitoring_logs2.py
wget https://raw.githubusercontent.com/broadinstitute/gatk-sv/master/scripts/cromwell/analyze_resource_acquisition.py
sed -i 's/^  job_id = call_info\[\x27jobId\x27\].split(\x27\/\x27)\[-1\]$/  job_id = call_info\[\x27jobId\x27\].split(\x27\/\x27)\[-1\] if \x27jobId\x27 in call_info else \x270\x27/' analyze_resource_acquisition.py
chmod a+x cromwell_monitoring_script2.sh get_cromwell_resource_usage2.sh analyze_monitoring_logs2.py analyze_resource_acquisition.py


# then you have copy the cromwell_monitoring_script2.sh somewhere on a Google bucket
gsutil cp cromwell_monitoring_script2.sh   gs://dsge-aoxing/cromwell/
gsutil cp get_cromwell_resource_usage2.sh  gs://dsge-aoxing/cromwell/
gsutil cp analyze_monitoring_logs2.py  gs://dsge-aoxing/cromwell/
gsutil cp analyze_resource_acquisition.py  gs://dsge-aoxing/cromwell/


# when you submit a job, include monitoring_script in the option file:
java -jar cromwell-55.jar submit mocha.wdl -i finngen.mocha.json -o options.google.json

In the options.google.json file you have to indicate the location of the script. This is an example of a configuration file:

{
  "final_workflow_outputs_dir": "gs://mccarroll-mocha/finngen",
  "use_relative_output_paths": true,
  "final_workflow_log_dir": "gs://mccarroll-mocha/cromwell/wf_logs",
  "final_call_logs_dir": "gs://mccarroll-mocha/cromwell/call_logs",
  "delete_intermediate_output_files": true,
  "monitoring_script": "gs://mccarroll-mocha/cromwell/cromwell_monitoring_script2.sh"
}



# by doing that, each task directory will contain a monitoring.log including all sort of statistics about the task
# there is then a script that collects all of these statistics and makes a table and a figure
get_cromwell_resource_usage2.sh -o $job_id.tsv gs://mccarroll-mocha/cromwell/cromwell-executions/mocha/$job_id
analyze_monitoring_logs2.py $job_id.tsv $job_id

# where job_id is something like a6fbd23b-8f8f-4d8d-b6c7-040e135ec0d4 and it generates a table and a plot, like this: 
curl -X GET http://localhost:8000/api/workflows/v1/$job_id/metadata > $job_id.json
analyze_resource_acquisition.py --override-warning $job_id.json $job_id

# will generate a plot


