
# From laptop (recommended) or VM instance. Command will just make a "bridge" between the machine and the cromwell server, and them the VM will be set by the cromwell pipeline you submitted.
gcloud compute ssh cromwell-runeberg --project finngen-refinery-dsgelab --zone europe-west1-b -- -NL 8000:localhost:8000   # then could asscess the cromwell server mode from http://localhost:8000



# Two options if the address 8000 already been used
gcloud compute ssh cromwell-runeberg --project finngen-refinery-dsgelab --zone europe-west1-b -- -NL 8081:localhost:8000
gcloud compute ssh cromwell-runeberg --project finngen-refinery-dsgelab --zone europe-west1-b -- -R  8000:localhost:8000



# Submit the pipelines
java -jar /home/vllorens/cromwell-48.jar submit -i mocha/B0102/input.json mocha/B0102/mocha.wdl



# Abort the pipeline
curl -X POST http://localhost:8000/api/workflows/v1/{id}/abort



# mount/unmount directory from google buscket to VM instance
gcsfuse --implicit-dirs  from-fg-datateam  from-fg-datateam   # mount
fusermount -u  /home/aoxliu/mCA/input/from-fg-datateam        # unmount

