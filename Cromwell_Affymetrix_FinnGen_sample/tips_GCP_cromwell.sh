
# From laptop (recommended) or VM instance. Command will just make a "bridge" between the machine and the cromwell server, and them the VM will be set by the cromwell pipeline you submitted.
gcloud compute ssh cromwell-runeberg --project finngen-refinery-dsgelab --zone europe-west1-b -- -NL 8000:localhost:8000   # then could asscess the cromwell server mode from http://localhost:8000



# Two options if the address 8000 already been used
gcloud compute ssh cromwell-runeberg --project finngen-refinery-dsgelab --zone europe-west1-b -- -NL 8081:localhost:8000
gcloud compute ssh cromwell-runeberg --project finngen-refinery-dsgelab --zone europe-west1-b -- -R  8000:localhost:8000


