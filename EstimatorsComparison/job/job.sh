#!/bin/bash

# Submit file passes these arguments: $(file) $(ClusterId) $(ProcId) $(job_name), I rename them here as well, just for clarity
file=${1}
cluster_id=${2}
proc_id=${3}
job_name=${4}
project_folder=${5}

# Usually I just need some version of ilcsoft with one additional lib
# source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/init_ilcsoft.sh
source /cvmfs/ilc.desy.de/key4hep/setup.sh &&
export MARLIN_DLL=/cvmfs/ilc.desy.de/key4hep/spackages/marlindd4hep/0.6.2/x86_64-centos7-gcc11.2.0-opt/sxr7sivtshqs5azecdxlhqdd2s47bytm/lib/libMarlinDD4hep.so:/afs/desy.de/user/d/dudarboh/analysis/iLCSoft/MarlinReco/lib/libMarlinReco.so

mkdir ${cluster_id}_${proc_id} && cd ${cluster_id}_${proc_id}

# cp ${1} .
# use copied file locally, not the one from cvmfs! Previously it caused a lot of jobs to crash. But now it seems fine...
#filename=$(basename ${1})

Marlin ${project_folder}/xml/steer.xml --global.LCIOInputFiles="${file}"
# # Marlin sometimes seg. faults after successful finish so don't do &&...
mv *.root ../../final/${job_name}_${cluster_id}_${proc_id}.root
cd .. && rm -r ${cluster_id}_${proc_id}
