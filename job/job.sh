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
export MARLIN_DLL=$MARLIN_DLL:${project_folder}/lib/libTrackLengthDebug.so && \

# This is if I need to recompile my own version of some pre-imported libraries of ilcsoft, like MarlinReco or MarlinUtils, etc.
# source /cvmfs/ilc.desy.de/key4hep/setup.sh && \
# export MARLIN_DLL=/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/MarlinDD4hep/v00-06/lib/libMarlinDD4hep.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/DDMarlinPandora/v00-12/lib/libDDMarlinPandora.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/MarlinReco/v01-32/lib/libMarlinReco.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/PandoraAnalysis/v02-00-01/lib/libPandoraAnalysis.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/LCFIVertex/v00-08/lib/libLCFIVertexProcessors.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/CEDViewer/v01-19/lib/libCEDViewer.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/Overlay/v00-22-03/lib/libOverlay.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/MarlinFastJet/v00-05-02/lib/libMarlinFastJet.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/LCTuple/v01-13/lib/libLCTuple.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/MarlinKinfit/v00-06/lib/libMarlinKinfit.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/MarlinTrkProcessors/v02-12/lib/libMarlinTrkProcessors.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/MarlinKinfitProcessors/v00-04-02/lib/libMarlinKinfitProcessors.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/ILDPerformance/v01-10/lib/libILDPerformance.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/Clupatra/v01-03/lib/libClupatra.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/Physsim/v00-04-01/lib/libPhyssim.so:/afs/desy.de/user/d/dudarboh/analysis/iLCSoft/LCFIPlus/lib/libLCFIPlus.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/FCalClusterer/v01-00-01/lib/libFCalClusterer.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/ForwardTracking/v01-14/lib/libForwardTracking.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/ConformalTracking/v01-11/lib/libConformalTracking.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/LICH/v00-01/lib/libLICH.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/Garlic/v03-01/lib/libGarlic.so && \


mkdir ${cluster_id}_${proc_id} && cd ${cluster_id}_${proc_id}

# cp ${1} .
# use copied file locally, not the one from cvmfs! Previously it caused a lot of jobs to crash. But now it seems fine...
#filename=$(basename ${1})

Marlin ${project_folder}/xml/steer.xml --global.LCIOInputFiles="${file}"
# # Marlin sometimes seg. faults after successful finish so don't do &&...
mv *.root ../../final/${job_name}_${cluster_id}_${proc_id}.root
cd .. && rm -r ${cluster_id}_${proc_id}
