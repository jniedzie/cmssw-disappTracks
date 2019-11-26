#!/bin/bash

# Set those variables:
config_path=/afs/cern.ch/work/j/jniedzie/private/disapp_tracks/friendTreeProducer/CMSSW_10_2_4_patch1/src/CharginoAnalysis/CharginoAnalyzer/python/ConfFile_cfg.py
cmssw_path=/afs/cern.ch/work/j/jniedzie/private/disapp_tracks/friendTreeProducer/CMSSW_10_2_4_patch1/src
output_path=/eos/cms/store/group/phys_exotica/xtracks/7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/Wino_500GeV10cm/treeProducerXtracks/friend_trees/

#----------------------------------------------------------------------------
export XRD_NETWORKSTACK=IPv42
export CMSSWVER="CMSSW_10_2_4_patch1"
export SCRAM_ARCH="slc6_amd64_gcc700"

cd $cmssw_path
cmsenv

source /afs/cern.ch/cms/cmsset_default.sh
eval `scramv1 runtime -sh`

mkdir -p $output_path/chunk_$1

echo "cmsRun $config_path fileNumber=$1"
cmsRun $config_path fileNumber=$1
