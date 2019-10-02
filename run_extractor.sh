#!/bin/bash

# Set those variables:
config_path=/afs/cern.ch/work/j/jniedzie/private/CMSSW_9_4_6_patch1/src/CharginoAnalysis/CharginoAnalyzer/python/ConfFile_cfg.py
cmssw_path=/afs/cern.ch/work/j/jniedzie/private/CMSSW_9_4_6_patch1/src
output_path=/afs/cern.ch/work/j/jniedzie/private/disapp_tracks/pionBackground/friendInfo/

#----------------------------------------------------------------------------
export XRD_NETWORKSTACK=IPv4
export CMSSWVER="CMSSW_9_4_6_patch1"
export SCRAM_ARCH="slc6_amd64_gcc630"

cd $cmssw_path
cmsenv

source /afs/cern.ch/cms/cmsset_default.sh
eval `scramv1 runtime -sh`

mkdir -p $output_path/chunk_$1

echo "cmsRun $config_path fileNumber=$1"
cmsRun $config_path fileNumber=$1
