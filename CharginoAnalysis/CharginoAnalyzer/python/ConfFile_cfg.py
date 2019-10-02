import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("CharginoAnalyzer")

options = VarParsing.VarParsing ()

options.register('fileNumber',
                 -1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "File number")
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("CharginoAnalysis.CharginoAnalyzer.CfiFile_cfi")

from CharginoAnalysis.CharginoAnalyzer.CfiFile_cfi import *

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#process.source = cms.Source("PoolSource",
#          fileNames = cms.untracked.vstring('root://eoscms.cern.ch///eos/cms/store/group/phys_susy/xtracks/V3/GEN-SIM-RAW-RECO/chargino300GeV_ctau10cm_GEN-SIM-RAW-RECO_{}.root'.format(options.fileNumber)),
#          secondaryFileNames=cms.untracked.vstring('root://eoscms.cern.ch///eos/cms/store/group/phys_susy/xtracks/V3/GEN-SIM-RAW/chargino300GeV_ctau10cm_GEN-SIM-RAW_{}.root'.format(options.fileNumber))
#                            )

# process.source = cms.Source("PoolSource",
#           fileNames = cms.untracked.vstring('root://eoscms.cern.ch///eos/cms/store/cmst3/user/avartak/XTracks/WJets_HT_200_400/SUS-RunIIFall17MiniAODv2-00077.root'),
#           secondaryFileNames=cms.untracked.vstring('root://eoscms.cern.ch///eos/cms/store/cmst3/user/avartak/XTracks/WJets_HT_200_400/SUS-RunIIFall17DRPremix-00082.root')
#                             )

process.source = cms.Source("PoolSource",
         fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/user/jniedzie/WToLNuJets_HT_200_400/crab_pickEvents/191001_121105/0000/pickevents_{}.root'.format(options.fileNumber)),
         secondaryFileNames=cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/user/jniedzie/WToLNuJets_HT_200_400/crab_pickEvents/191001_103329/0000/pickevents_{}.root'.format(options.fileNumber))
                           )

# process.source = cms.Source("PoolSource",
#          fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov//eos/cms/store/group/phys_exotica/xtracks/taggerStudy/pickedEvents/miniAOD/WToLNuJets_HT_200_400/crab_pickEvents_miniAOD/191002_080910/0000/pickevents_{}.root'.format(options.fileNumber)),
#          secondaryFileNames=cms.untracked.vstring('root://cmsxrootd.fnal.gov//eos/cms/store/group/phys_exotica/xtracks/taggerStudy/pickedEvents/HLTDIGI/WToLNuJets_HT_200_400/crab_pickEvents_HLTDIGI/191002_081009/0000/pickevents_{}.root'.format(options.fileNumber))
#                            )

process.TFileService = cms.Service("TFileService",
                                  fileName = cms.string("/afs/cern.ch/work/j/jniedzie/private/disapp_tracks/pionBackground/friendInfo/chunk_{}/tree_friend.root".format(options.fileNumber)))

#process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string("/afs/cern.ch/work/j/jniedzie/private/disapp_tracks/pionBackground/friendInfo/tree_friend.root"))


process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_ReReco_EOY17_v1', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v6', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v14', '')

process.p = cms.Path(process.CharginoAnalyzer)
