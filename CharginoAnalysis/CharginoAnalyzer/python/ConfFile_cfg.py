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

# process.source = cms.Source("PoolSource",
#          fileNames = cms.untracked.vstring('root://eoscms.cern.ch///eos/cms/store/group/phys_susy/xtracks/V3/GEN-SIM-RAW-RECO/chargino300GeV_ctau10cm_GEN-SIM-RAW-RECO_{}.root'.format(options.fileNumber)),
#          secondaryFileNames=cms.untracked.vstring('root://eoscms.cern.ch///eos/cms/store/group/phys_susy/xtracks/V3/GEN-SIM-RAW/chargino300GeV_ctau10cm_GEN-SIM-RAW_{}.root'.format(options.fileNumber))
#                            )

# process.source = cms.Source("PoolSource",
#           fileNames = cms.untracked.vstring('root://eoscms.cern.ch///eos/cms/store/cmst3/user/avartak/XTracks/WJets_HT_200_400/SUS-RunIIFall17MiniAODv2-00077.root'),
#           secondaryFileNames=cms.untracked.vstring('root://eoscms.cern.ch///eos/cms/store/cmst3/user/avartak/XTracks/WJets_HT_200_400/SUS-RunIIFall17DRPremix-00082.root')
#                             )

# process.source = cms.Source("PoolSource",
#          fileNames = cms.untracked.vstring('root://eoscms.cern.ch//eos/cms/store/group/phys_exotica/xtracks/taggerStudy/pickedEvents/miniAOD/WToLNuJets_HT_200_400/crab_pickEvents_miniAOD/191004_114459/0000/pickevents_{}.root'.format(options.fileNumber)),
#          secondaryFileNames=cms.untracked.vstring('root://eoscms.cern.ch//eos/cms/store/group/phys_exotica/xtracks/taggerStudy/pickedEvents/RECO/WToLNuJets_HT_200_400/crab_pickEvents_RECO/191004_140823/0000/pickevents_{}.root'.format(options.fileNumber))
#          # secondaryFileNames=cms.untracked.vstring('root://eoscms.cern.ch//eos/cms/store/group/phys_exotica/xtracks/taggerStudy/pickedEvents/HLTDIGI/WToLNuJets_HT_200_400/crab_pickEvents_HLTDIGI/191004_120843/0000/pickevents_{}.root'.format(options.fileNumber))
#                            )

# process.source = cms.Source("PoolSource",
#          fileNames = cms.untracked.vstring('root://eoscms.cern.ch//eos/cms/store/group/phys_exotica/xtracks/taggerStudy/pickedEvents/miniAODnoPU/WToLNuJets_HT_200_400/crab_pickEvents_miniAOD_noPU/191008_071346/0000/pickevents_{}.root'.format(options.fileNumber)),
#          secondaryFileNames=cms.untracked.vstring('root://eoscms.cern.ch//eos/cms/store/group/phys_exotica/xtracks/taggerStudy/pickedEvents/RECOnoPU/WToLNuJets_HT_200_400/crab_pickEvents_RECO_noPU/191008_071546/0000/pickevents_{}.root'.format(options.fileNumber))
#                            )

# process.source = cms.Source("PoolSource",
#          fileNames=cms.untracked.vstring('root://eoscms.cern.ch//eos/cms/store/group/phys_exotica/xtracks/taggerStudy/pickedEvents/MET_2018A/RAW-RECO/MET/crab_pickEvents_MET_2018A_RAW-RECO/191015_114245/0000/pickevents_{}.root'.format(options.fileNumber)),
#          secondaryFileNames = cms.untracked.vstring('root://eoscms.cern.ch//eos/cms/store/group/phys_exotica/xtracks/taggerStudy/pickedEvents/MET_2018A/miniAOD/MET/crab_pickEvents_MET_2018A_miniAOD/191015_121144/0000/pickevents_{}.root'.format(options.fileNumber))
#                            )

# process.source = cms.Source("PoolSource",
#          fileNames = cms.untracked.vstring('root://eoscms.cern.ch///eos/cms/store/group/phys_susy/xtracks/V3/GEN-SIM-RAW-RECO-PU/GEN-SIM-RAW-RECO-PU/chargino300GeV_ctau10cm_GEN-SIM-RAW-RECO_{}.root'.format(options.fileNumber)),
#          secondaryFileNames=cms.untracked.vstring('root://eoscms.cern.ch///eos/cms/store/group/phys_susy/xtracks/V3/GEN-SIM-RAW-PU/GEN-SIM-RAW-PU/chargino300GeV_ctau10cm_GEN-SIM-RAW_{}.root'.format(options.fileNumber))
#                            )

process.source = cms.Source("PoolSource",
         fileNames = cms.untracked.vstring('root://eoscms.cern.ch///eos/cms/store/group/phys_susy/xtracks/300GeV1cm/GEN-SIM-RAW-RECO-PU/chargino300GeV_ctau1cm_GEN-SIM-RAW-RECO_{}.root'.format(options.fileNumber)),
         secondaryFileNames=cms.untracked.vstring('root://eoscms.cern.ch///eos/cms/store/group/phys_susy/xtracks/300GeV1cm/GEN-SIM-RAW-PU/chargino300GeV_ctau1cm_GEN-SIM-RAW_{}.root'.format(options.fileNumber))
                           )

# process.source = cms.Source("PoolSource",
#          fileNames = cms.untracked.vstring('root://eoscms.cern.ch//eos/cms/store/group/phys_exotica/xtracks/taggerStudy/pickedEvents/miniAOD-ext/WToLNuJets_HT_200_400_ext/crab_pickEvents_miniAOD-ext/191202_102433/0000/pickevents_{}.root'.format(options.fileNumber)),
#          secondaryFileNames=cms.untracked.vstring('root://eoscms.cern.ch//eos/cms/store/group/phys_exotica/xtracks/taggerStudy/pickedEvents/RECO-ext/WToLNuJets_HT_200_400_ext/crab_pickEvents_RECO-ext/191202_103647/0000/pickevents_{}.root'.format(options.fileNumber))
#                            )

# process.source = cms.Source("PoolSource",
#          fileNames = cms.untracked.vstring('root://eoscms.cern.ch//eos/cms/store/group/phys_exotica/xtracks/private_reco/test.reco.aod.root'),
#          secondaryFileNames=cms.untracked.vstring('root://eoscms.cern.ch//eos/cms/store/group/phys_exotica/xtracks/private_reco/test.reco.root')
#                            )

process.TFileService = cms.Service("TFileService",
                                  fileName = cms.string("/eos/cms/store/group/phys_exotica/xtracks/7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/Wino_300GeV1cm/treeProducerXtracks/friend_trees/chunk_{}/tree_friend.root".format(options.fileNumber)))


# process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string("/eos/cms/store/group/phys_exotica/xtracks/7Sep2019/Calibrated-MC-ext-SR-2018-Hadded/WJetsToLNu_HT200to400_special/treeProducerXtracks/friend_trees_noHighPtHits/chunk_{}/tree_friend.root".format(options.fileNumber)))

# process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string("/eos/cms/store/group/phys_exotica/xtracks/private_reco/friend_trees_noHighPtHits.root"))


# process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string("/afs/cern.ch/work/j/jniedzie/private/disapp_tracks/friendTreeProducer/tree_friend.root"))


process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_ReReco_EOY17_v1', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v6', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v14', '')

process.p = cms.Path(process.CharginoAnalyzer)
