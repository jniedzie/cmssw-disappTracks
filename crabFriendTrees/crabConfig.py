from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'friendTreesProducer'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/j/jniedzie/private/disapp_tracks/friendTreeProducer/CMSSW_10_2_4_patch1/src/CharginoAnalysis/CharginoAnalyzer/python/ConfFile_cfg.py'

config.Data.inputDataset = '/MET/jniedzie-crab_pickEvents_MET_2018A_RAW-RECO-e5f6624f5ab7506b38712cf6e55471cc/USER'
config.Data.secondaryInputDataset = '/MET/jniedzie-crab_pickEvents_MET_2018A_miniAOD-e5f6624f5ab7506b38712cf6e55471cc/USER'


config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/group/phys_exotica/xtracks/taggerStudy/pickedEvents/MET_2018A/friendTrees'
config.Data.outputDatasetTag = 'friendTreesProducer'

config.Site.storageSite = 'T2_CH_CERN'