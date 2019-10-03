from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'friendTreesProducer'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/j/jniedzie/private/disapp_tracks/friendTreeProducer/CMSSW_9_4_11/src/CharginoAnalysis/CharginoAnalyzer/python/ConfFile_cfg.py'

config.Data.inputDataset = '/WToLNuJets_HT_200_400/jniedzie-crab_pickEvents_miniAOD-18783c0a07109245951450a1a4f55409/USER'
config.Data.secondaryInputDataset = '/WToLNuJets_HT_200_400/jniedzie-crab_pickEvents_HLTDIGI-424e4485a07f26f554e82f829d793003/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/group/phys_exotica/xtracks/taggerStudy/friend_info'
config.Data.publication = False
config.Data.outputDatasetTag = 'friendTreesProducer'

config.Site.storageSite = 'T2_CH_CSCS'