from CRABClient.UserUtilities import config
config = config()
config.General.transferOutputs = True
config.General.requestName = 'test-crab-job4'

config.JobType.psetName = 'experimental_forest_miniAOD_run2_MC.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['crab_test_out_PbPb.root']


config.Data.inputDataset = '/DiJet_pThat-15_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/HINPbPbSpring21MiniAOD-FixL1CaloGT_New_Release_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM'
config.Data.inputBlocks = ['/DiJet_pThat-15_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/HINPbPbSpring21MiniAOD-FixL1CaloGT_New_Release_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM#82022616-3127-422e-a0e0-24ba28e93eab'] 

config.Data.publication = False
# config.Data.runRange = '306773-306793'
# config.Data.totalUnits = -1
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.JobType.maxMemoryMB = 2000

config.section_("Site")
config.section_("User")

config.Data.outLFNDirBase = '/store/user/vavladim/'
# /store/user/vladimiv/test_crab_jobs/'
# config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/5TeV/ReReco/Cert_306546-306826_5TeV_EOY2017ReReco_Collisions17_JSON.txt'
config.Site.storageSite = 'T3_CH_CERNBOX'