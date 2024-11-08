from CRABClient.UserUtilities import config
jobTag = "PbPb_data_correct_jetID_HLTL1_v1Dataset"
config = config()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = jobTag

config.section_('JobType')
# config.JobType.psetName = '/afs/cern.ch/user/v/vavladim/private/CMSSW_12_5_0/src/HeavyIonsAnalysis/Configuration/test/data_run2_forest.py'
config.JobType.psetName = '/afs/cern.ch/user/v/vavladim/private/CMSSW_12_5_0/src/HeavyIonsAnalysis/Configuration/test/PbPb_data_trigger_experiments.py'
config.JobType.pluginName = 'Analysis'
# config.JobType.outputFiles = ['Run2_PbPb_data_ijet.root']
# config.JobType.maxMemoryMB = 2000
config.JobType.maxMemoryMB = 5000


config.section_('Data')
config.Data.inputDataset = '/HIHardProbes/HIRun2018A-PbPb18_MiniAODv1-v1/MINIAOD'
# config.Data.inputDataset = '/HIHardProbes/HIRun2018A-PbPb18_MiniAODv1-v2/MINIAOD' #<- V2 of dataset!
# config.Data.inputBlocks = ['/DiJet_pThat-15_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/HINPbPbSpring21MiniAOD-FixL1CaloGT_New_Release_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM#82022616-3127-422e-a0e0-24ba28e93eab'] 
config.Data.publication = False
# config.Data.runRange = '306773-306793'
# config.Data.totalUnits = -1
# config.Data.splitting = 'FileBased'
# config.Data.unitsPerJob = 1
config.Data.totalUnits = -1
config.Data.unitsPerJob = 30000
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/vavladim/' + config.General.requestName
# config.Data.outLFNDirBase = '/store/user/vavladim/' + config.General.requestName
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/PromptReco/Cert_326381-327564_HI_PromptReco_Collisions18_JSON.txt'

config.section_('Site')
# config.Site.storageSite = 'T2_IT_Rome'
# config.Site.blacklist = ['T2_US_Vanderbilt']
config.Site.storageSite = 'T2_CH_CERN'
print(config)
