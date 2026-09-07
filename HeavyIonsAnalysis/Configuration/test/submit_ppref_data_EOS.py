from CRABClient.UserUtilities import config

jobTag = "ppref_data_HLT_L1objects_energyFractions"
config = config()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = jobTag

config.section_('JobType')
config.JobType.psetName = '/afs/cern.ch/user/v/vavladim/private/CMSSW_12_5_0/src/HeavyIonsAnalysis/Configuration/test/ppref_data.py'
# config.JobType.psetName = '/afs/cern.ch/user/v/vavladim/private/CMSSW_12_5_0/src/HeavyIonsAnalysis/Configuration/test/ppRef_trigger_experiments.py'
config.JobType.pluginName = 'Analysis'
# ppRef_trigger_experiments.py
# config.JobType.outputFiles = ['Run2_pp_data_ijet.root']
# config.JobType.maxMemoryMB = 2000
config.JobType.maxMemoryMB = 5000

config.section_('Data')
config.Data.inputDataset = '/HighEGJet/Run2017G-UL2017_MiniAODv2-v1/MINIAOD'
config.Data.publication = False
# config.Data.runRange = '306773-306793'
# config.Data.totalUnits = -1
# config.Data.splitting = 'FileBased'
# config.Data.unitsPerJob = 1
config.Data.totalUnits = -1
config.Data.unitsPerJob = 50000
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/vavladim/' + config.General.requestName
# config.Data.outLFNDirBase ='/store/user/vavladim/' + config.General.requestName
# config.Data.outLFNDirBase = '/store/user/lcunquei/Run2_pp_data_CMT'
# config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/5TeV/ReReco/Cert_306546-306826_5TeV_EOY2017ReReco_Collisions17_JSON.txt'

config.section_('Site')
# config.Site.storageSite = 'T2_IT_Rome'
config.Site.storageSite = 'T2_CH_CERN'
print(config)
