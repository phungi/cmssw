from CRABClient.UserUtilities import config
jobTag = "pp_pythia_jetID_energyFraction_R4_NeutralCompensatory"
config = config()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = jobTag

config.section_('JobType')
config.JobType.psetName = '/afs/cern.ch/user/v/vavladim/private/CMSSW_12_5_0/src/HeavyIonsAnalysis/Configuration/test/ppref_MC.py'
config.JobType.pluginName = 'Analysis'
# config.JobType.outputFiles = [config.General.requestName + '.root']
# config.JobType.maxMemoryMB = 2000
config.JobType.maxMemoryMB = 5000


config.section_('Data')

config.Data.inputDataset = '/QCD_pThat-15_Dijet_TuneCP5_5p02TeV_pythia8/RunIIpp5Spring18MiniAOD-94X_mc2017_realistic_forppRef5TeV-v1/MINIAODSIM'
# Official herwig sample from Jelena
# config.Data.inputDataset = '/QCD_PtGT15_TuneCH3_5p02TeV_herwig7/RunIIpp5Spring18MiniAOD-94X_mc2017_realistic_forppRef5TeV-v2/MINIAODSIM'
# Lida's private sample below
# config.Data.inputDataset = '/Herwig_CH3_qcd_5TeV/mnguyen-Herwig_CH3_qcd_5TeV_MINI_v7-0574b4ed1d118382056a950d3f8bdf1b/USER'
# config.Data.inputBlocks = []
# config.Data.inputDBS = "phys03"
config.Data.publication = False
# config.Data.runRange = '306773-306793'
# config.Data.totalUnits = -1
# config.Data.splitting = 'FileBased'
# config.Data.unitsPerJob = 1
config.Data.totalUnits = -1
config.Data.unitsPerJob = 30000
config.Data.splitting = 'EventAwareLumiBased'
# config.Data.outLFNDirBase = '/store/user/vavladim/' + config.General.requestName
config.Data.outLFNDirBase = '/store/group/phys_heavyions/vavladim/' + config.General.requestName
# config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/5TeV/ReReco/Cert_306546-306826_5TeV_EOY2017ReReco_Collisions17_JSON.txt'

config.section_('Site')
# config.Site.storageSite = 'T2_IT_Rome'
config.Site.storageSite = 'T2_CH_CERN'
print(config)
