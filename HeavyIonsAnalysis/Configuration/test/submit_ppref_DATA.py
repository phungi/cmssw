from CRABClient.UserUtilities import config
config = config()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'Run2_pp_data_charged_qg'

config.section_('JobType')
config.JobType.psetName = '/afs/cern.ch/user/v/vavladim/private/CMSSW_12_5_0/src/HeavyIonsAnalysis/Configuration/test/ppref_data.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['Run2_pp_data_charged_qg.root']
# config.JobType.maxMemoryMB = 2000
config.JobType.maxMemoryMB = 3000


config.section_('Data')
config.Data.inputDataset = '/HighEGJet/Run2017G-UL2017_MiniAODv2-v1/MINIAOD'
#do it in blocks because apparently there is some memory leak somewhere (in charged or split matching) and for 40k+ events it exceeds 5GB RAM...
# config.Data.inputBlocks = ['/HighEGJet/Run2017G-UL2017_MiniAODv2-v1/MINIAOD#15ec4129-61e4-4a46-b85c-bda61009d068',
# 							'/HighEGJet/Run2017G-UL2017_MiniAODv2-v1/MINIAOD#170048c3-1a55-43a4-8e51-f5c22e421eb7',
# 							'/HighEGJet/Run2017G-UL2017_MiniAODv2-v1/MINIAOD#25cf044f-7a0d-4198-b06d-fccb462e018f',
# 							'/HighEGJet/Run2017G-UL2017_MiniAODv2-v1/MINIAOD#273e0a65-49bd-4bff-98ef-398342240b2b',
# 							'/HighEGJet/Run2017G-UL2017_MiniAODv2-v1/MINIAOD#48720ec6-cee3-4ad0-8e41-00452783c364',
# 							'/HighEGJet/Run2017G-UL2017_MiniAODv2-v1/MINIAOD#8cc2f72a-93e6-4d2c-9c5d-d49b99eacec9'] 
# config.Data.inputBlocks = ['/HighEGJet/Run2017G-UL2017_MiniAODv2-v1/MINIAOD#98993408-6630-4ba1-b3e7-79d646e237f0',
# 							'/HighEGJet/Run2017G-UL2017_MiniAODv2-v1/MINIAOD#a6c2d147-c06a-4118-b5ac-fdd8cbe925d7',
# 							'/HighEGJet/Run2017G-UL2017_MiniAODv2-v1/MINIAOD#aa57d15a-e4f9-4a93-af9f-0cbac474a690',
# 							'/HighEGJet/Run2017G-UL2017_MiniAODv2-v1/MINIAOD#aab36a39-3ebe-4684-9a2c-1dc4468ee147',
# 							'/HighEGJet/Run2017G-UL2017_MiniAODv2-v1/MINIAOD#b8107fae-9af0-4592-9e8a-11b6f804cac2',
# 							'/HighEGJet/Run2017G-UL2017_MiniAODv2-v1/MINIAOD#bee5b299-695d-4407-ad60-75884e7e0c29',
# 							'/HighEGJet/Run2017G-UL2017_MiniAODv2-v1/MINIAOD#c0df5069-f281-403a-9796-58e366ebc7e5',
# 							'/HighEGJet/Run2017G-UL2017_MiniAODv2-v1/MINIAOD#c2832573-6b24-410d-8d9b-950b87e9b55d'
# 							]
config.Data.publication = False
# config.Data.runRange = '306773-306793'
# config.Data.totalUnits = -1
# config.Data.splitting = 'FileBased'
# config.Data.unitsPerJob = 1
config.Data.totalUnits = -1
config.Data.unitsPerJob = 50000
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/vavladim/Run2_pp_data_charged_qg'
# config.Data.outLFNDirBase = '/store/user/lcunquei/Run2_pp_data_CMT'
# config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/5TeV/ReReco/Cert_306546-306826_5TeV_EOY2017ReReco_Collisions17_JSON.txt'

config.section_('Site')
# config.Site.storageSite = 'T2_IT_Rome'
config.Site.storageSite = 'T2_CH_CERN'
print(config)
