#!/usr/bin/env python
from CRABClient.UserUtilities import config

jobTag = "OO_pythia_30k_ChargedDown"
config = config()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = jobTag
config.General.workArea = 'crab_projects'

config.section_('JobType')
config.JobType.psetName = 'forest_miniAOD_run3_MC.py'

config.JobType.pluginName = 'Analysis'
# ppRef_trigger_experiments.py
# config.JobType.outputFiles = ['Run2_pp_data_ijet.root']
# config.JobType.maxMemoryMB = 2000
config.JobType.maxMemoryMB = 3000

config.section_('Data')
#embedded sample below
config.Data.inputDataset = '/QCD-dijet_pThat15-event-weighted_TuneCP5_5p36TeV_pythia8/HINOOSpring25MiniAOD-CustomTrack_150X_mcRun3_2025_forOO_realistic_v9-v2/MINIAODSIM'
#non-embedded sample below
# config.Data.inputDataset = '/Dijet_pThat-15to1200_TuneCP5_5p36TeV_pythia8/HINOOSpring25MiniAOD-NoPU_150X_mcRun3_2025_forOO_realistic_v9-v1/MINIAODSIM'
config.Data.publication = False
# config.Data.runRange = '306773-306793'
# config.Data.totalUnits = -1
# config.Data.splitting = 'FileBased'
# config.Data.unitsPerJob = 1
config.Data.totalUnits = -1
config.Data.unitsPerJob = 30000
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/vavladim/' + config.General.requestName
# config.Data.outLFNDirBase ='/store/user/vavladim/' + config.General.requestName
# config.Data.outLFNDirBase = '/store/user/lcunquei/Run2_pp_data_CMT'
# config.Data.lumiMask = '/afs/cern.ch/user/v/vavladim/public/recovery_OO_setup/CMSSW_15_0_11/src/HeavyIonsAnalysis/Configuration/test/Cert_Collisions2025OO_394153_394217_golden.json'
config.section_('Site')
# config.Site.storageSite = 'T2_IT_Rome'
config.Site.storageSite = 'T2_CH_CERN'
config.Site.blacklist        = ['T3_UK_ScotGrid_GLA','T2_DE_DESY']
print(config)



# (Uncomment if you need extra JDL flags)
# config.section_('Debug')
# config.Debug.extraJDL       = ['+CMS_ALLOW_OVERFLOW=False']