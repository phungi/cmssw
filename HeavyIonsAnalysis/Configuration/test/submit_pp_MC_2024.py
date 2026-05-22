#!/usr/bin/env python
from CRABClient.UserUtilities import config

jobTag = "2024_pp_pythia_2GeV_constit_cut"
config = config()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = jobTag
config.General.workArea = 'crab_projects'

config.section_('JobType')
config.JobType.psetName = 'forest_miniAOD_run3_ppref_MC.py'

config.JobType.pluginName = 'Analysis'
# ppRef_trigger_experiments.py
# config.JobType.outputFiles = ['Run2_pp_data_ijet.root']
# config.JobType.maxMemoryMB = 2000
config.JobType.maxMemoryMB = 3000

config.section_('Data')
# default sample for pp
config.Data.inputDataset = '/QCD_pThat-15to1200_TuneCP5_5p36TeV_pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v2/MINIAODSIM'
# herwig7 sample
# config.Data.inputDataset = '/QCD_Pt-15to1200_TuneCH3_Flat_5p36TeV_herwig7/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v2/MINIAODSIM'
config.Data.publication = False
# config.Data.runRange = '306773-306793'
# config.Data.totalUnits = -1
# config.Data.splitting = 'FileBased'
# config.Data.unitsPerJob = 1
config.Data.totalUnits = -1
config.Data.unitsPerJob = 100000
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/vavladim/' + config.General.requestName
# config.Data.outLFNDirBase ='/store/user/vavladim/' + config.General.requestName
# config.Data.outLFNDirBase = '/store/user/lcunquei/Run2_pp_data_CMT'
# config.Data.lumiMask = '/afs/cern.ch/user/v/vavladim/public/recovery_OO_setup/CMSSW_15_0_11/src/HeavyIonsAnalysis/Configuration/test/Cert_Collisions2025OO_394153_394217_golden.json'
config.section_('Site')
# config.Site.storageSite = 'T2_IT_Rome'
config.Site.storageSite = 'T2_CH_CERN'
# config.Site.blacklist        = ['T3_UK_ScotGrid_GLA','T2_DE_DESY']
print(config)



# (Uncomment if you need extra JDL flags)
# config.section_('Debug')
# config.Debug.extraJDL       = ['+CMS_ALLOW_OVERFLOW=False']