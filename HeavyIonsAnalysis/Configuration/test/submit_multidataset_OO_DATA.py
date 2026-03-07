#!/usr/bin/env python
from CRABClient.UserUtilities import config

jobTag = "OO_data_filtered"
config = config()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = jobTag
config.General.workArea = 'crab_projects'

config.section_('JobType')
config.JobType.psetName = 'forest_miniAOD_run3_DATA.py'

config.JobType.pluginName = 'Analysis'
# ppRef_trigger_experiments.py
# config.JobType.outputFiles = ['Run2_pp_data_ijet.root']
# config.JobType.maxMemoryMB = 2000
config.JobType.maxMemoryMB = 5000

config.section_('Data')
config.Data.inputDataset = '/IonPhysics0/OORun2025-PromptReco-v1/MINIAOD'
# config.Data.inputDataset = '/HIMinimumBias0/HIRun2024A-PromptReco-v1/MINIAOD'
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
# config.Data.lumiMask = '/afs/cern.ch/user/v/vavladim/public/recovery_OO_setup/CMSSW_15_0_11/src/HeavyIonsAnalysis/Configuration/test/Cert_Collisions2025OO_394153_394217_golden.json'
config.section_('Site')
# config.Site.storageSite = 'T2_IT_Rome'
config.Site.storageSite = 'T2_CH_CERN'
config.Site.blacklist        = ['T3_UK_ScotGrid_GLA','T2_DE_DESY']
print(config)



# (Uncomment if you need extra JDL flags)
# config.section_('Debug')
# config.Debug.extraJDL       = ['+CMS_ALLOW_OVERFLOW=False']

# if __name__ == '__main__':

#     from CRABAPI.RawCommand import crabCommand
    
#     dataset_list_2024PbPb = [
#     "/IonPhysics0/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics1/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics10/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics11/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics12/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics13/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics14/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics15/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics16/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics17/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics18/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics19/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics2/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics20/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics21/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics22/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics23/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics24/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics25/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics26/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics27/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics28/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics29/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics3/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics30/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics31/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics32/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics33/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics34/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics35/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics36/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics37/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics38/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics39/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics4/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics40/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics41/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics42/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics43/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics44/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics45/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics46/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics47/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics48/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics49/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics5/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics50/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics51/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics52/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics53/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics54/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics55/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics56/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics57/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics58/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics59/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics6/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics7/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics8/OORun2025-PromptReco-v1/MINIAOD",
#     "/IonPhysics9/OORun2025-PromptReco-v1/MINIAOD"
#     ]
#     for dataset in dataset_list_2024PbPb:
#         config.Data.inputDataset = dataset
#         config.General.requestName = 'OO_' + dataset.split('/')[1]
#         crabCommand('submit', config = config)
