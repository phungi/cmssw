#!/usr/bin/env python
from CRABClient.UserUtilities import config

jobTag = "2024_ppRef_ZB"
config = config()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = jobTag
config.General.workArea = 'crab_projects'

config.section_('JobType')
config.JobType.psetName = 'forest_miniAOD_run3_ppref_DATA.py'

config.JobType.pluginName = 'Analysis'
# ppRef_trigger_experiments.py
# config.JobType.outputFiles = ['Run2_pp_data_ijet.root']
# config.JobType.maxMemoryMB = 2000
config.JobType.maxMemoryMB = 3000

config.section_('Data')
# config.Data.inputDataset = '/IonPhysics0/OORun2025-PromptReco-v1/MINIAOD'
# config.Data.inputDataset = '/HIMinimumBias0/HIRun2024A-PromptReco-v1/MINIAOD'
config.Data.inputDataset = '/PPRefZeroBiasPlusForward0/Run2024J-PromptReco-v1/MINIAOD'

config.Data.publication = False
# config.Data.runRange = '306773-306793'
# config.Data.totalUnits = -1
# config.Data.splitting = 'FileBased'
# config.Data.unitsPerJob = 1
config.Data.totalUnits = -1
config.Data.unitsPerJob = 200000
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/vavladim/' + config.General.requestName
# config.Data.outLFNDirBase ='/store/user/vavladim/' + config.General.requestName
# config.Data.outLFNDirBase = '/store/user/lcunquei/Run2_pp_data_CMT'
config.Data.lumiMask = '/afs/cern.ch/user/v/vavladim/public/OO_setup/CMSSW_15_0_11/src/HeavyIonsAnalysis/Configuration/test/Cert_Collisions2024_ppref_387474_387721_golden.json'
config.section_('Site')
# config.Site.storageSite = 'T2_IT_Rome'
config.Site.storageSite = 'T2_CH_CERN'
config.Site.blacklist        = ['T3_UK_ScotGrid_GLA','T2_DE_DESY']
print(config)



# (Uncomment if you need extra JDL flags)
# config.section_('Debug')
# config.Debug.extraJDL       = ['+CMS_ALLOW_OVERFLOW=False']

if __name__ == '__main__':
    from CRABAPI.RawCommand import crabCommand
    dataset_list_2024Ppp = [
    '/PPRefHardProbes0/Run2024J-PromptReco-v1/MINIAOD',
	'/PPRefHardProbes1/Run2024J-PromptReco-v1/MINIAOD',
	'/PPRefHardProbes2/Run2024J-PromptReco-v1/MINIAOD',
	'/PPRefHardProbes3/Run2024J-PromptReco-v1/MINIAOD',
	'/PPRefHardProbes4/Run2024J-PromptReco-v1/MINIAOD',
    "/PPRefZeroBiasPlusForward0/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward1/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward10/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward11/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward12/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward13/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward14/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward15/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward16/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward17/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward18/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward19/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward2/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward20/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward21/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward22/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward23/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward24/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward3/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward4/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward5/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward6/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward7/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward8/Run2024J-PromptReco-v1/MINIAOD",
    "/PPRefZeroBiasPlusForward9/Run2024J-PromptReco-v1/MINIAOD"
    ]
    for dataset in dataset_list_2024Ppp:
        config.Data.inputDataset = dataset
        config.General.requestName = '2024_filtered_' + dataset.split('/')[1]
        crabCommand('submit', config = config)
