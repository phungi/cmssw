### HiForest Configuration
# Collisions: pp
# Type: data
# Input: miniAOD

import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process('HiForest', Run3)
process.options = cms.untracked.PSet()

#####################################################################################
# HiForest labelling info
#####################################################################################

process.load("HeavyIonsAnalysis.EventAnalysis.HiForestInfo_cfi")
process.HiForestInfo.info = cms.vstring("HiForest, miniAOD, 125X, ak4PFJetSequence_ppref_data_cff")

#####################################################################################
# Input source
#####################################################################################

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
        # '/HighEGJet/Run2017G-UL2017_MiniAODv2-v1/MINIAOD'
        '/store/data/Run2017G/HighEGJet/MINIAOD/UL2017_MiniAODv2-v1/00000/0555A405-79C7-0744-A78E-61B94060F6C7.root'
    )
)

# Number of events we want to process, -1 = all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

#####################################################################################
# Load Global Tag, Geometry, etc.
#####################################################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

# TODO: Global tag complete guess from the list. Probably wrong. But does not crash
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '125X_mcRun3_2022_realistic_v3', '')
process.HiForestInfo.GlobalTagLabel = process.GlobalTag.globaltag

# TODO: Old calibration here, might need to update
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
             tag = cms.string("JPcalib_MC94X_2017pp_v2"),
             connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")

         )
      ])

#####################################################################################
# Define tree output
#####################################################################################

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("Run2_pp_data_fix.root"))

#####################################################################################
# Additional Reconstruction and Analysis: Main Body
#####################################################################################

#############################
# Jets
#############################
process.load("HeavyIonsAnalysis.JetAnalysis.ak4PFJetSequence_ppref_data_cff")
#####################################################################################

############################
# Event Analysis
############################
# use data version to avoid PbPb MC
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.hiEvtAnalyzer.Vertex = cms.InputTag("offlineSlimmedPrimaryVertices")
process.hiEvtAnalyzer.doCentrality = cms.bool(False)
process.hiEvtAnalyzer.doEvtPlane = cms.bool(False)
process.hiEvtAnalyzer.doEvtPlaneFlat = cms.bool(False)
process.hiEvtAnalyzer.doMC = cms.bool(True) # general MC info
process.hiEvtAnalyzer.doHiMC = cms.bool(False) # HI specific MC info
process.hiEvtAnalyzer.doHFfilters = cms.bool(False) # Disable HF filters for ppRef

process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.l1object_cfi')

# TODO: Many of these triggers are not available in the test file
from HeavyIonsAnalysis.EventAnalysis.hltobject_cfi import trigger_list_data
process.hltobject.triggerNames = trigger_list_data

# Gen particles
process.load('HeavyIonsAnalysis.EventAnalysis.HiGenAnalyzer_cfi')

#####################################################################################

#########################
# Track Analyzer
#########################
# Tracks cause an error due to missing product: packedPFCandidateTrackChi2
#process.load('HeavyIonsAnalysis.TrackAnalysis.TrackAnalyzers_cff')

#####################################################################################

#####################
# photons
######################
process.load('HeavyIonsAnalysis.EGMAnalysis.ggHiNtuplizer_cfi')
process.ggHiNtuplizer.doGenParticles = cms.bool(True)
process.ggHiNtuplizer.doMuons = cms.bool(False) # unpackedMuons collection not found from file
process.ggHiNtuplizer.useValMapIso = cms.bool(False) # True here causes seg fault
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

####################################################################################

#########################
# Main analysis list
#########################

process.forest = cms.Path(
    process.hltobject +
    process.HiForestInfo +
    process.hltanalysis +
    process.hiEvtAnalyzer    
    # process.l1object +
    # process.HiGenParticleAna +
    # process.ggHiNtuplizer #+
#    process.trackSequencePP
)

#####################################################################################


addR4Jets = True
addR2Jets = True

if addR2Jets or addR4Jets :
    process.load("HeavyIonsAnalysis.JetAnalysis.extraJets_cff")
    from HeavyIonsAnalysis.JetAnalysis.clusterJetsFromMiniAOD_cff import setupPprefJets

    if addR2Jets :
        process.jetsR2 = cms.Sequence()
        setupPprefJets('ak2PF', process.jetsR2, process, isMC = 0, radius = 0.20, JECTag = 'AK2PF')
        process.ak2PFpatJetCorrFactors.levels = ['L2Relative', 'L3Absolute']
        process.ak2PFpatJetCorrFactors.primaryVertices = "offlineSlimmedPrimaryVertices"
        process.load("HeavyIonsAnalysis.JetAnalysis.candidateBtaggingMiniAOD_cff")
        process.ak2PFJetAnalyzer = process.ak4PFJetAnalyzer.clone(jetTag = "ak2PFpatJets", jetName = 'ak2PF')
        process.forest += process.extraJetsData * process.jetsR2 * process.ak2PFJetAnalyzer

    if addR4Jets :
        # Recluster using an alias "0" in order not to get mixed up with the default AK4 collections
        process.jetsR4 = cms.Sequence()
        setupPprefJets('ak04PF', process.jetsR4, process, isMC = 0, radius = 0.40, JECTag = 'AK4PF')
        process.ak04PFpatJetCorrFactors.levels = ['L2Relative', 'L3Absolute']
        process.ak04PFpatJetCorrFactors.primaryVertices = "offlineSlimmedPrimaryVertices"
        process.load("HeavyIonsAnalysis.JetAnalysis.candidateBtaggingMiniAOD_cff")
        process.ak4PFJetAnalyzer.jetTag = 'ak04PFpatJets'
        process.ak4PFJetAnalyzer.jetName = 'ak04PF'
        process.ak4PFJetAnalyzer.doSubEvent = False # Need to disable this, since there is some issue with the gen jet constituents. More debugging needed is want to use constituents.
        process.forest += process.extraJetsData * process.jetsR4 * process.ak4PFJetAnalyzer
