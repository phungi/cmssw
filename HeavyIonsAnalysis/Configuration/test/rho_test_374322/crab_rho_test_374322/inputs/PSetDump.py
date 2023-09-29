import FWCore.ParameterSet.Config as cms
from HeterogeneousCore.CUDACore.SwitchProducerCUDA import SwitchProducerCUDA

process = cms.Process("HiForest")

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring('file:/eos/cms/store/group/phys_heavyions/wangj/RECO2023/miniaod_PhysicsHIPhysicsRawPrime0_374289/reco_run374289_ls0228_streamPhysicsHIPhysicsRawPrime0_StorageManager.root')
)
process.CondDB = cms.PSet(
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        connectionTimeout = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    connect = cms.string('')
)

process.HFRecalParameterBlock = cms.PSet(
    HFdepthOneParameterA = cms.vdouble(
        0.004123, 0.00602, 0.008201, 0.010489, 0.013379,
        0.016997, 0.021464, 0.027371, 0.034195, 0.044807,
        0.058939, 0.125497
    ),
    HFdepthOneParameterB = cms.vdouble(
        -4e-06, -2e-06, 0.0, 4e-06, 1.5e-05,
        2.6e-05, 6.3e-05, 8.4e-05, 0.00016, 0.000107,
        0.000425, 0.000209
    ),
    HFdepthTwoParameterA = cms.vdouble(
        0.002861, 0.004168, 0.0064, 0.008388, 0.011601,
        0.014425, 0.018633, 0.023232, 0.028274, 0.035447,
        0.051579, 0.086593
    ),
    HFdepthTwoParameterB = cms.vdouble(
        -2e-06, -0.0, -7e-06, -6e-06, -2e-06,
        1e-06, 1.9e-05, 3.1e-05, 6.7e-05, 1.2e-05,
        0.000157, -3e-06
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

process.maxLuminosityBlocks = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.nanoDQMIO_perLSoutput = cms.PSet(
    MEsToSave = cms.untracked.vstring(
        'Hcal/DigiTask/Occupancy/depth/depth1',
        'Hcal/DigiTask/Occupancy/depth/depth2',
        'Hcal/DigiTask/Occupancy/depth/depth3',
        'Hcal/DigiTask/Occupancy/depth/depth4',
        'Hcal/DigiTask/Occupancy/depth/depth5',
        'Hcal/DigiTask/Occupancy/depth/depth6',
        'Hcal/DigiTask/Occupancy/depth/depth7',
        'Hcal/DigiTask/Occupancy/depth/depthHO',
        'Hcal/DigiTask/OccupancyCut/depth/depth1',
        'Hcal/DigiTask/OccupancyCut/depth/depth2',
        'Hcal/DigiTask/OccupancyCut/depth/depth3',
        'Hcal/DigiTask/OccupancyCut/depth/depth4',
        'Hcal/DigiTask/OccupancyCut/depth/depth5',
        'Hcal/DigiTask/OccupancyCut/depth/depth6',
        'Hcal/DigiTask/OccupancyCut/depth/depth7',
        'Hcal/DigiTask/OccupancyCut/depth/depthHO',
        'EcalBarrel/EBOccupancyTask/EBOT digi occupancy',
        'EcalEndcap/EEOccupancyTask/EEOT digi occupancy EE -',
        'EcalEndcap/EEOccupancyTask/EEOT digi occupancy EE +',
        'EcalBarrel/EBOccupancyTask/EBOT DCC entries',
        'EcalEndcap/EEOccupancyTask/EEOT DCC entries',
        'Ecal/EventInfo/processedEvents',
        'PixelPhase1/Tracks/charge_PXBarrel',
        'PixelPhase1/Tracks/charge_PXForward',
        'PixelPhase1/Tracks/PXBarrel/charge_PXLayer_1',
        'PixelPhase1/Tracks/PXBarrel/charge_PXLayer_2',
        'PixelPhase1/Tracks/PXBarrel/charge_PXLayer_3',
        'PixelPhase1/Tracks/PXBarrel/charge_PXLayer_4',
        'PixelPhase1/Tracks/PXForward/charge_PXDisk_+1',
        'PixelPhase1/Tracks/PXForward/charge_PXDisk_+2',
        'PixelPhase1/Tracks/PXForward/charge_PXDisk_+3',
        'PixelPhase1/Tracks/PXForward/charge_PXDisk_-1',
        'PixelPhase1/Tracks/PXForward/charge_PXDisk_-2',
        'PixelPhase1/Tracks/PXForward/charge_PXDisk_-3',
        'PixelPhase1/Tracks/PXBarrel/size_PXLayer_1',
        'PixelPhase1/Tracks/PXBarrel/size_PXLayer_2',
        'PixelPhase1/Tracks/PXBarrel/size_PXLayer_3',
        'PixelPhase1/Tracks/PXBarrel/size_PXLayer_4',
        'PixelPhase1/Tracks/PXForward/size_PXDisk_+1',
        'PixelPhase1/Tracks/PXForward/size_PXDisk_+2',
        'PixelPhase1/Tracks/PXForward/size_PXDisk_+3',
        'PixelPhase1/Tracks/PXForward/size_PXDisk_-1',
        'PixelPhase1/Tracks/PXForward/size_PXDisk_-2',
        'PixelPhase1/Tracks/PXForward/size_PXDisk_-3',
        'HLT/Vertexing/hltPixelVertices/hltPixelVertices/goodvtxNbr',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/mon_eta',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/mon_hits',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/mon_phi',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/mon_pt',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/mon_unMatched_eta',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/mon_unMatched_hits',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/mon_unMatched_phi',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/mon_unMatched_pt',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/ref_eta',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/ref_hits',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/ref_matched_eta',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/ref_matched_hits',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/ref_matched_phi',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/ref_matched_pt',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/ref_phi',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/ref_pt',
        'HLT/Tracking/pixelTracks/GeneralProperties/Chi2Prob_GenTk',
        'HLT/Tracking/pixelTracks/GeneralProperties/Chi2oNDFVsEta_ImpactPoint_GenTk',
        'HLT/Tracking/pixelTracks/GeneralProperties/DeltaZToPVZoom_GenTk',
        'HLT/Tracking/pixelTracks/GeneralProperties/DistanceOfClosestApproachToPVVsPhi_GenTk',
        'HLT/Tracking/pixelTracks/GeneralProperties/DistanceOfClosestApproachToPVZoom_GenTk',
        'HLT/Tracking/pixelTracks/GeneralProperties/NumberOfTracks_GenTk',
        'HLT/Tracking/tracks/GeneralProperties/Chi2Prob_GenTk',
        'HLT/Tracking/tracks/GeneralProperties/Chi2oNDFVsEta_ImpactPoint_GenTk',
        'HLT/Tracking/tracks/GeneralProperties/DeltaZToPVZoom_GenTk',
        'HLT/Tracking/tracks/GeneralProperties/DistanceOfClosestApproachToPVVsPhi_GenTk',
        'HLT/Tracking/tracks/GeneralProperties/DistanceOfClosestApproachToPVZoom_GenTk',
        'HLT/Tracking/tracks/GeneralProperties/NumberOfTracks_GenTk',
        'HLT/Tracking/tracks/LUMIanalysis/NumberEventsVsLUMI',
        'HLT/Tracking/tracks/PUmonitoring/NumberEventsVsGoodPVtx',
        'PixelPhase1/Tracks/num_clusters_ontrack_PXBarrel',
        'PixelPhase1/Tracks/num_clusters_ontrack_PXForward',
        'PixelPhase1/Tracks/clusterposition_zphi_ontrack',
        'PixelPhase1/Tracks/PXBarrel/clusterposition_zphi_ontrack_PXLayer_1',
        'PixelPhase1/Tracks/PXBarrel/clusterposition_zphi_ontrack_PXLayer_2',
        'PixelPhase1/Tracks/PXBarrel/clusterposition_zphi_ontrack_PXLayer_3',
        'PixelPhase1/Tracks/PXBarrel/clusterposition_zphi_ontrack_PXLayer_4',
        'PixelPhase1/Tracks/PXForward/clusterposition_xy_ontrack_PXDisk_+1',
        'PixelPhase1/Tracks/PXForward/clusterposition_xy_ontrack_PXDisk_+2',
        'PixelPhase1/Tracks/PXForward/clusterposition_xy_ontrack_PXDisk_+3',
        'PixelPhase1/Tracks/PXForward/clusterposition_xy_ontrack_PXDisk_-1',
        'PixelPhase1/Tracks/PXForward/clusterposition_xy_ontrack_PXDisk_-2',
        'PixelPhase1/Tracks/PXForward/clusterposition_xy_ontrack_PXDisk_-3',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_1/NormalizedHitResiduals_TEC__wheel__1',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_2/NormalizedHitResiduals_TEC__wheel__2',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_3/NormalizedHitResiduals_TEC__wheel__3',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_4/NormalizedHitResiduals_TEC__wheel__4',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_5/NormalizedHitResiduals_TEC__wheel__5',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_6/NormalizedHitResiduals_TEC__wheel__6',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_7/NormalizedHitResiduals_TEC__wheel__7',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_8/NormalizedHitResiduals_TEC__wheel__8',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_9/NormalizedHitResiduals_TEC__wheel__9',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_1/NormalizedHitResiduals_TEC__wheel__1',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_2/NormalizedHitResiduals_TEC__wheel__2',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_3/NormalizedHitResiduals_TEC__wheel__3',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_4/NormalizedHitResiduals_TEC__wheel__4',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_5/NormalizedHitResiduals_TEC__wheel__5',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_6/NormalizedHitResiduals_TEC__wheel__6',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_7/NormalizedHitResiduals_TEC__wheel__7',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_8/NormalizedHitResiduals_TEC__wheel__8',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_9/NormalizedHitResiduals_TEC__wheel__9',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_1/Summary_ClusterStoNCorr__OnTrack__TEC__PLUS__wheel__1',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_2/Summary_ClusterStoNCorr__OnTrack__TEC__PLUS__wheel__2',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_3/Summary_ClusterStoNCorr__OnTrack__TEC__PLUS__wheel__3',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_4/Summary_ClusterStoNCorr__OnTrack__TEC__PLUS__wheel__4',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_5/Summary_ClusterStoNCorr__OnTrack__TEC__PLUS__wheel__5',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_6/Summary_ClusterStoNCorr__OnTrack__TEC__PLUS__wheel__6',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_7/Summary_ClusterStoNCorr__OnTrack__TEC__PLUS__wheel__7',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_8/Summary_ClusterStoNCorr__OnTrack__TEC__PLUS__wheel__8',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_9/Summary_ClusterStoNCorr__OnTrack__TEC__PLUS__wheel__9',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_1/Summary_ClusterStoNCorr__OnTrack__TEC__MINUS__wheel__1',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_2/Summary_ClusterStoNCorr__OnTrack__TEC__MINUS__wheel__2',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_3/Summary_ClusterStoNCorr__OnTrack__TEC__MINUS__wheel__3',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_4/Summary_ClusterStoNCorr__OnTrack__TEC__MINUS__wheel__4',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_5/Summary_ClusterStoNCorr__OnTrack__TEC__MINUS__wheel__5',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_6/Summary_ClusterStoNCorr__OnTrack__TEC__MINUS__wheel__6',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_7/Summary_ClusterStoNCorr__OnTrack__TEC__MINUS__wheel__7',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_8/Summary_ClusterStoNCorr__OnTrack__TEC__MINUS__wheel__8',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_9/Summary_ClusterStoNCorr__OnTrack__TEC__MINUS__wheel__9',
        'SiStrip/MechanicalView/TIB/layer_1/NormalizedHitResiduals_TIB__Layer__1',
        'SiStrip/MechanicalView/TIB/layer_2/NormalizedHitResiduals_TIB__Layer__2',
        'SiStrip/MechanicalView/TIB/layer_3/NormalizedHitResiduals_TIB__Layer__3',
        'SiStrip/MechanicalView/TIB/layer_4/NormalizedHitResiduals_TIB__Layer__4',
        'SiStrip/MechanicalView/TIB/layer_1/Summary_ClusterStoNCorr__OnTrack__TIB__layer__1',
        'SiStrip/MechanicalView/TIB/layer_2/Summary_ClusterStoNCorr__OnTrack__TIB__layer__2',
        'SiStrip/MechanicalView/TIB/layer_3/Summary_ClusterStoNCorr__OnTrack__TIB__layer__3',
        'SiStrip/MechanicalView/TIB/layer_4/Summary_ClusterStoNCorr__OnTrack__TIB__layer__4',
        'SiStrip/MechanicalView/TID/PLUS/wheel_1/NormalizedHitResiduals_TID__wheel__1',
        'SiStrip/MechanicalView/TID/PLUS/wheel_2/NormalizedHitResiduals_TID__wheel__2',
        'SiStrip/MechanicalView/TID/PLUS/wheel_3/NormalizedHitResiduals_TID__wheel__3',
        'SiStrip/MechanicalView/TID/MINUS/wheel_1/NormalizedHitResiduals_TID__wheel__1',
        'SiStrip/MechanicalView/TID/MINUS/wheel_2/NormalizedHitResiduals_TID__wheel__2',
        'SiStrip/MechanicalView/TID/MINUS/wheel_3/NormalizedHitResiduals_TID__wheel__3',
        'SiStrip/MechanicalView/TID/PLUS/wheel_1/Summary_ClusterStoNCorr__OnTrack__TID__PLUS__wheel__1',
        'SiStrip/MechanicalView/TID/PLUS/wheel_2/Summary_ClusterStoNCorr__OnTrack__TID__PLUS__wheel__2',
        'SiStrip/MechanicalView/TID/PLUS/wheel_3/Summary_ClusterStoNCorr__OnTrack__TID__PLUS__wheel__3',
        'SiStrip/MechanicalView/TID/MINUS/wheel_1/Summary_ClusterStoNCorr__OnTrack__TID__MINUS__wheel__1',
        'SiStrip/MechanicalView/TID/MINUS/wheel_2/Summary_ClusterStoNCorr__OnTrack__TID__MINUS__wheel__2',
        'SiStrip/MechanicalView/TID/MINUS/wheel_3/Summary_ClusterStoNCorr__OnTrack__TID__MINUS__wheel__3',
        'SiStrip/MechanicalView/TOB/layer_1/NormalizedHitResiduals_TOB__Layer__1',
        'SiStrip/MechanicalView/TOB/layer_2/NormalizedHitResiduals_TOB__Layer__2',
        'SiStrip/MechanicalView/TOB/layer_3/NormalizedHitResiduals_TOB__Layer__3',
        'SiStrip/MechanicalView/TOB/layer_4/NormalizedHitResiduals_TOB__Layer__4',
        'SiStrip/MechanicalView/TOB/layer_5/NormalizedHitResiduals_TOB__Layer__5',
        'SiStrip/MechanicalView/TOB/layer_6/NormalizedHitResiduals_TOB__Layer__6',
        'SiStrip/MechanicalView/TOB/layer_1/Summary_ClusterStoNCorr__OnTrack__TOB__layer__1',
        'SiStrip/MechanicalView/TOB/layer_2/Summary_ClusterStoNCorr__OnTrack__TOB__layer__2',
        'SiStrip/MechanicalView/TOB/layer_3/Summary_ClusterStoNCorr__OnTrack__TOB__layer__3',
        'SiStrip/MechanicalView/TOB/layer_4/Summary_ClusterStoNCorr__OnTrack__TOB__layer__4',
        'SiStrip/MechanicalView/TOB/layer_5/Summary_ClusterStoNCorr__OnTrack__TOB__layer__5',
        'SiStrip/MechanicalView/TOB/layer_6/Summary_ClusterStoNCorr__OnTrack__TOB__layer__6',
        'SiStrip/MechanicalView/MainDiagonal Position',
        'SiStrip/MechanicalView/NumberOfClustersInPixel',
        'SiStrip/MechanicalView/NumberOfClustersInStrip',
        'Tracking/TrackParameters/generalTracks/LSanalysis/Chi2oNDF_lumiFlag_GenTk',
        'Tracking/TrackParameters/generalTracks/LSanalysis/NumberOfRecHitsPerTrack_lumiFlag_GenTk',
        'Tracking/TrackParameters/generalTracks/LSanalysis/NumberOfTracks_lumiFlag_GenTk',
        'Tracking/TrackParameters/highPurityTracks/pt_1/GeneralProperties/SIPDxyToPV_GenTk',
        'Tracking/TrackParameters/highPurityTracks/pt_1/GeneralProperties/SIPDzToPV_GenTk',
        'Tracking/TrackParameters/highPurityTracks/pt_1/GeneralProperties/SIP3DToPV_GenTk',
        'Tracking/TrackParameters/generalTracks/HitProperties/NumberOfMissingOuterRecHitsPerTrack_GenTk',
        'Tracking/TrackParameters/generalTracks/HitProperties/NumberMORecHitsPerTrackVsPt_GenTk',
        'OfflinePV/offlinePrimaryVertices/tagVtxProb',
        'OfflinePV/offlinePrimaryVertices/tagType',
        'OfflinePV/Resolution/PV/pull_x',
        'OfflinePV/Resolution/PV/pull_y',
        'OfflinePV/Resolution/PV/pull_z',
        'JetMET/Jet/Cleanedak4PFJetsCHS/CHFrac_highPt_Barrel',
        'JetMET/Jet/Cleanedak4PFJetsCHS/CHFrac_highPt_EndCap',
        'JetMET/Jet/Cleanedak4PFJetsCHS/CHFrac_mediumPt_Barrel',
        'JetMET/Jet/Cleanedak4PFJetsCHS/CHFrac_mediumPt_EndCap',
        'JetMET/Jet/Cleanedak4PFJetsCHS/CHFrac_lowPt_Barrel',
        'JetMET/Jet/Cleanedak4PFJetsCHS/CHFrac_lowPt_EndCap',
        'JetMET/Jet/Cleanedak4PFJetsCHS/ChMultiplicity_highPt_Barrel',
        'JetMET/Jet/Cleanedak4PFJetsCHS/ChMultiplicity_highPt_EndCap',
        'JetMET/Jet/Cleanedak4PFJetsCHS/ChMultiplicity_mediumPt_Barrel',
        'JetMET/Jet/Cleanedak4PFJetsCHS/ChMultiplicity_mediumPt_EndCap',
        'JetMET/Jet/Cleanedak4PFJetsCHS/ChMultiplicity_lowPt_Barrel',
        'JetMET/Jet/Cleanedak4PFJetsCHS/ChMultiplicity_lowPt_EndCap',
        'JetMET/Jet/Cleanedak4PFJetsCHS/Constituents',
        'JetMET/Jet/Cleanedak4PFJetsCHS/Eta',
        'JetMET/Jet/Cleanedak4PFJetsCHS/Eta_uncor',
        'JetMET/Jet/Cleanedak4PFJetsCHS/JetEnergyCorr',
        'JetMET/Jet/Cleanedak4PFJetsCHS/NJets',
        'JetMET/Jet/Cleanedak4PFJetsCHS/Phi',
        'JetMET/Jet/Cleanedak4PFJetsCHS/Phi_Barrel',
        'JetMET/Jet/Cleanedak4PFJetsCHS/Phi_EndCap',
        'JetMET/Jet/Cleanedak4PFJetsCHS/Pt',
        'JetMET/MET/pfMETT1/Cleaned/METSig',
        'JetMET/vertices'
    )
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

process.variableJTAPars = cms.PSet(
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5)
)

process.ZDC2018Pedestal_0 = cms.VPSet(
    cms.PSet(
        object = cms.untracked.string('hZDCM_EM0'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_EM1'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_EM2'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_EM3'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_EM4'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_EM5'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_EM6'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_EM7'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_EM8'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_EM9'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_EM10'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_EM11'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_EM12'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_EM13'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_EM14'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_EM15'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_HAD0'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_HAD1'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_HAD2'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_HAD3'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_HAD4'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_HAD5'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_HAD6'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_HAD7'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_HAD8'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_HAD9'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_HAD10'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_HAD11'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_HAD12'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_HAD13'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_HAD14'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_HAD15'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_EM0'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_EM1'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_EM2'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_EM3'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_EM4'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_EM5'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_EM6'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_EM7'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_EM8'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_EM9'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_EM10'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_EM11'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_EM12'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_EM13'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_EM14'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_EM15'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_HAD0'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_HAD1'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_HAD2'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_HAD3'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_HAD4'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_HAD5'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_HAD6'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_HAD7'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_HAD8'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_HAD9'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_HAD10'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_HAD11'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_HAD12'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_HAD13'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_HAD14'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_HAD15'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_RPD0'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_RPD1'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_RPD2'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_RPD3'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_RPD4'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_RPD5'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_RPD6'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_RPD7'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_RPD8'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_RPD9'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_RPD10'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_RPD11'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_RPD12'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_RPD13'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_RPD14'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCM_RPD15'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_RPD0'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_RPD1'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_RPD2'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_RPD3'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_RPD4'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_RPD5'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_RPD6'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_RPD7'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_RPD8'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_RPD9'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_RPD10'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_RPD11'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_RPD12'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_RPD13'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_RPD14'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    ),
    cms.PSet(
        object = cms.untracked.string('hZDCP_RPD15'),
        ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
    )
)

process.gainFill7435 = cms.VPSet(
    cms.PSet(
        calib = cms.untracked.double(0.11),
        object = cms.untracked.string('hZDCP_EM1')
    ),
    cms.PSet(
        calib = cms.untracked.double(0.11),
        object = cms.untracked.string('hZDCP_EM2')
    ),
    cms.PSet(
        calib = cms.untracked.double(0.08),
        object = cms.untracked.string('hZDCP_EM3')
    ),
    cms.PSet(
        calib = cms.untracked.double(0.1),
        object = cms.untracked.string('hZDCP_EM4')
    ),
    cms.PSet(
        calib = cms.untracked.double(0.06),
        object = cms.untracked.string('hZDCP_EM5')
    ),
    cms.PSet(
        calib = cms.untracked.double(1.0),
        object = cms.untracked.string('hZDCP_HAD1')
    ),
    cms.PSet(
        calib = cms.untracked.double(0.61),
        object = cms.untracked.string('hZDCP_HAD2')
    ),
    cms.PSet(
        calib = cms.untracked.double(0.46),
        object = cms.untracked.string('hZDCP_HAD3')
    ),
    cms.PSet(
        calib = cms.untracked.double(0.52),
        object = cms.untracked.string('hZDCP_HAD4')
    ),
    cms.PSet(
        calib = cms.untracked.double(0.12),
        object = cms.untracked.string('hZDCM_EM1')
    ),
    cms.PSet(
        calib = cms.untracked.double(0.12),
        object = cms.untracked.string('hZDCM_EM2')
    ),
    cms.PSet(
        calib = cms.untracked.double(0.08),
        object = cms.untracked.string('hZDCM_EM3')
    ),
    cms.PSet(
        calib = cms.untracked.double(0.1),
        object = cms.untracked.string('hZDCM_EM4')
    ),
    cms.PSet(
        calib = cms.untracked.double(0.11),
        object = cms.untracked.string('hZDCM_EM5')
    ),
    cms.PSet(
        calib = cms.untracked.double(1.0),
        object = cms.untracked.string('hZDCM_HAD1')
    ),
    cms.PSet(
        calib = cms.untracked.double(0.81),
        object = cms.untracked.string('hZDCM_HAD2')
    ),
    cms.PSet(
        calib = cms.untracked.double(0.46),
        object = cms.untracked.string('hZDCM_HAD3')
    ),
    cms.PSet(
        calib = cms.untracked.double(0.58),
        object = cms.untracked.string('hZDCM_HAD4')
    ),
    cms.PSet(
        calib = cms.untracked.double(2.67),
        object = cms.untracked.string('Pscale')
    ),
    cms.PSet(
        calib = cms.untracked.double(4.45),
        object = cms.untracked.string('Mscale')
    )
)

process.PFTowers = cms.EDProducer("PFTowers",
    mightGet = cms.optional.untracked.vstring,
    src = cms.InputTag("particleFlow"),
    useHF = cms.bool(True)
)


process.PackedPFTowers = cms.EDProducer("PackedPFTowers",
    mightGet = cms.optional.untracked.vstring,
    src = cms.InputTag("packedPFCandidates"),
    useHF = cms.bool(True)
)


process.QWzdcreco = cms.EDProducer("QWZDC2018RecHit",
    ZDCCalib = cms.VPSet(
        cms.PSet(
            calib = cms.untracked.double(0.11),
            object = cms.untracked.string('hZDCP_EM1')
        ),
        cms.PSet(
            calib = cms.untracked.double(0.11),
            object = cms.untracked.string('hZDCP_EM2')
        ),
        cms.PSet(
            calib = cms.untracked.double(0.08),
            object = cms.untracked.string('hZDCP_EM3')
        ),
        cms.PSet(
            calib = cms.untracked.double(0.1),
            object = cms.untracked.string('hZDCP_EM4')
        ),
        cms.PSet(
            calib = cms.untracked.double(0.06),
            object = cms.untracked.string('hZDCP_EM5')
        ),
        cms.PSet(
            calib = cms.untracked.double(1.0),
            object = cms.untracked.string('hZDCP_HAD1')
        ),
        cms.PSet(
            calib = cms.untracked.double(0.61),
            object = cms.untracked.string('hZDCP_HAD2')
        ),
        cms.PSet(
            calib = cms.untracked.double(0.46),
            object = cms.untracked.string('hZDCP_HAD3')
        ),
        cms.PSet(
            calib = cms.untracked.double(0.52),
            object = cms.untracked.string('hZDCP_HAD4')
        ),
        cms.PSet(
            calib = cms.untracked.double(0.12),
            object = cms.untracked.string('hZDCM_EM1')
        ),
        cms.PSet(
            calib = cms.untracked.double(0.12),
            object = cms.untracked.string('hZDCM_EM2')
        ),
        cms.PSet(
            calib = cms.untracked.double(0.08),
            object = cms.untracked.string('hZDCM_EM3')
        ),
        cms.PSet(
            calib = cms.untracked.double(0.1),
            object = cms.untracked.string('hZDCM_EM4')
        ),
        cms.PSet(
            calib = cms.untracked.double(0.11),
            object = cms.untracked.string('hZDCM_EM5')
        ),
        cms.PSet(
            calib = cms.untracked.double(1.0),
            object = cms.untracked.string('hZDCM_HAD1')
        ),
        cms.PSet(
            calib = cms.untracked.double(0.81),
            object = cms.untracked.string('hZDCM_HAD2')
        ),
        cms.PSet(
            calib = cms.untracked.double(0.46),
            object = cms.untracked.string('hZDCM_HAD3')
        ),
        cms.PSet(
            calib = cms.untracked.double(0.58),
            object = cms.untracked.string('hZDCM_HAD4')
        ),
        cms.PSet(
            calib = cms.untracked.double(2.67),
            object = cms.untracked.string('Pscale')
        ),
        cms.PSet(
            calib = cms.untracked.double(4.45),
            object = cms.untracked.string('Mscale')
        )
    ),
    srcDetId = cms.untracked.InputTag("zdcdigi","DetId"),
    srcfC = cms.untracked.InputTag("zdcdigi","regularfC")
)


process.ak4PFJetsForFlow = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(1),
    GhostArea = cms.double(0.01),
    Ghost_EtaMax = cms.double(5.0),
    Rho_EtaMax = cms.double(4.4),
    applyWeight = cms.bool(False),
    doAreaDiskApprox = cms.bool(False),
    doAreaFastjet = cms.bool(False),
    doFastJetNonUniform = cms.bool(True),
    doPUOffsetCorr = cms.bool(True),
    doPVCorrection = cms.bool(False),
    doRhoFastjet = cms.bool(False),
    inputEMin = cms.double(0.0),
    inputEtMin = cms.double(0.0),
    jetAlgorithm = cms.string('AntiKt'),
    jetPtMin = cms.double(15.0),
    jetType = cms.string('BasicJet'),
    maxBadEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxRecoveredHcalCells = cms.uint32(9999999),
    minSeed = cms.uint32(14327),
    nExclude = cms.uint32(2),
    nSigmaPU = cms.double(1.0),
    puCenters = cms.vdouble(
        -5, -4, -3, -2, -1,
        0, 1, 2, 3, 4,
        5
    ),
    puPtMin = cms.double(25),
    puWidth = cms.double(0.8),
    rParam = cms.double(0.4),
    radiusPU = cms.double(0.5),
    src = cms.InputTag("PackedPFTowers"),
    srcPVs = cms.InputTag(""),
    subtractorName = cms.string('MultipleAlgoIterator'),
    sumRecHits = cms.bool(False),
    useDeterministicSeed = cms.bool(True),
    voronoiRfact = cms.double(-0.9)
)


process.akCs0PFJets = cms.EDProducer("CSJetProducer",
    Active_Area_Repeats = cms.int32(1),
    GhostArea = cms.double(0.005),
    Ghost_EtaMax = cms.double(6.5),
    Rho_EtaMax = cms.double(4.5),
    applyWeight = cms.bool(False),
    csAlpha = cms.double(2.0),
    csRParam = cms.double(-1.0),
    doAreaDiskApprox = cms.bool(False),
    doAreaFastjet = cms.bool(True),
    doFastJetNonUniform = cms.bool(True),
    doPUOffsetCorr = cms.bool(False),
    doPVCorrection = cms.bool(False),
    doRhoFastjet = cms.bool(True),
    etaMap = cms.InputTag("hiPuRho","mapEtaEdges"),
    inputEMin = cms.double(0.0),
    inputEtMin = cms.double(0.0),
    jetAlgorithm = cms.string('AntiKt'),
    jetCollInstanceName = cms.string('pfParticlesCs'),
    jetPtMin = cms.double(0.0),
    jetType = cms.string('PFJet'),
    maxBadEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxRecoveredHcalCells = cms.uint32(9999999),
    minSeed = cms.uint32(14327),
    nExclude = cms.uint32(2),
    nSigmaPU = cms.double(1.0),
    puCenters = cms.vdouble(
        -5, -4, -3, -2, -1,
        0, 1, 2, 3, 4,
        5
    ),
    puWidth = cms.double(0.8),
    rParam = cms.double(0.4),
    radiusPU = cms.double(0.5),
    rho = cms.InputTag("hiPuRho","mapToRho"),
    rhoFlowFitParams = cms.InputTag("hiFJRhoFlowModulation","rhoFlowFitParams"),
    rhom = cms.InputTag("hiPuRho","mapToRhoM"),
    src = cms.InputTag("packedPFCandidates"),
    srcPVs = cms.InputTag(""),
    useDeterministicSeed = cms.bool(True),
    useExplicitGhosts = cms.bool(True),
    useModulatedRho = cms.bool(False),
    voronoiRfact = cms.double(-0.9),
    writeJetsWithConst = cms.bool(True)
)


process.akCs0PFpatJetCorrFactors = cms.EDProducer("JetCorrFactorsProducer",
    emf = cms.bool(False),
    extraJPTOffset = cms.string('L1FastJet'),
    flavorType = cms.string('J'),
    levels = cms.vstring(
        'L2Relative',
        'L2L3Residual'
    ),
    mightGet = cms.optional.untracked.vstring,
    payload = cms.string('AK4PF'),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    src = cms.InputTag("akCs0PFJets"),
    useNPV = cms.bool(False),
    useRho = cms.bool(False)
)


process.akCs0PFpatJets = cms.EDProducer("PATJetProducer",
    JetFlavourInfoSource = cms.InputTag("akCs0PFpatJetFlavourAssociation"),
    JetPartonMapSource = cms.InputTag("akCs0PFpatJetFlavourAssociationLegacy"),
    addAssociatedTracks = cms.bool(False),
    addBTagInfo = cms.bool(True),
    addDiscriminators = cms.bool(False),
    addEfficiencies = cms.bool(False),
    addGenJetMatch = cms.bool(True),
    addGenPartonMatch = cms.bool(True),
    addJetCharge = cms.bool(False),
    addJetCorrFactors = cms.bool(True),
    addJetFlavourInfo = cms.bool(True),
    addJetID = cms.bool(False),
    addPartonJetMatch = cms.bool(False),
    addResolutions = cms.bool(False),
    addTagInfos = cms.bool(False),
    discriminatorSources = cms.VInputTag(cms.InputTag("akCs0PFpfDeepCSVJetTags","probb"), cms.InputTag("akCs0PFpfDeepCSVJetTags","probc"), cms.InputTag("akCs0PFpfDeepCSVJetTags","probudsg"), cms.InputTag("akCs0PFpfDeepCSVJetTags","probbb"), cms.InputTag("akCs0PFpfJetProbabilityBJetTags")),
    efficiencies = cms.PSet(

    ),
    embedCaloTowers = cms.bool(False),
    embedGenJetMatch = cms.bool(True),
    embedGenPartonMatch = cms.bool(True),
    embedPFCandidates = cms.bool(False),
    genJetMatch = cms.InputTag("akCs0PFpatJetGenJetMatch"),
    genPartonMatch = cms.InputTag("akCs0PFpatJetPartonMatch"),
    getJetMCFlavour = cms.bool(True),
    jetChargeSource = cms.InputTag("patJetCharge"),
    jetCorrFactorsSource = cms.VInputTag("akCs0PFpatJetCorrFactors"),
    jetIDMap = cms.InputTag("ak4JetID"),
    jetSource = cms.InputTag("akCs0PFJets"),
    mightGet = cms.optional.untracked.vstring,
    partonJetSource = cms.InputTag("NOT_IMPLEMENTED"),
    resolutions = cms.PSet(

    ),
    tagInfoSources = cms.VInputTag(),
    trackAssociationSource = cms.InputTag(""),
    useLegacyJetMCFlavour = cms.bool(True),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.akCs0PFpfDeepCSVJetTags = cms.EDProducer("DeepFlavourJetTagsProducer",
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepCSV_PhaseI.json'),
    checkSVForDefaults = cms.bool(True),
    meanPadding = cms.bool(True),
    src = cms.InputTag("akCs0PFpfDeepCSVTagInfos"),
    toAdd = cms.PSet(

    )
)


process.akCs0PFpfDeepCSVTagInfos = cms.EDProducer("DeepNNTagInfoProducer",
    computer = cms.PSet(
        SoftLeptonFlip = cms.bool(False),
        charmCut = cms.double(1.5),
        correctVertexMass = cms.bool(True),
        minimumTrackWeight = cms.double(0.5),
        pseudoMultiplicityMin = cms.uint32(2),
        pseudoVertexV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        ),
        trackFlip = cms.bool(False),
        trackMultiplicityMin = cms.uint32(2),
        trackPairV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.03)
        ),
        trackPseudoSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(2.0),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(99999.9),
            sip3dSigMin = cms.double(-99999.9),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(-99999.9),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(99999.9),
            sip3dSigMin = cms.double(-99999.9),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSort = cms.string('sip2dSig'),
        useTrackWeights = cms.bool(True),
        vertexFlip = cms.bool(False)
    ),
    svTagInfos = cms.InputTag("akCs0PFpfSecondaryVertexTagInfos")
)


process.akCs0PFpfImpactParameterTagInfos = cms.EDProducer("CandIPProducer",
    candidates = cms.InputTag("packedPFCandidates"),
    computeGhostTrack = cms.bool(True),
    computeProbabilities = cms.bool(True),
    ghostTrackPriorDeltaR = cms.double(0.03),
    jetDirectionUsingGhostTrack = cms.bool(False),
    jetDirectionUsingTracks = cms.bool(False),
    jets = cms.InputTag("akCs0PFJets"),
    maxDeltaR = cms.double(0.4),
    maximumChiSquared = cms.double(5.0),
    maximumLongitudinalImpactParameter = cms.double(17.0),
    maximumTransverseImpactParameter = cms.double(0.2),
    minimumNumberOfHits = cms.int32(0),
    minimumNumberOfPixelHits = cms.int32(1),
    minimumTransverseMomentum = cms.double(1.0),
    primaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    useTrackQuality = cms.bool(False)
)


process.akCs0PFpfJetProbabilityBJetTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('candidateJetProbabilityComputer'),
    tagInfos = cms.VInputTag("akCs0PFpfImpactParameterTagInfos")
)


process.akCs0PFpfSecondaryVertexTagInfos = cms.EDProducer("CandSecondaryVertexProducer",
    beamSpotTag = cms.InputTag("offlineBeamSpot"),
    constraint = cms.string('BeamSpot'),
    extSVCollection = cms.InputTag("secondaryVertices"),
    extSVDeltaRToJet = cms.double(0.3),
    minimumTrackWeight = cms.double(0.5),
    trackIPTagInfos = cms.InputTag("akCs0PFpfImpactParameterTagInfos"),
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(99999.9),
        maxDistToAxis = cms.double(0.2),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(1),
        ptMin = cms.double(1.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSort = cms.string('sip3dSig'),
    useExternalSV = cms.bool(False),
    usePVError = cms.bool(True),
    vertexCuts = cms.PSet(
        distSig2dMax = cms.double(99999.9),
        distSig2dMin = cms.double(3.0),
        distSig3dMax = cms.double(99999.9),
        distSig3dMin = cms.double(-99999.9),
        distVal2dMax = cms.double(2.5),
        distVal2dMin = cms.double(0.01),
        distVal3dMax = cms.double(99999.9),
        distVal3dMin = cms.double(-99999.9),
        fracPV = cms.double(0.65),
        massMax = cms.double(6.5),
        maxDeltaRToJetAxis = cms.double(0.4),
        minimumTrackWeight = cms.double(0.5),
        multiplicityMin = cms.uint32(2),
        useTrackWeights = cms.bool(True),
        v0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        )
    ),
    vertexReco = cms.PSet(
        finder = cms.string('avr'),
        minweight = cms.double(0.5),
        primcut = cms.double(1.8),
        seccut = cms.double(6.0),
        smoothing = cms.bool(False),
        weightthreshold = cms.double(0.001)
    ),
    vertexSelection = cms.PSet(
        sortCriterium = cms.string('dist3dError')
    )
)


process.allPartons = cms.EDProducer("PartonSelector",
    src = cms.InputTag("hiSignalGenParticles"),
    withLeptons = cms.bool(False)
)


process.hiFJRhoFlowModulation = cms.EDProducer("HiFJRhoFlowModulationProducer",
    EvtPlane = cms.InputTag("hiEvtPlane"),
    doEvtPlane = cms.bool(False),
    doFreePlaneFit = cms.bool(False),
    doJettyExclusion = cms.bool(False),
    evtPlaneLevel = cms.int32(0),
    exclusionRadius = cms.double(0.4),
    firstFittedVn = cms.int32(2),
    jetTag = cms.InputTag("ak4PFJetsForFlow"),
    lastFittedVn = cms.int32(3),
    mightGet = cms.optional.untracked.vstring,
    minPfCandidatesPerEvent = cms.int32(100),
    pfCandSource = cms.InputTag("packedPFCandidates"),
    pfCandidateEtaCut = cms.double(1),
    pfCandidateMaxPtCut = cms.double(3),
    pfCandidateMinPtCut = cms.double(0.3)
)


process.hiFJRhoFlowModulationProducer = cms.EDProducer("HiFJRhoFlowModulationProducer",
    EvtPlane = cms.InputTag("hiEvtPlane"),
    doEvtPlane = cms.bool(False),
    doFreePlaneFit = cms.bool(False),
    doJettyExclusion = cms.bool(False),
    evtPlaneLevel = cms.int32(0),
    exclusionRadius = cms.double(0.4),
    firstFittedVn = cms.int32(2),
    jetTag = cms.InputTag("ak4PFJetsForFlow"),
    lastFittedVn = cms.int32(3),
    mightGet = cms.optional.untracked.vstring,
    minPfCandidatesPerEvent = cms.int32(100),
    pfCandSource = cms.InputTag("packedPFCandidates"),
    pfCandidateEtaCut = cms.double(1),
    pfCandidateMaxPtCut = cms.double(3),
    pfCandidateMinPtCut = cms.double(0.3)
)


process.hiPuRho = cms.EDProducer("HiPuRhoProducer",
    dropZeroTowers = cms.bool(True),
    medianWindowWidth = cms.int32(2),
    mightGet = cms.optional.untracked.vstring,
    minimumTowersFraction = cms.double(0.7),
    nSigmaPU = cms.double(1),
    puPtMin = cms.double(15),
    rParam = cms.double(0.3),
    radiusPU = cms.double(0.5),
    src = cms.InputTag("PackedPFTowers"),
    towSigmaCut = cms.double(5)
)


process.hiSignalGenParticles = cms.EDProducer("HiSignalParticleProducer",
    mightGet = cms.optional.untracked.vstring,
    src = cms.InputTag("prunedGenParticles")
)


process.pfDeepCSVJetTags = cms.EDProducer("DeepFlavourJetTagsProducer",
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepCSV_PhaseI.json'),
    checkSVForDefaults = cms.bool(True),
    meanPadding = cms.bool(True),
    src = cms.InputTag("pfDeepCSVTagInfos"),
    toAdd = cms.PSet(

    )
)


process.pfDeepCSVTagInfos = cms.EDProducer("DeepNNTagInfoProducer",
    computer = cms.PSet(
        SoftLeptonFlip = cms.bool(False),
        charmCut = cms.double(1.5),
        correctVertexMass = cms.bool(True),
        minimumTrackWeight = cms.double(0.5),
        pseudoMultiplicityMin = cms.uint32(2),
        pseudoVertexV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        ),
        trackFlip = cms.bool(False),
        trackMultiplicityMin = cms.uint32(2),
        trackPairV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.03)
        ),
        trackPseudoSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(2.0),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(99999.9),
            sip3dSigMin = cms.double(-99999.9),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(-99999.9),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(99999.9),
            sip3dSigMin = cms.double(-99999.9),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSort = cms.string('sip2dSig'),
        useTrackWeights = cms.bool(True),
        vertexFlip = cms.bool(False)
    ),
    svTagInfos = cms.InputTag("pfSecondaryVertexTagInfos")
)


process.pfImpactParameterTagInfos = cms.EDProducer("CandIPProducer",
    candidates = cms.InputTag("packedPFCandidates"),
    computeGhostTrack = cms.bool(True),
    computeProbabilities = cms.bool(True),
    ghostTrackPriorDeltaR = cms.double(0.03),
    jetDirectionUsingGhostTrack = cms.bool(False),
    jetDirectionUsingTracks = cms.bool(False),
    jets = cms.InputTag("slimmedJets"),
    maxDeltaR = cms.double(0.4),
    maximumChiSquared = cms.double(5.0),
    maximumLongitudinalImpactParameter = cms.double(17.0),
    maximumTransverseImpactParameter = cms.double(0.2),
    minimumNumberOfHits = cms.int32(0),
    minimumNumberOfPixelHits = cms.int32(1),
    minimumTransverseMomentum = cms.double(1.0),
    primaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    useTrackQuality = cms.bool(False)
)


process.pfJetProbabilityBJetTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('candidateJetProbabilityComputer'),
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfos"))
)


process.pfSecondaryVertexTagInfos = cms.EDProducer("CandSecondaryVertexProducer",
    beamSpotTag = cms.InputTag("offlineBeamSpot"),
    constraint = cms.string('BeamSpot'),
    extSVCollection = cms.InputTag("secondaryVertices"),
    extSVDeltaRToJet = cms.double(0.3),
    minimumTrackWeight = cms.double(0.5),
    trackIPTagInfos = cms.InputTag("pfImpactParameterTagInfos"),
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(99999.9),
        maxDistToAxis = cms.double(0.2),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(1),
        ptMin = cms.double(1.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSort = cms.string('sip3dSig'),
    useExternalSV = cms.bool(False),
    usePVError = cms.bool(True),
    vertexCuts = cms.PSet(
        distSig2dMax = cms.double(99999.9),
        distSig2dMin = cms.double(3.0),
        distSig3dMax = cms.double(99999.9),
        distSig3dMin = cms.double(-99999.9),
        distVal2dMax = cms.double(2.5),
        distVal2dMin = cms.double(0.01),
        distVal3dMax = cms.double(99999.9),
        distVal3dMin = cms.double(-99999.9),
        fracPV = cms.double(0.65),
        massMax = cms.double(6.5),
        maxDeltaRToJetAxis = cms.double(0.4),
        minimumTrackWeight = cms.double(0.5),
        multiplicityMin = cms.uint32(2),
        useTrackWeights = cms.bool(True),
        v0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        )
    ),
    vertexReco = cms.PSet(
        finder = cms.string('avr'),
        minweight = cms.double(0.5),
        primcut = cms.double(1.8),
        seccut = cms.double(6.0),
        smoothing = cms.bool(False),
        weightthreshold = cms.double(0.001)
    ),
    vertexSelection = cms.PSet(
        sortCriterium = cms.string('dist3dError')
    )
)


process.randomEngineStateProducer = cms.EDProducer("RandomEngineStateProducer")


process.siPixelRecHits = cms.EDProducer("SiPixelRecHitConverter",
    CPE = cms.string('PixelCPEGeneric'),
    VerboseLevel = cms.untracked.int32(0),
    src = cms.InputTag("siPixelClusters")
)


process.siPixelRecHitsPreSplittingCPU = cms.EDProducer("SiPixelRecHitSoAFromLegacyHIonPhase1",
    CPE = cms.string('PixelCPEFastHIonPhase1'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    convertToLegacy = cms.bool(True),
    mightGet = cms.optional.untracked.vstring,
    src = cms.InputTag("siPixelClustersPreSplitting")
)


process.siPixelRecHitsPreSplittingCUDA = cms.EDProducer("SiPixelRecHitCUDAHIonPhase1",
    CPE = cms.string('PixelCPEFastHIonPhase1'),
    beamSpot = cms.InputTag("offlineBeamSpotToCUDA"),
    mightGet = cms.optional.untracked.vstring,
    src = cms.InputTag("siPixelClustersPreSplittingCUDA")
)


process.unpackedMuons = cms.EDProducer("MuonUnpacker",
    addPropToMuonSt = cms.bool(False),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    cosmicPropagationHypothesis = cms.bool(False),
    fallbackToME1 = cms.bool(True),
    mightGet = cms.optional.untracked.vstring,
    muonSelectors = cms.vstring(
        'AllTrackerMuons',
        'TMOneStationTight'
    ),
    muons = cms.InputTag("slimmedMuons"),
    primaryVertices = cms.InputTag("unpackedTracksAndVertices"),
    propagatorAlong = cms.ESInputTag("","hltESPSteppingHelixPropagatorAlong"),
    propagatorAny = cms.ESInputTag("","SteppingHelixPropagatorAny"),
    propagatorOpposite = cms.ESInputTag("","hltESPSteppingHelixPropagatorOpposite"),
    tracks = cms.InputTag("unpackedTracksAndVertices"),
    useMB2InOverlap = cms.bool(True),
    useSimpleGeometry = cms.bool(True),
    useState = cms.string('atVertex'),
    useStation2 = cms.bool(True),
    useTrack = cms.string('tracker')
)


process.unpackedTracksAndVertices = cms.EDProducer("TrackAndVertexUnpacker",
    lostTrackNormChi2Map = cms.InputTag("lostTrackChi2"),
    lostTracks = cms.InputTag("lostTracks"),
    mightGet = cms.optional.untracked.vstring,
    packedPFCandidateNormChi2Map = cms.InputTag("packedPFCandidateTrackChi2"),
    packedPFCandidates = cms.InputTag("packedPFCandidates"),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    recoverTracks = cms.bool(False),
    secondaryVertices = cms.InputTag("slimmedSecondaryVertices")
)


process.zdcdigi = cms.EDProducer("QWZDC2018Producer2",
    Debug = cms.untracked.bool(False),
    HardCode = cms.untracked.bool(True),
    Pedestal = cms.VPSet(
        cms.PSet(
            object = cms.untracked.string('hZDCM_EM0'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_EM1'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_EM2'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_EM3'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_EM4'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_EM5'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_EM6'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_EM7'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_EM8'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_EM9'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_EM10'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_EM11'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_EM12'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_EM13'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_EM14'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_EM15'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_HAD0'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_HAD1'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_HAD2'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_HAD3'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_HAD4'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_HAD5'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_HAD6'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_HAD7'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_HAD8'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_HAD9'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_HAD10'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_HAD11'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_HAD12'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_HAD13'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_HAD14'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_HAD15'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_EM0'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_EM1'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_EM2'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_EM3'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_EM4'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_EM5'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_EM6'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_EM7'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_EM8'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_EM9'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_EM10'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_EM11'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_EM12'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_EM13'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_EM14'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_EM15'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_HAD0'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_HAD1'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_HAD2'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_HAD3'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_HAD4'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_HAD5'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_HAD6'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_HAD7'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_HAD8'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_HAD9'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_HAD10'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_HAD11'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_HAD12'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_HAD13'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_HAD14'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_HAD15'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_RPD0'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_RPD1'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_RPD2'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_RPD3'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_RPD4'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_RPD5'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_RPD6'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_RPD7'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_RPD8'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_RPD9'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_RPD10'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_RPD11'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_RPD12'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_RPD13'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_RPD14'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCM_RPD15'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_RPD0'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_RPD1'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_RPD2'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_RPD3'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_RPD4'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_RPD5'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_RPD6'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_RPD7'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_RPD8'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_RPD9'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_RPD10'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_RPD11'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_RPD12'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_RPD13'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_RPD14'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        ),
        cms.PSet(
            object = cms.untracked.string('hZDCP_RPD15'),
            ped = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0)
        )
    ),
    SOI = cms.untracked.int32(4),
    Src = cms.untracked.InputTag("hcalDigis","ZDC")
)


process.siPixelRecHitsPreSplitting = SwitchProducerCUDA(
    cpu = cms.EDProducer("SiPixelRecHitConverter",
        CPE = cms.string('PixelCPEGeneric'),
        VerboseLevel = cms.untracked.int32(0),
        src = cms.InputTag("siPixelClustersPreSplitting")
    )
)


process.siPixelRecHitsPreSplittingSoA = SwitchProducerCUDA(
    cpu = cms.EDAlias(
        siPixelRecHitsPreSplittingCPU = cms.VPSet(
            cms.PSet(
                type = cms.string('pixelTopologyHIonPhase1TrackingRecHitSoAHost')
            ),
            cms.PSet(
                type = cms.string('uintAsHostProduct')
            )
        )
    )
)


process.clusterCompatibilityFilter = cms.EDFilter("HIClusterCompatibilityFilter",
    cluscomSrc = cms.InputTag("hiClusterCompatibility"),
    clusterPars = cms.vdouble(0.0, 0.0045),
    clusterTrunc = cms.double(2.0),
    maxZ = cms.double(20.05),
    minZ = cms.double(-20.0),
    nhitsTrunc = cms.int32(150)
)


process.hltPixelClusterShapeFilter = cms.EDFilter("HLTPixelClusterShapeFilter",
    clusterPars = cms.vdouble(0, 0.0045),
    clusterTrunc = cms.double(2),
    inputTag = cms.InputTag("siPixelRecHits"),
    maxZ = cms.double(20.05),
    mightGet = cms.optional.untracked.vstring,
    minZ = cms.double(-20),
    nhitsTrunc = cms.int32(150),
    saveTags = cms.bool(True),
    zStep = cms.double(0.2)
)


process.phfCoincFilter1Th3 = cms.EDFilter("HiHFFilter",
    HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
    minnumtowers = cms.int32(1),
    threshold = cms.int32(3)
)


process.phfCoincFilter1Th4 = cms.EDFilter("HiHFFilter",
    HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
    minnumtowers = cms.int32(1),
    threshold = cms.int32(4)
)


process.phfCoincFilter1Th5 = cms.EDFilter("HiHFFilter",
    HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
    minnumtowers = cms.int32(1),
    threshold = cms.int32(5)
)


process.phfCoincFilter1Th6 = cms.EDFilter("HiHFFilter",
    HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
    minnumtowers = cms.int32(1),
    threshold = cms.int32(6)
)


process.phfCoincFilter2Th3 = cms.EDFilter("HiHFFilter",
    HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
    minnumtowers = cms.int32(2),
    threshold = cms.int32(3)
)


process.phfCoincFilter2Th4 = cms.EDFilter("HiHFFilter",
    HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
    minnumtowers = cms.int32(2),
    threshold = cms.int32(4)
)


process.phfCoincFilter2Th5 = cms.EDFilter("HiHFFilter",
    HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
    minnumtowers = cms.int32(2),
    threshold = cms.int32(5)
)


process.phfCoincFilter2Th6 = cms.EDFilter("HiHFFilter",
    HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
    minnumtowers = cms.int32(2),
    threshold = cms.int32(6)
)


process.phfCoincFilter3Th3 = cms.EDFilter("HiHFFilter",
    HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
    minnumtowers = cms.int32(3),
    threshold = cms.int32(3)
)


process.phfCoincFilter3Th4 = cms.EDFilter("HiHFFilter",
    HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
    minnumtowers = cms.int32(3),
    threshold = cms.int32(4)
)


process.phfCoincFilter3Th5 = cms.EDFilter("HiHFFilter",
    HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
    minnumtowers = cms.int32(3),
    threshold = cms.int32(5)
)


process.phfCoincFilter3Th6 = cms.EDFilter("HiHFFilter",
    HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
    minnumtowers = cms.int32(3),
    threshold = cms.int32(6)
)


process.phfCoincFilter4Th2 = cms.EDFilter("HiHFFilter",
    HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
    minnumtowers = cms.int32(4),
    threshold = cms.int32(2)
)


process.phfCoincFilter4Th3 = cms.EDFilter("HiHFFilter",
    HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
    minnumtowers = cms.int32(4),
    threshold = cms.int32(3)
)


process.phfCoincFilter4Th4 = cms.EDFilter("HiHFFilter",
    HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
    minnumtowers = cms.int32(4),
    threshold = cms.int32(4)
)


process.phfCoincFilter4Th5 = cms.EDFilter("HiHFFilter",
    HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
    minnumtowers = cms.int32(4),
    threshold = cms.int32(5)
)


process.phfCoincFilter4Th6 = cms.EDFilter("HiHFFilter",
    HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
    minnumtowers = cms.int32(4),
    threshold = cms.int32(6)
)


process.phfCoincFilter5Th3 = cms.EDFilter("HiHFFilter",
    HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
    minnumtowers = cms.int32(5),
    threshold = cms.int32(3)
)


process.phfCoincFilter5Th4 = cms.EDFilter("HiHFFilter",
    HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
    minnumtowers = cms.int32(5),
    threshold = cms.int32(4)
)


process.phfCoincFilter5Th5 = cms.EDFilter("HiHFFilter",
    HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
    minnumtowers = cms.int32(5),
    threshold = cms.int32(5)
)


process.phfCoincFilter5Th6 = cms.EDFilter("HiHFFilter",
    HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
    minnumtowers = cms.int32(5),
    threshold = cms.int32(6)
)


process.primaryVertexFilter = cms.EDFilter("VertexSelector",
    cut = cms.string('!isFake && abs(z) <= 25 && position.Rho <= 2'),
    filter = cms.bool(True),
    src = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.HiForestInfo = cms.EDAnalyzer("HiForestInfo",
    GlobalTagLabel = cms.string('124X_dataRun3_Prompt_v10'),
    HiForestVersion = cms.string(''),
    info = cms.vstring('HiForest, miniAOD, 132X, data')
)


process.PbPbTracks = cms.EDAnalyzer("TrackAnalyzer",
    beamSpotSrc = cms.untracked.InputTag("offlineBeamSpot"),
    chi2Map = cms.InputTag("packedPFCandidateTrackChi2"),
    chi2MapLost = cms.InputTag("lostTrackChi2"),
    doTrack = cms.untracked.bool(True),
    lostTracksSrc = cms.InputTag("lostTracks"),
    packedCandSrc = cms.InputTag("packedPFCandidates"),
    trackPtMin = cms.untracked.double(0.01),
    vertexSrc = cms.InputTag("unpackedTracksAndVertices")
)


process.akCs4PFJetAnalyzer = cms.EDAnalyzer("HiInclusiveJetAnalyzer",
    caloJetTag = cms.InputTag("slimmedCaloJets"),
    doCandidateBtagging = cms.untracked.bool(True),
    doGenSym = cms.untracked.bool(False),
    doGenTaus = cms.untracked.bool(False),
    doHiJetID = cms.untracked.bool(True),
    doJetConstituents = cms.untracked.bool(False),
    doLegacyBtagging = cms.untracked.bool(False),
    doStandardJetID = cms.untracked.bool(False),
    doSubEvent = cms.untracked.bool(False),
    doSubJets = cms.untracked.bool(False),
    doWTARecluster = cms.untracked.bool(True),
    eventInfoTag = cms.InputTag("generator"),
    fillGenJets = cms.untracked.bool(False),
    genjetTag = cms.InputTag("ak4HiGenJets"),
    hltTrgResults = cms.untracked.string('TriggerResults::HISIGNAL'),
    isMC = cms.untracked.bool(False),
    jetName = cms.untracked.string('akCs0PF'),
    jetPtMin = cms.double(5.0),
    jetTag = cms.InputTag("akCs0PFpatJets"),
    matchJets = cms.untracked.bool(False),
    matchTag = cms.untracked.InputTag("akPu4PFpatJets"),
    pfCandidateLabel = cms.untracked.InputTag("packedPFCandidates"),
    rParam = cms.double(0.4),
    towersSrc = cms.InputTag("towerMaker"),
    trackQuality = cms.untracked.string('highPurity'),
    trackTag = cms.InputTag("hiTracks"),
    useHepMC = cms.untracked.bool(False),
    useQuality = cms.untracked.bool(True),
    useRawPt = cms.untracked.bool(False),
    vtxTag = cms.InputTag("hiSelectedVertex")
)


process.akFlowPuCs4PFJetAnalyzer = cms.EDAnalyzer("HiInclusiveJetAnalyzer",
    caloJetTag = cms.InputTag("slimmedCaloJets"),
    doCandidateBtagging = cms.untracked.bool(True),
    doGenSym = cms.untracked.bool(False),
    doGenTaus = cms.untracked.bool(False),
    doHiJetID = cms.untracked.bool(False),
    doJetConstituents = cms.untracked.bool(False),
    doLegacyBtagging = cms.untracked.bool(False),
    doStandardJetID = cms.untracked.bool(False),
    doSubEvent = cms.untracked.bool(False),
    doSubJets = cms.untracked.bool(False),
    doWTARecluster = cms.untracked.bool(False),
    eventInfoTag = cms.InputTag("generator"),
    fillGenJets = cms.untracked.bool(False),
    genjetTag = cms.InputTag("ak4HiGenJets"),
    hltTrgResults = cms.untracked.string('TriggerResults::HISIGNAL'),
    isMC = cms.untracked.bool(False),
    jetName = cms.untracked.string('akCs4PF'),
    jetPtMin = cms.double(5.0),
    jetTag = cms.InputTag("slimmedJets"),
    matchJets = cms.untracked.bool(False),
    matchTag = cms.untracked.InputTag("akPu4PFpatJets"),
    pfCandidateLabel = cms.untracked.InputTag("packedPFCandidates"),
    rParam = cms.double(0.4),
    towersSrc = cms.InputTag("towerMaker"),
    trackQuality = cms.untracked.string('highPurity'),
    trackTag = cms.InputTag("hiTracks"),
    useHepMC = cms.untracked.bool(False),
    useQuality = cms.untracked.bool(True),
    useRawPt = cms.untracked.bool(False),
    vtxTag = cms.InputTag("hiSelectedVertex")
)


process.ggHiNtuplizer = cms.EDAnalyzer("ggHiNtuplizer",
    beamSpotSrc = cms.InputTag("offlineBeamSpot"),
    conversionsSrc = cms.InputTag("reducedEgamma","reducedConversions"),
    doEffectiveAreas = cms.bool(True),
    doElectrons = cms.bool(True),
    doGenParticles = cms.bool(False),
    doMuons = cms.bool(False),
    doPackedGenParticle = cms.bool(True),
    doPfIso = cms.bool(True),
    doPhoERegression = cms.bool(True),
    doPhotons = cms.bool(True),
    doRecHitsEB = cms.bool(False),
    doRecHitsEE = cms.bool(False),
    doSuperClusters = cms.bool(False),
    effAreasConfigFile = cms.FileInPath('HeavyIonsAnalysis/EGMAnalysis/data/EffectiveAreas_94X_v0'),
    electronSrc = cms.InputTag("slimmedElectrons"),
    genParticleSrc = cms.InputTag("packedGenParticles"),
    isPackedPFCandidate = cms.bool(True),
    isParticleGun = cms.bool(False),
    muonSrc = cms.InputTag("unpackedMuons"),
    particleFlowCollection = cms.InputTag("packedPFCandidates"),
    photonSrc = cms.InputTag("slimmedPhotons"),
    pileupSrc = cms.InputTag("slimmedAddPileupInfo"),
    recHitsEB = cms.untracked.InputTag("reducedEgamma","reducedEBRecHits"),
    recHitsEE = cms.untracked.InputTag("reducedEgamma","reducedEERecHits"),
    rhoSrc = cms.InputTag("fixedGridRhoFastjetAll"),
    signalGenParticleSrc = cms.InputTag("packedGenParticlesSignal"),
    superClusters = cms.InputTag("reducedEgamma","reducedSuperClusters"),
    useValMapIso = cms.bool(True),
    vertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.hiEvtAnalyzer = cms.EDAnalyzer("HiEvtAnalyzer",
    CentralityBinSrc = cms.InputTag("centralityBin","HFtowers"),
    CentralitySrc = cms.InputTag("hiCentrality"),
    EvtPlane = cms.InputTag("hiEvtPlane"),
    EvtPlaneFlat = cms.InputTag("hiEvtPlaneFlat"),
    HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
    HiMC = cms.InputTag("heavyIon"),
    Vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    doCentrality = cms.bool(True),
    doEvtPlane = cms.bool(True),
    doEvtPlaneFlat = cms.bool(True),
    doHFfilters = cms.bool(True),
    doHiMC = cms.bool(False),
    doMC = cms.bool(False),
    doVertex = cms.bool(True),
    evtPlaneLevel = cms.int32(0),
    pfCandidateSrc = cms.InputTag("packedPFCandidates"),
    useHepMC = cms.bool(False)
)


process.hltanalysis = cms.EDAnalyzer("TriggerAnalyzer",
    HLTProcessName = cms.string('HLT'),
    hltPSProvCfg = cms.PSet(
        stageL1Trigger = cms.uint32(2)
    ),
    hltdummybranches = cms.vstring( (
        'DST_Physics_v7',
        'AlCa_EcalEtaEBonlyForHI_v1',
        'AlCa_EcalEtaEEonlyForHI_v1',
        'AlCa_EcalPhiSymForHI_v1',
        'AlCa_EcalPi0EBonlyForHI_v1',
        'AlCa_EcalPi0EEonlyForHI_v1',
        'AlCa_RPCMuonNormalisationForHI_v1',
        'HLT_EcalCalibration_v4',
        'HLT_HICastorMinimumBias_part0_v1',
        'HLT_HICastorMinimumBias_part10_v1',
        'HLT_HICastorMinimumBias_part11_v1',
        'HLT_HICastorMinimumBias_part12_v1',
        'HLT_HICastorMinimumBias_part13_v1',
        'HLT_HICastorMinimumBias_part14_v1',
        'HLT_HICastorMinimumBias_part15_v1',
        'HLT_HICastorMinimumBias_part16_v1',
        'HLT_HICastorMinimumBias_part17_v1',
        'HLT_HICastorMinimumBias_part18_v1',
        'HLT_HICastorMinimumBias_part19_v1',
        'HLT_HICastorMinimumBias_part1_v1',
        'HLT_HICastorMinimumBias_part2_v1',
        'HLT_HICastorMinimumBias_part3_v1',
        'HLT_HICastorMinimumBias_part4_v1',
        'HLT_HICastorMinimumBias_part5_v1',
        'HLT_HICastorMinimumBias_part6_v1',
        'HLT_HICastorMinimumBias_part7_v1',
        'HLT_HICastorMinimumBias_part8_v1',
        'HLT_HICastorMinimumBias_part9_v1',
        'HLT_HICastor_HighJet_BptxAND_v1',
        'HLT_HICastor_HighJet_MBHF1AND_BptxAND_v1',
        'HLT_HICastor_HighJet_MBHF1OR_BptxAND_v1',
        'HLT_HICastor_HighJet_MBHF1OR_v1',
        'HLT_HICastor_HighJet_MBHF2AND_BptxAND_v1',
        'HLT_HICastor_HighJet_NotMBHF2AND_v1',
        'HLT_HICastor_HighJet_NotMBHF2OR_v1',
        'HLT_HICastor_HighJet_v1',
        'HLT_HICastor_MediumJet_BptxAND_v1',
        'HLT_HICastor_MediumJet_MBHF1OR_BptxAND_v1',
        'HLT_HICastor_MediumJet_MBHF1OR_v1',
        'HLT_HICastor_MediumJet_NotMBHF2AND_v1',
        'HLT_HICastor_MediumJet_NotMBHF2OR_v1',
        'HLT_HICastor_MediumJet_SingleEG5_MBHF1OR_BptxAND_v1',
        'HLT_HICastor_MediumJet_SingleEG5_MBHF1OR_v1',
        'HLT_HICastor_MediumJet_SingleMu0_MBHF1OR_BptxAND_v1',
        'HLT_HICastor_MediumJet_SingleMu0_MBHF1OR_v1',
        'HLT_HICastor_MediumJet_v1',
        'HLT_HICastor_Muon_BptxAND_v1',
        'HLT_HICastor_Muon_v1',
        'HLT_HICentrality30100_FirstCollisionAfterAbortGap_v1',
        'HLT_HICentralityTag20100_v1',
        'HLT_HICentralityTag30100_v1',
        'HLT_HICentralityTag50100_v1',
        'HLT_HICentralityVeto_Beamspot_v1',
        'HLT_HICentralityVeto_v1',
        'HLT_HICsAK4PFJet100Eta1p5_Beamspot_v1',
        'HLT_HICsAK4PFJet100Eta1p5_Centrality_30_100_v1',
        'HLT_HICsAK4PFJet100Eta1p5_Centrality_50_100_v1',
        'HLT_HICsAK4PFJet100Eta1p5_v1',
        'HLT_HICsAK4PFJet120Eta1p5_v1',
        'HLT_HICsAK4PFJet60Eta1p5_Centrality_30_100_v1',
        'HLT_HICsAK4PFJet60Eta1p5_Centrality_50_100_v1',
        'HLT_HICsAK4PFJet60Eta1p5_v1',
        'HLT_HICsAK4PFJet80Eta1p5_Centrality_30_100_v1',
        'HLT_HICsAK4PFJet80Eta1p5_Centrality_50_100_v1',
        'HLT_HICsAK4PFJet80Eta1p5_v1',
        'HLT_HIDmesonPPTrackingGlobal_Dpt15_v1',
        'HLT_HIDmesonPPTrackingGlobal_Dpt20_v1',
        'HLT_HIDmesonPPTrackingGlobal_Dpt30_v1',
        'HLT_HIDmesonPPTrackingGlobal_Dpt40_v1',
        'HLT_HIDmesonPPTrackingGlobal_Dpt50_v1',
        'HLT_HIDmesonPPTrackingGlobal_Dpt60_v1',
        'HLT_HIDoubleEle10GsfMass50_v1',
        'HLT_HIDoubleEle10Gsf_v1',
        'HLT_HIDoubleEle15GsfMass50_v1',
        'HLT_HIDoubleEle15Gsf_v1',
        'HLT_HIDsPPTrackingGlobal_Dpt20_v1',
        'HLT_HIDsPPTrackingGlobal_Dpt30_v1',
        'HLT_HIDsPPTrackingGlobal_Dpt40_v1',
        'HLT_HIDsPPTrackingGlobal_Dpt50_v1',
        'HLT_HIDsPPTrackingGlobal_Dpt60_v1',
        'HLT_HIEle10Gsf_PuAK4CaloJet100Eta2p1_v1',
        'HLT_HIEle10Gsf_PuAK4CaloJet40Eta2p1_v1',
        'HLT_HIEle10Gsf_PuAK4CaloJet60Eta2p1_v1',
        'HLT_HIEle10Gsf_PuAK4CaloJet80Eta2p1_v1',
        'HLT_HIEle10Gsf_v1',
        'HLT_HIEle15Ele10GsfMass50_v1',
        'HLT_HIEle15Ele10Gsf_v1',
        'HLT_HIEle15Gsf_PuAK4CaloJet100Eta2p1_v1',
        'HLT_HIEle15Gsf_PuAK4CaloJet40Eta2p1_v1',
        'HLT_HIEle15Gsf_PuAK4CaloJet60Eta2p1_v1',
        'HLT_HIEle15Gsf_PuAK4CaloJet80Eta2p1_v1',
        'HLT_HIEle15Gsf_v1',
        'HLT_HIEle20Gsf_PuAK4CaloJet100Eta2p1_v1',
        'HLT_HIEle20Gsf_PuAK4CaloJet40Eta2p1_v1',
        'HLT_HIEle20Gsf_PuAK4CaloJet60Eta2p1_v1',
        'HLT_HIEle20Gsf_PuAK4CaloJet80Eta2p1_v1',
        'HLT_HIEle20Gsf_v1',
        'HLT_HIEle30Gsf_v1',
        'HLT_HIEle40Gsf_v1',
        'HLT_HIEle50Gsf_v1',
        'HLT_HIFullTracks2018_HighPt18_v1',
        'HLT_HIFullTracks2018_HighPt24_v1',
        'HLT_HIFullTracks2018_HighPt34_v1',
        'HLT_HIFullTracks2018_HighPt45_v1',
        'HLT_HIFullTracks2018_HighPt56_v1',
        'HLT_HIFullTracks2018_HighPt60_v1',
        'HLT_HIFullTracks_Multiplicity020_HF1AND_v1',
        'HLT_HIFullTracks_Multiplicity020_HF1OR_v1',
        'HLT_HIFullTracks_Multiplicity020_HF2OR_v1',
        'HLT_HIFullTracks_Multiplicity020_v1',
        'HLT_HIFullTracks_Multiplicity2040_HF1AND_v1',
        'HLT_HIFullTracks_Multiplicity2040_HF1OR_v1',
        'HLT_HIFullTracks_Multiplicity2040_HF2OR_v1',
        'HLT_HIFullTracks_Multiplicity2040_v1',
        'HLT_HIFullTracks_Multiplicity335_HF1OR_v1',
        'HLT_HIFullTracks_Multiplicity4060_v1',
        'HLT_HIFullTracks_Multiplicity6080_v1',
        'HLT_HIFullTracks_Multiplicity80100_v1',
        'HLT_HIGEDPhoton10_Cent30_100_v1',
        'HLT_HIGEDPhoton10_Cent50_100_v1',
        'HLT_HIGEDPhoton10_EB_HECut_v1',
        'HLT_HIGEDPhoton10_EB_v1',
        'HLT_HIGEDPhoton10_HECut_v1',
        'HLT_HIGEDPhoton10_v1',
        'HLT_HIGEDPhoton20_Cent30_100_v1',
        'HLT_HIGEDPhoton20_Cent50_100_v1',
        'HLT_HIGEDPhoton20_EB_HECut_v1',
        'HLT_HIGEDPhoton20_EB_v1',
        'HLT_HIGEDPhoton20_HECut_v1',
        'HLT_HIGEDPhoton20_v1',
        'HLT_HIGEDPhoton30_Cent30_100_v1',
        'HLT_HIGEDPhoton30_Cent50_100_v1',
        'HLT_HIGEDPhoton30_EB_HECut_v1',
        'HLT_HIGEDPhoton30_EB_v1',
        'HLT_HIGEDPhoton30_HECut_v1',
        'HLT_HIGEDPhoton30_v1',
        'HLT_HIGEDPhoton40_Cent30_100_v1',
        'HLT_HIGEDPhoton40_Cent50_100_v1',
        'HLT_HIGEDPhoton40_EB_HECut_v1',
        'HLT_HIGEDPhoton40_EB_v1',
        'HLT_HIGEDPhoton40_HECut_v1',
        'HLT_HIGEDPhoton40_v1',
        'HLT_HIGEDPhoton50_EB_HECut_v1',
        'HLT_HIGEDPhoton50_EB_v1',
        'HLT_HIGEDPhoton50_HECut_v1',
        'HLT_HIGEDPhoton50_v1',
        'HLT_HIGEDPhoton60_EB_HECut_v1',
        'HLT_HIGEDPhoton60_EB_v1',
        'HLT_HIGEDPhoton60_HECut_v1',
        'HLT_HIGEDPhoton60_v1',
        'HLT_HIHcalNZS_v1',
        'HLT_HIHcalPhiSym_v1',
        'HLT_HIIslandPhoton10_Eta1p5_v1',
        'HLT_HIIslandPhoton10_Eta2p4_Cent30_100_v1',
        'HLT_HIIslandPhoton10_Eta2p4_Cent50_100_v1',
        'HLT_HIIslandPhoton10_Eta2p4_v1',
        'HLT_HIIslandPhoton10_Eta3p1_Cent30_100_v1',
        'HLT_HIIslandPhoton10_Eta3p1_Cent50_100_v1',
        'HLT_HIIslandPhoton10_Eta3p1_v1',
        'HLT_HIIslandPhoton20_Eta1p5_v1',
        'HLT_HIIslandPhoton20_Eta2p4_Cent30_100_v1',
        'HLT_HIIslandPhoton20_Eta2p4_Cent50_100_v1',
        'HLT_HIIslandPhoton20_Eta2p4_v1',
        'HLT_HIIslandPhoton20_Eta3p1_Cent30_100_v1',
        'HLT_HIIslandPhoton20_Eta3p1_Cent50_100_v1',
        'HLT_HIIslandPhoton20_Eta3p1_v1',
        'HLT_HIIslandPhoton30_Eta1p5_v1',
        'HLT_HIIslandPhoton30_Eta2p4_Cent30_100_v1',
        'HLT_HIIslandPhoton30_Eta2p4_Cent50_100_v1',
        'HLT_HIIslandPhoton30_Eta2p4_v1',
        'HLT_HIIslandPhoton30_Eta3p1_Cent30_100_v1',
        'HLT_HIIslandPhoton30_Eta3p1_Cent50_100_v1',
        'HLT_HIIslandPhoton30_Eta3p1_v1',
        'HLT_HIIslandPhoton40_Eta1p5_v1',
        'HLT_HIIslandPhoton40_Eta2p4_Cent30_100_v1',
        'HLT_HIIslandPhoton40_Eta2p4_Cent50_100_v1',
        'HLT_HIIslandPhoton40_Eta2p4_v1',
        'HLT_HIIslandPhoton40_Eta3p1_Cent30_100_v1',
        'HLT_HIIslandPhoton40_Eta3p1_Cent50_100_v1',
        'HLT_HIIslandPhoton40_Eta3p1_v1',
        'HLT_HIIslandPhoton50_Eta1p5_v1',
        'HLT_HIIslandPhoton50_Eta2p4_v1',
        'HLT_HIIslandPhoton50_Eta3p1_v1',
        'HLT_HIIslandPhoton60_Eta1p5_v1',
        'HLT_HIIslandPhoton60_Eta2p4_v1',
        'HLT_HIIslandPhoton60_Eta3p1_v1',
        'HLT_HIL1DoubleMu0_v1',
        'HLT_HIL1DoubleMu10_v1',
        'HLT_HIL1DoubleMuOpen_Centrality_30_100_v1',
        'HLT_HIL1DoubleMuOpen_Centrality_40_100_v1',
        'HLT_HIL1DoubleMuOpen_Centrality_50_100_v1',
        'HLT_HIL1DoubleMuOpen_MAXdR2p0_v1',
        'HLT_HIL1DoubleMuOpen_MAXdR3p5_v1',
        'HLT_HIL1DoubleMuOpen_OS_Centrality_30_100_v1',
        'HLT_HIL1DoubleMuOpen_OS_Centrality_40_100_v1',
        'HLT_HIL1DoubleMuOpen_OS_MAXdR2p0_v1',
        'HLT_HIL1DoubleMuOpen_OS_er1p6_v1',
        'HLT_HIL1DoubleMuOpen_OS_v1',
        'HLT_HIL1DoubleMuOpen_SS_v1',
        'HLT_HIL1DoubleMuOpen_er1p6_v1',
        'HLT_HIL1DoubleMuOpen_v1',
        'HLT_HIL1Mu3Eta2p5_Ele10Gsf_v1',
        'HLT_HIL1Mu3Eta2p5_Ele15Gsf_v1',
        'HLT_HIL1Mu3Eta2p5_Ele20Gsf_v1',
        'HLT_HIL1Mu3_Centrality_70_100_v1',
        'HLT_HIL1Mu3_Centrality_80_100_v1',
        'HLT_HIL1Mu3_v1',
        'HLT_HIL1Mu5Eta2p5_Ele10Gsf_v1',
        'HLT_HIL1Mu5Eta2p5_Ele15Gsf_v1',
        'HLT_HIL1Mu5Eta2p5_Ele20Gsf_v1',
        'HLT_HIL1Mu5_v1',
        'HLT_HIL1Mu7Eta2p5_Ele10Gsf_v1',
        'HLT_HIL1Mu7Eta2p5_Ele15Gsf_v1',
        'HLT_HIL1Mu7Eta2p5_Ele20Gsf_v1',
        'HLT_HIL1Mu7_v1',
        'HLT_HIL1MuOpen_Centrality_70_100_v1',
        'HLT_HIL1MuOpen_Centrality_80_100_v1',
        'HLT_HIL1NotBptxOR_v1',
        'HLT_HIL1UnpairedBunchBptxMinus_v1',
        'HLT_HIL1UnpairedBunchBptxPlus_v1',
        'HLT_HIL1_ETT10_ETTAsym50_MinimumBiasHF1_OR_BptxAND_v1',
        'HLT_HIL1_ETT60_ETTAsym65_MinimumBiasHF2_OR_PixelTracks10_v1',
        'HLT_HIL1_ETT65_ETTAsym80_MinimumBiasHF2_OR_PixelTracks10_v1',
        'HLT_HIL1_ETT8_ETTAsym50_MinimumBiasHF1_OR_BptxAND_v1',
        'HLT_HIL1_ZDC_AND_OR_MinimumBiasHF1_AND_BptxAND_v1',
        'HLT_HIL1_ZDC_AND_OR_MinimumBiasHF2_AND_BptxAND_v1',
        'HLT_HIL2DoubleMuOpen_v1',
        'HLT_HIL2DoubleMu_L1DoubleMuOpen_OS_MAXdR2p0_v1',
        'HLT_HIL2Mu12_v1',
        'HLT_HIL2Mu15_v1',
        'HLT_HIL2Mu20_v1',
        'HLT_HIL2Mu3_NHitQ10_v1',
        'HLT_HIL2Mu3_NHitQ15_tagging_v1',
        'HLT_HIL2Mu3_NHitQ15_v1',
        'HLT_HIL2Mu3_v1',
        'HLT_HIL2Mu5_NHitQ10_v1',
        'HLT_HIL2Mu5_NHitQ15_tagging_v1',
        'HLT_HIL2Mu5_NHitQ15_v1',
        'HLT_HIL2Mu5_v1',
        'HLT_HIL2Mu7_NHitQ10_v1',
        'HLT_HIL2Mu7_NHitQ15_tagging_v1',
        'HLT_HIL2Mu7_NHitQ15_v1',
        'HLT_HIL2Mu7_v1',
        'HLT_HIL2_L1DoubleMu10_v1',
        'HLT_HIL3DoubleMuOpen_JpsiPsi_v1',
        'HLT_HIL3DoubleMuOpen_M60120_v1',
        'HLT_HIL3DoubleMuOpen_Upsi_v1',
        'HLT_HIL3DoubleMuOpen_v1',
        'HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1',
        'HLT_HIL3Mu0_L2Mu0_v1',
        'HLT_HIL3Mu12_v1',
        'HLT_HIL3Mu15_v1',
        'HLT_HIL3Mu20_v1',
        'HLT_HIL3Mu2p5NHitQ10_L2Mu2_M7toinf_v1',
        'HLT_HIL3Mu2p5NHitQ10_L2Mu2_v1',
        'HLT_HIL3Mu2p5_L1DoubleMu0_v1',
        'HLT_HIL3Mu3Eta2p5_PuAK4CaloJet100Eta2p1_FilterDr_v1',
        'HLT_HIL3Mu3Eta2p5_PuAK4CaloJet100Eta2p1_v1',
        'HLT_HIL3Mu3Eta2p5_PuAK4CaloJet40Eta2p1_FilterDr_v1',
        'HLT_HIL3Mu3Eta2p5_PuAK4CaloJet40Eta2p1_v1',
        'HLT_HIL3Mu3Eta2p5_PuAK4CaloJet60Eta2p1_FilterDr_v1',
        'HLT_HIL3Mu3Eta2p5_PuAK4CaloJet60Eta2p1_v1',
        'HLT_HIL3Mu3Eta2p5_PuAK4CaloJet80Eta2p1_FilterDr_v1',
        'HLT_HIL3Mu3Eta2p5_PuAK4CaloJet80Eta2p1_v1',
        'HLT_HIL3Mu3NHitQ10_L1DoubleMuOpen_v1',
        'HLT_HIL3Mu3_EG10HECut_v1',
        'HLT_HIL3Mu3_EG15HECut_v1',
        'HLT_HIL3Mu3_EG20HECut_v1',
        'HLT_HIL3Mu3_EG30HECut_v1',
        'HLT_HIL3Mu3_L1DoubleMu0_v1',
        'HLT_HIL3Mu3_L1DoubleMuOpen_OS_v1',
        'HLT_HIL3Mu3_L1DoubleMuOpen_v1',
        'HLT_HIL3Mu3_L1TripleMuOpen_v1',
        'HLT_HIL3Mu3_NHitQ10_tagging_v1',
        'HLT_HIL3Mu3_NHitQ10_v1',
        'HLT_HIL3Mu3_v1',
        'HLT_HIL3Mu5Eta2p5_PuAK4CaloJet100Eta2p1_FilterDr_v1',
        'HLT_HIL3Mu5Eta2p5_PuAK4CaloJet100Eta2p1_v1',
        'HLT_HIL3Mu5Eta2p5_PuAK4CaloJet40Eta2p1_FilterDr_v1',
        'HLT_HIL3Mu5Eta2p5_PuAK4CaloJet40Eta2p1_v1',
        'HLT_HIL3Mu5Eta2p5_PuAK4CaloJet60Eta2p1_FilterDr_v1',
        'HLT_HIL3Mu5Eta2p5_PuAK4CaloJet60Eta2p1_v1',
        'HLT_HIL3Mu5Eta2p5_PuAK4CaloJet80Eta2p1_FilterDr_v1',
        'HLT_HIL3Mu5Eta2p5_PuAK4CaloJet80Eta2p1_v1',
        'HLT_HIL3Mu5_EG10HECut_v1',
        'HLT_HIL3Mu5_EG15HECut_v1',
        'HLT_HIL3Mu5_EG20HECut_v1',
        'HLT_HIL3Mu5_EG30HECut_v1',
        'HLT_HIL3Mu5_NHitQ10_tagging_v1',
        'HLT_HIL3Mu5_NHitQ10_v1',
        'HLT_HIL3Mu5_v1',
        'HLT_HIL3Mu7_EG10HECut_v1',
        'HLT_HIL3Mu7_EG15HECut_v1',
        'HLT_HIL3Mu7_EG20HECut_v1',
        'HLT_HIL3Mu7_EG30HECut_v1',
        'HLT_HIL3Mu7_NHitQ10_tagging_v1',
        'HLT_HIL3Mu7_NHitQ10_v1',
        'HLT_HIL3Mu7_v1',
        'HLT_HIL3_L1DoubleMu10_v1',
        'HLT_HILcPPTrackingGlobal_Dpt20_v1',
        'HLT_HILcPPTrackingGlobal_Dpt30_v1',
        'HLT_HILcPPTrackingGlobal_Dpt40_v1',
        'HLT_HILcPPTrackingGlobal_Dpt50_v1',
        'HLT_HILcPPTrackingGlobal_Dpt60_v1',
        'HLT_HIMinimumBiasHFOR_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part0_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part10_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part11_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part12_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part13_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part14_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part15_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part16_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part17_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part18_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part19_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part1_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part20_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part21_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part22_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part23_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part2_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part3_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part4_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part5_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part6_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part7_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part8_v1',
        'HLT_HIMinimumBiasRF_SinglePixelTrack_part9_v1',
        'HLT_HIMinimumBiasRF_part0_v1',
        'HLT_HIMinimumBiasRF_part10_v1',
        'HLT_HIMinimumBiasRF_part11_v1',
        'HLT_HIMinimumBiasRF_part12_v1',
        'HLT_HIMinimumBiasRF_part13_v1',
        'HLT_HIMinimumBiasRF_part14_v1',
        'HLT_HIMinimumBiasRF_part15_v1',
        'HLT_HIMinimumBiasRF_part16_v1',
        'HLT_HIMinimumBiasRF_part17_v1',
        'HLT_HIMinimumBiasRF_part18_v1',
        'HLT_HIMinimumBiasRF_part19_v1',
        'HLT_HIMinimumBiasRF_part1_v1',
        'HLT_HIMinimumBiasRF_part20_v1',
        'HLT_HIMinimumBiasRF_part21_v1',
        'HLT_HIMinimumBiasRF_part22_v1',
        'HLT_HIMinimumBiasRF_part23_v1',
        'HLT_HIMinimumBiasRF_part2_v1',
        'HLT_HIMinimumBiasRF_part3_v1',
        'HLT_HIMinimumBiasRF_part4_v1',
        'HLT_HIMinimumBiasRF_part5_v1',
        'HLT_HIMinimumBiasRF_part6_v1',
        'HLT_HIMinimumBiasRF_part7_v1',
        'HLT_HIMinimumBiasRF_part8_v1',
        'HLT_HIMinimumBiasRF_part9_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part0_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part10_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part11_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part12_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part13_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part14_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part15_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part16_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part17_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part18_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part19_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part1_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part2_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part3_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part4_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part5_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part6_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part7_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part8_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part9_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part0_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part10_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part11_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part12_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part13_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part14_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part15_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part16_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part17_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part18_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part19_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part1_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part2_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part3_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part4_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part5_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part6_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part7_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part8_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part9_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_part0_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_part10_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_part11_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_part12_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_part13_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_part14_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_part15_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_part16_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_part17_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_part18_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_part19_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_part1_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_part2_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_part3_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_part4_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_part5_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_part6_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_part7_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_part8_v1',
        'HLT_HIMinimumBias_SinglePixelTrack_part9_v1',
        'HLT_HIMinimumBias_part0_v1',
        'HLT_HIMinimumBias_part10_v1',
        'HLT_HIMinimumBias_part10_v6',
        'HLT_HIMinimumBias_part11_v1',
        'HLT_HIMinimumBias_part11_v6',
        'HLT_HIMinimumBias_part12_v1',
        'HLT_HIMinimumBias_part12_v7',
        'HLT_HIMinimumBias_part13_v1',
        'HLT_HIMinimumBias_part13_v7',
        'HLT_HIMinimumBias_part14_v1',
        'HLT_HIMinimumBias_part14_v8',
        'HLT_HIMinimumBias_part15_v1',
        'HLT_HIMinimumBias_part15_v8',
        'HLT_HIMinimumBias_part16_v1',
        'HLT_HIMinimumBias_part16_v9',
        'HLT_HIMinimumBias_part17_v1',
        'HLT_HIMinimumBias_part17_v9',
        'HLT_HIMinimumBias_part18_v1',
        'HLT_HIMinimumBias_part18_v10',
        'HLT_HIMinimumBias_part18_v11',
        'HLT_HIMinimumBias_part19_v1',
        'HLT_HIMinimumBias_part19_v10',
        'HLT_HIMinimumBias_part19_v11',
        'HLT_HIMinimumBias_part1_v1',
        'HLT_HIMinimumBias_part20_v1',
        'HLT_HIMinimumBias_part21_v1',
        'HLT_HIMinimumBias_part22_v1',
        'HLT_HIMinimumBias_part23_v1',
        'HLT_HIMinimumBias_part24_v1',
        'HLT_HIMinimumBias_part25_v1',
        'HLT_HIMinimumBias_part2_v1',
        'HLT_HIMinimumBias_part2_v2',
        'HLT_HIMinimumBias_part3_v1',
        'HLT_HIMinimumBias_part3_v2',
        'HLT_HIMinimumBias_part4_v1',
        'HLT_HIMinimumBias_part4_v3',
        'HLT_HIMinimumBias_part5_v1',
        'HLT_HIMinimumBias_part5_v3',
        'HLT_HIMinimumBias_part6_v1',
        'HLT_HIMinimumBias_part6_v4',
        'HLT_HIMinimumBias_part7_v1',
        'HLT_HIMinimumBias_part7_v4',
        'HLT_HIMinimumBias_part8_v1',
        'HLT_HIMinimumBias_part8_v5',
        'HLT_HIMinimumBias_part9_v1',
        'HLT_HIMinimumBias_part9_v5',
        'HLT_HIPhysicsForZS_v1',
        'HLT_HIPhysics_v1',
        'HLT_HIPuAK4CaloJet100Eta2p4_CSVv2WP0p75_v1',
        'HLT_HIPuAK4CaloJet100Eta2p4_CSVv2WP0p8_v1',
        'HLT_HIPuAK4CaloJet100Eta2p4_DeepCSV0p4_v1',
        'HLT_HIPuAK4CaloJet100Eta2p4_DeepCSV0p71_v1',
        'HLT_HIPuAK4CaloJet100Eta5p1_Centrality_30_100_v1',
        'HLT_HIPuAK4CaloJet100Eta5p1_Centrality_50_100_v1',
        'HLT_HIPuAK4CaloJet100Eta5p1_v1',
        'HLT_HIPuAK4CaloJet100Fwd_v1',
        'HLT_HIPuAK4CaloJet100_35_Eta0p7_v1',
        'HLT_HIPuAK4CaloJet100_35_Eta1p1_v1',
        'HLT_HIPuAK4CaloJet120Eta5p1_v1',
        'HLT_HIPuAK4CaloJet120Fwd_v1',
        'HLT_HIPuAK4CaloJet40Eta5p1_Centrality_30_100_v1',
        'HLT_HIPuAK4CaloJet40Eta5p1_Centrality_50_100_v1',
        'HLT_HIPuAK4CaloJet40Eta5p1_v1',
        'HLT_HIPuAK4CaloJet40Fwd_v1',
        'HLT_HIPuAK4CaloJet60Eta2p4_CSVv2WP0p75_v1',
        'HLT_HIPuAK4CaloJet60Eta2p4_CSVv2WP0p8_v1',
        'HLT_HIPuAK4CaloJet60Eta2p4_DeepCSV0p4_v1',
        'HLT_HIPuAK4CaloJet60Eta2p4_DeepCSV0p71_v1',
        'HLT_HIPuAK4CaloJet60Eta5p1_Centrality_30_100_v1',
        'HLT_HIPuAK4CaloJet60Eta5p1_Centrality_50_100_v1',
        'HLT_HIPuAK4CaloJet60Eta5p1_v1',
        'HLT_HIPuAK4CaloJet60Fwd_v1',
        'HLT_HIPuAK4CaloJet80Eta2p4_CSVv2WP0p75_v1',
        'HLT_HIPuAK4CaloJet80Eta2p4_CSVv2WP0p8_v1',
        'HLT_HIPuAK4CaloJet80Eta2p4_DeepCSV0p4_v1',
        'HLT_HIPuAK4CaloJet80Eta2p4_DeepCSV0p71_v1',
        'HLT_HIPuAK4CaloJet80Eta5p1_Centrality_30_100_v1',
        'HLT_HIPuAK4CaloJet80Eta5p1_Centrality_50_100_v1',
        'HLT_HIPuAK4CaloJet80Eta5p1_v1',
        'HLT_HIPuAK4CaloJet80Fwd_v1',
        'HLT_HIPuAK4CaloJet80_35_Eta0p7_v1',
        'HLT_HIPuAK4CaloJet80_35_Eta1p1_v1',
        'HLT_HIPuAK4CaloJet80_45_45_Eta2p1_v1',
        'HLT_HIRandom_v1',
        'HLT_HIUPC_DoubleEG2_BptxAND_MaxPixelTrack_v1',
        'HLT_HIUPC_DoubleEG2_BptxAND_SinglePixelTrack_MaxPixelTrack_v1',
        'HLT_HIUPC_DoubleEG2_NotMBHF2AND_SinglePixelTrack_MaxPixelTrack_v1',
        'HLT_HIUPC_DoubleEG2_NotMBHF2AND_SinglePixelTrack_v1',
        'HLT_HIUPC_DoubleEG2_NotMBHF2AND_v1',
        'HLT_HIUPC_DoubleEG2_NotMBHF2OR_SinglePixelTrack_MaxPixelTrack_v1',
        'HLT_HIUPC_DoubleEG2_NotMBHF2OR_SinglePixelTrack_v1',
        'HLT_HIUPC_DoubleEG2_NotMBHF2OR_v1',
        'HLT_HIUPC_DoubleEG5_BptxAND_MaxPixelTrack_v1',
        'HLT_HIUPC_DoubleEG5_BptxAND_SinglePixelTrack_MaxPixelTrack_v1',
        'HLT_HIUPC_DoubleEG5_NotMBHF2AND_SinglePixelTrack_MaxPixelTrack_v1',
        'HLT_HIUPC_DoubleEG5_NotMBHF2AND_SinglePixelTrack_v1',
        'HLT_HIUPC_DoubleEG5_NotMBHF2AND_v1',
        'HLT_HIUPC_DoubleEG5_NotMBHF2OR_SinglePixelTrack_MaxPixelTrack_v1',
        'HLT_HIUPC_DoubleEG5_NotMBHF2OR_SinglePixelTrack_v1',
        'HLT_HIUPC_DoubleEG5_NotMBHF2OR_v1',
        'HLT_HIUPC_DoubleMu0_BptxAND_MaxPixelTrack_v1',
        'HLT_HIUPC_DoubleMu0_NotMBHF2AND_MaxPixelTrack_v1',
        'HLT_HIUPC_DoubleMu0_NotMBHF2AND_v1',
        'HLT_HIUPC_DoubleMu0_NotMBHF2OR_MaxPixelTrack_v1',
        'HLT_HIUPC_DoubleMu0_NotMBHF2OR_v1',
        'HLT_HIUPC_DoubleMuOpen_BptxAND_MaxPixelTrack_v1',
        'HLT_HIUPC_DoubleMuOpen_NotMBHF2OR_MaxPixelTrack_v1',
        'HLT_HIUPC_DoubleMuOpen_NotMBHF2OR_v1',
        'HLT_HIUPC_ETT5_Asym50_NotMBHF2OR_SinglePixelTrack_v1',
        'HLT_HIUPC_ETT5_Asym50_NotMBHF2OR_v1',
        'HLT_HIUPC_Mu8_Mu13_MaxPixelTrack_v1',
        'HLT_HIUPC_Mu8_Mu13_v1',
        'HLT_HIUPC_NotMBHF2AND_SinglePixelTrack_MaxPixelTrack_v1',
        'HLT_HIUPC_NotMBHF2AND_SinglePixelTrack_v1',
        'HLT_HIUPC_NotMBHF2AND_v1',
        'HLT_HIUPC_NotMBHF2OR_BptxAND_SinglePixelTrack_MaxPixelTrack_v1',
        'HLT_HIUPC_NotMBHF2OR_BptxAND_SinglePixelTrack_v1',
        'HLT_HIUPC_SingleEG3_BptxAND_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleEG3_BptxAND_SinglePixelTrack_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleEG3_NotMBHF2AND_SinglePixelTrack_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleEG3_NotMBHF2AND_SinglePixelTrack_v1',
        'HLT_HIUPC_SingleEG3_NotMBHF2AND_v1',
        'HLT_HIUPC_SingleEG3_NotMBHF2OR_SinglePixelTrack_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleEG3_NotMBHF2OR_SinglePixelTrack_v1',
        'HLT_HIUPC_SingleEG3_NotMBHF2OR_v1',
        'HLT_HIUPC_SingleEG5_BptxAND_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleEG5_BptxAND_SinglePixelTrack_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleEG5_NotMBHF2AND_SinglePixelTrack_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleEG5_NotMBHF2AND_SinglePixelTrack_v1',
        'HLT_HIUPC_SingleEG5_NotMBHF2AND_v1',
        'HLT_HIUPC_SingleEG5_NotMBHF2OR_SinglePixelTrack_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleEG5_NotMBHF2OR_SinglePixelTrack_v1',
        'HLT_HIUPC_SingleEG5_NotMBHF2OR_v1',
        'HLT_HIUPC_SingleMu0_BptxAND_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleMu0_NotMBHF2AND_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleMu0_NotMBHF2AND_v1',
        'HLT_HIUPC_SingleMu0_NotMBHF2OR_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleMu0_NotMBHF2OR_v1',
        'HLT_HIUPC_SingleMu3_BptxAND_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleMu3_NotMBHF2OR_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleMu3_NotMBHF2OR_v1',
        'HLT_HIUPC_SingleMuOpen_BptxAND_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleMuOpen_NotMBHF2AND_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleMuOpen_NotMBHF2AND_v1',
        'HLT_HIUPC_SingleMuOpen_NotMBHF2OR_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleMuOpen_NotMBHF2OR_v1',
        'HLT_HIUPC_ZeroBias_MaxPixelCluster_v1',
        'HLT_HIUPC_ZeroBias_SinglePixelTrack_v1',
        'HLT_HIZeroBias_FirstCollisionAfterAbortGap_v1',
        'HLT_HIZeroBias_v1',
        'HLT_HI_CastorMuon_v1',
        'HLT_HMinimumBias_part12_v1',
        'HLT_HMinimumBias_part18_v1',
        'HLT_HcalCalibration_v5',
        'HLT_Mu8_Mu13_MaxPixelTrack_v1',
        'HLT_Mu8_Mu13_v1',
        'HLT_ZDCAND_MBHF1AND_BptxAND_v1',
        'HLT_ZDCAND_MBHF1OR_BptxAND_v1',
        'HLT_ZDCAND_MBHF2AND_BptxAND_v1',
        'HLT_ZDCAND_MBHF2OR_BptxAND_v1',
        'HLT_ZDCM_BptxAND_v1',
        'HLT_ZDCM_ZDCP_BptxAND_v1',
        'HLT_ZDCM_ZDCP_MBHF1AND_BptxAND_v1',
        'HLT_ZDCM_v1',
        'HLT_ZDCOR_MBHF1OR_BptxAND_v1',
        'HLT_ZDCOR_MBHF2OR_BptxAND_v1',
        'HLT_ZDCP_BptxAND_v1',
        'HLT_ZDCP_v1'
     ) ),
    hltresults = cms.InputTag("TriggerResults","","HLT"),
    l1dummybranches = cms.vstring( (
        'L1_AlwaysTrue',
        'L1_BPTX_AND_Ref1_VME',
        'L1_BPTX_AND_Ref3_VME',
        'L1_BPTX_AND_Ref4_VME',
        'L1_BPTX_BeamGas_B1_VME',
        'L1_BPTX_BeamGas_B2_VME',
        'L1_BPTX_BeamGas_Ref1_VME',
        'L1_BPTX_BeamGas_Ref2_VME',
        'L1_BPTX_NotOR_VME',
        'L1_BPTX_OR_Ref3_VME',
        'L1_BPTX_OR_Ref4_VME',
        'L1_BPTX_RefAND_VME',
        'L1_BptxMinus',
        'L1_BptxMinus_NotBptxPlus',
        'L1_BptxOR',
        'L1_BptxPlus',
        'L1_BptxPlus_NotBptxMinus',
        'L1_BptxXOR',
        'L1_Castor1',
        'L1_CastorHighJet',
        'L1_CastorHighJet_BptxAND',
        'L1_CastorHighJet_MinimumBiasHF1_OR_BptxAND',
        'L1_CastorHighJet_NotMinimumBiasHF2_AND_BptxAND',
        'L1_CastorHighJet_NotMinimumBiasHF2_OR_BptxAND',
        'L1_CastorHighJet_OR_MinimumBiasHF1_AND_BptxAND',
        'L1_CastorHighJet_OR_MinimumBiasHF2_AND_BptxAND',
        'L1_CastorMediumJet',
        'L1_CastorMediumJet_BptxAND',
        'L1_CastorMediumJet_MinimumBiasHF1_OR_BptxAND',
        'L1_CastorMediumJet_NotMinimumBiasHF2_AND_BptxAND',
        'L1_CastorMediumJet_NotMinimumBiasHF2_OR_BptxAND',
        'L1_CastorMediumJet_SingleEG5_MinimumBiasHF1_OR_BptxAND',
        'L1_CastorMediumJet_SingleMu0_MinimumBiasHF1_OR_BptxAND',
        'L1_CastorMuon',
        'L1_CastorMuon_BptxAND',
        'L1_Centrality_20_100_MinimumBiasHF1_AND_BptxAND',
        'L1_Centrality_30_100',
        'L1_Centrality_30_100_MinimumBiasHF1_AND_BptxAND',
        'L1_Centrality_50_100',
        'L1_Centrality_Saturation',
        'L1_DoubleEG10_BptxAND',
        'L1_DoubleEG2',
        'L1_DoubleEG2_BptxAND',
        'L1_DoubleEG2_NotMinimumBiasHF2_AND_BptxAND',
        'L1_DoubleEG2_NotMinimumBiasHF2_OR_BptxAND',
        'L1_DoubleEG5',
        'L1_DoubleEG5_BptxAND',
        'L1_DoubleEG5_NotMinimumBiasHF2_AND_BptxAND',
        'L1_DoubleEG5_NotMinimumBiasHF2_OR_BptxAND',
        'L1_DoubleEG8_BptxAND',
        'L1_DoubleJet16And12_MidEta2p7_BptxAND',
        'L1_DoubleJet16And12_MidEta2p7_Centrality_30_100_BptxAND',
        'L1_DoubleJet16And12_MidEta2p7_Centrality_50_100_BptxAND',
        'L1_DoubleJet16And8_MidEta2p7_BptxAND',
        'L1_DoubleJet16And8_MidEta2p7_Centrality_30_100_BptxAND',
        'L1_DoubleJet16And8_MidEta2p7_Centrality_50_100_BptxAND',
        'L1_DoubleJet20And12_MidEta2p7_BptxAND',
        'L1_DoubleJet20And12_MidEta2p7_Centrality_30_100_BptxAND',
        'L1_DoubleJet20And12_MidEta2p7_Centrality_50_100_BptxAND',
        'L1_DoubleJet20And8_MidEta2p7_BptxAND',
        'L1_DoubleJet20And8_MidEta2p7_Centrality_30_100_BptxAND',
        'L1_DoubleJet20And8_MidEta2p7_Centrality_50_100_BptxAND',
        'L1_DoubleJet28And16_MidEta2p7_BptxAND',
        'L1_DoubleJet28And16_MidEta2p7_Centrality_30_100_BptxAND',
        'L1_DoubleJet28And16_MidEta2p7_Centrality_50_100_BptxAND',
        'L1_DoubleMu0',
        'L1_DoubleMu0_BptxAND',
        'L1_DoubleMu0_Centrality_10_100_MinimumBiasHF1_AND_BptxAND',
        'L1_DoubleMu0_Centrality_30_100_MinimumBiasHF1_AND_BptxAND',
        'L1_DoubleMu0_Centrality_50_100_MinimumBiasHF1_AND_BptxAND',
        'L1_DoubleMu0_Mass_Min1',
        'L1_DoubleMu0_MinimumBiasHF1_AND_BptxAND',
        'L1_DoubleMu0_NotMinimumBiasHF2_AND_BptxAND',
        'L1_DoubleMu0_NotMinimumBiasHF2_OR_BptxAND',
        'L1_DoubleMu0_SQ',
        'L1_DoubleMu0_SQ_OS',
        'L1_DoubleMu10_BptxAND',
        'L1_DoubleMuOpen',
        'L1_DoubleMuOpen_BptxAND',
        'L1_DoubleMuOpen_Centrality_10_100_BptxAND',
        'L1_DoubleMuOpen_Centrality_30_100_BptxAND',
        'L1_DoubleMuOpen_Centrality_40_100_BptxAND',
        'L1_DoubleMuOpen_Centrality_50_100_BptxAND',
        'L1_DoubleMuOpen_MaxDr2p0_BptxAND',
        'L1_DoubleMuOpen_MaxDr2p0_OS_BptxAND',
        'L1_DoubleMuOpen_MaxDr3p5',
        'L1_DoubleMuOpen_MaxDr3p5_BptxAND',
        'L1_DoubleMuOpen_NotMinimumBiasHF2_AND_BptxAND',
        'L1_DoubleMuOpen_NotMinimumBiasHF2_OR_BptxAND',
        'L1_DoubleMuOpen_OS',
        'L1_DoubleMuOpen_OS_BptxAND',
        'L1_DoubleMuOpen_SS',
        'L1_DoubleMuOpen_SS_BptxAND',
        'L1_ETMHF100',
        'L1_ETMHF100_HTT60er',
        'L1_ETMHF120',
        'L1_ETMHF120_HTT60er',
        'L1_ETT10_ETTAsym50_MinimumBiasHF1_OR_BptxAND',
        'L1_ETT10_ETTAsym55_MinimumBiasHF1_OR_BptxAND',
        'L1_ETT10_ETTAsym60_MinimumBiasHF1_OR_BptxAND',
        'L1_ETT10_ETTAsym65_MinimumBiasHF1_OR_BptxAND',
        'L1_ETT10_ETTAsym70_MinimumBiasHF1_OR_BptxAND',
        'L1_ETT10_ETTAsym80_MinimumBiasHF1_OR_BptxAND',
        'L1_ETT1200',
        'L1_ETT1600',
        'L1_ETT2000',
        'L1_ETT35_NotETT80_BptxAND',
        'L1_ETT40_NotETT95_BptxAND',
        'L1_ETT45_NotETT110_BptxAND',
        'L1_ETT5',
        'L1_ETT50_ETTAsym40_MinimumBiasHF2_OR_BptxAND',
        'L1_ETT50_ETTAsym40_NotMinimumBiasHF2_OR_BptxAND',
        'L1_ETT50_ETTAsym50_MinimumBiasHF2_OR_BptxAND',
        'L1_ETT50_ETTAsym50_NotMinimumBiasHF2_OR_BptxAND',
        'L1_ETT50_ETTAsym55_MinimumBiasHF2_OR_BptxAND',
        'L1_ETT50_ETTAsym60_MinimumBiasHF2_OR_BptxAND',
        'L1_ETT50_ETTAsym60_NotMinimumBiasHF2_OR_BptxAND',
        'L1_ETT50_ETTAsym65_MinimumBiasHF2_OR_BptxAND',
        'L1_ETT50_ETTAsym70_MinimumBiasHF2_OR_BptxAND',
        'L1_ETT50_ETTAsym70_NotMinimumBiasHF2_OR_BptxAND',
        'L1_ETT50_ETTAsym80_MinimumBiasHF2_OR_BptxAND',
        'L1_ETT50_ETTAsym80_NotMinimumBiasHF2_OR_BptxAND',
        'L1_ETT50_NotETT120_BptxAND',
        'L1_ETT55_NotETT130_BptxAND',
        'L1_ETT5_BptxAND',
        'L1_ETT5_ETTAsym40_BptxAND',
        'L1_ETT5_ETTAsym40_MinimumBiasHF1_OR_BptxAND',
        'L1_ETT5_ETTAsym40_NotMinimumBiasHF2_OR_BptxAND',
        'L1_ETT5_ETTAsym50_BptxAND',
        'L1_ETT5_ETTAsym50_MinimumBiasHF1_OR_BptxAND',
        'L1_ETT5_ETTAsym50_NotMinimumBiasHF2_OR_BptxAND',
        'L1_ETT5_ETTAsym55_MinimumBiasHF1_OR_BptxAND',
        'L1_ETT5_ETTAsym60_BptxAND',
        'L1_ETT5_ETTAsym60_MinimumBiasHF1_OR_BptxAND',
        'L1_ETT5_ETTAsym60_NotMinimumBiasHF2_OR_BptxAND',
        'L1_ETT5_ETTAsym65_MinimumBiasHF1_OR_BptxAND',
        'L1_ETT5_ETTAsym70_BptxAND',
        'L1_ETT5_ETTAsym70_MinimumBiasHF1_OR_BptxAND',
        'L1_ETT5_ETTAsym70_NotMinimumBiasHF2_OR_BptxAND',
        'L1_ETT5_ETTAsym80_BptxAND',
        'L1_ETT5_ETTAsym80_MinimumBiasHF1_OR_BptxAND',
        'L1_ETT5_ETTAsym80_NotMinimumBiasHF2_OR_BptxAND',
        'L1_ETT5_MinimumBiasHF1_OR_BptxAND',
        'L1_ETT5_NotETT30_BptxAND',
        'L1_ETT5_NotMinimumBiasHF2_OR',
        'L1_ETT60_ETTAsym60_MinimumBiasHF2_OR_BptxAND',
        'L1_ETT60_ETTAsym65_MinimumBiasHF2_OR_BptxAND',
        'L1_ETT65_ETTAsym70_MinimumBiasHF2_OR_BptxAND',
        'L1_ETT65_ETTAsym80_MinimumBiasHF2_OR_BptxAND',
        'L1_ETT8_ETTAsym50_MinimumBiasHF1_OR_BptxAND',
        'L1_ETT8_ETTAsym55_MinimumBiasHF1_OR_BptxAND',
        'L1_ETT8_ETTAsym60_MinimumBiasHF1_OR_BptxAND',
        'L1_ETT8_ETTAsym65_MinimumBiasHF1_OR_BptxAND',
        'L1_ETT8_ETTAsym70_MinimumBiasHF1_OR_BptxAND',
        'L1_ETT8_ETTAsym80_MinimumBiasHF1_OR_BptxAND',
        'L1_ETTAsym40',
        'L1_ETTAsym40_BptxAND',
        'L1_ETTAsym40_MinimumBiasHF1_OR_BptxAND',
        'L1_ETTAsym40_NotMinimumBiasHF2_OR_BptxAND',
        'L1_ETTAsym50',
        'L1_ETTAsym50_BptxAND',
        'L1_ETTAsym50_MinimumBiasHF1_OR_BptxAND',
        'L1_ETTAsym50_NotMinimumBiasHF2_OR_BptxAND',
        'L1_ETTAsym60',
        'L1_ETTAsym60_BptxAND',
        'L1_ETTAsym60_MinimumBiasHF1_OR_BptxAND',
        'L1_ETTAsym60_NotMinimumBiasHF2_OR_BptxAND',
        'L1_ETTAsym70',
        'L1_ETTAsym70_BptxAND',
        'L1_ETTAsym70_MinimumBiasHF1_OR_BptxAND',
        'L1_ETTAsym70_NotMinimumBiasHF2_OR_BptxAND',
        'L1_ETTAsym80',
        'L1_ETTAsym80_BptxAND',
        'L1_ETTAsym80_MinimumBiasHF1_OR_BptxAND',
        'L1_ETTAsym80_NotMinimumBiasHF2_OR_BptxAND',
        'L1_FirstBunchAfterTrain',
        'L1_FirstBunchBeforeTrain',
        'L1_FirstBunchInTrain',
        'L1_FirstCollisionInOrbit',
        'L1_FirstCollisionInOrbit_Centrality30_100_BptxAND',
        'L1_FirstCollisionInTrain',
        'L1_HCAL_LaserMon_Trig',
        'L1_HCAL_LaserMon_Veto',
        'L1_HTT120er',
        'L1_HTT200er',
        'L1_HTT280er',
        'L1_HTT360er',
        'L1_HTT450er',
        'L1_IsolatedBunch',
        'L1_LastBunchInTrain',
        'L1_LastCollisionInTrain',
        'L1_MinimumBiasHF0_AND_BptxAND',
        'L1_MinimumBiasHF0_OR_BptxAND',
        'L1_MinimumBiasHF1_AND',
        'L1_MinimumBiasHF1_AND_BptxAND',
        'L1_MinimumBiasHF1_AND_OR_ETT10_BptxAND',
        'L1_MinimumBiasHF1_OR',
        'L1_MinimumBiasHF1_OR_BptxAND',
        'L1_MinimumBiasHF1_XOR_BptxAND',
        'L1_MinimumBiasHF2_AND',
        'L1_MinimumBiasHF2_AND_BptxAND',
        'L1_MinimumBiasHF2_OR',
        'L1_MinimumBiasHF2_OR_BptxAND',
        'L1_NotBptxOR',
        'L1_NotETT100_MinimumBiasHF1_AND_BptxAND',
        'L1_NotETT110_MinimumBiasHF1_OR_BptxAND',
        'L1_NotETT110_MinimumBiasHF2_OR_BptxAND',
        'L1_NotETT150_MinimumBiasHF1_AND_BptxAND',
        'L1_NotETT150_MinimumBiasHF1_OR_BptxAND',
        'L1_NotETT150_MinimumBiasHF2_OR_BptxAND',
        'L1_NotETT200_MinimumBiasHF1_AND_BptxAND',
        'L1_NotETT20_MinimumBiasHF1_AND_BptxAND',
        'L1_NotETT20_MinimumBiasHF1_OR_BptxAND',
        'L1_NotETT20_MinimumBiasHF2_OR_BptxAND',
        'L1_NotETT80_MinimumBiasHF1_AND_BptxAND',
        'L1_NotETT80_MinimumBiasHF1_OR_BptxAND',
        'L1_NotETT80_MinimumBiasHF2_OR_BptxAND',
        'L1_NotETT95_MinimumBiasHF1_AND_BptxAND',
        'L1_NotETT95_MinimumBiasHF1_OR_BptxAND',
        'L1_NotETT95_MinimumBiasHF2_OR_BptxAND',
        'L1_NotMinimumBiasHF0_AND_BptxAND',
        'L1_NotMinimumBiasHF0_AND_BptxAND_TOTEM_1',
        'L1_NotMinimumBiasHF0_AND_BptxAND_TOTEM_2',
        'L1_NotMinimumBiasHF0_AND_BptxAND_TOTEM_4',
        'L1_NotMinimumBiasHF0_OR_BptxAND',
        'L1_NotMinimumBiasHF0_OR_BptxAND_TOTEM_1',
        'L1_NotMinimumBiasHF0_OR_BptxAND_TOTEM_2',
        'L1_NotMinimumBiasHF0_OR_BptxAND_TOTEM_4',
        'L1_NotMinimumBiasHF1_AND',
        'L1_NotMinimumBiasHF1_OR',
        'L1_NotMinimumBiasHF1_OR_BptxAND',
        'L1_NotMinimumBiasHF2_AND',
        'L1_NotMinimumBiasHF2_AND_BptxAND',
        'L1_NotMinimumBiasHF2_OR',
        'L1_NotMinimumBiasHF2_OR_BptxAND',
        'L1_SecondBunchInTrain',
        'L1_SecondLastBunchInTrain',
        'L1_SingleEG10er2p5',
        'L1_SingleEG12_BptxAND',
        'L1_SingleEG12_SingleJet28_MidEta2p7_BptxAND',
        'L1_SingleEG12_SingleJet28_MidEta2p7_MinDr0p4_BptxAND',
        'L1_SingleEG12_SingleJet32_MidEta2p7_BptxAND',
        'L1_SingleEG12_SingleJet40_MidEta2p7_BptxAND',
        'L1_SingleEG12_SingleJet44_MidEta2p7_BptxAND',
        'L1_SingleEG12_SingleJet44_MidEta2p7_MinDr0p4_BptxAND',
        'L1_SingleEG12_SingleJet56_MidEta2p7_BptxAND',
        'L1_SingleEG12_SingleJet56_MidEta2p7_MinDr0p4_BptxAND',
        'L1_SingleEG12_SingleJet60_MidEta2p7_BptxAND',
        'L1_SingleEG12_SingleJet60_MidEta2p7_MinDr0p4_BptxAND',
        'L1_SingleEG15_BptxAND',
        'L1_SingleEG15_Centrality_30_100_BptxAND',
        'L1_SingleEG15_Centrality_50_100_BptxAND',
        'L1_SingleEG15_SingleJet28_MidEta2p7_BptxAND',
        'L1_SingleEG15_SingleJet28_MidEta2p7_MinDr0p4_BptxAND',
        'L1_SingleEG15_SingleJet44_MidEta2p7_BptxAND',
        'L1_SingleEG15_SingleJet44_MidEta2p7_MinDr0p4_BptxAND',
        'L1_SingleEG15_SingleJet56_MidEta2p7_BptxAND',
        'L1_SingleEG15_SingleJet56_MidEta2p7_MinDr0p4_BptxAND',
        'L1_SingleEG15_SingleJet60_MidEta2p7_BptxAND',
        'L1_SingleEG15_SingleJet60_MidEta2p7_MinDr0p4_BptxAND',
        'L1_SingleEG15er2p5',
        'L1_SingleEG21_BptxAND',
        'L1_SingleEG21_Centrality_30_100_BptxAND',
        'L1_SingleEG21_Centrality_50_100_BptxAND',
        'L1_SingleEG26er2p5',
        'L1_SingleEG3',
        'L1_SingleEG30_BptxAND',
        'L1_SingleEG3_BptxAND',
        'L1_SingleEG3_Centrality_30_100_BptxAND',
        'L1_SingleEG3_Centrality_50_100_BptxAND',
        'L1_SingleEG3_NotMinimumBiasHF2_AND_BptxAND',
        'L1_SingleEG3_NotMinimumBiasHF2_OR_BptxAND',
        'L1_SingleEG5',
        'L1_SingleEG50',
        'L1_SingleEG5_BptxAND',
        'L1_SingleEG5_NotMinimumBiasHF2_AND_BptxAND',
        'L1_SingleEG5_NotMinimumBiasHF2_OR_BptxAND',
        'L1_SingleEG5_SingleJet28_MidEta2p7_BptxAND',
        'L1_SingleEG5_SingleJet32_MidEta2p7_BptxAND',
        'L1_SingleEG5_SingleJet40_MidEta2p7_BptxAND',
        'L1_SingleEG7_BptxAND',
        'L1_SingleEG7_Centrality_30_100_BptxAND',
        'L1_SingleEG7_Centrality_50_100_BptxAND',
        'L1_SingleEG7_SingleJet28_MidEta2p7_BptxAND',
        'L1_SingleEG7_SingleJet28_MidEta2p7_MinDr0p4_BptxAND',
        'L1_SingleEG7_SingleJet32_MidEta2p7_BptxAND',
        'L1_SingleEG7_SingleJet40_MidEta2p7_BptxAND',
        'L1_SingleEG7_SingleJet44_MidEta2p7_BptxAND',
        'L1_SingleEG7_SingleJet44_MidEta2p7_MinDr0p4_BptxAND',
        'L1_SingleEG7_SingleJet56_MidEta2p7_BptxAND',
        'L1_SingleEG7_SingleJet56_MidEta2p7_MinDr0p4_BptxAND',
        'L1_SingleEG7_SingleJet60_MidEta2p7_BptxAND',
        'L1_SingleEG7_SingleJet60_MidEta2p7_MinDr0p4_BptxAND',
        'L1_SingleEG8er2p5',
        'L1_SingleIsoEG12_BptxAND',
        'L1_SingleIsoEG15_BptxAND',
        'L1_SingleIsoEG21_BptxAND',
        'L1_SingleIsoEG3_BptxAND',
        'L1_SingleIsoEG7_BptxAND',
        'L1_SingleJet120',
        'L1_SingleJet120_FWD3p0',
        'L1_SingleJet120er2p5',
        'L1_SingleJet16_BptxAND',
        'L1_SingleJet16_Centrality_30_100_BptxAND',
        'L1_SingleJet16_Centrality_50_100_BptxAND',
        'L1_SingleJet16_FWD_BptxAND',
        'L1_SingleJet16_FWD_Centrality_30_100_BptxAND',
        'L1_SingleJet16_FWD_Centrality_50_100_BptxAND',
        'L1_SingleJet180er2p5',
        'L1_SingleJet200',
        'L1_SingleJet24_BptxAND',
        'L1_SingleJet24_Centrality_30_100_BptxAND',
        'L1_SingleJet24_Centrality_50_100_BptxAND',
        'L1_SingleJet28_BptxAND',
        'L1_SingleJet28_Centrality_30_100_BptxAND',
        'L1_SingleJet28_Centrality_50_100_BptxAND',
        'L1_SingleJet28_FWD_BptxAND',
        'L1_SingleJet28_FWD_Centrality_30_100_BptxAND',
        'L1_SingleJet28_FWD_Centrality_50_100_BptxAND',
        'L1_SingleJet32_BptxAND',
        'L1_SingleJet32_Centrality_30_100_BptxAND',
        'L1_SingleJet32_Centrality_50_100_BptxAND',
        'L1_SingleJet35',
        'L1_SingleJet35_FWD3p0',
        'L1_SingleJet35er2p5',
        'L1_SingleJet36_BptxAND',
        'L1_SingleJet36_Centrality_30_100_BptxAND',
        'L1_SingleJet36_Centrality_50_100_BptxAND',
        'L1_SingleJet36_FWD_BptxAND',
        'L1_SingleJet36_FWD_Centrality_30_100_BptxAND',
        'L1_SingleJet36_FWD_Centrality_50_100_BptxAND',
        'L1_SingleJet40_BptxAND',
        'L1_SingleJet40_Centrality_30_100_BptxAND',
        'L1_SingleJet40_Centrality_50_100_BptxAND',
        'L1_SingleJet44_BptxAND',
        'L1_SingleJet44_Centrality_30_100_BptxAND',
        'L1_SingleJet44_Centrality_50_100_BptxAND',
        'L1_SingleJet44_FWD_BptxAND',
        'L1_SingleJet44_FWD_Centrality_30_100_BptxAND',
        'L1_SingleJet44_FWD_Centrality_50_100_BptxAND',
        'L1_SingleJet48_BptxAND',
        'L1_SingleJet48_Centrality_30_100_BptxAND',
        'L1_SingleJet48_Centrality_50_100_BptxAND',
        'L1_SingleJet56_BptxAND',
        'L1_SingleJet56_Centrality_30_100_BptxAND',
        'L1_SingleJet56_Centrality_50_100_BptxAND',
        'L1_SingleJet56_FWD_BptxAND',
        'L1_SingleJet56_FWD_Centrality_30_100_BptxAND',
        'L1_SingleJet56_FWD_Centrality_50_100_BptxAND',
        'L1_SingleJet60',
        'L1_SingleJet60_BptxAND',
        'L1_SingleJet60_Centrality_30_100_BptxAND',
        'L1_SingleJet60_Centrality_50_100_BptxAND',
        'L1_SingleJet60_FWD3p0',
        'L1_SingleJet60er2p5',
        'L1_SingleJet64_BptxAND',
        'L1_SingleJet64_Centrality_30_100_BptxAND',
        'L1_SingleJet64_Centrality_50_100_BptxAND',
        'L1_SingleJet64_FWD_BptxAND',
        'L1_SingleJet64_FWD_Centrality_30_100_BptxAND',
        'L1_SingleJet64_FWD_Centrality_50_100_BptxAND',
        'L1_SingleJet72_BptxAND',
        'L1_SingleJet8',
        'L1_SingleJet80_BptxAND',
        'L1_SingleJet8_BptxAND',
        'L1_SingleJet8_Centrality_30_100_BptxAND',
        'L1_SingleJet8_Centrality_50_100_BptxAND',
        'L1_SingleJet8_FWD_BptxAND',
        'L1_SingleJet8_FWD_Centrality_30_100_BptxAND',
        'L1_SingleJet8_FWD_Centrality_50_100_BptxAND',
        'L1_SingleJet90',
        'L1_SingleJet90_FWD3p0',
        'L1_SingleJet90er2p5',
        'L1_SingleMu0',
        'L1_SingleMu0_BMTF',
        'L1_SingleMu0_BptxAND',
        'L1_SingleMu0_DQ',
        'L1_SingleMu0_EMTF',
        'L1_SingleMu0_NotMinimumBiasHF2_AND_BptxAND',
        'L1_SingleMu0_NotMinimumBiasHF2_OR_BptxAND',
        'L1_SingleMu0_OMTF',
        'L1_SingleMu12',
        'L1_SingleMu12_BptxAND',
        'L1_SingleMu12_DQ_BMTF',
        'L1_SingleMu12_DQ_EMTF',
        'L1_SingleMu12_DQ_OMTF',
        'L1_SingleMu12_MinimumBiasHF1_AND_BptxAND',
        'L1_SingleMu12_SingleEG7_BptxAND',
        'L1_SingleMu16',
        'L1_SingleMu16_BptxAND',
        'L1_SingleMu16_MinimumBiasHF1_AND_BptxAND',
        'L1_SingleMu22',
        'L1_SingleMu22_BMTF',
        'L1_SingleMu22_EMTF',
        'L1_SingleMu22_OMTF',
        'L1_SingleMu3',
        'L1_SingleMu3Open_BptxAND',
        'L1_SingleMu3_BptxAND',
        'L1_SingleMu3_Centrality_70_100_BptxAND',
        'L1_SingleMu3_Centrality_80_100_BptxAND',
        'L1_SingleMu3_MinimumBiasHF1_AND_BptxAND',
        'L1_SingleMu3_NotMinimumBiasHF2_OR_BptxAND',
        'L1_SingleMu3_SingleEG12_BptxAND',
        'L1_SingleMu3_SingleEG15_BptxAND',
        'L1_SingleMu3_SingleEG20_BptxAND',
        'L1_SingleMu3_SingleEG30_BptxAND',
        'L1_SingleMu3_SingleJet28_MidEta2p7_BptxAND',
        'L1_SingleMu3_SingleJet32_MidEta2p7_BptxAND',
        'L1_SingleMu3_SingleJet40_MidEta2p7_BptxAND',
        'L1_SingleMu5',
        'L1_SingleMu5_BptxAND',
        'L1_SingleMu5_MinimumBiasHF1_AND_BptxAND',
        'L1_SingleMu5_SingleEG10_BptxAND',
        'L1_SingleMu5_SingleEG12_BptxAND',
        'L1_SingleMu5_SingleEG15_BptxAND',
        'L1_SingleMu5_SingleEG20_BptxAND',
        'L1_SingleMu7',
        'L1_SingleMu7_BptxAND',
        'L1_SingleMu7_DQ',
        'L1_SingleMu7_MinimumBiasHF1_AND_BptxAND',
        'L1_SingleMu7_SingleEG10_BptxAND',
        'L1_SingleMu7_SingleEG12_BptxAND',
        'L1_SingleMu7_SingleEG15_BptxAND',
        'L1_SingleMu7_SingleEG7_BptxAND',
        'L1_SingleMuCosmics',
        'L1_SingleMuCosmics_BMTF',
        'L1_SingleMuCosmics_EMTF',
        'L1_SingleMuCosmics_OMTF',
        'L1_SingleMuOpen',
        'L1_SingleMuOpen_BptxAND',
        'L1_SingleMuOpen_Centrality_70_100_BptxAND',
        'L1_SingleMuOpen_Centrality_70_100_MinimumBiasHF1_AND_BptxAND',
        'L1_SingleMuOpen_Centrality_80_100_BptxAND',
        'L1_SingleMuOpen_Centrality_80_100_MinimumBiasHF1_AND_BptxAND',
        'L1_SingleMuOpen_NotMinimumBiasHF2_AND_BptxAND',
        'L1_SingleMuOpen_NotMinimumBiasHF2_OR_BptxAND',
        'L1_SingleMuOpen_SingleEG15_BptxAND',
        'L1_SingleMuOpen_SingleJet28_MidEta2p7_BptxAND',
        'L1_SingleMuOpen_SingleJet44_MidEta2p7_BptxAND',
        'L1_SingleMuOpen_SingleJet56_MidEta2p7_BptxAND',
        'L1_SingleMuOpen_SingleJet64_MidEta2p7_BptxAND',
        'L1_TOTEM_1',
        'L1_TOTEM_2',
        'L1_TOTEM_3',
        'L1_TOTEM_4',
        'L1_UnpairedBunchBptxMinus',
        'L1_UnpairedBunchBptxPlus',
        'L1_ZDCM',
        'L1_ZDCM_BptxAND',
        'L1_ZDCM_ZDCP_BptxAND',
        'L1_ZDCP',
        'L1_ZDCP_BptxAND',
        'L1_ZDC_AND_OR_MinimumBiasHF1_AND_BptxAND',
        'L1_ZDC_AND_OR_MinimumBiasHF1_OR_BptxAND',
        'L1_ZDC_AND_OR_MinimumBiasHF2_AND_BptxAND',
        'L1_ZDC_AND_OR_MinimumBiasHF2_OR_BptxAND',
        'L1_ZDC_OR_OR_MinimumBiasHF1_OR_BptxAND',
        'L1_ZDC_OR_OR_MinimumBiasHF2_OR_BptxAND',
        'L1_ZeroBias',
        'L1_ZeroBias_copy'
     ) ),
    l1results = cms.InputTag("gtStage2Digis")
)


process.hltobject = cms.EDAnalyzer("TriggerObjectAnalyzer",
    processName = cms.string('HLT'),
    treeName = cms.string('JetTriggers'),
    triggerObjects = cms.InputTag("slimmedPatTrigger","","PAT"),
    triggerResults = cms.InputTag("TriggerResults","","HLT")
)


process.inclusiveJetAnalyzer = cms.EDAnalyzer("HiInclusiveJetAnalyzer",
    caloJetTag = cms.InputTag("slimmedCaloJets"),
    doCandidateBtagging = cms.untracked.bool(True),
    doGenSym = cms.untracked.bool(False),
    doGenTaus = cms.untracked.bool(False),
    doHiJetID = cms.untracked.bool(False),
    doJetConstituents = cms.untracked.bool(False),
    doLegacyBtagging = cms.untracked.bool(False),
    doStandardJetID = cms.untracked.bool(False),
    doSubEvent = cms.untracked.bool(False),
    doSubJets = cms.untracked.bool(False),
    doWTARecluster = cms.untracked.bool(False),
    eventInfoTag = cms.InputTag("generator"),
    fillGenJets = cms.untracked.bool(False),
    genjetTag = cms.InputTag("ak4HiGenJets"),
    isMC = cms.untracked.bool(False),
    jetPtMin = cms.double(5.0),
    jetTag = cms.InputTag("ak4PFJets"),
    matchJets = cms.untracked.bool(False),
    matchTag = cms.untracked.InputTag("akPu4PFpatJets"),
    pfCandidateLabel = cms.untracked.InputTag("packedPFCandidates"),
    rParam = cms.double(0.4),
    towersSrc = cms.InputTag("towerMaker"),
    trackQuality = cms.untracked.string('highPurity'),
    trackTag = cms.InputTag("hiTracks"),
    useHepMC = cms.untracked.bool(False),
    useQuality = cms.untracked.bool(True),
    useRawPt = cms.untracked.bool(False),
    vtxTag = cms.InputTag("hiSelectedVertex")
)


process.l1object = cms.EDAnalyzer("L1UpgradeFlatTreeProducer",
    doEg = cms.bool(True),
    doJet = cms.bool(True),
    doMuon = cms.bool(True),
    doSum = cms.bool(True),
    doTau = cms.bool(True),
    egToken = cms.untracked.InputTag("caloStage2Digis","EGamma"),
    jetToken = cms.untracked.InputTag("caloStage2Digis","Jet"),
    maxL1Upgrade = cms.uint32(60),
    muonToken = cms.untracked.InputTag("gmtStage2Digis","Muon"),
    sumToken = cms.untracked.InputTag("caloStage2Digis","EtSum"),
    tauTokens = cms.untracked.VInputTag(cms.InputTag("caloStage2Digis","Tau"))
)


process.muonAnalyzer = cms.EDAnalyzer("MuonAnalyzer",
    doGen = cms.bool(False),
    doReco = cms.untracked.bool(True),
    genparticle = cms.InputTag("packedGenParticles"),
    muonSrc = cms.InputTag("unpackedMuons"),
    simtrack = cms.InputTag("mergedtruth","MergedTrackTruth"),
    vertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.particleFlowAnalyser = cms.EDAnalyzer("ParticleFlowAnalyser",
    absEtaMax = cms.double(5.0),
    pfCandidateSrc = cms.InputTag("packedPFCandidates"),
    ptMin = cms.double(5.0)
)


process.ppTracks = cms.EDAnalyzer("TrackAnalyzer",
    beamSpotSrc = cms.untracked.InputTag("offlineBeamSpot"),
    chi2Map = cms.InputTag("packedPFCandidateTrackChi2"),
    chi2MapLost = cms.InputTag("lostTrackChi2"),
    doTrack = cms.untracked.bool(True),
    lostTracksSrc = cms.InputTag("lostTracks"),
    packedCandSrc = cms.InputTag("packedPFCandidates"),
    trackPtMin = cms.untracked.double(0.01),
    vertexSrc = cms.InputTag("unpackedTracksAndVertices")
)


process.randomConeAnalysisR2 = cms.EDAnalyzer("RandomConeAnalysis",
    etaMap = cms.InputTag("hiPuRho","mapEtaEdges"),
    maxFlowChi2Prob = cms.double(0.95),
    minFlowChi2Prob = cms.double(0.05),
    pfCandSource = cms.InputTag("packedPFCandidates"),
    rho = cms.InputTag("hiPuRho","mapToRho"),
    rhoFlowFitParams = cms.InputTag("hiFJRhoFlowModulationProducer","rhoFlowFitParams"),
    rngConeRadius = cms.double(0.2),
    useModulatedRho = cms.bool(True)
)


process.randomConeAnalysisR4 = cms.EDAnalyzer("RandomConeAnalysis",
    etaMap = cms.InputTag("hiPuRho","mapEtaEdges"),
    maxFlowChi2Prob = cms.double(0.95),
    minFlowChi2Prob = cms.double(0.05),
    pfCandSource = cms.InputTag("packedPFCandidates"),
    rho = cms.InputTag("hiPuRho","mapToRho"),
    rhoFlowFitParams = cms.InputTag("hiFJRhoFlowModulationProducer","rhoFlowFitParams"),
    rngConeRadius = cms.double(0.4),
    useModulatedRho = cms.bool(True)
)


process.rhoAnalysis = cms.EDAnalyzer("RhoAnalysis",
    etaMap = cms.InputTag("hiPuRho","mapEtaEdges"),
    rho = cms.InputTag("hiPuRho","mapToRho")
)


process.skimanalysis = cms.EDAnalyzer("FilterAnalyzer",
    hltresults = cms.InputTag("TriggerResults","","HiForest"),
    superFilters = cms.vstring('')
)


process.trackAnalyzer = cms.EDAnalyzer("TrackAnalyzer",
    beamSpotSrc = cms.untracked.InputTag("offlineBeamSpot"),
    chi2Map = cms.InputTag("packedPFCandidateTrackChi2"),
    chi2MapLost = cms.InputTag("lostTrackChi2"),
    doTrack = cms.untracked.bool(True),
    lostTracksSrc = cms.InputTag("lostTracks"),
    packedCandSrc = cms.InputTag("packedPFCandidates"),
    trackPtMin = cms.untracked.double(0.01),
    vertexSrc = cms.InputTag("unpackedTracksAndVertices")
)


process.zdcanalyzer = cms.EDAnalyzer("ZDCTreeProducer",
    calZDCDigi = cms.bool(False),
    doZDCDigi = cms.bool(True),
    doZDCRecHit = cms.bool(False),
    nZdcTs = cms.int32(10),
    verbose = cms.bool(False),
    zdcDigiSrc = cms.InputTag("hcalDigis","ZDC"),
    zdcRecHitSrc = cms.InputTag("QWzdcreco")
)


process.DQMStore = cms.Service("DQMStore",
    MEsToSave = cms.untracked.vstring(
        'Hcal/DigiTask/Occupancy/depth/depth1',
        'Hcal/DigiTask/Occupancy/depth/depth2',
        'Hcal/DigiTask/Occupancy/depth/depth3',
        'Hcal/DigiTask/Occupancy/depth/depth4',
        'Hcal/DigiTask/Occupancy/depth/depth5',
        'Hcal/DigiTask/Occupancy/depth/depth6',
        'Hcal/DigiTask/Occupancy/depth/depth7',
        'Hcal/DigiTask/Occupancy/depth/depthHO',
        'Hcal/DigiTask/OccupancyCut/depth/depth1',
        'Hcal/DigiTask/OccupancyCut/depth/depth2',
        'Hcal/DigiTask/OccupancyCut/depth/depth3',
        'Hcal/DigiTask/OccupancyCut/depth/depth4',
        'Hcal/DigiTask/OccupancyCut/depth/depth5',
        'Hcal/DigiTask/OccupancyCut/depth/depth6',
        'Hcal/DigiTask/OccupancyCut/depth/depth7',
        'Hcal/DigiTask/OccupancyCut/depth/depthHO',
        'EcalBarrel/EBOccupancyTask/EBOT digi occupancy',
        'EcalEndcap/EEOccupancyTask/EEOT digi occupancy EE -',
        'EcalEndcap/EEOccupancyTask/EEOT digi occupancy EE +',
        'EcalBarrel/EBOccupancyTask/EBOT DCC entries',
        'EcalEndcap/EEOccupancyTask/EEOT DCC entries',
        'Ecal/EventInfo/processedEvents',
        'PixelPhase1/Tracks/charge_PXBarrel',
        'PixelPhase1/Tracks/charge_PXForward',
        'PixelPhase1/Tracks/PXBarrel/charge_PXLayer_1',
        'PixelPhase1/Tracks/PXBarrel/charge_PXLayer_2',
        'PixelPhase1/Tracks/PXBarrel/charge_PXLayer_3',
        'PixelPhase1/Tracks/PXBarrel/charge_PXLayer_4',
        'PixelPhase1/Tracks/PXForward/charge_PXDisk_+1',
        'PixelPhase1/Tracks/PXForward/charge_PXDisk_+2',
        'PixelPhase1/Tracks/PXForward/charge_PXDisk_+3',
        'PixelPhase1/Tracks/PXForward/charge_PXDisk_-1',
        'PixelPhase1/Tracks/PXForward/charge_PXDisk_-2',
        'PixelPhase1/Tracks/PXForward/charge_PXDisk_-3',
        'PixelPhase1/Tracks/PXBarrel/size_PXLayer_1',
        'PixelPhase1/Tracks/PXBarrel/size_PXLayer_2',
        'PixelPhase1/Tracks/PXBarrel/size_PXLayer_3',
        'PixelPhase1/Tracks/PXBarrel/size_PXLayer_4',
        'PixelPhase1/Tracks/PXForward/size_PXDisk_+1',
        'PixelPhase1/Tracks/PXForward/size_PXDisk_+2',
        'PixelPhase1/Tracks/PXForward/size_PXDisk_+3',
        'PixelPhase1/Tracks/PXForward/size_PXDisk_-1',
        'PixelPhase1/Tracks/PXForward/size_PXDisk_-2',
        'PixelPhase1/Tracks/PXForward/size_PXDisk_-3',
        'HLT/Vertexing/hltPixelVertices/hltPixelVertices/goodvtxNbr',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/mon_eta',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/mon_hits',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/mon_phi',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/mon_pt',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/mon_unMatched_eta',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/mon_unMatched_hits',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/mon_unMatched_phi',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/mon_unMatched_pt',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/ref_eta',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/ref_hits',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/ref_matched_eta',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/ref_matched_hits',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/ref_matched_phi',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/ref_matched_pt',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/ref_phi',
        'HLT/Tracking/ValidationWRTOffline/hltMergedWrtHighPurityPV/ref_pt',
        'HLT/Tracking/pixelTracks/GeneralProperties/Chi2Prob_GenTk',
        'HLT/Tracking/pixelTracks/GeneralProperties/Chi2oNDFVsEta_ImpactPoint_GenTk',
        'HLT/Tracking/pixelTracks/GeneralProperties/DeltaZToPVZoom_GenTk',
        'HLT/Tracking/pixelTracks/GeneralProperties/DistanceOfClosestApproachToPVVsPhi_GenTk',
        'HLT/Tracking/pixelTracks/GeneralProperties/DistanceOfClosestApproachToPVZoom_GenTk',
        'HLT/Tracking/pixelTracks/GeneralProperties/NumberOfTracks_GenTk',
        'HLT/Tracking/tracks/GeneralProperties/Chi2Prob_GenTk',
        'HLT/Tracking/tracks/GeneralProperties/Chi2oNDFVsEta_ImpactPoint_GenTk',
        'HLT/Tracking/tracks/GeneralProperties/DeltaZToPVZoom_GenTk',
        'HLT/Tracking/tracks/GeneralProperties/DistanceOfClosestApproachToPVVsPhi_GenTk',
        'HLT/Tracking/tracks/GeneralProperties/DistanceOfClosestApproachToPVZoom_GenTk',
        'HLT/Tracking/tracks/GeneralProperties/NumberOfTracks_GenTk',
        'HLT/Tracking/tracks/LUMIanalysis/NumberEventsVsLUMI',
        'HLT/Tracking/tracks/PUmonitoring/NumberEventsVsGoodPVtx',
        'PixelPhase1/Tracks/num_clusters_ontrack_PXBarrel',
        'PixelPhase1/Tracks/num_clusters_ontrack_PXForward',
        'PixelPhase1/Tracks/clusterposition_zphi_ontrack',
        'PixelPhase1/Tracks/PXBarrel/clusterposition_zphi_ontrack_PXLayer_1',
        'PixelPhase1/Tracks/PXBarrel/clusterposition_zphi_ontrack_PXLayer_2',
        'PixelPhase1/Tracks/PXBarrel/clusterposition_zphi_ontrack_PXLayer_3',
        'PixelPhase1/Tracks/PXBarrel/clusterposition_zphi_ontrack_PXLayer_4',
        'PixelPhase1/Tracks/PXForward/clusterposition_xy_ontrack_PXDisk_+1',
        'PixelPhase1/Tracks/PXForward/clusterposition_xy_ontrack_PXDisk_+2',
        'PixelPhase1/Tracks/PXForward/clusterposition_xy_ontrack_PXDisk_+3',
        'PixelPhase1/Tracks/PXForward/clusterposition_xy_ontrack_PXDisk_-1',
        'PixelPhase1/Tracks/PXForward/clusterposition_xy_ontrack_PXDisk_-2',
        'PixelPhase1/Tracks/PXForward/clusterposition_xy_ontrack_PXDisk_-3',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_1/NormalizedHitResiduals_TEC__wheel__1',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_2/NormalizedHitResiduals_TEC__wheel__2',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_3/NormalizedHitResiduals_TEC__wheel__3',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_4/NormalizedHitResiduals_TEC__wheel__4',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_5/NormalizedHitResiduals_TEC__wheel__5',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_6/NormalizedHitResiduals_TEC__wheel__6',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_7/NormalizedHitResiduals_TEC__wheel__7',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_8/NormalizedHitResiduals_TEC__wheel__8',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_9/NormalizedHitResiduals_TEC__wheel__9',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_1/NormalizedHitResiduals_TEC__wheel__1',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_2/NormalizedHitResiduals_TEC__wheel__2',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_3/NormalizedHitResiduals_TEC__wheel__3',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_4/NormalizedHitResiduals_TEC__wheel__4',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_5/NormalizedHitResiduals_TEC__wheel__5',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_6/NormalizedHitResiduals_TEC__wheel__6',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_7/NormalizedHitResiduals_TEC__wheel__7',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_8/NormalizedHitResiduals_TEC__wheel__8',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_9/NormalizedHitResiduals_TEC__wheel__9',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_1/Summary_ClusterStoNCorr__OnTrack__TEC__PLUS__wheel__1',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_2/Summary_ClusterStoNCorr__OnTrack__TEC__PLUS__wheel__2',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_3/Summary_ClusterStoNCorr__OnTrack__TEC__PLUS__wheel__3',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_4/Summary_ClusterStoNCorr__OnTrack__TEC__PLUS__wheel__4',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_5/Summary_ClusterStoNCorr__OnTrack__TEC__PLUS__wheel__5',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_6/Summary_ClusterStoNCorr__OnTrack__TEC__PLUS__wheel__6',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_7/Summary_ClusterStoNCorr__OnTrack__TEC__PLUS__wheel__7',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_8/Summary_ClusterStoNCorr__OnTrack__TEC__PLUS__wheel__8',
        'SiStrip/MechanicalView/TEC/PLUS/wheel_9/Summary_ClusterStoNCorr__OnTrack__TEC__PLUS__wheel__9',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_1/Summary_ClusterStoNCorr__OnTrack__TEC__MINUS__wheel__1',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_2/Summary_ClusterStoNCorr__OnTrack__TEC__MINUS__wheel__2',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_3/Summary_ClusterStoNCorr__OnTrack__TEC__MINUS__wheel__3',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_4/Summary_ClusterStoNCorr__OnTrack__TEC__MINUS__wheel__4',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_5/Summary_ClusterStoNCorr__OnTrack__TEC__MINUS__wheel__5',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_6/Summary_ClusterStoNCorr__OnTrack__TEC__MINUS__wheel__6',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_7/Summary_ClusterStoNCorr__OnTrack__TEC__MINUS__wheel__7',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_8/Summary_ClusterStoNCorr__OnTrack__TEC__MINUS__wheel__8',
        'SiStrip/MechanicalView/TEC/MINUS/wheel_9/Summary_ClusterStoNCorr__OnTrack__TEC__MINUS__wheel__9',
        'SiStrip/MechanicalView/TIB/layer_1/NormalizedHitResiduals_TIB__Layer__1',
        'SiStrip/MechanicalView/TIB/layer_2/NormalizedHitResiduals_TIB__Layer__2',
        'SiStrip/MechanicalView/TIB/layer_3/NormalizedHitResiduals_TIB__Layer__3',
        'SiStrip/MechanicalView/TIB/layer_4/NormalizedHitResiduals_TIB__Layer__4',
        'SiStrip/MechanicalView/TIB/layer_1/Summary_ClusterStoNCorr__OnTrack__TIB__layer__1',
        'SiStrip/MechanicalView/TIB/layer_2/Summary_ClusterStoNCorr__OnTrack__TIB__layer__2',
        'SiStrip/MechanicalView/TIB/layer_3/Summary_ClusterStoNCorr__OnTrack__TIB__layer__3',
        'SiStrip/MechanicalView/TIB/layer_4/Summary_ClusterStoNCorr__OnTrack__TIB__layer__4',
        'SiStrip/MechanicalView/TID/PLUS/wheel_1/NormalizedHitResiduals_TID__wheel__1',
        'SiStrip/MechanicalView/TID/PLUS/wheel_2/NormalizedHitResiduals_TID__wheel__2',
        'SiStrip/MechanicalView/TID/PLUS/wheel_3/NormalizedHitResiduals_TID__wheel__3',
        'SiStrip/MechanicalView/TID/MINUS/wheel_1/NormalizedHitResiduals_TID__wheel__1',
        'SiStrip/MechanicalView/TID/MINUS/wheel_2/NormalizedHitResiduals_TID__wheel__2',
        'SiStrip/MechanicalView/TID/MINUS/wheel_3/NormalizedHitResiduals_TID__wheel__3',
        'SiStrip/MechanicalView/TID/PLUS/wheel_1/Summary_ClusterStoNCorr__OnTrack__TID__PLUS__wheel__1',
        'SiStrip/MechanicalView/TID/PLUS/wheel_2/Summary_ClusterStoNCorr__OnTrack__TID__PLUS__wheel__2',
        'SiStrip/MechanicalView/TID/PLUS/wheel_3/Summary_ClusterStoNCorr__OnTrack__TID__PLUS__wheel__3',
        'SiStrip/MechanicalView/TID/MINUS/wheel_1/Summary_ClusterStoNCorr__OnTrack__TID__MINUS__wheel__1',
        'SiStrip/MechanicalView/TID/MINUS/wheel_2/Summary_ClusterStoNCorr__OnTrack__TID__MINUS__wheel__2',
        'SiStrip/MechanicalView/TID/MINUS/wheel_3/Summary_ClusterStoNCorr__OnTrack__TID__MINUS__wheel__3',
        'SiStrip/MechanicalView/TOB/layer_1/NormalizedHitResiduals_TOB__Layer__1',
        'SiStrip/MechanicalView/TOB/layer_2/NormalizedHitResiduals_TOB__Layer__2',
        'SiStrip/MechanicalView/TOB/layer_3/NormalizedHitResiduals_TOB__Layer__3',
        'SiStrip/MechanicalView/TOB/layer_4/NormalizedHitResiduals_TOB__Layer__4',
        'SiStrip/MechanicalView/TOB/layer_5/NormalizedHitResiduals_TOB__Layer__5',
        'SiStrip/MechanicalView/TOB/layer_6/NormalizedHitResiduals_TOB__Layer__6',
        'SiStrip/MechanicalView/TOB/layer_1/Summary_ClusterStoNCorr__OnTrack__TOB__layer__1',
        'SiStrip/MechanicalView/TOB/layer_2/Summary_ClusterStoNCorr__OnTrack__TOB__layer__2',
        'SiStrip/MechanicalView/TOB/layer_3/Summary_ClusterStoNCorr__OnTrack__TOB__layer__3',
        'SiStrip/MechanicalView/TOB/layer_4/Summary_ClusterStoNCorr__OnTrack__TOB__layer__4',
        'SiStrip/MechanicalView/TOB/layer_5/Summary_ClusterStoNCorr__OnTrack__TOB__layer__5',
        'SiStrip/MechanicalView/TOB/layer_6/Summary_ClusterStoNCorr__OnTrack__TOB__layer__6',
        'SiStrip/MechanicalView/MainDiagonal Position',
        'SiStrip/MechanicalView/NumberOfClustersInPixel',
        'SiStrip/MechanicalView/NumberOfClustersInStrip',
        'Tracking/TrackParameters/generalTracks/LSanalysis/Chi2oNDF_lumiFlag_GenTk',
        'Tracking/TrackParameters/generalTracks/LSanalysis/NumberOfRecHitsPerTrack_lumiFlag_GenTk',
        'Tracking/TrackParameters/generalTracks/LSanalysis/NumberOfTracks_lumiFlag_GenTk',
        'Tracking/TrackParameters/highPurityTracks/pt_1/GeneralProperties/SIPDxyToPV_GenTk',
        'Tracking/TrackParameters/highPurityTracks/pt_1/GeneralProperties/SIPDzToPV_GenTk',
        'Tracking/TrackParameters/highPurityTracks/pt_1/GeneralProperties/SIP3DToPV_GenTk',
        'Tracking/TrackParameters/generalTracks/HitProperties/NumberOfMissingOuterRecHitsPerTrack_GenTk',
        'Tracking/TrackParameters/generalTracks/HitProperties/NumberMORecHitsPerTrackVsPt_GenTk',
        'OfflinePV/offlinePrimaryVertices/tagVtxProb',
        'OfflinePV/offlinePrimaryVertices/tagType',
        'OfflinePV/Resolution/PV/pull_x',
        'OfflinePV/Resolution/PV/pull_y',
        'OfflinePV/Resolution/PV/pull_z',
        'JetMET/Jet/Cleanedak4PFJetsCHS/CHFrac_highPt_Barrel',
        'JetMET/Jet/Cleanedak4PFJetsCHS/CHFrac_highPt_EndCap',
        'JetMET/Jet/Cleanedak4PFJetsCHS/CHFrac_mediumPt_Barrel',
        'JetMET/Jet/Cleanedak4PFJetsCHS/CHFrac_mediumPt_EndCap',
        'JetMET/Jet/Cleanedak4PFJetsCHS/CHFrac_lowPt_Barrel',
        'JetMET/Jet/Cleanedak4PFJetsCHS/CHFrac_lowPt_EndCap',
        'JetMET/Jet/Cleanedak4PFJetsCHS/ChMultiplicity_highPt_Barrel',
        'JetMET/Jet/Cleanedak4PFJetsCHS/ChMultiplicity_highPt_EndCap',
        'JetMET/Jet/Cleanedak4PFJetsCHS/ChMultiplicity_mediumPt_Barrel',
        'JetMET/Jet/Cleanedak4PFJetsCHS/ChMultiplicity_mediumPt_EndCap',
        'JetMET/Jet/Cleanedak4PFJetsCHS/ChMultiplicity_lowPt_Barrel',
        'JetMET/Jet/Cleanedak4PFJetsCHS/ChMultiplicity_lowPt_EndCap',
        'JetMET/Jet/Cleanedak4PFJetsCHS/Constituents',
        'JetMET/Jet/Cleanedak4PFJetsCHS/Eta',
        'JetMET/Jet/Cleanedak4PFJetsCHS/Eta_uncor',
        'JetMET/Jet/Cleanedak4PFJetsCHS/JetEnergyCorr',
        'JetMET/Jet/Cleanedak4PFJetsCHS/NJets',
        'JetMET/Jet/Cleanedak4PFJetsCHS/Phi',
        'JetMET/Jet/Cleanedak4PFJetsCHS/Phi_Barrel',
        'JetMET/Jet/Cleanedak4PFJetsCHS/Phi_EndCap',
        'JetMET/Jet/Cleanedak4PFJetsCHS/Pt',
        'JetMET/MET/pfMETT1/Cleaned/METSig',
        'JetMET/vertices'
    ),
    assertLegacySafe = cms.untracked.bool(False),
    enableMultiThread = cms.untracked.bool(True),
    onlineMode = cms.untracked.bool(False),
    saveByLumi = cms.untracked.bool(False),
    trackME = cms.untracked.string(''),
    verbose = cms.untracked.int32(0)
)


process.MessageLogger = cms.Service("MessageLogger",
    cerr = cms.untracked.PSet(
        FwkReport = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            reportEvery = cms.untracked.int32(1)
        ),
        FwkSummary = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            reportEvery = cms.untracked.int32(1)
        ),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000)
        ),
        enable = cms.untracked.bool(True),
        enableStatistics = cms.untracked.bool(False),
        lineLength = cms.optional.untracked.int32,
        noLineBreaks = cms.optional.untracked.bool,
        noTimeStamps = cms.untracked.bool(False),
        resetStatistics = cms.untracked.bool(False),
        statisticsThreshold = cms.untracked.string('WARNING'),
        threshold = cms.untracked.string('INFO'),
        allowAnyLabel_=cms.optional.untracked.PSetTemplate(
            limit = cms.optional.untracked.int32,
            reportEvery = cms.untracked.int32(1),
            timespan = cms.optional.untracked.int32
        )
    ),
    cout = cms.untracked.PSet(
        enable = cms.untracked.bool(False),
        enableStatistics = cms.untracked.bool(False),
        lineLength = cms.optional.untracked.int32,
        noLineBreaks = cms.optional.untracked.bool,
        noTimeStamps = cms.optional.untracked.bool,
        resetStatistics = cms.untracked.bool(False),
        statisticsThreshold = cms.optional.untracked.string,
        threshold = cms.optional.untracked.string,
        allowAnyLabel_=cms.optional.untracked.PSetTemplate(
            limit = cms.optional.untracked.int32,
            reportEvery = cms.untracked.int32(1),
            timespan = cms.optional.untracked.int32
        )
    ),
    debugModules = cms.untracked.vstring(),
    default = cms.untracked.PSet(
        limit = cms.optional.untracked.int32,
        lineLength = cms.untracked.int32(80),
        noLineBreaks = cms.untracked.bool(False),
        noTimeStamps = cms.untracked.bool(False),
        reportEvery = cms.untracked.int32(1),
        statisticsThreshold = cms.untracked.string('INFO'),
        threshold = cms.untracked.string('INFO'),
        timespan = cms.optional.untracked.int32,
        allowAnyLabel_=cms.optional.untracked.PSetTemplate(
            limit = cms.optional.untracked.int32,
            reportEvery = cms.untracked.int32(1),
            timespan = cms.optional.untracked.int32
        )
    ),
    files = cms.untracked.PSet(
        allowAnyLabel_=cms.optional.untracked.PSetTemplate(
            enableStatistics = cms.untracked.bool(False),
            extension = cms.optional.untracked.string,
            filename = cms.optional.untracked.string,
            lineLength = cms.optional.untracked.int32,
            noLineBreaks = cms.optional.untracked.bool,
            noTimeStamps = cms.optional.untracked.bool,
            output = cms.optional.untracked.string,
            resetStatistics = cms.untracked.bool(False),
            statisticsThreshold = cms.optional.untracked.string,
            threshold = cms.optional.untracked.string,
            allowAnyLabel_=cms.optional.untracked.PSetTemplate(
                limit = cms.optional.untracked.int32,
                reportEvery = cms.untracked.int32(1),
                timespan = cms.optional.untracked.int32
            )
        )
    ),
    suppressDebug = cms.untracked.vstring(),
    suppressFwkInfo = cms.untracked.vstring(),
    suppressInfo = cms.untracked.vstring(),
    suppressWarning = cms.untracked.vstring(),
    allowAnyLabel_=cms.optional.untracked.PSetTemplate(
        limit = cms.optional.untracked.int32,
        reportEvery = cms.untracked.int32(1),
        timespan = cms.optional.untracked.int32
    )
)


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    CTPPSFastRecHits = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(1357987)
    ),
    LHCTransport = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(87654321)
    ),
    MuonSimHits = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(987346)
    ),
    RPSiDetDigitizer = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(137137)
    ),
    RPixDetDigitizer = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(137137)
    ),
    VtxSmeared = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(98765432)
    ),
    ecalPreshowerRecHit = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(6541321)
    ),
    ecalRecHit = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(654321)
    ),
    externalLHEProducer = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(234567)
    ),
    famosPileUp = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(918273)
    ),
    fastSimProducer = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(13579)
    ),
    fastTrackerRecHits = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(24680)
    ),
    g4SimHits = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(11)
    ),
    generator = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(123456789)
    ),
    hbhereco = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(541321)
    ),
    hfreco = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(541321)
    ),
    hiSignal = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(123456789)
    ),
    hiSignalG4SimHits = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(11)
    ),
    hiSignalLHCTransport = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(88776655)
    ),
    horeco = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(541321)
    ),
    l1ParamMuons = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(6453209)
    ),
    mix = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(12345)
    ),
    mixData = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(12345)
    ),
    mixGenPU = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(918273)
    ),
    mixRecoTracks = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(918273)
    ),
    mixSimCaloHits = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(918273)
    ),
    paramMuons = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(54525)
    ),
    saveFileName = cms.untracked.string(''),
    simBeamSpotFilter = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(87654321)
    ),
    simMuonCSCDigis = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(11223344)
    ),
    simMuonDTDigis = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(1234567)
    ),
    simMuonGEMDigis = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(1234567)
    ),
    simMuonRPCDigis = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(1234567)
    ),
    simSiStripDigiSimLink = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(1234567)
    )
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('Run3_data_rho_test_new_file.root')
)


process.CSCGeometryESModule = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    debugV = cms.untracked.bool(False),
    fromDD4hep = cms.bool(False),
    fromDDD = cms.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useGangedStripsInME1a = cms.bool(False),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.CaloGeometryBuilder = cms.ESProducer("CaloGeometryBuilder",
    SelectedCalos = cms.vstring(
        'HCAL',
        'ZDC',
        'EcalBarrel',
        'EcalEndcap',
        'EcalPreshower',
        'TOWER'
    )
)


process.CaloTopologyBuilder = cms.ESProducer("CaloTopologyBuilder")


process.CaloTowerGeometryFromDBEP = cms.ESProducer("CaloTowerGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.CaloTowerTopologyEP = cms.ESProducer("CaloTowerTopologyEP")


process.CastorDbProducer = cms.ESProducer("CastorDbProducer",
    appendToDataLabel = cms.string('')
)


process.CastorGeometryFromDBEP = cms.ESProducer("CastorGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.DDCompactViewESProducer = cms.ESProducer("DDCompactViewESProducer",
    appendToDataLabel = cms.string('')
)


process.DDSpecParRegistryESProducer = cms.ESProducer("DDSpecParRegistryESProducer",
    appendToDataLabel = cms.string('')
)


process.DDVectorRegistryESProducer = cms.ESProducer("DDVectorRegistryESProducer",
    appendToDataLabel = cms.string('')
)


process.DTGeometryESModule = cms.ESProducer("DTGeometryESModule",
    DDDetector = cms.ESInputTag("",""),
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    attribute = cms.string('MuStructure'),
    fromDD4hep = cms.bool(False),
    fromDDD = cms.bool(False),
    value = cms.string('MuonBarrelDT')
)


process.EcalBarrelGeometryFromDBEP = cms.ESProducer("EcalBarrelGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalElectronicsMappingBuilder = cms.ESProducer("EcalElectronicsMappingBuilder")


process.EcalEndcapGeometryFromDBEP = cms.ESProducer("EcalEndcapGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalLaserCorrectionService = cms.ESProducer("EcalLaserCorrectionService",
    maxExtrapolationTimeInSec = cms.uint32(0)
)


process.EcalLaserCorrectionServiceMC = cms.ESProducer("EcalLaserCorrectionServiceMC",
    appendToDataLabel = cms.string('')
)


process.EcalPreshowerGeometryFromDBEP = cms.ESProducer("EcalPreshowerGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalTrigTowerConstituentsMapBuilder = cms.ESProducer("EcalTrigTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/EcalMapping/data/EndCap_TTMap.txt')
)


process.GEMGeometryESModule = cms.ESProducer("GEMGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDD4hep = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.GlobalTrackingGeometryESProducer = cms.ESProducer("GlobalTrackingGeometryESProducer")


process.HcalAlignmentEP = cms.ESProducer("HcalAlignmentEP")


process.HcalGeometryFromDBEP = cms.ESProducer("HcalGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.MuonDetLayerGeometryESProducer = cms.ESProducer("MuonDetLayerGeometryESProducer")


process.ParabolicParametrizedMagneticFieldProducer = cms.ESProducer("AutoParametrizedMagneticFieldProducer",
    label = cms.untracked.string('ParabolicMf'),
    valueOverride = cms.int32(18268),
    version = cms.string('Parabolic')
)


process.RPCGeometryESModule = cms.ESProducer("RPCGeometryESModule",
    fromDD4hep = cms.untracked.bool(False),
    fromDDD = cms.untracked.bool(False)
)


process.SiStripRecHitMatcherESProducer = cms.ESProducer("SiStripRecHitMatcherESProducer",
    ComponentName = cms.string('StandardMatcher'),
    NSigmaInside = cms.double(3.0),
    PreFilter = cms.bool(False)
)


process.StripCPEfromTrackAngleESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('StripCPEfromTrackAngle'),
    ComponentType = cms.string('StripCPEfromTrackAngle'),
    parameters = cms.PSet(
        mLC_P0 = cms.double(-0.326),
        mLC_P1 = cms.double(0.618),
        mLC_P2 = cms.double(0.3),
        mTEC_P0 = cms.double(-1.885),
        mTEC_P1 = cms.double(0.471),
        mTIB_P0 = cms.double(-0.742),
        mTIB_P1 = cms.double(0.202),
        mTID_P0 = cms.double(-1.427),
        mTID_P1 = cms.double(0.433),
        mTOB_P0 = cms.double(-1.026),
        mTOB_P1 = cms.double(0.253),
        maxChgOneMIP = cms.double(6000.0),
        useLegacyError = cms.bool(False)
    )
)


process.TrackerRecoGeometryESProducer = cms.ESProducer("TrackerRecoGeometryESProducer",
    usePhase2Stacks = cms.bool(False)
)


process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName = cms.string('TransientTrackBuilder')
)


process.VolumeBasedMagneticFieldESProducer = cms.ESProducer("DD4hep_VolumeBasedMagneticFieldESProducerFromDB",
    debugBuilder = cms.untracked.bool(False),
    label = cms.untracked.string(''),
    valueOverride = cms.int32(18268)
)


process.ZdcGeometryFromDBEP = cms.ESProducer("ZdcGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.caloSimulationParameters = cms.ESProducer("CaloSimParametersESModule",
    appendToDataLabel = cms.string(''),
    fromDD4hep = cms.bool(True)
)


process.candidateJetProbabilityComputer = cms.ESProducer("CandidateJetProbabilityESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(0.3),
    impactParameterType = cms.int32(0),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    maximumDecayLength = cms.double(5.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5),
    minimumProbability = cms.double(0.005),
    trackIpSign = cms.int32(1),
    trackQualityClass = cms.string('any'),
    useVariableJTA = cms.bool(False)
)


process.ctppsBeamParametersFromLHCInfoESSource = cms.ESProducer("CTPPSBeamParametersFromLHCInfoESSource",
    appendToDataLabel = cms.string(''),
    beamDivX45 = cms.double(0.1),
    beamDivX56 = cms.double(0.1),
    beamDivY45 = cms.double(0.1),
    beamDivY56 = cms.double(0.1),
    lhcInfoLabel = cms.string(''),
    vtxOffsetX45 = cms.double(0.01),
    vtxOffsetX56 = cms.double(0.01),
    vtxOffsetY45 = cms.double(0.01),
    vtxOffsetY56 = cms.double(0.01),
    vtxOffsetZ45 = cms.double(0.01),
    vtxOffsetZ56 = cms.double(0.01),
    vtxStddevX = cms.double(0.02),
    vtxStddevY = cms.double(0.02),
    vtxStddevZ = cms.double(0.02)
)


process.ctppsInterpolatedOpticalFunctionsESSource = cms.ESProducer("CTPPSInterpolatedOpticalFunctionsESSource",
    appendToDataLabel = cms.string(''),
    lhcInfoLabel = cms.string(''),
    opticsLabel = cms.string('')
)


process.ecalSimulationParametersEB = cms.ESProducer("EcalSimParametersESModule",
    appendToDataLabel = cms.string(''),
    fromDD4hep = cms.bool(True),
    name = cms.string('EcalHitsEB')
)


process.ecalSimulationParametersEE = cms.ESProducer("EcalSimParametersESModule",
    appendToDataLabel = cms.string(''),
    fromDD4hep = cms.bool(True),
    name = cms.string('EcalHitsEE')
)


process.ecalSimulationParametersES = cms.ESProducer("EcalSimParametersESModule",
    appendToDataLabel = cms.string(''),
    fromDD4hep = cms.bool(True),
    name = cms.string('EcalHitsES')
)


process.fakeForIdealAlignment = cms.ESProducer("FakeAlignmentProducer",
    appendToDataLabel = cms.string('fakeForIdeal')
)


process.hcalDDDRecConstants = cms.ESProducer("HcalDDDRecConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalDDDSimConstants = cms.ESProducer("HcalDDDSimConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalSimulationConstants = cms.ESProducer("HcalSimulationConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalSimulationParameters = cms.ESProducer("HcalSimParametersESModule",
    appendToDataLabel = cms.string(''),
    fromDD4hep = cms.bool(True)
)


process.hcalTopologyIdeal = cms.ESProducer("HcalTopologyIdealEP",
    Exclude = cms.untracked.string(''),
    MergePosition = cms.untracked.bool(False),
    appendToDataLabel = cms.string('')
)


process.hcal_db_producer = cms.ESProducer("HcalDbProducer",
    dump = cms.untracked.vstring(''),
    file = cms.untracked.string('')
)


process.idealForDigiCSCGeometry = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    debugV = cms.untracked.bool(False),
    fromDD4hep = cms.bool(False),
    fromDDD = cms.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useGangedStripsInME1a = cms.bool(False),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.idealForDigiDTGeometry = cms.ESProducer("DTGeometryESModule",
    DDDetector = cms.ESInputTag("",""),
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    attribute = cms.string('MuStructure'),
    fromDD4hep = cms.bool(False),
    fromDDD = cms.bool(False),
    value = cms.string('MuonBarrelDT')
)


process.idealForDigiTrackerGeometry = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.multipleScatteringParametrisationMakerESProducer = cms.ESProducer("MultipleScatteringParametrisationMakerESProducer",
    appendToDataLabel = cms.string('')
)


process.muonGeometryConstants = cms.ESProducer("MuonGeometryConstantsESModule",
    appendToDataLabel = cms.string(''),
    fromDD4hep = cms.bool(True)
)


process.muonOffsetESProducer = cms.ESProducer("MuonOffsetESProducer",
    appendToDataLabel = cms.string(''),
    fromDD4hep = cms.bool(True),
    names = cms.vstring(
        'MuonCommonNumbering',
        'MuonBarrel',
        'MuonEndcap',
        'MuonBarrelWheels',
        'MuonBarrelStation1',
        'MuonBarrelStation2',
        'MuonBarrelStation3',
        'MuonBarrelStation4',
        'MuonBarrelSuperLayer',
        'MuonBarrelLayer',
        'MuonBarrelWire',
        'MuonRpcPlane1I',
        'MuonRpcPlane1O',
        'MuonRpcPlane2I',
        'MuonRpcPlane2O',
        'MuonRpcPlane3S',
        'MuonRpcPlane4',
        'MuonRpcChamberLeft',
        'MuonRpcChamberMiddle',
        'MuonRpcChamberRight',
        'MuonRpcEndcap1',
        'MuonRpcEndcap2',
        'MuonRpcEndcap3',
        'MuonRpcEndcap4',
        'MuonRpcEndcapSector',
        'MuonRpcEndcapChamberB1',
        'MuonRpcEndcapChamberB2',
        'MuonRpcEndcapChamberB3',
        'MuonRpcEndcapChamberC1',
        'MuonRpcEndcapChamberC2',
        'MuonRpcEndcapChamberC3',
        'MuonRpcEndcapChamberE1',
        'MuonRpcEndcapChamberE2',
        'MuonRpcEndcapChamberE3',
        'MuonRpcEndcapChamberF1',
        'MuonRpcEndcapChamberF2',
        'MuonRpcEndcapChamberF3',
        'MuonRpcEndcapChamberG1',
        'MuonRpcEndcapChamberH1',
        'MuonEndcapStation1',
        'MuonEndcapStation2',
        'MuonEndcapStation3',
        'MuonEndcapStation4',
        'MuonEndcapSubrings',
        'MuonEndcapSectors',
        'MuonEndcapLayers',
        'MuonEndcapRing1',
        'MuonEndcapRing2',
        'MuonEndcapRing3',
        'MuonEndcapRingA',
        'MuonGEMEndcap',
        'MuonGEMEndcap2',
        'MuonGEMSector',
        'MuonGEMChamber'
    )
)


process.siPixelQualityESProducer = cms.ESProducer("SiPixelQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(
        cms.PSet(
            record = cms.string('SiPixelQualityFromDbRcd'),
            tag = cms.string('')
        ),
        cms.PSet(
            record = cms.string('SiPixelDetVOffRcd'),
            tag = cms.string('')
        )
    ),
    appendToDataLabel = cms.string(''),
    siPixelQualityLabel = cms.string(''),
    siPixelQualityLabel_RawToDigi = cms.string('')
)


process.siStripBackPlaneCorrectionDepESProducer = cms.ESProducer("SiStripBackPlaneCorrectionDepESProducer",
    BackPlaneCorrectionDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    BackPlaneCorrectionPeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    )
)


process.siStripGainESProducer = cms.ESProducer("SiStripGainESProducer",
    APVGain = cms.VPSet(
        cms.PSet(
            Label = cms.untracked.string(''),
            NormalizationFactor = cms.untracked.double(1.0),
            Record = cms.string('SiStripApvGainRcd')
        ),
        cms.PSet(
            Label = cms.untracked.string(''),
            NormalizationFactor = cms.untracked.double(1.0),
            Record = cms.string('SiStripApvGain2Rcd')
        )
    ),
    AutomaticNormalization = cms.bool(False),
    appendToDataLabel = cms.string(''),
    printDebug = cms.untracked.bool(False)
)


process.siStripLorentzAngleDepESProducer = cms.ESProducer("SiStripLorentzAngleDepESProducer",
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    ),
    LorentzAngleDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripLorentzAngleRcd')
    ),
    LorentzAnglePeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripLorentzAngleRcd')
    )
)


process.siStripQualityESProducer = cms.ESProducer("SiStripQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(
        cms.PSet(
            record = cms.string('SiStripDetVOffRcd'),
            tag = cms.string('')
        ),
        cms.PSet(
            record = cms.string('SiStripDetCablingRcd'),
            tag = cms.string('')
        ),
        cms.PSet(
            record = cms.string('RunInfoRcd'),
            tag = cms.string('')
        ),
        cms.PSet(
            record = cms.string('SiStripBadChannelRcd'),
            tag = cms.string('')
        ),
        cms.PSet(
            record = cms.string('SiStripBadFiberRcd'),
            tag = cms.string('')
        ),
        cms.PSet(
            record = cms.string('SiStripBadModuleRcd'),
            tag = cms.string('')
        ),
        cms.PSet(
            record = cms.string('SiStripBadStripRcd'),
            tag = cms.string('')
        )
    ),
    PrintDebugOutput = cms.bool(False),
    ReduceGranularity = cms.bool(False),
    ThresholdForReducedGranularity = cms.double(0.3),
    UseEmptyRunInfo = cms.bool(False),
    appendToDataLabel = cms.string('')
)


process.sistripconn = cms.ESProducer("SiStripConnectivity")


process.stripCPEESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('stripCPE'),
    ComponentType = cms.string('SimpleStripCPE'),
    parameters = cms.PSet(

    )
)


process.trackerGeometryDB = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDDD = cms.bool(False)
)


process.trackerNumberingGeometryDB = cms.ESProducer("TrackerGeometricDetESModule",
    appendToDataLabel = cms.string(''),
    fromDD4hep = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.trackerTopology = cms.ESProducer("TrackerTopologyEP",
    appendToDataLabel = cms.string('')
)


process.BTagRecord = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('JetTagComputerRecord')
)


process.DDDetectorESProducerFromDB = cms.ESSource("DDDetectorESProducer",
    appendToDataLabel = cms.string(''),
    fromDB = cms.bool(True),
    label = cms.string('Extended'),
    rootDDName = cms.string('cms:OCMS')
)


process.GlobalTag = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        connectionTimeout = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    DumpStat = cms.untracked.bool(False),
    ReconnectEachRun = cms.untracked.bool(False),
    RefreshAlways = cms.untracked.bool(False),
    RefreshEachRun = cms.untracked.bool(False),
    RefreshOpenIOVs = cms.untracked.bool(False),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    globaltag = cms.string('124X_dataRun3_Prompt_v10'),
    pfnPostfix = cms.untracked.string(''),
    pfnPrefix = cms.untracked.string(''),
    snapshotTime = cms.string(''),
    toGet = cms.VPSet()
)


process.HcalTimeSlewEP = cms.ESSource("HcalTimeSlewEP",
    appendToDataLabel = cms.string('HBHE'),
    timeSlewParametersM2 = cms.VPSet(
        cms.PSet(
            slope = cms.double(-3.178648),
            tmax = cms.double(16.0),
            tzero = cms.double(23.960177)
        ),
        cms.PSet(
            slope = cms.double(-1.5610227),
            tmax = cms.double(10.0),
            tzero = cms.double(11.977461)
        ),
        cms.PSet(
            slope = cms.double(-1.075824),
            tmax = cms.double(6.25),
            tzero = cms.double(9.109694)
        )
    ),
    timeSlewParametersM3 = cms.VPSet(
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(12.2999),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-2.19142),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(0.0),
            tspar2_siPM = cms.double(0.0)
        ),
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(15.5),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-3.2),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(32.0),
            tspar2_siPM = cms.double(0.0)
        ),
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(12.2999),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-2.19142),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(0.0),
            tspar2_siPM = cms.double(0.0)
        ),
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(12.2999),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-2.19142),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(0.0),
            tspar2_siPM = cms.double(0.0)
        )
    )
)


process.HepPDTESSource = cms.ESSource("HepPDTESSource",
    pdtFileName = cms.FileInPath('SimGeneral/HepPDTESSource/data/pythiaparticle.tbl')
)


process.eegeom = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('EcalMappingRcd')
)


process.es_hardcode = cms.ESSource("HcalHardcodeCalibrations",
    GainWidthsForTrigPrims = cms.bool(False),
    HBRecalibration = cms.bool(False),
    HBmeanenergies = cms.FileInPath('CalibCalorimetry/HcalPlugins/data/meanenergiesHB.txt'),
    HBreCalibCutoff = cms.double(100.0),
    HERecalibration = cms.bool(False),
    HEmeanenergies = cms.FileInPath('CalibCalorimetry/HcalPlugins/data/meanenergiesHE.txt'),
    HEreCalibCutoff = cms.double(100.0),
    HFRecalParameterBlock = cms.PSet(
        HFdepthOneParameterA = cms.vdouble(
            0.004123, 0.00602, 0.008201, 0.010489, 0.013379,
            0.016997, 0.021464, 0.027371, 0.034195, 0.044807,
            0.058939, 0.125497
        ),
        HFdepthOneParameterB = cms.vdouble(
            -4e-06, -2e-06, 0.0, 4e-06, 1.5e-05,
            2.6e-05, 6.3e-05, 8.4e-05, 0.00016, 0.000107,
            0.000425, 0.000209
        ),
        HFdepthTwoParameterA = cms.vdouble(
            0.002861, 0.004168, 0.0064, 0.008388, 0.011601,
            0.014425, 0.018633, 0.023232, 0.028274, 0.035447,
            0.051579, 0.086593
        ),
        HFdepthTwoParameterB = cms.vdouble(
            -2e-06, -0.0, -7e-06, -6e-06, -2e-06,
            1e-06, 1.9e-05, 3.1e-05, 6.7e-05, 1.2e-05,
            0.000157, -3e-06
        )
    ),
    HFRecalibration = cms.bool(False),
    SiPMCharacteristics = cms.VPSet(
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(36000)
        ),
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(2500)
        ),
        cms.PSet(
            crosstalk = cms.double(0.17),
            nonlin1 = cms.double(1.00985),
            nonlin2 = cms.double(7.84089e-06),
            nonlin3 = cms.double(2.86282e-10),
            pixels = cms.int32(27370)
        ),
        cms.PSet(
            crosstalk = cms.double(0.196),
            nonlin1 = cms.double(1.00546),
            nonlin2 = cms.double(6.40239e-06),
            nonlin3 = cms.double(1.27011e-10),
            pixels = cms.int32(38018)
        ),
        cms.PSet(
            crosstalk = cms.double(0.17),
            nonlin1 = cms.double(1.00985),
            nonlin2 = cms.double(7.84089e-06),
            nonlin3 = cms.double(2.86282e-10),
            pixels = cms.int32(27370)
        ),
        cms.PSet(
            crosstalk = cms.double(0.196),
            nonlin1 = cms.double(1.00546),
            nonlin2 = cms.double(6.40239e-06),
            nonlin3 = cms.double(1.27011e-10),
            pixels = cms.int32(38018)
        ),
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(0)
        )
    ),
    hb = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.19),
        gainWidth = cms.vdouble(0.0),
        mcShape = cms.int32(125),
        noiseCorrelation = cms.vdouble(0.0),
        noiseThreshold = cms.double(0.0),
        pedestal = cms.double(3.285),
        pedestalWidth = cms.double(0.809),
        photoelectronsToAnalog = cms.double(0.3305),
        qieOffset = cms.vdouble(-0.49, 1.8, 7.2, 37.9),
        qieSlope = cms.vdouble(0.912, 0.917, 0.922, 0.923),
        qieType = cms.int32(0),
        recoShape = cms.int32(105),
        seedThreshold = cms.double(0.1),
        zsThreshold = cms.int32(8)
    ),
    hbUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.01, 0.015),
        doRadiationDamage = cms.bool(True),
        gain = cms.vdouble(0.0006252),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(206),
        noiseCorrelation = cms.vdouble(0.26, 0.254),
        noiseThreshold = cms.double(0.0),
        pedestal = cms.double(17.3),
        pedestalWidth = cms.double(1.5),
        photoelectronsToAnalog = cms.double(40.0),
        qieOffset = cms.vdouble(0.0, 0.0, 0.0, 0.0),
        qieSlope = cms.vdouble(0.05376, 0.05376, 0.05376, 0.05376),
        qieType = cms.int32(2),
        radiationDamage = cms.PSet(
            depVsNeutrons = cms.vdouble(5.543e-10, 8.012e-10),
            depVsTemp = cms.double(0.0631),
            intlumiOffset = cms.double(150),
            intlumiToNeutrons = cms.double(367000000.0),
            temperatureBase = cms.double(20),
            temperatureNew = cms.double(-5)
        ),
        recoShape = cms.int32(208),
        seedThreshold = cms.double(0.1),
        zsThreshold = cms.int32(16)
    ),
    he = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.23),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(125),
        noiseCorrelation = cms.vdouble(0.0),
        noiseThreshold = cms.double(0.0),
        pedestal = cms.double(3.163),
        pedestalWidth = cms.double(0.9698),
        photoelectronsToAnalog = cms.double(0.3305),
        qieOffset = cms.vdouble(-0.38, 2.0, 7.6, 39.6),
        qieSlope = cms.vdouble(0.912, 0.916, 0.92, 0.922),
        qieType = cms.int32(0),
        recoShape = cms.int32(105),
        seedThreshold = cms.double(0.1),
        zsThreshold = cms.int32(9)
    ),
    heUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.01, 0.015),
        doRadiationDamage = cms.bool(True),
        gain = cms.vdouble(0.0006252),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(206),
        noiseCorrelation = cms.vdouble(0.26, 0.254),
        noiseThreshold = cms.double(0.0),
        pedestal = cms.double(17.3),
        pedestalWidth = cms.double(1.5),
        photoelectronsToAnalog = cms.double(40.0),
        qieOffset = cms.vdouble(0.0, 0.0, 0.0, 0.0),
        qieSlope = cms.vdouble(0.05376, 0.05376, 0.05376, 0.05376),
        qieType = cms.int32(2),
        radiationDamage = cms.PSet(
            depVsNeutrons = cms.vdouble(5.543e-10, 8.012e-10),
            depVsTemp = cms.double(0.0631),
            intlumiOffset = cms.double(75),
            intlumiToNeutrons = cms.double(29200000.0),
            temperatureBase = cms.double(20),
            temperatureNew = cms.double(5)
        ),
        recoShape = cms.int32(208),
        seedThreshold = cms.double(0.1),
        zsThreshold = cms.int32(16)
    ),
    hf = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.14, 0.135),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(301),
        noiseCorrelation = cms.vdouble(0.0),
        noiseThreshold = cms.double(0.0),
        pedestal = cms.double(9.354),
        pedestalWidth = cms.double(2.516),
        photoelectronsToAnalog = cms.double(0.0),
        qieOffset = cms.vdouble(-0.87, 1.4, 7.8, -29.6),
        qieSlope = cms.vdouble(0.359, 0.358, 0.36, 0.367),
        qieType = cms.int32(0),
        recoShape = cms.int32(301),
        seedThreshold = cms.double(0.1),
        zsThreshold = cms.int32(-9999)
    ),
    hfUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.14, 0.135),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(301),
        noiseCorrelation = cms.vdouble(0.0),
        noiseThreshold = cms.double(0.0),
        pedestal = cms.double(13.33),
        pedestalWidth = cms.double(3.33),
        photoelectronsToAnalog = cms.double(0.0),
        qieOffset = cms.vdouble(0.0697, -0.7405, 12.38, -671.9),
        qieSlope = cms.vdouble(0.297, 0.298, 0.298, 0.313),
        qieType = cms.int32(1),
        recoShape = cms.int32(301),
        seedThreshold = cms.double(0.1),
        zsThreshold = cms.int32(-9999)
    ),
    ho = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.006, 0.0087),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(201),
        noiseCorrelation = cms.vdouble(0.0),
        noiseThreshold = cms.double(0.0),
        pedestal = cms.double(12.06),
        pedestalWidth = cms.double(0.6285),
        photoelectronsToAnalog = cms.double(4.0),
        qieOffset = cms.vdouble(-0.44, 1.4, 7.1, 38.5),
        qieSlope = cms.vdouble(0.907, 0.915, 0.92, 0.921),
        qieType = cms.int32(0),
        recoShape = cms.int32(201),
        seedThreshold = cms.double(0.1),
        zsThreshold = cms.int32(24)
    ),
    iLumi = cms.double(-1.0),
    killHE = cms.bool(False),
    testHEPlan1 = cms.bool(False),
    testHFQIE10 = cms.bool(False),
    toGet = cms.untracked.vstring('GainWidths'),
    useHBUpgrade = cms.bool(True),
    useHEUpgrade = cms.bool(True),
    useHFUpgrade = cms.bool(True),
    useHOUpgrade = cms.bool(True),
    useIeta18depth1 = cms.bool(False),
    useLayer0Weight = cms.bool(True)
)


process.prefer("es_hardcode")

process.siPixelRecHitsPreSplittingTask = cms.Task(process.siPixelRecHitsPreSplitting)


process.trackSequencePP = cms.Sequence(process.unpackedTracksAndVertices+process.ppTracks)


process.trackSequencePbPb = cms.Sequence(process.unpackedTracksAndVertices+process.PbPbTracks)


process.muonSequence = cms.Sequence(process.unpackedTracksAndVertices+process.unpackedMuons+process.muonAnalyzer)


process.extraFlowJetsData = cms.Sequence(process.PackedPFTowers+process.hiPuRho+process.ak4PFJetsForFlow+process.hiFJRhoFlowModulation)


process.extraFlowJetsMC = cms.Sequence(process.PackedPFTowers+process.hiPuRho+process.hiSignalGenParticles+process.allPartons+process.ak4PFJetsForFlow+process.hiFJRhoFlowModulation)


process.extraJetsData = cms.Sequence(process.PackedPFTowers+process.hiPuRho)


process.extraJetsMC = cms.Sequence(process.PackedPFTowers+process.hiPuRho+process.hiSignalGenParticles+process.allPartons)


process.candidateBtagging = cms.Sequence(process.pfImpactParameterTagInfos+process.pfSecondaryVertexTagInfos+process.pfJetProbabilityBJetTags+process.pfDeepCSVTagInfos+process.pfDeepCSVJetTags)


process.jetsR4 = cms.Sequence(process.akCs0PFJets+process.akCs0PFpatJetCorrFactors+process.akCs0PFpfImpactParameterTagInfos+process.akCs0PFpfSecondaryVertexTagInfos+process.akCs0PFpfDeepCSVTagInfos+process.akCs0PFpfDeepCSVJetTags+process.akCs0PFpfJetProbabilityBJetTags+process.akCs0PFpatJets)


process.forest = cms.Path(process.HiForestInfo+process.hiEvtAnalyzer+process.hltanalysis+process.extraJetsData+process.jetsR4+process.akCs4PFJetAnalyzer+process.hiFJRhoFlowModulationProducer+process.rhoAnalysis+process.randomConeAnalysisR4+process.randomConeAnalysisR2)


process.pclusterCompatibilityFilter = cms.Path(process.clusterCompatibilityFilter)


process.pprimaryVertexFilter = cms.Path(process.primaryVertexFilter)


process.pphfCoincFilter4Th2 = cms.Path(process.phfCoincFilter4Th2)


process.pphfCoincFilter1Th3 = cms.Path(process.phfCoincFilter1Th3)


process.pphfCoincFilter2Th3 = cms.Path(process.phfCoincFilter2Th3)


process.pphfCoincFilter3Th3 = cms.Path(process.phfCoincFilter3Th3)


process.pphfCoincFilter4Th3 = cms.Path(process.phfCoincFilter4Th3)


process.pphfCoincFilter5Th3 = cms.Path(process.phfCoincFilter5Th3)


process.pphfCoincFilter1Th4 = cms.Path(process.phfCoincFilter1Th4)


process.pphfCoincFilter2Th4 = cms.Path(process.phfCoincFilter2Th4)


process.pphfCoincFilter3Th4 = cms.Path(process.phfCoincFilter3Th4)


process.pphfCoincFilter4Th4 = cms.Path(process.phfCoincFilter4Th4)


process.pphfCoincFilter5Th4 = cms.Path(process.phfCoincFilter5Th4)


process.pphfCoincFilter1Th5 = cms.Path(process.phfCoincFilter1Th5)


process.pphfCoincFilter2Th5 = cms.Path(process.phfCoincFilter2Th5)


process.pphfCoincFilter3Th5 = cms.Path(process.phfCoincFilter3Th5)


process.pphfCoincFilter4Th5 = cms.Path(process.phfCoincFilter4Th5)


process.pphfCoincFilter5Th5 = cms.Path(process.phfCoincFilter5Th5)


process.pphfCoincFilter1Th6 = cms.Path(process.phfCoincFilter1Th6)


process.pphfCoincFilter2Th6 = cms.Path(process.phfCoincFilter2Th6)


process.pphfCoincFilter3Th6 = cms.Path(process.phfCoincFilter3Th6)


process.pphfCoincFilter4Th6 = cms.Path(process.phfCoincFilter4Th6)


process.pphfCoincFilter5Th6 = cms.Path(process.phfCoincFilter5Th6)


process.pAna = cms.EndPath(process.skimanalysis)


