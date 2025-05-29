import FWCore.ParameterSet.Config as cms

HiMergedGenParticleAna = cms.EDAnalyzer(
    'HiMergedGenAnalyser',
    doVertex = cms.untracked.bool(False),
    etaMax = cms.untracked.double(2.5),
    ptMin = cms.untracked.double(2),
    chargedOnly = cms.untracked.bool(False),
    stableOnly = cms.untracked.bool(False),
    src = cms.untracked.InputTag("generator"),
    prunedGenParticlesSrc = cms.InputTag("prunedGenParticles"),
    packedGenParticlesSignalSrc = cms.InputTag("packedGenParticlesSignal"),
    # genParticleSrc = cms.InputTag("packedGenParticles"),
    # signalGenParticleSrc = cms.InputTag("packedGenParticlesSignal"),
    genHIsrc = cms.untracked.InputTag("heavyIon"),
    doParticles = cms.untracked.bool(True),
    doHI = cms.untracked.bool(False)  ## Relevant info (its AOD counterpart is edm::GenHIEvent "heavyIon") is missing currently.
    )
