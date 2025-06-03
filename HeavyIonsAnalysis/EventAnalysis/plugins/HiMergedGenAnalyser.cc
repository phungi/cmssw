// -*- C++ -*-
//
// Package:    HiMergedGenAnalyser
// Class:      HiMergedGenAnalyser
//
/**\class HiMergedGenAnalyser HiMergedGenAnalyser.cc

   Description: Analyzer that studies (HI) gen event info in miniAOD

   Implementation:
   This analyzer is copied from its AOD counterpart https://github.com/CmsHI/cmssw/blob/2c806f88506f7ef732b725142ae85750a31dc646/HeavyIonsAnalysis/EventAnalysis/src/HiEvtAnalyzer.cc and adapted for gen info in miniAOD
*/

// system include files
#include <memory>
#include <string>
#include <vector>

// user include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/HiGenData/interface/GenHIEvent.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "HepMC/GenEvent.h"
#include "HepMC/HeavyIon.h"

// root include file
#include "TFile.h"
#include "TTree.h"

using namespace std;

static const Int_t ETABINS = 3;  // Fix also in branch string

//
// class decleration
//

struct GenStruct {
  Int_t event;
  Float_t b;
  Float_t npart;
  Float_t ncoll;
  Float_t nhard;
  Float_t phi0;
  Float_t scale;

  Int_t n[ETABINS];
  Float_t ptav[ETABINS];

  Int_t mult;
  std::vector<Float_t> pt;
  std::vector<Float_t> eta;
  std::vector<Float_t> phi;
  std::vector<Int_t> pdg;
  std::vector<Int_t> chg;
  std::vector<Int_t> sube;
  std::vector<Int_t> sta;
  std::vector<Int_t> matchingID;
  std::vector<Int_t> nMothers;
  std::vector<std::vector<Int_t>> motherIndex;
  std::vector<std::vector<Int_t>> custom_motherIndex;
  std::vector<std::vector<Int_t>> custom_daughterIndex;
  // Float_t tge;
  // std::vector<std::vector<Int_t>> haha;
  std::vector<Int_t> nDaughters;
  std::vector<std::vector<Int_t>> daughterIndex;
  // std::vector<Bool_t> isFromHardScatter;

  Float_t vx;
  Float_t vy;
  Float_t vz;
  Float_t vr;
};

class HiMergedGenAnalyser : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
public:
  explicit HiMergedGenAnalyser(const edm::ParameterSet&);
  ~HiMergedGenAnalyser() override;

private:
  HiMergedGenAnalyser(){}
  edm::InputTag prunedInTag;
  edm::InputTag packedInTag;
  void beginRun(const edm::Run& run, const edm::EventSetup& iSetup) override;
  void endRun(const edm::Run& run, const edm::EventSetup& iSetup) override;
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  Bool_t passesSelection(const reco::Candidate* p);
  vector<int> getMotherIdx(std::vector<const reco::Candidate*> parts, const reco::Candidate *);
  vector<int> getDaughterIdx(std::vector<const reco::Candidate*> parts, const reco::Candidate *);
  void getMotherIdxCustom(vector<int> &, std::vector<const reco::Candidate*> parts, const reco::Candidate *);
  void getDaughterIdxCustom(vector<int> &, std::vector<const reco::Candidate*> parts, const reco::Candidate *);



  // ----------member data ---------------------------

  edm::EDGetTokenT<edm::SimVertexContainer> g4Label;

  TTree* hydjetTree_;
  GenStruct hev_;

  Bool_t doVertex_;
  Bool_t useHepMCProduct_;
  Bool_t doHI_;
  Bool_t doParticles_;
  std::vector<int> motherDaughterPDGsToSave_;

  Double_t etaMax_;
  Double_t ptMin_;
  Bool_t chargedOnly_;
  Bool_t stableOnly_;
  // edm::InputTag packedInTag;
  // edm::InputTag prunedInTag;

  // edm::InputTag prunedGenParticlesSrc_;
  // edm::InputTag packedGenParticlesSignalSrc_;
  edm::EDGetTokenT<edm::HepMCProduct> src_;
  // edm::EDGetTokenT<reco::CandidateView>
  edm::EDGetTokenT<reco::CandidateView> prunedGenParticlesSrc_;
  edm::EDGetTokenT<reco::CandidateView> packedGenParticlesSignalSrc_;
  edm::EDGetTokenT<edm::GenHIEvent> genHIsrc_;
  edm::ESGetToken<HepPDT::ParticleDataTable, PDTRecord> tok_pdt_;
  edm::Service<TFileService> f;
};
//
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HiMergedGenAnalyser::HiMergedGenAnalyser(const edm::ParameterSet& iConfig):
  prunedInTag(iConfig.getParameter<edm::InputTag>("prunedGenParticlesSrc")),
  packedInTag(iConfig.getParameter<edm::InputTag>("packedGenParticlesSignalSrc")){
  //now do what ever initialization is needed
  useHepMCProduct_ = iConfig.getUntrackedParameter<Bool_t>("useHepMCProduct", false);
  doHI_ = iConfig.getUntrackedParameter<Bool_t>("doHI", true);
  
  doVertex_ = iConfig.getUntrackedParameter<Bool_t>("doVertex", false);
  etaMax_ = iConfig.getUntrackedParameter<Double_t>("etaMax", 2);
  ptMin_ = iConfig.getUntrackedParameter<Double_t>("ptMin", 0);
  chargedOnly_ = iConfig.getUntrackedParameter<Bool_t>("chargedOnly", false);
  stableOnly_ = iConfig.getUntrackedParameter<Bool_t>("stableOnly", false);
  if (useHepMCProduct_) {
    src_ = consumes<edm::HepMCProduct>(iConfig.getUntrackedParameter<edm::InputTag>("src", edm::InputTag("generator")));
  } else {
    //these are only decays of interest, stable + unstable, with some selections on the daughters (e.g. missing gammas in final state below 10 GeV etc.)
    // prunedGenParticlesSrc_ = consumes<reco::CandidateView>(prunedGenParticlesSrc_);
        // consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("prunedGenParticlesSrc"));
    //these contain only stable particles which go to the jet clustering algorithm
    // packedGenParticlesSignalSrc_ = consumes<reco::CandidateView>(signalGenParticleSrc_);
        // consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("signalGenParticleSrc"));
  }
  prunedGenParticlesSrc_ = consumes<reco::CandidateView>(prunedInTag);
  packedGenParticlesSignalSrc_ = consumes<reco::CandidateView>(packedInTag);
  if (doHI_) {
    genHIsrc_ =
        consumes<edm::GenHIEvent>(iConfig.getUntrackedParameter<edm::InputTag>("genHiSrc", edm::InputTag("heavyIon")));
  }
  tok_pdt_ = esConsumes<HepPDT::ParticleDataTable, PDTRecord>();
  doParticles_ = iConfig.getUntrackedParameter<Bool_t>("doParticles", true);
  vector<int> defaultPDGs;
  motherDaughterPDGsToSave_ = iConfig.getUntrackedParameter<std::vector<int>>("motherDaughterPDGsToSave", defaultPDGs);

  if (doVertex_) {
    g4Label = consumes<edm::SimVertexContainer>(iConfig.getUntrackedParameter<std::string>("ModuleLabel", "g4SimHits"));
  }
}

HiMergedGenAnalyser::~HiMergedGenAnalyser() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

Bool_t HiMergedGenAnalyser::passesSelection(const reco::Candidate* p){
  if (stableOnly_ && p->status() != 1) return false;
  if (p->pt() < ptMin_) return false;
  // if (fabs(p->eta()) > etaMax_) return false;
  if (chargedOnly_ && p->charge() == 0) return false;
  return true;
}
void HiMergedGenAnalyser::getMotherIdxCustom(vector<int> &motherArr, std::vector<const reco::Candidate*> parts, const reco::Candidate *pin){
  // std::cout << "we are in the custom function " << std::endl;
  unsigned int nMo=pin->numberOfMothers();
  // std::cout << "Number of mothers " << nMo << std::endl;
  if ( !nMo or nMo==0 ) return;
  for(unsigned int i{0}; i < nMo; ++i){
    const reco::Candidate* mo = pin->mother(i);
    Bool_t keep = passesSelection(mo);
    // std::cout << "keeping this mother particle!" << std::endl;
    if(keep){
      // std::cout << "Go in loop of size " << parts->size() << std::endl;
      // double min_dr = 99999;
      double closest_pdgid = 0;
      // std::cout << pin->numberOfMothers() << " mothers " << std::endl;
      // std::cout << "collection in function " << parts.size() << " mother count " << nMo << std::endl;
      for (UInt_t j = 0; j < parts.size(); ++j) {
        const reco::Candidate * p = (parts)[j];
        // unsigned int nDa = p->numberOfDaughters();
        for (unsigned int idx = 0; idx < p->numberOfDaughters(); idx++) {
          // double rij = fabs(p->daughter(idx)->eta() - pin->eta())*fabs(p->daughter(idx)->eta() - pin->eta()) + fabs(p->daughter(idx)->phi() - pin->phi())*fabs(p->daughter(idx)->phi() - pin->phi());
          // if(rij < min_dr){
          //   min_dr = rij;
          //   closest_pdgid = p->daughter(idx)->pdgId();
          // }
          // don't use the dR dpT method, just check if pointers are the same
          // if (fabs(p->daughter(idx)->pt() - pin->pt()) < 0.001 && fabs(p->daughter(idx)->eta() - pin->eta()) < 0.001 && fabs(p->daughter(idx)->phi() - pin->phi()) < 0.001){
          // }
          if(p->daughter(idx) == pin){
            // std::cout << "coming from particle " << pin->pdgId() << " saving idx " << j << std::endl;
            if( find(motherArr.begin(), motherArr.end(), j) == motherArr.end() ){
              motherArr.push_back(j);
              break;
            }
          }
        }
      }
    }
    //if we don't keep the daughter, record that we skipped a particle 
    else{
      motherArr.push_back(-999);
    }
    getMotherIdxCustom(motherArr, parts, mo);
  }
}

void HiMergedGenAnalyser::getDaughterIdxCustom(vector<int> &daughterArr, std::vector<const reco::Candidate*> parts, const reco::Candidate *pin){
  // std::cout << "we are in the custom function " << std::endl;
  unsigned int nDa=pin->numberOfDaughters();
  // std::cout << "Number of mothers " << nMo << std::endl;
  if ( !nDa or nDa==0 ){
    daughterArr.push_back(-999);
    return;
  }
  for(unsigned int i{0}; i < 1; ++i){
    const reco::Candidate* da = pin->daughter(i);
    bool keep = passesSelection(da);
    if(keep){
      // std::cout << "Go in loop of size " << parts->size() << std::endl;
      // double min_dr = 99999;
      // double closest_pdgid = 0;
      // std::cout << pin->numberOfMothers() << " mothers " << std::endl;
      // std::cout << "collection in function " << parts.size() << " daughter count " << nDa << std::endl;
      for (UInt_t j = 0; j < parts.size(); ++j) {
        const reco::Candidate * p = (parts)[j];
        // unsigned int nMo = p->numberOfMothers();
        for (unsigned int idx = 0; idx < p->numberOfMothers(); idx++){
          // double rij = fabs(p->daughter(idx)->eta() - pin->eta())*fabs(p->daughter(idx)->eta() - pin->eta()) + fabs(p->daughter(idx)->phi() - pin->phi())*fabs(p->daughter(idx)->phi() - pin->phi());
          // if(rij < min_dr){
          //   min_dr = rij;
          //   closest_pdgid = p->daughter(idx)->pdgId();
          // }
          // don't use the dR dpT method, just check if pointers are the same
          // if (fabs(p->daughter(idx)->pt() - pin->pt()) < 0.001 && fabs(p->daughter(idx)->eta() - pin->eta()) < 0.001 && fabs(p->daughter(idx)->phi() - pin->phi()) < 0.001){
          // }
          if(p->mother(idx) == pin){
            // std::cout << "coming from particle " << pin->pdgId() << " with pt=" << pin->pt() << " saving idx " << j << " with id " << p->pdgId() << "and pt=" << p->pt() << std::endl;
            if( find(daughterArr.begin(), daughterArr.end(), j) == daughterArr.end() ){
              daughterArr.push_back(j);
            }
          }
        }
      }
    }
    //if we don't keep the daughter, record that we skipped a particle 
    else{
      daughterArr.push_back(-999);
    }
    getDaughterIdxCustom(daughterArr, parts, da);
  }
}

vector<int> HiMergedGenAnalyser::getMotherIdx(std::vector<const reco::Candidate*> parts, const reco::Candidate *pin) {
  vector<int> motherArr;
  if (!motherDaughterPDGsToSave_.empty()) {
    for (UInt_t i = 0; i < parts.size(); ++i) {
      const reco::Candidate* p = (parts)[i];
      if (stableOnly_ && p->status() != 1)
        continue;
      if (p->pt() < ptMin_)
        continue;
      if (chargedOnly_ && p->charge() == 0)
        continue;
      bool saveFlag = false;
      for (unsigned int ipdg = 0; ipdg < motherDaughterPDGsToSave_.size(); ipdg++) {
        if (p->pdgId() == motherDaughterPDGsToSave_.at(ipdg))
          saveFlag = true;
      }
      if (!motherDaughterPDGsToSave_.empty() && saveFlag != true)
        continue;  //save all particles in vector unless vector is empty, then save all particles
      if (p->status() == 3)
        continue;  //don't match to the initial collision particles
      for (unsigned int idx = 0; idx < p->numberOfDaughters(); idx++) {
        //if (p->daughter(idx)->pt()*p->daughter(idx)->eta()*p->daughter(idx)->phi() == pin->pt()*pin->eta()*pin->phi()) motherArr.push_back(i);
        if (fabs(p->daughter(idx)->pt() - pin->pt()) < 0.001 && fabs(p->daughter(idx)->eta() - pin->eta()) < 0.001 &&
            fabs(p->daughter(idx)->phi() - pin->phi()) < 0.001)
          motherArr.push_back(i);
      }
    }
  }
  if (motherArr.empty())
    motherArr.push_back(-999);
  return motherArr;
}

// //----------------------------------------------------------
vector<int> HiMergedGenAnalyser::getDaughterIdx(std::vector<const reco::Candidate*> parts, const reco::Candidate *pin) {
  vector<int> daughterArr;
  if (!motherDaughterPDGsToSave_.empty()) {
    for (UInt_t i = 0; i < parts.size(); ++i) {
      const reco::Candidate* p = (parts)[i];
      if (stableOnly_ && p->status() != 1)
        continue;
      if (p->pt() < ptMin_)
        continue;
      if (chargedOnly_ && p->charge() == 0)
        continue;
      bool saveFlag = false;
      for (unsigned int ipdg = 0; ipdg < motherDaughterPDGsToSave_.size(); ipdg++) {
        if (p->pdgId() == motherDaughterPDGsToSave_.at(ipdg))
          saveFlag = true;
      }
      if (!motherDaughterPDGsToSave_.empty() && saveFlag != true)
        continue;  //save all particles in vector unless vector is empty, then save all particles
      if (p->status() == 3)
        continue;  //don't match to the initial collision particles
      for (unsigned int idx = 0; idx < p->numberOfMothers(); idx++) {
        //if (p->mother(idx)->pt()*p->mother(idx)->eta()*p->mother(idx)->phi() == pin->pt()*pin->eta()*pin->phi()) daughterArr.push_back(i);
        if (fabs(p->mother(idx)->pt() - pin->pt()) < 0.001 && fabs(p->mother(idx)->eta() - pin->eta()) < 0.001 &&
            fabs(p->mother(idx)->phi() - pin->phi()) < 0.001)
          daughterArr.push_back(i);
      }
    }
  }
  if (daughterArr.empty())
    daughterArr.push_back(-999);
  return daughterArr;
}

// ------------ method called to for each event  ------------
void HiMergedGenAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace HepMC;

  const HepPDT::ParticleDataTable* pdt = &iSetup.getData(tok_pdt_);

  hev_.pt.clear();
  hev_.eta.clear();
  hev_.phi.clear();
  hev_.pdg.clear();
  hev_.chg.clear();
  hev_.sube.clear();
  hev_.sta.clear();
  hev_.matchingID.clear();
  hev_.nMothers.clear();
  hev_.motherIndex.clear();
  hev_.custom_motherIndex.clear();
  hev_.custom_daughterIndex.clear();
  hev_.nDaughters.clear();
  hev_.daughterIndex.clear();
  // hev_.haha.clear();
  // hev_.isFromHardScatter.clear();

  hev_.event = iEvent.id().event();
  for (Int_t ieta = 0; ieta < ETABINS; ++ieta){
    hev_.n[ieta] = 0;
    hev_.ptav[ieta] = 0;
  }
  hev_.mult = 0;

  Double_t phi0 = 0;
  Double_t b = -1;
  Double_t scale = -1;
  Int_t npart = -1;
  Int_t ncoll = -1;
  Int_t nhard = -1;
  Double_t vx = -99;
  Double_t vy = -99;
  Double_t vz = -99;
  Double_t vr = -99;
  const GenEvent* evt;
  // if (useHepMCProduct_) {
  //   Handle<edm::HepMCProduct> mc;
  //   iEvent.getByToken(src_, mc);
  //   evt = mc->GetEvent();
  //   scale = evt->event_scale();

  //   const HeavyIon* hi = evt->heavy_ion();
  //   if (hi) {
  //     b = hi->impact_parameter();
  //     npart = hi->Npart_proj() + hi->Npart_targ();
  //     ncoll = hi->Ncoll();
  //     nhard = hi->Ncoll_hard();
  //     phi0 = hi->event_plane_angle();
  //   }

  //   HepMC::GenEvent::particle_const_iterator begin = evt->particles_begin();
  //   HepMC::GenEvent::particle_const_iterator end = evt->particles_end();
  //   int nparticles = -1;
  //   for (HepMC::GenEvent::particle_const_iterator it = begin; it != end; ++it) {
  //     nparticles++;
  //     if ((*it)->momentum().perp() < ptMin_)
  //       continue;
  //     if (fabs((*it)->momentum().eta()) > etaMax_)
  //       continue;
  //     Int_t pdg_id = (*it)->pdg_id();
  //     Float_t eta = (*it)->momentum().eta();
  //     Float_t phi = (*it)->momentum().phi();
  //     Float_t pt = (*it)->momentum().perp();
  //     const ParticleData* part = pdt->particle(pdg_id);
  //     Int_t charge = static_cast<Int_t>(part->charge());
  //     if (chargedOnly_ && charge == 0)
  //       continue;

  //     hev_.pt.push_back(pt);
  //     hev_.eta.push_back(eta);
  //     hev_.phi.push_back(phi);
  //     hev_.pdg.push_back(pdg_id);
  //     hev_.chg.push_back(charge);
  //     hev_.sta.push_back((*it)->status());
  //     hev_.matchingID.push_back(nparticles);

  //     eta = fabs(eta);
  //     Int_t etabin = 0;
  //     if (eta > 0.5)
  //       etabin = 1;
  //     if (eta > 1.)
  //       etabin = 2;
  //     if (eta < 2.) {
  //       hev_.ptav[etabin] += pt;
  //       ++(hev_.n[etabin]);
  //     }
  //     ++(hev_.mult);
  //   }
  // }
  // else {

    edm::Handle<reco::CandidateView> parts;
    iEvent.getByToken(packedGenParticlesSignalSrc_, parts);
    edm::Handle<reco::CandidateView> pruned;
    iEvent.getByToken(prunedGenParticlesSrc_, pruned);

    std::vector<const reco::Candidate*> combined;
    for(auto iter=parts->begin(); iter!=parts->end(); ++iter){
      combined.push_back(&*iter);
    }
    for(auto iter=pruned->begin(); iter!=pruned->end(); ++iter){
      //skip stable pruned particles
      if(iter->status()==1) continue;
      combined.push_back(&*iter);
    }

    // std::cout << "Combined container size " << combined.size() << std::endl;
    //clear out the unwanted particles 
    std::vector<const reco::Candidate*> particles;
    for (UInt_t i = 0; i < combined.size(); ++i){
      const reco::Candidate* p = (combined)[i];
      Bool_t keep = passesSelection(p);
      if(keep) particles.push_back(p);
    }
    // std::cout << "Reduced to container size " << particles.size() << std::endl;
    for (UInt_t i = 0; i < particles.size(); ++i){
      // const reco::GenParticle& p = (*parts)[i];
      const reco::Candidate* p = (particles)[i];
      // if (stableOnly_ && p->status() != 1)
      //   continue;
      // if (p->pt() < ptMin_)
      //   continue;
      // if (fabs(p->eta()) > etaMax_)
      //   continue;
      // if (chargedOnly_ && p->charge() == 0)
      //   continue;
      unsigned int daN = p->numberOfDaughters();
      // std::cout << "Pushing back combined particle: pt=" << p->pt() << " eta=" << p->eta() << " phi=" << p->phi() << std::endl;
      hev_.pt.push_back(p->pt());
      hev_.eta.push_back(p->eta());
      hev_.phi.push_back(p->phi());
      hev_.pdg.push_back(p->pdgId());
      hev_.chg.push_back(p->charge());

      unsigned int nDa = p->numberOfDaughters();
      for(unsigned int i{0}; i < nDa; ++i){
        // std::cout << "particle of " << p->pdgId() << " with daughter " << p->daughter(i)->pdgId() << " in pos " << i << " with pt=" << p->daughter(i)->pt() << std::endl; 
      }
      //these are final state particles anyway, they don't come from the hard scatter
      // hev_.isFromHardScatter.push_back(0);
      // collisionId_ is not kept in pat::PackedGenParticle, use "packedGenParticlesSignal" (added by https://github.com/cms-sw/cmssw/pull/32668/) to tag particles from signal process
      // if (hasSignalPackedGen) {
      //   int tmpSube = 1;
      //   for (auto pSig = signalPackedGenParticles->begin(); pSig != signalPackedGenParticles->end(); ++pSig) {
      //     if (&(*pSig) == &(*parts)[i]){
      //       tmpSube = 0;
      //       break;
      //     }
      //   }
      //   hev_.sube.push_back(tmpSube);
      // } else {
      //   hev_.sube.push_back(-999);
      // }
      hev_.sta.push_back(p->status());
      hev_.matchingID.push_back(i);
      hev_.nMothers.push_back(p->numberOfMothers());
      vector<int> tempMothers = getMotherIdx(particles, p);
      hev_.motherIndex.push_back(tempMothers);

      vector<int> tempcustomMothers = {};
      getMotherIdxCustom(tempcustomMothers, particles, p);
      hev_.custom_motherIndex.push_back(tempcustomMothers);

      vector<int> tempcustomDaughters = {};
      getDaughterIdxCustom(tempcustomDaughters, particles, p);
      hev_.custom_daughterIndex.push_back(tempcustomDaughters);

      hev_.nDaughters.push_back(p->numberOfDaughters());
      vector<int> tempDaughters = getDaughterIdx(particles, p);
      hev_.daughterIndex.push_back(tempDaughters);
      Double_t eta = fabs(p->eta());

      Int_t etabin = 0;
      if (eta > 0.5)
        etabin = 1;
      if (eta > 1.)
        etabin = 2;
      if (eta < 2.) {
        hev_.ptav[etabin] += p->pt();
        ++(hev_.n[etabin]);
      }
      ++(hev_.mult);
    }

    if (doHI_) {
      edm::Handle<edm::GenHIEvent> higen;
      iEvent.getByToken(genHIsrc_, higen);

      b = higen->b();
      npart = higen->Npart();
      ncoll = higen->Ncoll();
      nhard = higen->Nhard();
      phi0 = higen->evtPlane();
    }
  // }

  if (doVertex_) {
    edm::Handle<edm::SimVertexContainer> simVertices;
    iEvent.getByToken(g4Label, simVertices);

    if (!simVertices.isValid())
      throw cms::Exception("FatalError") << "No vertices found\n";

    edm::SimVertexContainer::const_iterator it = simVertices->begin();
    if (it != simVertices->end()) {
      SimVertex vertex = (*it);
      vx = vertex.position().x();
      vy = vertex.position().y();
      vz = vertex.position().z();
      vr = vertex.position().rho();
    }
  }

  for (Int_t i = 0; i < 3; ++i) {
    hev_.ptav[i] = hev_.ptav[i] / hev_.n[i];
  }

  hev_.b = b;
  hev_.scale = scale;
  hev_.npart = npart;
  hev_.ncoll = ncoll;
  hev_.nhard = nhard;
  hev_.phi0 = phi0;
  hev_.vx = vx;
  hev_.vy = vy;
  hev_.vz = vz;
  hev_.vr = vr;

  hydjetTree_->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void HiMergedGenAnalyser::beginRun(const edm::Run& run, const edm::EventSetup& iSetup) {}

// ------------ method called once each job just after finishing event loop  ------------
void HiMergedGenAnalyser::endRun(const edm::Run& run, const edm::EventSetup& iSetup) {}

void HiMergedGenAnalyser::beginJob() {
  hydjetTree_ = f->make<TTree>("hi", "Tree of Hi gen Event");
  hydjetTree_->Branch("event", &hev_.event, "event/I");
  hydjetTree_->Branch("b", &hev_.b, "b/F");
  hydjetTree_->Branch("npart", &hev_.npart, "npart/F");
  hydjetTree_->Branch("ncoll", &hev_.ncoll, "ncoll/F");
  hydjetTree_->Branch("nhard", &hev_.nhard, "nhard/F");
  hydjetTree_->Branch("phi0", &hev_.phi0, "phi0/F");
  hydjetTree_->Branch("scale", &hev_.scale, "scale/F");

  hydjetTree_->Branch("n", hev_.n, "n[3]/I");
  hydjetTree_->Branch("ptav", hev_.ptav, "ptav[3]/F");

  if (doParticles_) {
    hydjetTree_->Branch("mult", &hev_.mult, "mult/I");
    hydjetTree_->Branch("pt", &hev_.pt);
    hydjetTree_->Branch("eta", &hev_.eta);
    hydjetTree_->Branch("phi", &hev_.phi);
    hydjetTree_->Branch("pdg", &hev_.pdg);
    hydjetTree_->Branch("chg", &hev_.chg);
    hydjetTree_->Branch("matchingID", &hev_.matchingID);
    hydjetTree_->Branch("nMothers", &hev_.nMothers);
    hydjetTree_->Branch("motherIdx", &hev_.motherIndex);
    hydjetTree_->Branch("custom_motherIdx", &hev_.custom_motherIndex);
    hydjetTree_->Branch("custom_daughterIdx", &hev_.custom_daughterIndex);
    hydjetTree_->Branch("nDaughters", &hev_.nDaughters);
    hydjetTree_->Branch("daughterIdx", &hev_.daughterIndex);
    if (!stableOnly_) {
      hydjetTree_->Branch("sta", &hev_.sta);
    }
    hydjetTree_->Branch("sube", &hev_.sube);

    hydjetTree_->Branch("vx", &hev_.vx, "vx/F");
    hydjetTree_->Branch("vy", &hev_.vy, "vy/F");
    hydjetTree_->Branch("vz", &hev_.vz, "vz/F");
    hydjetTree_->Branch("vr", &hev_.vr, "vr/F");
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void HiMergedGenAnalyser::endJob() {}

//define this as a plug-in
DEFINE_FWK_MODULE(HiMergedGenAnalyser);
