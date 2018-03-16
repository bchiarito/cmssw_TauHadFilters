// -*- C++ -*-
//
// Package:    TauHadFilters/ZtoTauHadRecoSelector
// Class:      ZtoTauHadRecoSelector
// 
/**\class ZtoTauHadRecoSelector ZtoTauHadRecoSelector.cc TauHadFilters/ZtoTauHadRecoSelector/plugins/ZtoTauHadRecoSelector.cc

 Description: [one line class summary]

*/
//
// Original Author:  Brandon Chiarito
//         Created:  Wed, 05 Jul 2017 22:31:58 GMT
//
//

// cpp includes
#include <memory>
#include <math.h>
#include <algorithm>
#include <iostream>
// ROOT includes
#include "TVector2.h"
#include "TH1.h"
#include "TTree.h"
// cmssw framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
// pat includes
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
// edm utilities includes
#include "DataFormats/Math/interface/deltaR.h"

// namesspace defaults
using std::vector;
using std::string;

//
// class declaration
//
class ZtoTauHadRecoSelector : public edm::stream::EDFilter<> {
   public:
      explicit ZtoTauHadRecoSelector(const edm::ParameterSet&);
      ~ZtoTauHadRecoSelector();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
      virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // Configuration Parameters

      // EDM Collection lables
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      edm::EDGetTokenT<pat::TauCollection> tauToken_;
      edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
      edm::EDGetTokenT<pat::JetCollection> jetToken_;
      edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
      edm::EDGetTokenT<pat::METCollection> metToken_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfcandsToken_;
      edm::EDGetTokenT<reco::BeamSpot> beamToken_;
      edm::EDGetTokenT<double> rhoToken_;
};

//
// constructors
//
ZtoTauHadRecoSelector::ZtoTauHadRecoSelector(const edm::ParameterSet& iConfig)
{
  vtxToken_ = consumes<vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));
  muonToken_ = consumes<std::vector<pat::Muon>>(edm::InputTag("slimmedMuons"));
  electronToken_ = consumes<std::vector<pat::Electron>>(edm::InputTag("slimmedElectrons"));
  tauToken_ = consumes<std::vector<pat::Tau>>(edm::InputTag("slimmedTaus"));
  photonToken_ = consumes<std::vector<pat::Photon>>(edm::InputTag("slimmedPhotons"));
  jetToken_ = consumes<std::vector<pat::Jet>>(edm::InputTag("slimmedJets"));
  fatjetToken_ = consumes<std::vector<pat::Jet>>(edm::InputTag("slimmedJetsAK8"));
  metToken_ = consumes<std::vector<pat::MET>>(edm::InputTag("slimmedMETs"));
  pfcandsToken_ = consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates"));
  beamToken_ = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  rhoToken_ = consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralChargedPileUp"));
}

//
// destructor
//
ZtoTauHadRecoSelector::~ZtoTauHadRecoSelector()
{
}

// ------------ method called on each new Event  ------------
bool
ZtoTauHadRecoSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // get event content
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  const reco::Vertex &PV = vertices->front();

  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);
 
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);

  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);
  pat::MET MET = (*mets)[0];

  edm::Handle<double> rho;
  iEvent.getByToken(rhoToken_, rho);

  // debug: verify both accessors for muon tag are equivilant
  for (const pat::Muon &muon : *muons) {
    if ( muon.isGlobalMuon() != muon.muonID("AllGlobalMuons") )
      std::cout << "disagree global: " << muon.isGlobalMuon() << ", " << muon.muonID("AllGlobalMuons") << std::endl;
    if ( muon.isTrackerMuon() != muon.muonID("AllTrackerMuons") )
      std::cout << "diagree tracker: " << muon.isTrackerMuon() << ", " << muon.muonID("AllTrackerMuons") << std::endl; 
  }

  // muons
  vector<const pat::Muon *> passedMuons;
  for (const pat::Muon &muon : *muons) {
    if (muon.pt() > 19.0 &&
        fabs(muon.eta()) < 2.1 &&
        (muon.chargedHadronIso() + muon.neutralHadronIso() + muon.photonIso())/muon.pt() - 0.5 * (*rho) < 0.1 &&
        muon.muonBestTrack()->dz(PV.position()) < 0.2 &&
        abs(muon.muonBestTrack()->dxy(PV.position())) < 0.045 &&
        muon.isMediumMuon()
        )
      passedMuons.push_back(&muon);
  }

  // extra lepton veto
  bool skip = false;
  bool weakMuon = false;
  for (const pat::Muon &muon : *muons) {
    if (muon.pt() > 10.0 &&
        fabs(muon.eta()) < 2.4 &&
        (muon.chargedHadronIso() + muon.neutralHadronIso() + muon.photonIso())/muon.pt() - 0.5 * (*rho) < 0.3 &&
        muon.muonBestTrack()->dz(PV.position()) < 0.2 &&
        abs(muon.muonBestTrack()->dxy(PV.position())) < 0.045 &&
        muon.isMediumMuon()
        ) {

      if (muon.pt() > 19.0 &&
          fabs(muon.eta()) < 2.1 &&
          (muon.chargedHadronIso() + muon.neutralHadronIso() + muon.photonIso())/muon.pt() - 0.5 * (*rho) < 0.1 &&
          muon.muonBestTrack()->dz(PV.position()) < 0.2 &&
          abs(muon.muonBestTrack()->dxy(PV.position())) < 0.045 &&
          muon.isMediumMuon()
          ) {
        if (skip) weakMuon = true;
        if (!skip) skip = true;
      } else {
        weakMuon = true;
      }
    }
  }
  bool extraMuon = false;
  if (passedMuons.size() > 0 && weakMuon) extraMuon = true;
  bool extraElectron = false;
  for (const pat::Electron &electron : *electrons) {
    if (electron.pt() > 10.0 &&
        fabs(electron.eta()) < 2.5 &&
        electron.gsfTrack()->dz(PV.position()) < 0.2 &&
        abs(electron.gsfTrack()->dxy(PV.position())) < 0.045 &&
        electron.passConversionVeto() &&
        electron.gsfTrack()->lost() <= 1)
      extraElectron = true;
  }
  bool diMuon = false;
  for (const pat::Muon &muon1 : *muons) {
    for (const pat::Muon &muon2 : *muons) {
      double dR = deltaR(muon1.eta(), muon1.phi(), muon2.eta(), muon2.phi());
      if (dR < 0.15) continue;
      if (muon1.pt() > 15 &&
          fabs(muon1.eta()) < 2.4 &&
          (muon1.chargedHadronIso() + muon1.neutralHadronIso() + muon1.photonIso())/muon1.pt() - 0.5 * (*rho) < 0.3 &&
          muon1.muonBestTrack()->dz(PV.position()) < 0.2 &&
          abs(muon1.muonBestTrack()->dxy(PV.position())) < 0.045 &&
          muon1.isPFMuon() &&
          muon1.isGlobalMuon() &&
          muon1.isTrackerMuon() &&
          muon2.pt() > 15 &&
          fabs(muon2.eta()) < 2.4 &&
          (muon2.chargedHadronIso() + muon2.neutralHadronIso() + muon2.photonIso())/muon2.pt() - 0.5 * (*rho) < 0.3 &&
          muon2.muonBestTrack()->dz(PV.position()) < 0.2 &&
          abs(muon2.muonBestTrack()->dxy(PV.position())) < 0.045 &&
          muon2.isPFMuon() &&
          muon2.isGlobalMuon() &&
          muon2.isTrackerMuon() &&
          muon1.charge() * muon2.charge() < 0
          )
        diMuon = true;
    }
  }  

  // tau cands and btag veto
  bool atLeastOneBTag = false;
  vector<const pat::Jet *> tauCands;
  vector<int> tauCandsCharge;
  for (const pat::Jet &jet : *jets) {
    if (jet.pt() > 20 && fabs(jet.eta()) < 2.4 && jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.890)
      atLeastOneBTag = true;
    bool noNearbyGlobalMuon = true;
    for (const pat::Muon &muon : *muons) {
      if (muon.pt() < 5) continue; 
      if (!muon.isGlobalMuon()) continue;
      double deltaR = reco::deltaR(muon.eta(),muon.phi(),jet.eta(),jet.phi());
      if (deltaR < 0.4) noNearbyGlobalMuon = false;
    }
    if (jet.pt() > 20.0 &&
        fabs(jet.eta()) < 2.3 &&
        noNearbyGlobalMuon
        ) {
      double max_pt = 0;
      int leading_track_charge = 0;
      for (unsigned int i = 0; i < jet.numberOfDaughters(); i++) {
        if (jet.daughter(i)->pdgId() == 22 || jet.daughter(i)->pdgId() == 111) continue;
        if (jet.daughter(i)->pt() > max_pt) {
          max_pt = jet.daughter(i)->pt();
          leading_track_charge = jet.daughter(i)->charge();
        }
      }
      if (max_pt > 5.0)
        tauCands.push_back(&jet);
        tauCandsCharge.push_back(leading_track_charge);
    }
  }

  // choose single tau-muon system
  int muon_index = -1;
  int tau_index = -1;
  for (unsigned int i = 0; i < passedMuons.size(); i++) {
    for (unsigned int j = 0; j < tauCands.size(); j++) {
      const pat::Muon & muon = *passedMuons[i];
      const pat::Jet & tau = *tauCands[j];
      int tau_charge = tauCandsCharge[j];
      double dR = reco::deltaR(muon.eta(),muon.phi(),tau.eta(),tau.phi());
      int charge = muon.charge() * tau_charge;
      if (dR < 0.5) continue;
      if (charge > 0) continue;
      if (muon_index == -1) {
        muon_index = i;
        tau_index = j; }
      else if (passedMuons[i]->pt() + tauCands[j]->pt() > passedMuons[muon_index]->pt() + tauCands[tau_index]->pt()) {
        muon_index = i;
        tau_index = j; }
    }
  }
  
  // filter decision
  bool passTauMuonPair = false;
  bool passDR = false;
  bool passOppSign = false;
  bool passMT = false;
  bool passPzeta = false;
  bool passExtraLep = false;
  bool passBTag = false;
  
  // at least one muon and one tau cand
  if (muon_index != -1 && tau_index != -1) passTauMuonPair = true;

  // tau-muon system: dR, opp sign, MT, Pzeta
  if (muon_index != -1 && tau_index != -1) {
    passTauMuonPair = true;
    const pat::Muon & theMuon = *passedMuons[muon_index];
    const pat::Jet & theTau = *tauCands[tau_index];
    int tau_charge = tauCandsCharge[tau_index];

    double dPhi = theMuon.phi() - MET.phi();
    double MT = sqrt(2 * theMuon.pt() * MET.pt() * (1 - cos(dPhi)));    
    TVector2 pTmuon;
    pTmuon.SetMagPhi(theMuon.pt(), theMuon.phi());  
    TVector2 pTtau;
    pTtau.SetMagPhi(theTau.pt(), theTau.phi());  
    TVector2 pTmet;
    pTmet.SetMagPhi(MET.pt(), MET.phi());  
    TVector2 zeta;
    zeta.SetMagPhi(1.0, (theMuon.phi()-theTau.phi())/2.0 + theTau.phi());
    double PzetaAll = zeta * (pTmuon + pTtau + pTmet);
    double PzetaVis = zeta * (pTmuon + pTtau);
    double Pzeta = PzetaAll - 0.85 * PzetaVis;

    if (reco::deltaR(theMuon.eta(), theMuon.phi(), theTau.eta(), theTau.phi()) > 0.5) passDR = true;
    if (theMuon.charge() * tau_charge < 0) passOppSign = true;
    if (MT < 40) passMT = true;
    if (Pzeta > -25) passPzeta = true;
  }

  // extra lepton veto
  if (!extraMuon && !extraElectron && !diMuon) passExtraLep = true;

  // btag veto
  if (!atLeastOneBTag) passBTag = true;

  // final filter decision
  bool passAll = passTauMuonPair && passDR && passOppSign && passMT && passPzeta && passExtraLep && passBTag;

  return passAll;
}

void
ZtoTauHadRecoSelector::beginStream(edm::StreamID)
{
}

void
ZtoTauHadRecoSelector::endStream() {
}

void
ZtoTauHadRecoSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
 
void
ZtoTauHadRecoSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
 
void
ZtoTauHadRecoSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

void
ZtoTauHadRecoSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

void
ZtoTauHadRecoSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  descriptions.add("ZtoTauHadRecoSelectorFilter", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZtoTauHadRecoSelector);
