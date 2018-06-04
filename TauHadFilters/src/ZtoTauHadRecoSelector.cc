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
#include "TLorentzVector.h"
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
// trigger inlcudes
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
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
using std::cout;
using std::endl;

//
// class declaration
//
class ZtoTauHadRecoSelector : public edm::EDFilter {
   public:
      explicit ZtoTauHadRecoSelector(const edm::ParameterSet&);
      ~ZtoTauHadRecoSelector();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;

      virtual void beginJob() override;
      virtual void endJob() override;
      virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // Configuration Parameters
      bool cfg_tauObjs;
      bool cfg_dumpCutflow;
      double cfg_isoCut;

      // EDM Collection lables
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
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

      // cutflow print
      int count_events;
      int count_foundMuonTrigger;
      int count_passMuonTrigger;
      int count_passMuon;
      int count_passMuon_barrel;
      int count_passMuon_endcap;
      int count_passMuon_and_trigger;
      int count_passMuon_and_trigger_barrel;
      int count_passMuon_and_trigger_endcap;
      int count_passTauMuonPair;
      int count_passDR_n1;
      int count_passMT_n1;
      int count_passPzeta_n1;
      int count_passExtraLep_n1;
      int count_BTag_n1;
      int count_passMuonTrigger_n1;
      int count_passAll;
};

//
// constructors
//
ZtoTauHadRecoSelector::ZtoTauHadRecoSelector(const edm::ParameterSet& iConfig) :
  cfg_tauObjs(iConfig.getUntrackedParameter<bool>("tauObjs")),
  cfg_dumpCutflow(iConfig.getUntrackedParameter<bool>("dumpCutflow")),
  cfg_isoCut(iConfig.getUntrackedParameter<double>("isoCut"))
{
  triggerBits_ = consumes<edm::TriggerResults>( edm::InputTag("TriggerResults","","HLT") );
  triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>( edm::InputTag("selectedPatTrigger") );
  triggerPrescales_ = consumes<pat::PackedTriggerPrescales>( edm::InputTag("patTrigger") );
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

  // cutflow print
  count_events = 0;
  count_foundMuonTrigger = 0;
  count_passMuonTrigger = 0;
  count_passMuon = 0;
  count_passMuon_barrel = 0;
  count_passMuon_endcap = 0;
  count_passMuon_and_trigger = 0;
  count_passMuon_and_trigger_barrel = 0;
  count_passMuon_and_trigger_endcap = 0;
  count_passTauMuonPair = 0;
  count_passDR_n1 = 0;
  count_passMT_n1 = 0;
  count_passPzeta_n1 = 0;
  count_passExtraLep_n1 = 0;
  count_BTag_n1 = 0;
  count_passMuonTrigger_n1 = 0;
  count_passAll = 0;
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
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  const reco::Vertex &PV = vertices->front();

  edm::Handle<pat::TauCollection> taus;
  iEvent.getByToken(tauToken_, taus);

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

  // trigger 
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  string trigger_muon = "HLT_IsoMu20";
  bool found_muon_trigger = false;
  bool bit_muon = false;
  string name_muon = "";
  for (unsigned int i = 0, n = triggerBits->size(); i < n; i++)
  {
     string triggerName = names.triggerName(i);

     std::size_t pos = triggerName.find(trigger_muon);
     if ( pos != std::string::npos ) {
       found_muon_trigger = true;
       bit_muon = triggerBits->accept(i);
       name_muon = triggerName;
     }
  }
  bool passMuonTrigger = bit_muon;

  // muons
  vector<const pat::Muon *> passedMuons;
  for (const pat::Muon &muon : *muons) {
    if (muon.pt() > 22.0 &&
        fabs(muon.eta()) < 2.1 &&
        //(muon.chargedHadronIso() + muon.neutralHadronIso() + muon.photonIso())/muon.pt() - 0.5 * (*rho) < 0.1 &&
        (muon.chargedHadronIso() + muon.neutralHadronIso() + muon.photonIso())/muon.pt() < cfg_isoCut &&
        muon.muonBestTrack()->dz(PV.position()) < 0.2 &&
        abs(muon.muonBestTrack()->dxy(PV.position())) < 0.045 &&
        muon.isMediumMuon() )
      passedMuons.push_back(&muon);
  }

  // extra lepton veto
  bool skipped_one_passing_muon = false;
  bool extraMuon = false;
  for (const pat::Muon &muon : *muons) {
    if (muon.pt() > 10.0 &&
        fabs(muon.eta()) < 2.4 &&
        //(muon.chargedHadronIso() + muon.neutralHadronIso() + muon.photonIso())/muon.pt() - 0.5 * (*rho) < 0.3 &&
        (muon.chargedHadronIso() + muon.neutralHadronIso() + muon.photonIso())/muon.pt() < 0.3 &&
        muon.muonBestTrack()->dz(PV.position()) < 0.2 &&
        abs(muon.muonBestTrack()->dxy(PV.position())) < 0.045 &&
        muon.isMediumMuon() ) {
      if (muon.pt() > 19.0 &&
          fabs(muon.eta()) < 2.1 &&
          //(muon.chargedHadronIso() + muon.neutralHadronIso() + muon.photonIso())/muon.pt() - 0.5 * (*rho) < 0.1 &&
          (muon.chargedHadronIso() + muon.neutralHadronIso() + muon.photonIso())/muon.pt() < 0.1 &&
          muon.muonBestTrack()->dz(PV.position()) < 0.2 &&
          abs(muon.muonBestTrack()->dxy(PV.position())) < 0.045 &&
          muon.isMediumMuon() ) {
        if (skipped_one_passing_muon) extraMuon = true;
        if (!skipped_one_passing_muon) skipped_one_passing_muon = true;
      } else {
        extraMuon = true;
      }
    }
  }
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
      if (dR <= 0.15) continue;
      if (muon1.pt() > 15 &&
          fabs(muon1.eta()) < 2.4 &&
          //(muon1.chargedHadronIso() + muon1.neutralHadronIso() + muon1.photonIso())/muon1.pt() - 0.5 * (*rho) < 0.3 &&
          (muon1.chargedHadronIso() + muon1.neutralHadronIso() + muon1.photonIso())/muon1.pt() < 0.3 &&
          muon1.muonBestTrack()->dz(PV.position()) < 0.2 &&
          abs(muon1.muonBestTrack()->dxy(PV.position())) < 0.045 &&
          muon1.isPFMuon() &&
          muon1.isGlobalMuon() &&
          muon1.isTrackerMuon() &&
          muon2.pt() > 15 &&
          fabs(muon2.eta()) < 2.4 &&
          //(muon2.chargedHadronIso() + muon2.neutralHadronIso() + muon2.photonIso())/muon2.pt() - 0.5 * (*rho) < 0.3 &&
          (muon2.chargedHadronIso() + muon2.neutralHadronIso() + muon2.photonIso())/muon2.pt() < 0.3 &&
          muon2.muonBestTrack()->dz(PV.position()) < 0.2 &&
          abs(muon2.muonBestTrack()->dxy(PV.position())) < 0.045 &&
          muon2.isPFMuon() &&
          muon2.isGlobalMuon() &&
          muon2.isTrackerMuon() &&
          muon1.charge() * muon2.charge() < 0)
        diMuon = true;
    }
  }  

  // tau cands from jets and btag veto
  bool atLeastOneBTag = false;
  vector<const pat::Jet *> tauJetCands;
  vector<int> tauJetCandsCharge;
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
        noNearbyGlobalMuon) {
      double leading_track_pt = 0;
      int leading_track_charge = 0;
      for (unsigned int i = 0; i < jet.numberOfDaughters(); i++) {
        if (jet.daughter(i)->pdgId() == 22 || jet.daughter(i)->pdgId() == 111) continue;
        if (jet.daughter(i)->pt() > leading_track_pt) {
          leading_track_pt = jet.daughter(i)->pt();
          leading_track_charge = jet.daughter(i)->charge();
        }
      }
      if (leading_track_pt > 5.0)
        tauJetCands.push_back(&jet);
        tauJetCandsCharge.push_back(leading_track_charge);
    }
  }

  // tau cands from tau collection
  vector<const pat::Tau *> tauCands;
  for (const pat::Tau &tau : *taus) {
    bool noNearbyGlobalMuon = true;
    for (const pat::Muon &muon : *muons) {
      if (muon.pt() < 5) continue; 
      if (!muon.isGlobalMuon()) continue;
      double deltaR = reco::deltaR(muon.eta(),muon.phi(),tau.eta(),tau.phi());
      if (deltaR < 0.4) noNearbyGlobalMuon = false;
    }
    double leading_track_pt = (tau.leadChargedHadrCand())->pt();

    if (tau.pt() > 20.0 &&
        fabs(tau.eta()) < 2.3 &&
        noNearbyGlobalMuon &&
        tau.tauID("againstMuonTight3") == 1 &&
        tau.tauID("againstElectronVLooseMVA6") == 1 &&
        leading_track_pt > 5.0)
          tauCands.push_back(&tau);
  }

  // choose single taujet-muon system
  int muon_jet_index = -1;
  int tau_jet_index = -1;
  for (unsigned int i = 0; i < passedMuons.size(); i++) {
    for (unsigned int j = 0; j < tauJetCands.size(); j++) {
      const pat::Muon & muon = *passedMuons[i];
      const pat::Jet & tau = *tauJetCands[j];
      int tau_charge = tauJetCandsCharge[j];
      double dR = reco::deltaR(muon.eta(),muon.phi(),tau.eta(),tau.phi());
      int charge = muon.charge() * tau_charge;
      if (dR < 0.5) continue;
      if (charge > 0) continue;
      if (muon_jet_index == -1) {
        muon_jet_index = i;
        tau_jet_index = j; }
      else if (passedMuons[i]->pt() + tauJetCands[j]->pt() > passedMuons[muon_jet_index]->pt() + tauJetCands[tau_jet_index]->pt()) {
        muon_jet_index = i;
        tau_jet_index = j; }
    }
  }

  // choose single tau-muon system
  int muon_index = -1;
  int tau_index = -1;
  for (unsigned int i = 0; i < passedMuons.size(); i++) {
    for (unsigned int j = 0; j < tauCands.size(); j++) {
      const pat::Muon & muon = *passedMuons[i];
      const pat::Tau & tau = *tauCands[j];
      int tau_charge = tau.charge();
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
  bool passMT = false;
  bool passPzeta = false;
  bool passExtraLep = false;
  bool passBTag = false;

  if (cfg_tauObjs) {  
    // at least one muon and one tau cand
    if (muon_index != -1 && tau_index != -1) passTauMuonPair = true;

    // tau-muon system: dR, opp sign, MT, Pzeta
    if (passTauMuonPair) {
      const pat::Muon & theMuon = *passedMuons[muon_index];
      const pat::Tau & theTau = *tauCands[tau_index];

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
      if (MT < 40) passMT = true;
      if (Pzeta > -25) passPzeta = true;
    }
  } else {
    // at least one muon and one tau cand
    if (muon_jet_index != -1 && tau_jet_index != -1) passTauMuonPair = true;

    // tau-muon system: dR, opp sign, MT, Pzeta
    if (passTauMuonPair) {
      const pat::Muon & theMuon = *passedMuons[muon_jet_index];
      const pat::Jet & theTau = *tauJetCands[tau_jet_index];

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
      if (MT < 40) passMT = true;
      if (Pzeta > -25) passPzeta = true;
    }
  }

  // extra lepton veto
  if (!extraMuon && !extraElectron && !diMuon) passExtraLep = true;

  // btag veto
  if (!atLeastOneBTag) passBTag = true;

  // barrel vs endcap muon if pass
  bool muonIsBarrel = false;
  bool muonIsEndcap = false;
  if (passedMuons.size() > 0 ) {

    const pat::Muon * poorestIsoMuon = passedMuons[0];
    double lowestIso = 2.0;
    for (const pat::Muon &muon : *muons) {
      double iso = (muon.chargedHadronIso() + muon.neutralHadronIso() + muon.photonIso())/muon.pt();
      if (iso < lowestIso) lowestIso = iso; poorestIsoMuon = &muon;
    }

    double eta = poorestIsoMuon->eta();
    if (eta < 1.479) muonIsBarrel = true; 
    else muonIsEndcap = true;
  }

  // final filter decision
  bool passAll = passMuonTrigger && passTauMuonPair && passDR && passMT && passPzeta && passExtraLep && passBTag;

  // cutflow print
  count_events++;
  if(found_muon_trigger) count_foundMuonTrigger++;
  if(passMuonTrigger) count_passMuonTrigger++;
  if(passedMuons.size() > 0) count_passMuon++;
  if(passedMuons.size() > 0 && muonIsBarrel) count_passMuon_barrel++;
  if(passedMuons.size() > 0 && muonIsEndcap) count_passMuon_endcap++;
  if(passedMuons.size() > 0 && passMuonTrigger) count_passMuon_and_trigger++;
  if(passedMuons.size() > 0 && passMuonTrigger && muonIsBarrel) count_passMuon_and_trigger_barrel++;
  if(passedMuons.size() > 0 && passMuonTrigger && muonIsEndcap) count_passMuon_and_trigger_endcap++;
  if(passTauMuonPair) count_passTauMuonPair++;
  if(passAll) count_passAll++;
  if(                   passTauMuonPair && passDR && passMT && passPzeta && passExtraLep && passBTag) count_passMuonTrigger_n1++;
  if(passMuonTrigger && passTauMuonPair &&           passMT && passPzeta && passExtraLep && passBTag) count_passDR_n1++;
  if(passMuonTrigger && passTauMuonPair && passDR &&           passPzeta && passExtraLep && passBTag) count_passMT_n1++;
  if(passMuonTrigger && passTauMuonPair && passDR && passMT &&              passExtraLep && passBTag) count_passPzeta_n1++;
  if(passMuonTrigger && passTauMuonPair && passDR && passMT && passPzeta &&                 passBTag) count_passExtraLep_n1++;
  if(passMuonTrigger && passTauMuonPair && passDR && passMT && passPzeta && passExtraLep            ) count_BTag_n1++;

  return passAll;
}

void
ZtoTauHadRecoSelector::beginJob()
{
}

void
ZtoTauHadRecoSelector::endJob()
{
  if (cfg_dumpCutflow) {
    // cutflow print
    std::cout << "" << std::endl;
    std::cout << "CUTFLOW" << std::endl;
    std::cout << "count_events " << count_events << std::endl;
    std::cout << "count_foundMuonTrigger " << count_foundMuonTrigger << std::endl;
    std::cout << "count_passMuonTrigger " << count_passMuonTrigger << std::endl;
    std::cout << "count_passMuon " << count_passMuon << std::endl;
    std::cout << "count_passMuon_barrel " << count_passMuon_barrel << std::endl;
    std::cout << "count_passMuon_endcap " << count_passMuon_endcap << std::endl;
    std::cout << "count_passMuon_and_trigger " << count_passMuon_and_trigger << std::endl;
    std::cout << "count_passMuon_and_trigger_barrel " << count_passMuon_and_trigger_barrel << std::endl;
    std::cout << "count_passMuon_and_trigger_endcap " << count_passMuon_and_trigger_endcap << std::endl;
    std::cout << "count_passTauMuonPair " << count_passTauMuonPair << std::endl;
    std::cout << "count_passDR_n1 " << count_passDR_n1 << std::endl;
    std::cout << "count_passMT_n1 " << count_passMT_n1 << std::endl;
    std::cout << "count_passPzeta_n1 " << count_passPzeta_n1 << std::endl;
    std::cout << "count_passExtraLep_n1 " << count_passExtraLep_n1 << std::endl;
    std::cout << "count_BTag_n1 " << count_BTag_n1 << std::endl;
    std::cout << "count_passMuonTrigger_n1 " << count_passMuonTrigger_n1 << std::endl;
    std::cout << "count_passAll " << count_passAll << std::endl;
  }
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
  desc.addUntracked<bool>("tauObjs", false);
  desc.addUntracked<bool>("dumpCutflow", true);
  desc.addUntracked<double>("isoCut", 0.1);
  descriptions.add("ZtoTauHadRecoSelectorFilter", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZtoTauHadRecoSelector);
