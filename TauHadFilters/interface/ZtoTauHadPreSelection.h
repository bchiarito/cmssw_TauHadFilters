#ifndef ZTOTAUHADPRESLECTION_H
#define ZTOTAUHADPRESLECTION_H

// cpp includes
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
// ROOT includes
#include "TVector2.h"
#include "TLorentzVector.h"
// cmssw framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
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

using std::string;
using std::vector;

namespace TauHadFilters
{
  struct PreSelectionResult
  {
    // objects
    bool foundTagMuon;
    int nTagMuons;
    const pat::Muon * tagMuon;
    double tagMuon_iso;
    bool foundProbeTau;
    int nProbeTaus;
    const pat::Jet * probeTau;
    // trigger
    string foundMuonTrigger;
    bool passMuonTrigger;
    // cut variables on pair
    double DR;
    double MT;
    double Pzeta;
    // descision on cut variables
    bool passDR;
    bool passMT;
    bool passPzeta;
    // event wide vetos
    bool passExtraMuonVeto;
    bool passDiMuonVeto;
    bool passExtraElectronVeto;
    bool passBtagVeto;
    double btagDiscriminant;
    // summary bools
    bool passTagAndProbeSelection;
    bool passPreSelection;
  };

  typedef struct PreSelectionResult PreSelectionResult;

  struct PreSelectionResult computePreSelectionResult(const edm::Event& iEvent, edm::Handle<edm::TriggerResults>& triggerBits, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, edm::Handle<pat::PackedTriggerPrescales> triggerPrescales, edm::Handle<reco::VertexCollection> vertices, edm::Handle<pat::TauCollection> taus, edm::Handle<pat::MuonCollection> muons, edm::Handle<pat::ElectronCollection> electrons, edm::Handle<pat::JetCollection> jets, edm::Handle<pat::METCollection> mets, edm::Handle<double> rho)
  {
    struct PreSelectionResult result;

    const reco::Vertex &PV = vertices->front();
    pat::MET MET = (*mets)[0];

    // trigger 
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    string trigger_muon = "HLT_IsoMu24";
    bool bit_muon = false;
    string name_muon_trigger = "";
    for (unsigned int i = 0, n = triggerBits->size(); i < n; i++)
    {
       string triggerName = names.triggerName(i);

       std::size_t pos = triggerName.find(trigger_muon);
       if ( pos != std::string::npos ) {
         bit_muon = triggerBits->accept(i);
         name_muon_trigger = triggerName;
       }
    }
    bool passMuonTrigger = bit_muon;

    // muons
    vector<const pat::Muon *> passedMuons;
    for (const pat::Muon &muon : *muons) {
      if (muon.pt() > 25.0 &&
          fabs(muon.eta()) < 2.1 &&
          (muon.chargedHadronIso() + muon.neutralHadronIso() + muon.photonIso())/muon.pt() - 0.5 * (*rho) < 0.1 &&
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
          (muon.chargedHadronIso() + muon.neutralHadronIso() + muon.photonIso())/muon.pt() - 0.5 * (*rho) < 0.3 &&
          muon.muonBestTrack()->dz(PV.position()) < 0.2 &&
          abs(muon.muonBestTrack()->dxy(PV.position())) < 0.045 &&
          muon.isMediumMuon() ) {
        if (muon.pt() > 19.0 &&
            fabs(muon.eta()) < 2.1 &&
            (muon.chargedHadronIso() + muon.neutralHadronIso() + muon.photonIso())/muon.pt() - 0.5 * (*rho) < 0.1 &&
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
        double DR = deltaR(muon1.eta(), muon1.phi(), muon2.eta(), muon2.phi());
        if (DR <= 0.15) continue;
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
            muon1.charge() * muon2.charge() < 0)
          diMuon = true;
      }
    }  

    // tau cands from jets and btag veto
    bool atLeastOneBTag = false;
    double highestBTagDiscriminant = -1.0;
    vector<const pat::Jet *> tauJetCands;
    vector<int> tauJetCandsCharge;
    for (const pat::Jet &jet : *jets) {
      if (jet.pt() > 20 && fabs(jet.eta()) < 2.4) {
        double btagdiscriminant = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        if (btagdiscriminant > highestBTagDiscriminant) highestBTagDiscriminant = btagdiscriminant;
      }
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

    // choose single taujet-muon system
    int muon_jet_index = -1;
    int tau_jet_index = -1;
    for (unsigned int i = 0; i < passedMuons.size(); i++) {
      for (unsigned int j = 0; j < tauJetCands.size(); j++) {
        const pat::Muon & muon = *passedMuons[i];
        const pat::Jet & tau = *tauJetCands[j];
        int tau_charge = tauJetCandsCharge[j];
        double DR = reco::deltaR(muon.eta(),muon.phi(),tau.eta(),tau.phi());
        int charge = muon.charge() * tau_charge;
        if (DR < 0.5) continue;
        if (charge > 0) continue;
        if (muon_jet_index == -1) {
          muon_jet_index = i;
          tau_jet_index = j; }
        else if (passedMuons[i]->pt() + tauJetCands[j]->pt() > passedMuons[muon_jet_index]->pt() + tauJetCands[tau_jet_index]->pt()) {
          muon_jet_index = i;
          tau_jet_index = j; }
      }
    }

    // filter decision
    bool passTauMuonPair = false;
    bool passDR = false;
    bool passMT = false;
    bool passPzeta = false;
    bool passExtraLep = false;
    bool passBTag = false;

    // at least one muon and one tau cand
    if (muon_jet_index != -1 && tau_jet_index != -1) passTauMuonPair = true;

    double MT = -100;
    double Pzeta = -100;
    double cut_DR = -1;
    // tau-muon system: DR, opp sign, MT, Pzeta
    if (passTauMuonPair) {
      const pat::Muon & theMuon = *passedMuons[muon_jet_index];
      const pat::Jet & theTau = *tauJetCands[tau_jet_index];

      double dPhi = theMuon.phi() - MET.phi();
      MT = sqrt(2 * theMuon.pt() * MET.pt() * (1 - cos(dPhi)));    
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
      Pzeta = PzetaAll - 0.85 * PzetaVis;
      cut_DR = reco::deltaR(theMuon.eta(), theMuon.phi(), theTau.eta(), theTau.phi());

      if (cut_DR > 0.5) passDR = true;
      if (MT < 40) passMT = true;
      if (Pzeta > -25) passPzeta = true;
    }

    // extra lepton veto
    if (!extraMuon && !extraElectron && !diMuon) passExtraLep = true;

    // btag veto
    if (!atLeastOneBTag) passBTag = true;

    // final filter decision
    bool passAll = passMuonTrigger && passTauMuonPair && passDR && passMT && passPzeta && passExtraLep && passBTag;

    // fill return struct
    result.nTagMuons = passedMuons.size();
    result.nProbeTaus = tauJetCands.size();
    result.foundTagMuon = (muon_jet_index != -1);
    result.foundProbeTau = (tau_jet_index != -1);
    result.tagMuon = NULL;
    result.probeTau = NULL;
    if (result.foundTagMuon) result.tagMuon = &(*passedMuons[muon_jet_index]);
    if (result.foundTagMuon) result.tagMuon_iso = (result.tagMuon->chargedHadronIso() + result.tagMuon->neutralHadronIso() + result.tagMuon->photonIso())/result.tagMuon->pt() - 0.5 * (*rho) < 0.1;
    if (result.foundProbeTau) result.probeTau = &(*tauJetCands[tau_jet_index]);
    // trigger
    result.foundMuonTrigger = name_muon_trigger;
    result.passMuonTrigger = passMuonTrigger;
    // cut variables on pair
    result.DR = 0;
    result.MT = 0;
    result.Pzeta = 0;
    result.passDR = false;
    result.passMT = false;
    result.passPzeta = false;
    if (result.foundTagMuon && result.foundProbeTau) {
      result.DR = cut_DR;
      result.MT = MT;
      result.Pzeta = Pzeta;
      result.passDR = passDR;
      result.passMT = passMT;
      result.passPzeta = passPzeta;
    }
    // event wide vetos
    result.passExtraMuonVeto = !extraMuon;
    result.passDiMuonVeto = !diMuon;
    result.passExtraElectronVeto = !extraElectron;
    result.passBtagVeto = !atLeastOneBTag;
    result.btagDiscriminant = highestBTagDiscriminant;
    // summary bools
    result.passTagAndProbeSelection = passTauMuonPair;
    result.passPreSelection = passAll;

    return result;
  }
}

#endif
