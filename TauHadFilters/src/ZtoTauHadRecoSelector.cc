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

#include "MyTauHadFilters/TauHadFilters/interface/ZtoTauHadPreSelection.h"

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
      bool cfg_dumpCutflow;
      bool cfg_tnpSelectionOnly;

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

      // cutflow counts
      int cutflow_total;
      int cutflow_foundTrigger;
      int cutflow_passTrigger;
      int cutflow_passMuon;
      int cutflow_passTau;
      int cutflow_passMuonAndTrigger;
      int cutflow_passMuonAndTau;
      int cutflow_passPreSelection;
      int cutflow_N1_DR;
      int cutflow_N1_MT;
      int cutflow_N1_Pzeta;
      int cutflow_N1_DiMuonVeto;
      int cutflow_N1_ExtraMuonVeto;
      int cutflow_N1_ExtraElectronVeto;
      int cutflow_N1_BtagVeto; 
};

//
// constructors
//
ZtoTauHadRecoSelector::ZtoTauHadRecoSelector(const edm::ParameterSet& iConfig) :
  cfg_dumpCutflow(iConfig.getUntrackedParameter<bool>("dumpCutflow")),
  cfg_tnpSelectionOnly(iConfig.getUntrackedParameter<bool>("tnpSelectionOnly"))
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

  cutflow_total = 0;
  cutflow_foundTrigger = 0;
  cutflow_passTrigger = 0;
  cutflow_passMuon = 0;
  cutflow_passTau = 0;
  cutflow_passMuonAndTrigger = 0;
  cutflow_passMuonAndTau = 0;
  cutflow_passPreSelection = 0;
  cutflow_N1_DR = 0;
  cutflow_N1_MT = 0;
  cutflow_N1_Pzeta = 0;
  cutflow_N1_DiMuonVeto = 0;
  cutflow_N1_ExtraMuonVeto = 0;
  cutflow_N1_ExtraElectronVeto = 0;
  cutflow_N1_BtagVeto = 0;
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

  // use preselection function
  TauHadFilters::PreSelectionResult result;
  result = TauHadFilters::computePreSelectionResult(iEvent, triggerBits, triggerObjects, triggerPrescales, vertices, taus, muons, electrons, jets, mets, rho);

  if (result.passPreSelection) cutflow_total += 1;
  if (result.foundMuonTrigger != "") cutflow_foundTrigger += 1;
  if (result.passMuonTrigger) cutflow_passTrigger += 1;
  if (result.foundTagMuon) cutflow_passMuon += 1;
  if (result.foundProbeTau) cutflow_passTau += 1;
  if (result.foundTagMuon && result.passMuonTrigger) cutflow_passMuonAndTrigger += 1;
  if (result.foundTagMuon && result.foundProbeTau) cutflow_passMuonAndTau += 1;
  if (result.passPreSelection) cutflow_passPreSelection += 1;
  if (result.passMuonTrigger &&
      result.foundTagMuon && result.foundProbeTau && 
                      result.passMT && result.passPzeta && 
      result.passDiMuonVeto && result.passExtraElectronVeto && result.passExtraMuonVeto && 
      result.passBtagVeto)
      cutflow_N1_DR += 1;
  if (result.passMuonTrigger &&
      result.foundTagMuon && result.foundProbeTau && 
      result.passDR &&                 result.passPzeta && 
      result.passDiMuonVeto && result.passExtraElectronVeto && result.passExtraMuonVeto && 
      result.passBtagVeto)
      cutflow_N1_MT += 1;
  if (result.passMuonTrigger &&
      result.foundTagMuon && result.foundProbeTau && 
      result.passDR && result.passMT &&                 
      result.passDiMuonVeto && result.passExtraElectronVeto && result.passExtraMuonVeto && 
      result.passBtagVeto)
      cutflow_N1_Pzeta += 1;
  if (result.passMuonTrigger &&
      result.foundTagMuon && result.foundProbeTau && 
      result.passDR && result.passMT && result.passPzeta && 
                               result.passExtraElectronVeto && result.passExtraMuonVeto && 
      result.passBtagVeto)
      cutflow_N1_DiMuonVeto += 1;
  if (result.passMuonTrigger &&
      result.foundTagMuon && result.foundProbeTau && 
      result.passDR && result.passMT && result.passPzeta && 
      result.passDiMuonVeto && result.passExtraElectronVeto &&
      result.passBtagVeto)
      cutflow_N1_ExtraMuonVeto += 1;
  if (result.passMuonTrigger &&
      result.foundTagMuon && result.foundProbeTau && 
      result.passDR && result.passMT && result.passPzeta && 
      result.passDiMuonVeto &&                                 result.passExtraMuonVeto && 
      result.passBtagVeto)
      cutflow_N1_ExtraElectronVeto += 1;
  if (result.passMuonTrigger &&
      result.foundTagMuon && result.foundProbeTau && 
      result.passDR && result.passMT && result.passPzeta && 
      result.passDiMuonVeto && result.passExtraElectronVeto && result.passExtraMuonVeto
      )
      cutflow_N1_BtagVeto += 1;

  if (!cfg_tnpSelectionOnly) return result.passPreSelection;
  else return result.passTagAndProbeSelection;
}

void
ZtoTauHadRecoSelector::beginJob()
{
}

void
ZtoTauHadRecoSelector::endJob()
{
  if (!cfg_tnpSelectionOnly)
    cout << "ran filter selecting muon and tau events only" << endl;
  else
    cout << "ran filter selecting full preselection" << endl;
  if (cfg_dumpCutflow) {
    // cutflow print
    cout << "total " << cutflow_total << endl;
    cout << "foundTrigger " << cutflow_foundTrigger << endl;
    cout << "passTrigger " << cutflow_passTrigger << endl;
    cout << "passMuon " << cutflow_passMuon << endl;
    cout << "passTau " << cutflow_passTau << endl;
    cout << "passMuonAndTrigger " << cutflow_passMuonAndTrigger << endl;
    cout << "passMuonAndTau " << cutflow_passMuonAndTau << endl;
    cout << "passPreSelection " << cutflow_passPreSelection << endl;
    cout << "N-1 DR " << cutflow_N1_DR << endl;
    cout << "N-1 MT " << cutflow_N1_MT << endl;
    cout << "N-1 Pzeta " << cutflow_N1_Pzeta << endl;
    cout << "N-1 DiMuonVeto " << cutflow_N1_DiMuonVeto << endl;
    cout << "N-1 ExtraMuonVeto " << cutflow_N1_ExtraMuonVeto << endl;
    cout << "N-1 ExtraElectronVeto " << cutflow_N1_ExtraElectronVeto << endl;
    cout << "N-1 BtagVeto " << cutflow_N1_BtagVeto << endl;
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
  desc.addUntracked<bool>("dumpCutflow", true);
  desc.addUntracked<bool>("tnpSelectionOnly", false);
  descriptions.add("ZtoTauHadRecoSelectorFilter", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZtoTauHadRecoSelector);
