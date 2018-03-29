// -*- C++ -*-
//
// Package:    TauHadFilters/ZtoTauHadTreeMaker
// Class:      ZtoTauHadTreeMaker
// 
/**\class ZtoTauHadTreeMaker ZtoTauHadTreeMaker.cc TauHadFilters/ZtoTauHadTreeMaker/plugins/ZtoTauHadTreeMaker.cc

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
class ZtoTauHadTreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources> {
   public:
      explicit ZtoTauHadTreeMaker(const edm::ParameterSet&);
      ~ZtoTauHadTreeMaker();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      vector<string> getDecay(const reco::Candidate&, int flag=0);
      bool isAncestorOfZ(const reco::Candidate *);
      bool notTerminalTau(const reco::Candidate *);
      const reco::Candidate * tauDaughter(const reco::Candidate *);

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
      edm::EDGetTokenT<vector<reco::GenParticle>> genToken_;
      edm::EDGetTokenT<double> rhoToken_;

      // tree and branches
      TTree* tree_;
      // event wide
      int bevent;
      int brun;
      int blumi;
      double bdecayType;
      // cutflow
      int bpassMuon;
      int bpassTauCand;
      int bpassTauAndMuon;
      int bpassTauMuonPair;
      int bpassDR;
      int bpassOppSign;
      int bpassMT;
      int bpassPzeta;
      int bpassExtraLep;
      int bpassDiMuon;
      int bpassExtraElec;
      int bpassExtraMuon;
      int bpassBTag;
      int bpassAll;
      // gencutflow
      int bgenpassTauCand;
      int bgenpassMuon;
      int bgenpassTauAndMuon;
      int bgenpassTauMuonPair;
      // genparticles
      int bnGenZ;
      vector<Double_t> bGenZ_pt;
      vector<Double_t> bGenZ_phi;
      vector<Double_t> bGenZ_eta;
      vector<Double_t> bGenZ_m;
      vector<Double_t> bGenZ_status;
      int bGenTau_nCon;
      vector<Double_t> bGenTau_con_pt;
      vector<Double_t> bGenTau_con_dr;
      double bGenTau_con_maxdr;
      int bnGenMuon;
      vector<Double_t> bGenMuon_pt;
      vector<Double_t> bGenMuon_phi;
      vector<Double_t> bGenMuon_eta;
      vector<Double_t> bGenMuon_m;
      int bnGenTau;
      vector<Double_t> bGenTau_pt;
      vector<Double_t> bGenTau_phi;
      vector<Double_t> bGenTau_eta;
      vector<Double_t> bGenTau_m;
      vector<Double_t> bGenTau_status;
      int bnGenHadTau;
      vector<Double_t> bGenHadTau_pt;
      vector<Double_t> bGenHadTau_phi;
      vector<Double_t> bGenHadTau_eta;
      vector<Double_t> bGenHadTau_m;
      vector<Double_t> bGenHadTau_status;
      // muons
      int bnMuon;
      vector<Double_t> bMuon_genTau_dR;
      // tau cands
      int bnTauCand;
      vector<Double_t> bTauCand_genTau_dR;
};

//
// constructors
//
ZtoTauHadTreeMaker::ZtoTauHadTreeMaker(const edm::ParameterSet& iConfig)
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
  genToken_ = consumes<vector<reco::GenParticle>>(edm::InputTag("prunedGenParticles"));
  rhoToken_ = consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralChargedPileUp"));

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree","ZtoTauHadTreeMakerTree");
  // event wide
  tree_->Branch("eventNum",&bevent,"eventNum/I");
  tree_->Branch("runNum",&brun,"runNum/I");
  tree_->Branch("lumiNum",&blumi,"lumiNum/I");
  tree_->Branch("decayType",&bdecayType,"decayType/D");
  // cutflow
  tree_->Branch("passMuon",&bpassMuon,"passMuon/I");
  tree_->Branch("passTauCand",&bpassTauCand,"passTauCand/I");
  tree_->Branch("passTauAndMuon",&bpassTauAndMuon,"passTauAndMuon/I");
  tree_->Branch("passTauMuonPair",&bpassTauMuonPair,"passTauMuonPair/I");
  tree_->Branch("passDR",&bpassDR,"passDR/I");
  tree_->Branch("passOppSign",&bpassOppSign,"passOppSign/I");
  tree_->Branch("passMT",&bpassMT,"passMT/I");
  tree_->Branch("passPzeta",&bpassPzeta,"passPzeta/I");
  tree_->Branch("passExtraLep",&bpassExtraLep,"passExtraLep/I");
  tree_->Branch("passDiMuon",&bpassDiMuon,"passDiMuon/I");
  tree_->Branch("passExtraElec",&bpassExtraElec,"passExtraElec/I");
  tree_->Branch("passExtraMuon",&bpassExtraMuon,"passExtraMuon/I");
  tree_->Branch("passBTag",&bpassBTag,"passBTag/I");
  tree_->Branch("passAll",&bpassAll,"passAll/I");
  // gen cutflow
  tree_->Branch("genpassMuon",&bgenpassMuon,"genpassMuon/I");
  tree_->Branch("genpassTauCand",&bgenpassTauCand,"genpassTauCand/I");
  tree_->Branch("genpassTauAndMuon",&bgenpassTauAndMuon,"genpassTauAndMuon/I");
  tree_->Branch("genpassTauMuonPair",&bgenpassTauMuonPair,"genpassTauMuonPair/I");
  // genparticles
  tree_->Branch("nGenZ",&bnGenZ,"nGenZ/I"); 
  tree_->Branch("GenZ_pt",&bGenZ_pt); 
  tree_->Branch("GenZ_m",&bGenZ_m);
  tree_->Branch("GenZ_phi",&bGenZ_phi);
  tree_->Branch("GenZ_eta",&bGenZ_eta);
  tree_->Branch("GenZ_status",&bGenZ_status);
  tree_->Branch("nGenMuon",&bnGenMuon,"nGenMuon/I"); 
  tree_->Branch("GenMuon_pt",&bGenMuon_pt); 
  tree_->Branch("GenMuon_m",&bGenMuon_m);
  tree_->Branch("GenMuon_phi",&bGenMuon_phi);
  tree_->Branch("GenMuon_eta",&bGenMuon_eta);
  tree_->Branch("nGenTau",&bnGenTau,"nGenTau/I"); 
  tree_->Branch("GenTau_pt",&bGenTau_pt); 
  tree_->Branch("GenTau_m",&bGenTau_m);
  tree_->Branch("GenTau_phi",&bGenTau_phi);
  tree_->Branch("GenTau_eta",&bGenTau_eta);
  tree_->Branch("GenTau_status",&bGenTau_status);
  tree_->Branch("GenTau_nCon",&bGenTau_nCon,"GenTau_nCon/I");
  tree_->Branch("GenTau_con_pt",&bGenTau_con_pt);
  tree_->Branch("GenTau_con_dr",&bGenTau_con_dr);
  tree_->Branch("GenTau_con_maxdr",&bGenTau_con_maxdr,"GenTau_con_maxdr/D");
  tree_->Branch("nGenHadTau",&bnGenHadTau,"nGenHadTau/I"); 
  tree_->Branch("GenHadTau_pt",&bGenHadTau_pt); 
  tree_->Branch("GenHadTau_m",&bGenHadTau_m);
  tree_->Branch("GenHadTau_phi",&bGenHadTau_phi);
  tree_->Branch("GenHadTau_eta",&bGenHadTau_eta);
  tree_->Branch("GenHadTau_status",&bGenHadTau_status);
  // muons
  tree_->Branch("nMuon",&bnMuon,"nMuon/I");
  tree_->Branch("Muon_genTau_dR",&bMuon_genTau_dR);
  // tau cands
  tree_->Branch("nTauCand",&bnTauCand,"nTauCand/I");
  tree_->Branch("TauCand_genTau_dR",&bTauCand_genTau_dR);
}

//
// destructor
//
ZtoTauHadTreeMaker::~ZtoTauHadTreeMaker()
{
}

// ------------ method called on each new Event  ------------
void
ZtoTauHadTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  edm::Handle<vector<reco::GenParticle>> genparticles;
  iEvent.getByToken(genToken_, genparticles);

  edm::Handle<double> rho;
  iEvent.getByToken(rhoToken_, rho);

  // initialize branches
  bevent = iEvent.id().event();
  brun = iEvent.id().run();
  blumi = iEvent.id().luminosityBlock();
  bdecayType = -1;
  bnGenHadTau = 0;

  bpassMuon = false;
  bpassTauCand = false;
  bpassTauAndMuon = false;
  bpassTauMuonPair = false;
  bpassDR = false;
  bpassOppSign = false;
  bpassMT = false;
  bpassPzeta = false;
  bpassExtraLep = false;
  bpassBTag = false;
  bpassAll = false;

  bgenpassTauCand = false;
  bgenpassTauMuonPair = false;
  bgenpassMuon = false;
  bgenpassTauAndMuon = false;

  bGenZ_pt.clear();
  bGenZ_m.clear();
  bGenZ_phi.clear();
  bGenZ_eta.clear();
  bGenZ_status.clear();
  bGenTau_con_pt.clear();
  bGenTau_con_dr.clear();
  bGenMuon_pt.clear();
  bGenMuon_m.clear();
  bGenMuon_phi.clear();
  bGenMuon_eta.clear();
  bGenTau_pt.clear();
  bGenTau_m.clear();
  bGenTau_phi.clear();
  bGenTau_eta.clear();
  bGenTau_status.clear();
  bGenHadTau_pt.clear();
  bGenHadTau_m.clear();
  bGenHadTau_phi.clear();
  bGenHadTau_eta.clear();
  bGenHadTau_status.clear();

  bMuon_genTau_dR.clear();
  bTauCand_genTau_dR.clear();

  // decay type
  vector<string> leptons;
  for (unsigned int i = 0; i < genparticles->size(); i++) {
    const reco::GenParticle & genparticle = (*genparticles)[i];
    if (genparticle.status() != 21 && genparticle.status() != 22) continue;
    leptons = getDecay(genparticle);
  }
  if (leptons.size() == 0) 
    bdecayType = 0;
  if (leptons.size() == 1) 
    bdecayType = -1;
  if (leptons.size() == 2) 
  {
    if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end()) bdecayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end()) bdecayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end()) bdecayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end()) bdecayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end()) bdecayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end()) bdecayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end()) bdecayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end()) bdecayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end()) bdecayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end()) bdecayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end()) bdecayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end()) bdecayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end()) bdecayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end()) bdecayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end()) bdecayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end()) bdecayType = 3;

    else if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end()) bdecayType = 4;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end()) bdecayType = 4;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end()) bdecayType = 4;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end()) bdecayType = 4;
    else if (std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end()) bdecayType = 4;
    else if (std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end()) bdecayType = 4;
    else if (std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end()) bdecayType = 4;
    else if (std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end()) bdecayType = 4;

    else if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) bdecayType = 5.2;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) bdecayType = 5.1;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) bdecayType = 5.4;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) bdecayType = 5.3;
    else if (std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end()) bdecayType = 5.2;
    else if (std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end()) bdecayType = 5.1;
    else if (std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end()) bdecayType = 5.4;
    else if (std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end()) bdecayType = 5.3;

    else if (std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end()) bdecayType = 6;

    else if (std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) bdecayType = 7;

    else if (std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) bdecayType = 8;
    else if (std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end()) bdecayType = 8;

    else if (std::find(leptons.begin(), leptons.end(), "e+") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "e-") != leptons.end()) bdecayType = 1;

    else if (std::find(leptons.begin(), leptons.end(), "mu+") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "mu-") != leptons.end()) bdecayType = 2;

    else
      bdecayType = 9;
  }
  if (leptons.size() > 2)
  {
    bdecayType = 10;
  }

  // apply kinemtics cuts to hadronic tau
  bool found_one_hadtau = false;
  std::vector<const reco::GenParticle *> hadronicTaus;
  for (unsigned int i = 0; i < genparticles->size(); i++) {
    const reco::GenParticle &genparticle = (*genparticles)[i];  
    if (abs(genparticle.pdgId()) != 15) continue;
    if (!isAncestorOfZ(&genparticle)) continue;
    leptons = getDecay(genparticle);
    if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() ||
        std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() ||
        std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() ||
        std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() ||
        std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end() ||
        std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end() ||
        std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end() ||
        std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end()
      ) {
      hadronicTaus.push_back(&genparticle);
      if (found_one_hadtau) continue;

      // hadronic tau decay daughters
      const reco::Candidate * terminaltau = & genparticle;
      while (notTerminalTau(terminaltau)) {
        terminaltau = tauDaughter(terminaltau);
      }
      bGenTau_nCon = terminaltau->numberOfDaughters();
      for (unsigned int i = 0; i < terminaltau->numberOfDaughters(); i++) {
        const reco::Candidate * con = terminaltau->daughter(i);
        bGenTau_con_pt.push_back(con->pt());
        bGenTau_con_dr.push_back(deltaR(terminaltau->eta(), terminaltau->phi(), con->eta(), con->phi()));
      }
      double maxdR = -1.0;
      for (unsigned int i = 0; i < terminaltau->numberOfDaughters(); i++) {
        for (unsigned int j = i+1; j < terminaltau->numberOfDaughters(); j++) {
          const reco::Candidate * con1 = terminaltau->daughter(i);
          const reco::Candidate * con2 = terminaltau->daughter(j);
          double dR = deltaR(con1->eta(), con1->phi(), con2->eta(), con2->phi());
          if (dR > maxdR) maxdR = dR;
        }
      }
      bGenTau_con_maxdr = maxdR;
      found_one_hadtau = true;
    }
  }
  //bnGenHadTau = hadronicTaus.size();

  // gen particles
  vector<const reco::GenParticle*> gen_taus;
  vector<const reco::GenParticle*> gen_had_taus;
  vector<const reco::GenParticle*> gen_muons;
  vector<const reco::GenParticle*> gen_zs;
  for (unsigned int i = 0; i < genparticles->size(); i++) {
    const reco::GenParticle &genparticle = (*genparticles)[i];
    if (genparticle.pdgId() == 23 && genparticle.status() == 62) {
      gen_zs.push_back(&genparticle);
      bGenZ_pt.push_back(genparticle.pt());
      bGenZ_m.push_back(genparticle.mass());
      bGenZ_phi.push_back(genparticle.phi());
      bGenZ_eta.push_back(genparticle.eta());
      bGenZ_status.push_back(genparticle.status());
    }
    if (abs(genparticle.pdgId()) == 13) {
      gen_muons.push_back(&genparticle);
      bGenMuon_pt.push_back(genparticle.pt());
      bGenMuon_m.push_back(genparticle.mass());
      bGenMuon_phi.push_back(genparticle.phi());
      bGenMuon_eta.push_back(genparticle.eta());
    }
    if (abs(genparticle.pdgId()) == 15) {
      gen_taus.push_back(&genparticle);
      bGenTau_pt.push_back(genparticle.pt());
      bGenTau_m.push_back(genparticle.mass());
      bGenTau_phi.push_back(genparticle.phi());
      bGenTau_eta.push_back(genparticle.eta());
      bGenTau_status.push_back(genparticle.status());
    }
    if (abs(genparticle.pdgId()) == 15) {
      if (!isAncestorOfZ(&genparticle)) continue;
      std::vector<string> leptons = getDecay(genparticle);
      if (! (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() ||
          std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() ||
          std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() ||
          std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() ||
          std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end() ||
          std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end() ||
          std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end() ||
          std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end() )
        ) continue;
      gen_had_taus.push_back(&genparticle);
      bGenHadTau_pt.push_back(genparticle.pt());
      bGenHadTau_m.push_back(genparticle.mass());
      bGenHadTau_phi.push_back(genparticle.phi());
      bGenHadTau_eta.push_back(genparticle.eta());
      bGenHadTau_status.push_back(genparticle.status());
    }
  }  
  bnGenZ = gen_zs.size();
  bnGenMuon = gen_muons.size();
  bnGenTau = gen_taus.size();
  bnGenHadTau = gen_had_taus.size();

  // gen cutflow 
  for (const reco::GenParticle* gen_muon : gen_muons) {
    if (gen_muon->pt() > 19 && fabs(gen_muon->eta()) < 2.1) bgenpassMuon = true;
  }
  for (const reco::GenParticle* gen_tau : gen_taus) {
    if (gen_tau->pt() > 20 && fabs(gen_tau->eta()) < 2.3) bgenpassTauCand = true;
  }
  if (bgenpassMuon && bgenpassTauCand) bgenpassTauAndMuon = true; 
  for (const reco::GenParticle* gen_muon : gen_muons) {
    for (const reco::GenParticle* gen_tau : gen_taus) {
      if (!(gen_muon->pt() > 19 && fabs(gen_muon->eta()) < 2.1)) continue;
      if (!(gen_tau->pt() > 20 && fabs(gen_tau->eta()) < 2.3)) continue;
      if (gen_muon->charge() * gen_tau->charge() < 0 &&
          reco::deltaR(gen_muon->eta(),gen_muon->phi(),gen_tau->eta(),gen_tau->phi()) > 0.5) bgenpassTauMuonPair = true;
    }
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
  // muon branches
  bnMuon = passedMuons.size();
  for(const pat::Muon * muon : passedMuons) {
    double closest_gentau_dr = 99.9;
    for (unsigned int i = 0; i < genparticles->size(); i++) {
      const reco::GenParticle & genparticle = (*genparticles)[i];
      if (abs(genparticle.pdgId()) != 15) continue;
      double deltaR = reco::deltaR(genparticle.eta(),genparticle.phi(),muon->eta(),muon->phi());
      if (deltaR < closest_gentau_dr) closest_gentau_dr = deltaR;
    }
    bMuon_genTau_dR.push_back(closest_gentau_dr);
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
  // tau cand branches
  bnTauCand = tauCands.size();
  for (const pat::Jet* tau : tauCands) {
    double closest_gentau_dr = 99.9;
    for (unsigned int i = 0; i < genparticles->size(); i++) {
      const reco::GenParticle & genparticle = (*genparticles)[i];
      if (abs(genparticle.pdgId()) != 15) continue;
      double deltaR = reco::deltaR(genparticle.eta(),genparticle.phi(),tau->eta(),tau->phi());
      if (deltaR < closest_gentau_dr) closest_gentau_dr = deltaR;
    }
    bTauCand_genTau_dR.push_back(closest_gentau_dr);
  }

  // find singel candidate tau-muon system
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
  
  // filter descision
  
  // at least one muon and one tau cand
  if (bnMuon > 0) bpassMuon = true;
  if (bnTauCand > 0) bpassTauCand = true;
  if (bnMuon > 0 && bnTauCand > 0) bpassTauAndMuon = true;
  if (muon_index != -1 && tau_index != -1) bpassTauMuonPair = true;

  // tau-muon system: dR, opp sign, MT, Pzeta
  if (muon_index != -1 && tau_index != -1) {
    bpassTauMuonPair = true;
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

    if (reco::deltaR(theMuon.eta(), theMuon.phi(), theTau.eta(), theTau.phi()) > 0.5) bpassDR = true;
    if (theMuon.charge() * tau_charge < 0) bpassOppSign = true;
    if (MT < 40) bpassMT = true;
    if (Pzeta > -25) bpassPzeta = true;
  }

  // extra lepton veto
  bpassExtraMuon = !extraMuon;
  bpassExtraElec = !extraElectron;
  bpassDiMuon = !diMuon;
  if (!extraMuon && !extraElectron && !diMuon) bpassExtraLep = true;

  // btag veto
  if (!atLeastOneBTag) bpassBTag = true;

  // filter decision
  bpassAll = bpassTauMuonPair && bpassDR && bpassOppSign && bpassMT && bpassPzeta && bpassExtraLep && bpassBTag;

  // fill tree
  tree_->Fill();

  // debug: verify both accessors for muon tag are equivilant
  for (const pat::Muon &muon : *muons) {
    if ( muon.isGlobalMuon() != muon.muonID("AllGlobalMuons") )
      std::cout << "disagree global: " << muon.isGlobalMuon() << ", " << muon.muonID("AllGlobalMuons") << std::endl;
    if ( muon.isTrackerMuon() != muon.muonID("AllTrackerMuons") )
      std::cout << "diagree tracker: " << muon.isTrackerMuon() << ", " << muon.muonID("AllTrackerMuons") << std::endl; 
  }

}

bool
ZtoTauHadTreeMaker::isAncestorOfZ(const reco::Candidate * particle)
{
  if(particle->pdgId() == 23) return true;
  for(size_t i=0; i<particle->numberOfMothers(); i++)
  {
    if(isAncestorOfZ(particle->mother(i))) return true;
  }
  return false;
}

bool
ZtoTauHadTreeMaker::notTerminalTau(const reco::Candidate * particle)
{
  for (unsigned int i = 0; i < particle->numberOfDaughters(); i++) {
    const reco::Candidate * d = particle->daughter(i);
    if (abs(d->pdgId()) == 15) return true;  
  }
  return false;
}

const reco::Candidate *
ZtoTauHadTreeMaker::tauDaughter(const reco::Candidate * particle)
{
  for (unsigned int i = 0; i < particle->numberOfDaughters(); i++) {
    const reco::Candidate * d = particle->daughter(i);
    if (abs(d->pdgId()) == 15) return d;
  }
  return NULL;
}

vector<string>
ZtoTauHadTreeMaker::getDecay(const reco::Candidate & genparticle, int flag)
{
  vector<string> products;

  if (flag == 1) { // parent is tau
    if (genparticle.pdgId() == 111) { products.push_back("npion"); return products; }
    if (abs(genparticle.pdgId()) == 211) { products.push_back("cpion"); return products; }
  }

  // ignore quarks and gluons and their daughters, unless status 21
  if (genparticle.status() != 21 && genparticle.pdgId() >= 1 && genparticle.pdgId() <= 8) { return products; }
  if (genparticle.status() != 21 && genparticle.pdgId() == 21) { return products; }

  if (genparticle.pdgId() == 11) { products.push_back("e-"); return products; }
  if (genparticle.pdgId() == -11) { products.push_back("e+"); return products; }
  if (genparticle.pdgId() == 13) { products.push_back("mu-"); return products; }
  if (genparticle.pdgId() == -13) { products.push_back("mu+"); return products; }

  if (abs(genparticle.pdgId()) == 15) {
    vector<string> tau_products; 
    for (unsigned int j = 0; j < genparticle.numberOfDaughters(); j++) {
      const reco::Candidate* daughter = genparticle.daughter(j);
      vector<string> daughter_products = getDecay(*daughter, 1);
      tau_products.insert(tau_products.end(), daughter_products.begin(), daughter_products.end());
    }
  
    if (tau_products.size() == 0) {
      if (genparticle.pdgId() == 15) { products.push_back("tau-unid"); return products; }
      if (genparticle.pdgId() == -15) { products.push_back("tau+unid"); return products; }
    }

    if (tau_products.size() == 1 && tau_products[0] == "e-") { products.push_back("tau-e"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "e+") { products.push_back("tau+e"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "mu-") { products.push_back("tau-mu"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "mu+") { products.push_back("tau+mu"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau-e") { products.push_back("tau-e"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau+e") { products.push_back("tau+e"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau-mu") { products.push_back("tau-mu"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau+mu") { products.push_back("tau+mu"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau+had10") { products.push_back("tau+had10"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau+had1") { products.push_back("tau+had1"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau+had30") { products.push_back("tau+had30"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau+had3") { products.push_back("tau+had3"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau-had10") { products.push_back("tau-had10"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau-had1") { products.push_back("tau-had1"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau-had30") { products.push_back("tau-had30"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau-had3") { products.push_back("tau-had3"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau+unid") { products.push_back("tau+unid"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau-unid") { products.push_back("tau-unid"); return products; }

    for (string prod : tau_products) {
      if (prod != "cpion" && prod != "npion") {  
        if (genparticle.pdgId() == 15) { products.push_back("tau-unid"); return products; }
        if (genparticle.pdgId() == -15) { products.push_back("tau+unid"); return products; }
      }
    }

    bool neutral = false;
    if (std::find(tau_products.begin(), tau_products.end(), "npion") != tau_products.end() ) neutral = true;

    int pion_count = 0;
    while( std::find(tau_products.begin(), tau_products.end(), "cpion") != tau_products.end() ) 
    {
      tau_products.erase( std::find(tau_products.begin(), tau_products.end(), "cpion") );
      pion_count += 1;
    }

    int charge = 0;
    if (genparticle.pdgId() == 15) charge = -1;
    if (genparticle.pdgId() == -15) charge = 1;

    if      (charge>0 && neutral  && pion_count == 1) products.push_back("tau+had10");
    else if (charge>0 && !neutral && pion_count == 1) products.push_back("tau+had1");
    else if (charge>0 && neutral  && pion_count == 3) products.push_back("tau+had30");
    else if (charge>0 && !neutral && pion_count == 3) products.push_back("tau+had3");
    else if (charge<0 && neutral  && pion_count == 1) products.push_back("tau-had10");
    else if (charge<0 && !neutral && pion_count == 1) products.push_back("tau-had1");
    else if (charge<0 && neutral  && pion_count == 3) products.push_back("tau-had30");
    else if (charge<0 && !neutral && pion_count == 3) products.push_back("tau-had3");
    else if (charge>0)                                products.push_back("tau+unid");
    else if (charge<0)                                products.push_back("tau-unid");

    return products;
    
  }

  for (unsigned int j = 0; j < genparticle.numberOfDaughters(); j++) {
    const reco::Candidate* daughter = genparticle.daughter(j);
    vector<string> daughter_products = getDecay(*daughter);
    products.insert(products.end(), daughter_products.begin(), daughter_products.end());
  }
  
  return products;
}


void
ZtoTauHadTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  descriptions.add("ZtoTauHadTreeMakerFilter", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZtoTauHadTreeMaker);
