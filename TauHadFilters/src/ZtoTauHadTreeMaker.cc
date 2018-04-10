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
      bool cfg_recoOnly;
      bool cfg_tauObjs;

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
      // cutflow
      int bpassMuon;
      int bpassTauCand;
      int bpassTauMuonPair;
      int bpassDR;
      int bpassMT;
      int bpassPzeta;
      int bpassExtraLep;
      int bpassBTag;
      int bpassAll;
      // genparticles
      double bdecayType;
      int bnGenHadTau;
      vector<Double_t> bGenTau_con_pt;
      vector<Double_t> bGenTau_con_dr;
      double bGenTau_con_maxdr;
      int bGenTau_nCon;
      // reco
      int bnTaucands;
      vector<Double_t> btaucand_pt;
      vector<Double_t> btaucand_eta;
      vector<Double_t> btaucand_phi;
      vector<Double_t> btaucand_e;
      int bnMuoncands;
      vector<Double_t> bmuoncand_pt;
      vector<Double_t> bmuoncand_eta;
      vector<Double_t> bmuoncand_phi;
      vector<Double_t> bmuoncand_e;
      int bmuonIndex;
      int btauIndex;
      double bmVis;
};

//
// constructors
//
ZtoTauHadTreeMaker::ZtoTauHadTreeMaker(const edm::ParameterSet& iConfig) :
  cfg_recoOnly(iConfig.getUntrackedParameter<bool>("recoOnly")),
  cfg_tauObjs(iConfig.getUntrackedParameter<bool>("tauObjs"))
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
  tree_->Branch("passTauMuonPair",&bpassTauMuonPair,"passTauMuonPair/I");
  tree_->Branch("passDR",&bpassDR,"passDR/I");
  tree_->Branch("passMT",&bpassMT,"passMT/I");
  tree_->Branch("passPzeta",&bpassPzeta,"passPzeta/I");
  tree_->Branch("passExtraLep",&bpassExtraLep,"passExtraLep/I");
  tree_->Branch("passBTag",&bpassBTag,"passBTag/I");
  tree_->Branch("passAll",&bpassAll,"passAll/I");
  // tau cands
  tree_->Branch("nTaucands",&bnTaucands,"nTaucands/I"); 
  tree_->Branch("taucand_pt",&btaucand_pt); 
  tree_->Branch("taucand_eta",&btaucand_eta); 
  tree_->Branch("taucand_phi",&btaucand_phi); 
  tree_->Branch("taucand_e",&btaucand_e); 
  // muon cands
  tree_->Branch("nMuoncands",&bnMuoncands,"nMuoncands/I"); 
  tree_->Branch("muoncand_pt",&bmuoncand_pt); 
  tree_->Branch("muoncand_eta",&bmuoncand_eta); 
  tree_->Branch("muoncand_phi",&bmuoncand_phi); 
  tree_->Branch("muoncand_e",&bmuoncand_e);
  // reco
  tree_->Branch("mVis",&bmVis,"mVis/D");
  tree_->Branch("muonIndex",&bmuonIndex,"muonIndex/I");
  tree_->Branch("tauIndex",&btauIndex,"tauIndex/I");
  // genparticles
  tree_->Branch("GenTau_nCon",&bGenTau_nCon,"GenTau_nCon/I");
  tree_->Branch("GenTau_con_pt",&bGenTau_con_pt);
  tree_->Branch("GenTau_con_dr",&bGenTau_con_dr);
  tree_->Branch("GenTau_con_maxdr",&bGenTau_con_maxdr,"GenTau_con_maxdr/D");
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

  edm::Handle<vector<reco::GenParticle>> genparticles;
  iEvent.getByToken(genToken_, genparticles);

  edm::Handle<double> rho;
  iEvent.getByToken(rhoToken_, rho);

  // initialize branches
  bevent = iEvent.id().event();
  brun = iEvent.id().run();
  blumi = iEvent.id().luminosityBlock();

  bpassMuon = false;
  bpassTauCand = false;
  bpassTauMuonPair = false;
  bpassDR = false;
  bpassMT = false;
  bpassPzeta = false;
  bpassExtraLep = false;
  bpassBTag = false;
  bpassAll = false;

  bnGenHadTau = 0;

  bdecayType = -1;
  bmVis = -1;
  bmuonIndex = -1;
  btauIndex = -1;
  bnMuoncands = -1;
  bmuoncand_pt.clear();
  bmuoncand_eta.clear();
  bmuoncand_phi.clear();
  bmuoncand_e.clear();
  bnTaucands = -1;
  btaucand_pt.clear();
  btaucand_eta.clear();
  btaucand_phi.clear();
  btaucand_e.clear();

  if (!cfg_recoOnly) {
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
  
  } // end cfg_recoOnly check

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
        //(muon.chargedHadronIso() + muon.neutralHadronIso() + muon.photonIso())/muon.pt() - 0.5 * (*rho) < 0.1 &&
        (muon.chargedHadronIso() + muon.neutralHadronIso() + muon.photonIso())/muon.pt() < 0.1 &&
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

  // final filter decision
  bool passAll = passTauMuonPair && passDR && passMT && passPzeta && passExtraLep && passBTag;

  // cutflow branches
  bpassMuon = (passedMuons.size() > 0);
  if (!cfg_tauObjs)  bpassTauCand = (tauJetCands.size() > 0);
  else bpassTauCand = (tauCands.size() > 0);
  bpassTauMuonPair = passTauMuonPair;
  bpassDR = passDR;
  bpassMT = passMT;
  bpassPzeta = passPzeta;
  bpassExtraLep = passExtraLep;
  bpassBTag = passBTag;
  bpassAll = passAll;

  // visible mass reconstruction
  if (passAll && cfg_tauObjs) {
      const pat::Muon & theMuon = *passedMuons[muon_index];
      const pat::Tau & theTau = *tauCands[tau_index];
      TLorentzVector muon_mom; muon_mom.SetPtEtaPhiM(theMuon.pt(), theMuon.eta(), theMuon.phi(), theMuon.mass());
      TLorentzVector tau_mom; tau_mom.SetPtEtaPhiM(theTau.pt(), theTau.eta(), theTau.phi(), theTau.mass());
      TLorentzVector Z_visible;
      Z_visible = muon_mom + tau_mom;
      bmVis = Z_visible.M();
      // reco branches
      bnMuoncands = passedMuons.size();
      bnTaucands = tauCands.size();
      for(const pat::Muon* mu : passedMuons) {
        bmuoncand_pt.push_back(mu->pt());
      }
      for(const pat::Tau* tau : tauCands) {
        btaucand_pt.push_back(tau->pt());
      }
      bmuonIndex = muon_index;
      btauIndex = tau_index;
  }
  if (passAll && !cfg_tauObjs) {
      const pat::Muon & theMuon = *passedMuons[muon_jet_index];
      const pat::Jet & theTau = *tauJetCands[tau_jet_index];
      TLorentzVector muon_mom; muon_mom.SetPtEtaPhiM(theMuon.pt(), theMuon.eta(), theMuon.phi(), theMuon.mass());
      TLorentzVector tau_mom; tau_mom.SetPtEtaPhiM(theTau.pt(), theTau.eta(), theTau.phi(), theTau.mass());
      TLorentzVector Z_visible;
      Z_visible = muon_mom + tau_mom;
      bmVis = Z_visible.M();
      // reco branches
      bnMuoncands = passedMuons.size();
      bnTaucands = tauJetCands.size();
      for(const pat::Muon* mu : passedMuons) {
        bmuoncand_pt.push_back(mu->pt());
      }
      for(const pat::Jet* tau : tauJetCands) {
        btaucand_pt.push_back(tau->pt());
      }
      bmuonIndex = muon_index;
      btauIndex = tau_jet_index;
  }


  // fill tree
  tree_->Fill();
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
  desc.addUntracked<bool>("recoOnly", false);
  desc.addUntracked<bool>("tauObjs", false);
  descriptions.add("ZtoTauHadTreeMakerFilter", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZtoTauHadTreeMaker);
