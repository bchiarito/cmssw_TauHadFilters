// -*- C++ -*-
//
// Package:    TauHadFilters/ZtoTauHadTruthSelector
// Class:      ZtoTauHadTruthSelector
// 
/**\class ZtoTauHadTruthSelector ZtoTauHadTruthSelector.cc TauHadFilters/ZtoTauHadTruthSelector/plugins/ZtoTauHadTruthSelector.cc

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
class ZtoTauHadTruthSelector : public edm::stream::EDFilter<> {
   public:
      explicit ZtoTauHadTruthSelector(const edm::ParameterSet&);
      ~ZtoTauHadTruthSelector();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      vector<string> getDecay(const reco::Candidate&, int flag=0);
      bool isAncestorOfZ(const reco::Candidate *);
      bool notTerminalTau(const reco::Candidate *);
      const reco::Candidate* tauDaughter(const reco::Candidate *);

      virtual void beginStream(edm::StreamID) override;
      virtual void endStream() override;
      virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // Configuration Parameters
      std::vector<double> cfg_filterByTruthDecayType;
      double cfg_ptcut;
      double cfg_etacut;

      // EDM Collection lables
      edm::EDGetTokenT<vector<reco::GenParticle>> genToken_;
};

//
// constructors
//
ZtoTauHadTruthSelector::ZtoTauHadTruthSelector(const edm::ParameterSet& iConfig) :
  cfg_filterByTruthDecayType(iConfig.getUntrackedParameter<std::vector<double>>("filterByTruthDecayType")),
  cfg_ptcut(iConfig.getUntrackedParameter<double>("ptMin")),
  cfg_etacut(iConfig.getUntrackedParameter<double>("absEtaMax"))
{
  genToken_ = consumes<vector<reco::GenParticle>>(edm::InputTag("prunedGenParticles"));
}

//
// destructor
//
ZtoTauHadTruthSelector::~ZtoTauHadTruthSelector()
{
}

// ------------ method called on each new Event  ------------
bool
ZtoTauHadTruthSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // get event content
  edm::Handle<vector<reco::GenParticle>> genparticles;
  iEvent.getByToken(genToken_, genparticles);

  // decay type
  double decayType = -10.0;
  vector<string> leptons;
  for (unsigned int i = 0; i < genparticles->size(); i++) {
    const reco::GenParticle & genparticle = (*genparticles)[i];
    if (genparticle.status() != 21 && genparticle.status() != 22) continue;
    leptons = getDecay(genparticle);
  }
  if (leptons.size() == 0) 
    decayType = 0;
  if (leptons.size() == 1) 
    decayType = -1;
  if (leptons.size() == 2) 
  {
    if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end()) decayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end()) decayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end()) decayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end()) decayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end()) decayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end()) decayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end()) decayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end()) decayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end()) decayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end()) decayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end()) decayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end()) decayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end()) decayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end()) decayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end()) decayType = 3;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end()) decayType = 3;

    else if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end()) decayType = 4;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end()) decayType = 4;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end()) decayType = 4;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end()) decayType = 4;
    else if (std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end()) decayType = 4;
    else if (std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end()) decayType = 4;
    else if (std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end()) decayType = 4;
    else if (std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end()) decayType = 4;

    else if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) decayType = 5.2;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) decayType = 5.1;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) decayType = 5.4;
    else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) decayType = 5.3;
    else if (std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end()) decayType = 5.2;
    else if (std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end()) decayType = 5.1;
    else if (std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end()) decayType = 5.4;
    else if (std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end()) decayType = 5.3;

    else if (std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end()) decayType = 6;

    else if (std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) decayType = 7;

    else if (std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) decayType = 8;
    else if (std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end()) decayType = 8;

    else if (std::find(leptons.begin(), leptons.end(), "e+") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "e-") != leptons.end()) decayType = 1;

    else if (std::find(leptons.begin(), leptons.end(), "mu+") != leptons.end() &&
        std::find(leptons.begin(), leptons.end(), "mu-") != leptons.end()) decayType = 2;

    else
      decayType = 9;
  }
  if (leptons.size() > 2)
  {
    decayType = 10;
  }
  bool passCfgDecayType = false;
  for (double type : cfg_filterByTruthDecayType) if (type == decayType) passCfgDecayType = true;

  // apply kinemtics cuts to hadronic tau
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
    }
  }
  bool hadTauPassKinematics = false;
  if (hadronicTaus.size() == 0)
  {
    hadTauPassKinematics = true;
  }
  else {
    for (unsigned int i=0; i<hadronicTaus.size(); i++)
    {
      hadTauPassKinematics = hadTauPassKinematics || (hadronicTaus[i]->pt() > cfg_ptcut && fabs(hadronicTaus[i]->eta())<cfg_etacut);
    }
  }

  if (cfg_filterByTruthDecayType.size() == 0) return true;
  else return (passCfgDecayType && hadTauPassKinematics);
}

bool
ZtoTauHadTruthSelector::isAncestorOfZ(const reco::Candidate * particle)
{
  if(particle->pdgId() == 23) return true;
  for(size_t i=0; i<particle->numberOfMothers(); i++)
  {
    if(isAncestorOfZ(particle->mother(i))) return true;
  }
  return false;
}

bool
ZtoTauHadTruthSelector::notTerminalTau(const reco::Candidate * particle)
{
  for (unsigned int i = 0; i < particle->numberOfDaughters(); i++) {
    const reco::Candidate * d = particle->daughter(i);
    if (abs(d->pdgId()) == 15) return true;  
  }
  return false;
}

const reco::Candidate *
ZtoTauHadTruthSelector::tauDaughter(const reco::Candidate * particle)
{
  for (unsigned int i = 0; i < particle->numberOfDaughters(); i++) {
    const reco::Candidate * d = particle->daughter(i);
    if (abs(d->pdgId()) == 15) return d;
  }
  return NULL;
}

vector<string>
ZtoTauHadTruthSelector::getDecay(const reco::Candidate & genparticle, int flag)
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
ZtoTauHadTruthSelector::beginStream(edm::StreamID)
{
}

void
ZtoTauHadTruthSelector::endStream() {
}

void
ZtoTauHadTruthSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
 
void
ZtoTauHadTruthSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
 
void
ZtoTauHadTruthSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

void
ZtoTauHadTruthSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

void
ZtoTauHadTruthSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  std::vector<double> empty;
  desc.addUntracked<vector<double>>("filterByTruthDecayType",empty);
  desc.addUntracked<double>("ptMin", 0.0);
  desc.addUntracked<double>("absEtaMax", 99999.9);
  descriptions.add("ZtoTauHadTruthSelectorFilter", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZtoTauHadTruthSelector);
