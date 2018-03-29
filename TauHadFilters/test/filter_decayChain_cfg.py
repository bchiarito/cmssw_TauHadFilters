import FWCore.ParameterSet.Config as cms

process = cms.Process("USER")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
  #'file:/cms/chiarito/samples/ztoll/miniaod/merged10k_30prong_pt40_eta25_MINIAOD.root'
  '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/00271851-3CC8-E611-AFB9-002590D0AFD8.root' # 10-50
  #'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/100000/00099D43-77ED-E611-8889-5065F381E1A1.root' # 50+

  ))

process.filt = cms.EDFilter('ZtoTauHadTruthSelector',
  filterByTruthDecayType = cms.untracked.vdouble(0.0),
  #ptMin = cms.untracked.double(40.0),
  #absEtaMax = cms.untracked.double(2.5), 
  )

# Draws gen particle decay chain
process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
                                   src = cms.InputTag("prunedGenParticles"),
                                   printStatus = cms.untracked.bool(True),
                                   printP4 = cms.untracked.bool(False),
                                   printPtEtaPhi = cms.untracked.bool(False),
                                   printVertex = cms.untracked.bool(False),
                                   printIndex = cms.untracked.bool(False),
                                   #status = cms.untracked.vint32( 3 )
)

process.p = cms.Path(process.filt * process.printTree)
