import FWCore.ParameterSet.Config as cms

process = cms.Process("USER")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
  #'file:/cms/chiarito/samples/ztoll/miniaod/merged10k_30prong_pt40_eta25_MINIAOD.root'
  #'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/00271851-3CC8-E611-AFB9-002590D0AFD8.root' # 10-50
  '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/100000/00099D43-77ED-E611-8889-5065F381E1A1.root' # 50+
  ))
process.TFileService = cms.Service("TFileService", fileName = cms.string( "output_tree_tauobjs.root" ) )

process.tree = cms.EDAnalyzer('ZtoTauHadCutflowMaker',
  tauObjs = cms.untracked.bool(False)
  )

process.p = cms.Path(process.tree)