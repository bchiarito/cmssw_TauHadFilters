import FWCore.ParameterSet.Config as cms

process = cms.Process("USER")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
  'file:/cms/chiarito/samples/ztoll/miniaod/merged10k_30prong_pt40_eta25_MINIAOD.root'
  ))
process.TFileService = cms.Service("TFileService", fileName = cms.string( "output_filter_tree.root" ) )

process.filt = cms.EDFilter('ZtoTauHadTruthSelector',
  filterByTruthDecayType = cms.untracked.vdouble(5.4),
  ptMin = cms.untracked.double(40.0),
  absEtaMax = cms.untracked.double(2.5), 
  )

process.tree = cms.EDAnalyzer('ZtoTauHadTreeMaker',
  )

process.p = cms.Path(process.filt * process.tree)
