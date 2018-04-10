import FWCore.ParameterSet.Config as cms

process = cms.Process("USER")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(3) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
  'file:/cms/chiarito/samples/ztoll/miniaod/DYJetsToLL_M-50_136kevents_02A210D6-F5C3-E611-B570-008CFA197BD4.root'
  ))

process.filt = cms.EDFilter('ZtoTauHadTruthSelector',
  filterByTruthDecayType = cms.untracked.vdouble(),
  )

process.p = cms.Path(process.filt)
