import FWCore.ParameterSet.Config as cms

process = cms.Process("demo")

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
    )

process.load("FWCore.MessageService.MessageLogger_cfi")

process.ZCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
                                    decay = cms.string("selectedPatMuons@+ selectedPatMuons@-"),
                                    cut   = cms.string("81. < mass < 101.")
                                    )

process.CompCandDump = cms.EDAnalyzer("ZcandHisto",
                                      ZCandCollection  = cms.untracked.string('ZCandidate'),
                                      pfCandCollection = cms.untracked.string('pfNoElectronPFlow')
                                      )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1)
                                        )



#process.load("AWhitbeck.JetHisto.testSample_cff")
#process.load("AWhitbeck.JetHisto.ZHbbSample_cff")
process.load("AWhitbeck.JetHisto.ZjetsSample_cff")

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('histo.root')
                                   )

process.p = cms.Path(process.ZCandidate
                     + process.CompCandDump)
