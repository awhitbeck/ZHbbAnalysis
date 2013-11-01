import FWCore.ParameterSet.Config as cms

process = cms.Process("demo")

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
    )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options   = cms.untracked.PSet(
    SkipEvent   = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
    )

##################################
## FOR SELECTING EVENTS WHICH HAVE
##        two or more muons
##        muon pt>10
##################################
process.minPtMuons = cms.EDFilter("PtMinCandViewSelector",
                                  src = cms.InputTag("selectedPatMuons"),
                                  ptMin = cms.double(10)
                                  )

process.atLeastTwoMuons = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("minPtMuons"),
                                       minNumber = cms.uint32(2)
                                       )

## BUILD Z CANDIDATES
process.ZCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
                                    decay = cms.string("minPtMuons@+ minPtMuons@-"),
                                    cut   = cms.string("81. < mass < 101.")
                                    )

process.minPtZCand = cms.EDFilter("PtMinCandViewSelector",
                                  src = cms.InputTag("ZCandidate"),
                                  ptMin = cms.double(120)
                                  )

##################################
## FOR SELECTING EVENTS WHICH HAVE
##        two or more jets
##        jet pt>20
##################################
process.minPtJets = cms.EDFilter("PtMinCandViewSelector",
                                 src = cms.InputTag("selectedPatJetsPFlowLoose"),
                                 ptMin = cms.double(10)
                                 )

process.atLeastTwoJets = cms.EDFilter("CandViewCountFilter",
                                      src = cms.InputTag("minPtJets"),
                                      minNumber = cms.uint32(2)
                                      )


## PERFORM TELESCOPING AND FILL TREE
process.TreeFiller = cms.EDAnalyzer("ZHbbTreeFiller",
                                    ZCandCollection  = cms.untracked.string('ZCandidate'),
                                    # USED FOR SIGNAL
                                    #pfCandCollection = cms.untracked.string('pfNoElectronPFlow')
                                    # USED FOR BACKGROUND
                                    pfCandCollection = cms.untracked.string('selectedPatJetsPFlow')
                                    # MUST FIX!!! PAT-TUPLES NOT CREATED WITH CONSISTENT VERSIONS
                                    # OF CODE... NEED TO REMAKE SOMETHING
                                    )

## CONFIGURE TFILESERVICE

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("ZHbbAnalysisTree.root"),
      closeFileFast = cms.untracked.bool(True)
  )

##  MAXIMUM NUMBER OF EVENTS TO PROCESS
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

##  LOAD DATAFILES
#process.load("AWhitbeck.JetHisto.testSample_cff")
#process.load("AWhitbeck.JetHisto.ZHbbSample_cff")
process.load("AWhitbeck.JetHisto.ZjetsSample_cff")


##  DEFINE PATH
process.p = cms.Path(process.minPtMuons*process.atLeastTwoMuons*process.ZCandidate*process.minPtZCand
                     +process.minPtJets*process.atLeastTwoJets
                     +process.TreeFiller)


##  OUPUT CONFIGURATION
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('test.root'),
                               #save only events passing the full path
                               SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               outputCommands = cms.untracked.vstring('keep *_*Candidate_*_*',
                                                                      'keep *_selectedPat*_*_*')
                               )

process.outpath = cms.EndPath(process.out)
