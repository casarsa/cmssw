import FWCore.ParameterSet.Config as cms

# AM-based pattern recognition default sequence
TTPatternsFromStub = cms.EDProducer("TrackFindingAMProducer",
   TTInputStubs       = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
   TTPatternName      = cms.string("AML1Patterns"),
   inputBankFile      = cms.string('/afs/cern.ch/work/s/sviret/public/PatternBanks/BE_5D/Eta6_Phi8/ss32/620_SLHC7_MUBANK_lowmidhig_sec24_ss32_cov40.pbk'),
   threshold          = cms.int32(5)
)


## Trackfit default sequence

doRetinaFit = True

TTTracksFromPattern = ( cms.EDProducer("TrackFitRetinaProducer",
                                       TTInputStubs       = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
                                       TTInputPatterns    = cms.InputTag("MergePROutput", "AML1Patterns"),
                                       TTTrackName        = cms.string("AML1Tracks"),
                                       verboseLevel       = cms.untracked.int32(1),
                                       fitPerTriggerTower = cms.untracked.bool(False),
                                       removeDuplicates   = cms.untracked.int32(1)
                                       )
                        if doRetinaFit else 
                        cms.EDProducer("TrackFitHoughProducer",
                                       TTInputStubs       = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
                                       TTInputPatterns    = cms.InputTag("MergePROutput", "AML1Patterns"),
                                       TTTrackName        = cms.string("AML1Tracks"),
                                       ) 
                        )


# AM output merging sequence
MergePROutput = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTStubsFromPixelDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("TTPatternsFromStub", "AML1Patterns")),
   TTFiltClustersName  = cms.string("ClusInPattern"),
   TTFiltStubsName     = cms.string("StubInPattern"),
   TTPatternsName      = cms.string("AML1Patterns")                         
)

MergeFITOutput = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTStubsFromPixelDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("TTTracksFromPattern", "AML1Tracks")),
   TTFiltClustersName  = cms.string("ClusInTrack"),
   TTFiltStubsName     = cms.string("StubInTrack"),
   TTPatternsName      = cms.string("AML1Tracks")                         
)
