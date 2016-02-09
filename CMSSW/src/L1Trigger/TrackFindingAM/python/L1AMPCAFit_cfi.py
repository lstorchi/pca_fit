import FWCore.ParameterSet.Config as cms

TTTracksFromTC = ( cms.EDProducer("TrackFitTCProducer",
                                       TTInputStubs       = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
                                       TTInputPatterns    = cms.InputTag("MergePROutput", "AML1Patterns"),
                                       TTTrackName        = cms.string("AML1TracksFromTC"),
                                       )
)


MergeFITOutputTC = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTStubsFromPixelDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("TTTracksFromTC", "AML1TracksFromTC")),                               
   TTFiltClustersName  = cms.string("ClusInTrack"),
   TTFiltStubsName     = cms.string("StubInTrack"),
   TTPatternsName      = cms.string("AML1TracksFromTC")                         
)

TTTracksFromPattern = ( cms.EDProducer("TrackFitPCAProducer",
                                       TTInputStubs       = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
                                       TTInputPatterns    = cms.InputTag("MergePROutput", "AML1Tracks"),
                                       TCBuilderTTracks   = cms.InputTag("MergeFITOutputTC", "AML1TracksFromTC"),
                                       TTTrackName        = cms.string("AML1TracksFromPCA"),
                                       )
)

MergeFITOutput = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTStubsFromPixelDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("TTTracksFromPattern", "AML1TracksFromPCA")),                               
   TTFiltClustersName  = cms.string("ClusInTrack"),
   TTFiltStubsName     = cms.string("StubInTrack"),
   TTPatternsName      = cms.string("AML1TracksFromPCA")                         
)
