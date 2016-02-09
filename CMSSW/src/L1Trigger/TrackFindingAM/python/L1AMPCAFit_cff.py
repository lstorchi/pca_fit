import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.MagneticField_cff import *
from L1Trigger.TrackFindingAM.L1AMPCAFit_cfi import *
#from SimTracker.TrackTriggerAssociation.TTStubAssociation_cfi import *
import SimTracker.TrackTriggerAssociation.TTStubAssociation_cfi 
from SimTracker.TrackTriggerAssociation.TTClusterAssociation_cfi import *
import SimTracker.TrackTriggerAssociation.TTTrackAssociation_cfi 



TTTrackAssociatorFromPixelDigisTC = SimTracker.TrackTriggerAssociation.TTTrackAssociation_cfi.TTTrackAssociatorFromPixelDigis.clone()
TTTrackAssociatorFromPixelDigisTC.TTTracks = cms.VInputTag( cms.InputTag("MergeFITOutputTC", "AML1TracksFromTC"))

TTTrackAssociatorFromPixelDigisPCA = SimTracker.TrackTriggerAssociation.TTTrackAssociation_cfi.TTTrackAssociatorFromPixelDigis.clone()
TTTrackAssociatorFromPixelDigisPCA.TTTracks = cms.VInputTag( cms.InputTag("MergeFITOutput", "AML1TracksFromPCA"))

TTStubAssociatorFromPixelDigisTC = SimTracker.TrackTriggerAssociation.TTStubAssociation_cfi.TTStubAssociatorFromPixelDigis.clone()
TTStubAssociatorFromPixelDigisTC.TTStubs = cms.VInputTag( cms.InputTag("MergeFITOutputTC", "StubInTrack"))
TTStubAssociatorFromPixelDigisTC.TTClusterTruth = cms.VInputTag( cms.InputTag("TTClusterAssociatorFromPixelDigis","ClusterAccepted"))

TTStubAssociatorFromPixelDigisPCA = SimTracker.TrackTriggerAssociation.TTStubAssociation_cfi.TTStubAssociatorFromPixelDigis.clone()
TTStubAssociatorFromPixelDigisPCA.TTStubs = cms.VInputTag( cms.InputTag("MergeFITOutput", "StubInTrack"))
TTStubAssociatorFromPixelDigisPCA.TTClusterTruth = cms.VInputTag( cms.InputTag("TTClusterAssociatorFromPixelDigis","ClusterAccepted"))

# The sequence. Note that we call the Merge plugins because the filtered containers are created
# here. We just merge one branch...

TTTracksFromTCwStubs   = cms.Sequence(TTTracksFromTC
                                      *MergeFITOutputTC 
                                      *TTStubAssociatorFromPixelDigisTC
                                      *TTTrackAssociatorFromPixelDigisTC
)

TTTracksFromPCAwStubs   = cms.Sequence(TTTracksFromPattern
                                            *MergeFITOutput 
                                            *TTStubAssociatorFromPixelDigisPCA
                                            *TTTrackAssociatorFromPixelDigisPCA
)
TTTracksFromPatternswStubs   = cms.Sequence(TTTracksFromTCwStubs*TTTracksFromPCAwStubs)
