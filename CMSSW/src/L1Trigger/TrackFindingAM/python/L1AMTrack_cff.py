import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.MagneticField_cff import *
from L1Trigger.TrackFindingAM.L1AMTrack_cfi import *
from SimTracker.TrackTriggerAssociation.TTStubAssociation_cfi import * 
from SimTracker.TrackTriggerAssociation.TTClusterAssociation_cfi import *
from SimTracker.TrackTriggerAssociation.TTTrackAssociation_cfi import *
import FWCore.ParameterSet.Config as cms


############################################
# STEP 1: AM-based pattern recognition 
############################################

# The simple sequence, creates only a pattern container
# used in principle only for multi-bank PR
#
# Indeed in this case we create the filtered stub container after merging all
# the pattern container (merging the filtered stub container is not possible
# due to persistency loss) 

TTPatternsFromStubs   = cms.Sequence(TTPatternsFromStub)


# The complete sequence, creates the pattern container and the 
# container of filtered stubs/clusters, with corresponding
# associators containers

TTPatternsFromStubswStubs   = cms.Sequence(TTPatternsFromStub*MergePROutput*TTStubAssociatorFromPixelDigis)

############################################
# STEP 2: TC builder 
############################################

# The simple sequence, creates only a track container
# used in principle only for debugging purposes
#



TTTCsFromPatterns  = cms.Sequence(TTTCsFromPattern)


# The sequence. Note that we call the Merge plugins because the filtered containers are created
# here. We just merge one branch...

TTTCsFromPatternswStubs   = cms.Sequence(TTTCsFromPattern*MergeTCOutput*MergeTCOutputb*TTStubAssociatorFromPixelDigis*TTTrackAssociatorFromPixelDigis)

############################################
# STEP 3: PCA fit 
############################################

# The simple sequence, creates only a track container
# used in principle only for debugging purposes
#

#TTTrackAssociatorFromPixelDigis.TTTracks = cms.VInputTag( cms.InputTag("MergeFITOutput", "AML1Tracks"))

TTTracksFromTCs  = cms.Sequence(TTTracksINFNFromTC)
#TTTracksFromTCs  = cms.Sequence(TTTracksTAMUFromTC)


# The sequence. Note that we call the Merge plugins because the filtered containers are created
# here. We just merge one branch...

#TTTracksFromTCswStubs   = cms.Sequence(TTTracksINFNFromTC*MergeFITOutput*TTStubAssociatorFromPixelDigis*TTTrackAssociatorFromPixelDigis)
TTTracksFromTCswStubs   = cms.Sequence(TTTracksINFNFromTC*MergeFITOutput*MergeFITOutputb*TTStubAssociatorFromPixelDigis*TTTrackAssociatorFromPixelDigis)

############################################
# STEP 4: MERGE PR outputs
############################################

# This sequence is used mainly the multi-bank merging process, please note that the filtered cluster container is
# not associated due to the lack of simPixelDigis in official samples

TTStubAssociatorFromPixelDigis.TTStubs        = cms.VInputTag( cms.InputTag("MergePROutput", "StubInPattern"))
TTStubAssociatorFromPixelDigis.TTClusterTruth = cms.VInputTag( cms.InputTag("TTClusterAssociatorFromPixelDigis","ClusterAccepted"))

MergePROutputs  = cms.Sequence(MergePROutput*TTStubAssociatorFromPixelDigis)


