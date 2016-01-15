/*! \class   TrackFitPCAProducer
 *
 * Interface plugin for the PCA fitter 
 *
 *  Copy and paste from plugin for the TC builder by S Viret / G Baulieu / G Galbit
 *  \author L Storchi / A Modak / SR Chowdhury
 *  \date   2016
 *
 */


#ifndef TRACK_FITTER_AM_H
#define TRACK_FITTER_AM_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"

#include "L1Trigger/TrackFindingAM/interface/CMSPatternLayer.h"
#include "L1Trigger/TrackFindingAM/interface/PatternFinder.h"
#include "L1Trigger/TrackFindingAM/interface/SectorTree.h"
#include "L1Trigger/TrackFindingAM/interface/Hit.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/shared_ptr.hpp>
#include <memory>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

//#ifndef __APPLE__
//BOOST_CLASS_EXPORT_IMPLEMENT(CMSPatternLayer)
//#endif

class TrackFitPCAProducer : public edm::EDProducer
{
  public:
    /// Constructor
    explicit TrackFitPCAProducer( const edm::ParameterSet& iConfig );

    /// Destructor;
    ~TrackFitPCAProducer();

  private:
  
  /// Data members
  double                       mMagneticField;
  unsigned int                 nSectors;
  unsigned int                 nWedges;
  std::string                  nBKName;
  int                          nThresh;
  const StackedTrackerGeometry *theStackedTracker;
  edm::InputTag                TTStubsInputTag;
  edm::InputTag                TTPatternsInputTag;
  std::string                  TTTrackOutputTag;

  /// Mandatory methods
  virtual void beginRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void endRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void produce( edm::Event& iEvent, const edm::EventSetup& iSetup );

}; /// Close class

/*! \brief   Implementation of methods
 */

/// Constructors
TrackFitPCAProducer::TrackFitPCAProducer( const edm::ParameterSet& iConfig )
{
  TTStubsInputTag          = iConfig.getParameter< edm::InputTag >( "TTInputStubs" );
  TTPatternsInputTag       = iConfig.getParameter< edm::InputTag >( "TTInputPatterns" );
  TTTrackOutputTag         = iConfig.getParameter< std::string >( "TTTrackName" );       // Container of C++ TCs
  //TTTrackBinaryOutputTag   = iConfig.getParameter< std::string >( "TTTrackBinaryName" ); // Container of bit-wise TCs

  produces< std::vector< TTTrack< Ref_PixelDigi_ > > >( TTTrackOutputTag );
  // produces< std::vector< TTTrack< Ref_PixelDigi_ > > >( TTTrackBinaryOutputTag );
}

/// Destructor
TrackFitPCAProducer::~TrackFitPCAProducer() {}

/// Begin run
void TrackFitPCAProducer::beginRun( const edm::Run& run, const edm::EventSetup& iSetup )
{
  /// Get the geometry references
  edm::ESHandle< StackedTrackerGeometry > StackedTrackerGeomHandle;
  iSetup.get< StackedTrackerGeometryRecord >().get( StackedTrackerGeomHandle );
  theStackedTracker = StackedTrackerGeomHandle.product();

  /// Get magnetic field
  edm::ESHandle<MagneticField> magneticFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticFieldHandle);
  const MagneticField* theMagneticField = magneticFieldHandle.product();
  double mMagneticFieldStrength = theMagneticField->inTesla(GlobalPoint(0,0,0)).z();
  mMagneticField = (floor(mMagneticFieldStrength*10.0 + 0.5))/10.0;


}

/// End run
void TrackFitPCAProducer::endRun( const edm::Run& run, const edm::EventSetup& iSetup ) {}

/// Implement the producer
void TrackFitPCAProducer::produce( edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  /// Prepare output
  /// The temporary collection is used to store tracks
  /// before removal of duplicates
  std::auto_ptr< std::vector< TTTrack< Ref_PixelDigi_ > > > TTTracksForOutput( new std::vector< TTTrack< Ref_PixelDigi_ > > );

  //std::auto_ptr< std::vector< TTTrack< Ref_PixelDigi_ > > > TTTracksBinForOutput( new std::vector< TTTrack< Ref_PixelDigi_ > > );

  /// Get the Stubs already stored away
  edm::Handle< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > > > TTStubHandle;
  edm::Handle< std::vector< TTTrack< Ref_PixelDigi_ > > > TTPatternHandle;

  iEvent.getByLabel( TTStubsInputTag, TTStubHandle );       // Get the stubs
  iEvent.getByLabel( TTPatternsInputTag, TTPatternHandle ); // Get the matched roads

  /// STEP 0
  /// Prepare output
  TTTracksForOutput->clear();
  //TTTracksBinForOutput->clear();

  int layer  = 0;
  int ladder = 0;
  int module = 0;
  int nbLayers = 0;

  unsigned int j     = 0;
  unsigned int tkCnt = 0;
  std::map< unsigned int , edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > > stubMap;

  // The fitter
  TCBuilder* TCB = new TCBuilder(nbLayers);
  PCATrackFitter* pcafitter = new PCATrackFitter(nbLayers);

  //
  // TCs are build pattern per pattern
  //

  std::vector<Hit*> m_hits;
  for(unsigned int i=0;i<m_hits.size();i++) delete m_hits[i];
  std::vector<Track*> tracks, tcb_tracks;
  for(unsigned int i=0;i<tracks.size();i++) delete tracks[i];

  edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >::const_iterator inputIter;
  edmNew::DetSet< TTStub< Ref_PixelDigi_ > >::const_iterator stubIter;

  std::vector< TTTrack< Ref_PixelDigi_ > >::const_iterator iterTTTrack;
  
  /// Go on only if there are Patterns from PixelDigis
  if ( TTPatternHandle->size() > 0 )
  {
    /// Loop over Patterns
 
    for ( iterTTTrack = TTPatternHandle->begin(); iterTTTrack != TTPatternHandle->end(); ++iterTTTrack )
    {
      edm::Ptr< TTTrack< Ref_PixelDigi_ > > tempTrackPtr( TTPatternHandle, tkCnt++ );

      // Initialize everything for each new road

      j = 0;
      m_hits.clear();
      tracks.clear();
      stubMap.clear();
      tcb_tracks.clear();

      /// Get everything relevant
      unsigned int seedSector = tempTrackPtr->getSector();
      nbLayers                = tempTrackPtr->getWedge();

      // Get the stubs in the road

      std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_  > >, 
        TTStub< Ref_PixelDigi_  > > > trackStubs = tempTrackPtr->getStubRefs();

      // And loop over them
      for(unsigned int i=0;i<trackStubs.size();i++)
      {
	++j;

	edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, 
          TTStub< Ref_PixelDigi_ > > tempStubRef = trackStubs.at(i);

	stubMap.insert( std::make_pair( j, tempStubRef ) );

	/// Calculate average coordinates col/row for inner/outer Cluster
	/// These are already corrected for being at the center of each pixel
	MeasurementPoint mp0 = tempStubRef->getClusterRef(0)->findAverageLocalCoordinates();
	GlobalPoint posStub  = theStackedTracker->findGlobalPosition( &(*tempStubRef) );
	
	StackedTrackerDetId detIdStub( tempStubRef->getDetId() );
	
	const GeomDetUnit* det0 = theStackedTracker->idToDetUnit( detIdStub, 0 );
	const GeomDetUnit* det1 = theStackedTracker->idToDetUnit( detIdStub, 1 );
	
	/// Find pixel pitch and topology related information
	const PixelGeomDetUnit* pix0 = dynamic_cast< const PixelGeomDetUnit* >( det0 );
	const PixelGeomDetUnit* pix1 = dynamic_cast< const PixelGeomDetUnit* >( det1 );
	const PixelTopology* top0    = dynamic_cast< const PixelTopology* >( &(pix0->specificTopology()) );
	const PixelTopology* top1    = dynamic_cast< const PixelTopology* >( &(pix1->specificTopology()) );
	
	/// Find the z-segment
	int cols0   = top0->ncolumns();
	int cols1   = top1->ncolumns();
	int ratio   = cols0/cols1; /// This assumes the ratio is integer!
	int segment = floor( mp0.y() / ratio );
	int strip   =  mp0.x();
     
	// Here we rearrange the number in order to be compatible with the TC builder
  
	// First of all we 
	// order the stubs per layers
	// and count the number of layers touched    

	// Layers are numbered as follows
	// Barrel      : 0,1,2,8,9,10
	// Disk z+/- PS: 3,4,5,6,7
	// Disk z+/- 2S: 11,12,13,14,15

	if ( detIdStub.isBarrel() )
	{
	  layer  = detIdStub.iLayer()-1;

	  if (layer>2) layer+=5;

	  ladder = detIdStub.iPhi()-1;
	  module = detIdStub.iZ()-1;

	  //	  cout << layer << " / " << detIdStub.iLayer()+4 << endl;
	}
	else if ( detIdStub.isEndcap() )
	{
	  //layer  = 10+detIdStub.iZ()+abs((int)(detIdStub.iSide())-2)*7;
	  layer = detIdStub.iZ()+2;
	  
	  if (ratio==1) layer+=8;

	  //	  cout << layer << " / " << ratio << " / " << 10+detIdStub.iZ()+abs((int)(detIdStub.iSide())-2)*7 << endl;

	  ladder = detIdStub.iRing()-1;
	  module = detIdStub.iPhi()-1;
	}
	
	Hit* h = new Hit(layer,ladder, module, segment, strip, 
			 j, -1, 0, 0, 0, 0, 
			 posStub.x(), posStub.y(), posStub.z(), 0, 0, 0, 
			 tempStubRef->getTriggerDisplacement()-tempStubRef->getTriggerOffset());

	m_hits.push_back(h);

      } /// End of loop over track stubs

      // Do the TC
      TCB->setSectorID(tempTrackPtr->getSector());
      TCB->fit(m_hits);

      // After this we have a first estimation of track parameters
      // here add the PCA fitter 
      
      //Recover it...
      tcb_tracks = TCB->getTracks();
      TCB->clean();

      // Store the tracks (no duplicate cleaning yet)
      
      // run the PCA fitter 
      pcafitter->setTracks (tcb_tracks);
      pcafitter->fit();
      tracks = pcafitter->getTracks();

      std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, 
        TTStub< Ref_PixelDigi_ > > > tempVec;

      for(unsigned int tt=0;tt<tracks.size();tt++)
      {	
	tempVec.clear();

	//Stubs used for the fit 
	vector<int> stubs = tracks[tt]->getStubs();
	for(unsigned int sti=0;sti<stubs.size();sti++) 
          tempVec.push_back( stubMap[ stubs[sti] ] );

	TTTrack< Ref_PixelDigi_ > tempTrack( tempVec );

	double pz = tracks[tt]->getCurve()/(tan(2*atan(exp(-tracks[tt]->getEta0()))));
	GlobalPoint POCA(0.,0.,tracks[tt]->getZ0());
	GlobalVector mom(tracks[tt]->getCurve()*cos(tracks[tt]->getPhi0()),
			 tracks[tt]->getCurve()*sin(tracks[tt]->getPhi0()),
			 pz);
		
	tempTrack.setSector( seedSector );
	tempTrack.setWedge( -1 );
	tempTrack.setMomentum( mom , 5);
	tempTrack.setPOCA( POCA , 5);
	
	TTTracksForOutput->push_back( tempTrack );
	
	delete tracks[tt];
      } // End of loop over TCs
    } // End of loop over patterns
  }

  delete(TCB);    
  delete(pcafitter);
  /// Put in the event content
  iEvent.put( TTTracksForOutput, TTTrackOutputTag);
}

// DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(TrackFitPCAProducer);

#endif

