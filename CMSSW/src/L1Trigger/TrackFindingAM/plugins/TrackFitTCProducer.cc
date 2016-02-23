/*! \class   TrackFitTCProducer
 *
 *  \author S Viret / G Baulieu / G Galbit
 *  \date   2015, Mar 10
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

#include "L1Trigger/TrackFindingAM/interface/PCATrackFitter.h"

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

class TrackFitTCProducer : public edm::EDProducer
{
  public:
    /// Constructor
    explicit TrackFitTCProducer( const edm::ParameterSet& iConfig );

    /// Destructor;
    ~TrackFitTCProducer();

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
  std::string                  TTTrackBinaryOutputTag;

  /// Mandatory methods
  virtual void beginRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void endRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void produce( edm::Event& iEvent, const edm::EventSetup& iSetup );

}; /// Close class

/*! \brief   Implementation of methods
 */

/// Constructors
TrackFitTCProducer::TrackFitTCProducer( const edm::ParameterSet& iConfig )
{
  TTStubsInputTag          = iConfig.getParameter< edm::InputTag >( "TTInputStubs" );
  TTPatternsInputTag       = iConfig.getParameter< edm::InputTag >( "TTInputPatterns" );
  TTTrackOutputTag         = iConfig.getParameter< std::string >( "TTTrackName" );
  TTTrackBinaryOutputTag   = iConfig.getParameter< std::string >( "TTTrackBinaryName" );

  produces< std::vector< TTTrack< Ref_PixelDigi_ > > >( TTTrackOutputTag );
  produces< std::vector< TTTrack< Ref_PixelDigi_ > > >( TTTrackBinaryOutputTag );
}

/// Destructor
TrackFitTCProducer::~TrackFitTCProducer() {}

/// Begin run
void TrackFitTCProducer::beginRun( const edm::Run& run, const edm::EventSetup& iSetup )
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
void TrackFitTCProducer::endRun( const edm::Run& run, const edm::EventSetup& iSetup ) {}

/// Implement the producer
void TrackFitTCProducer::produce( edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  /// Prepare output
  /// The temporary collection is used to store tracks
  /// before removal of duplicates
  std::auto_ptr< std::vector< TTTrack< Ref_PixelDigi_ > > > TTTracksForOutput( new std::vector< TTTrack< Ref_PixelDigi_ > > );
  std::auto_ptr< std::vector< TTTrack< Ref_PixelDigi_ > > > TTTracksBinForOutput( new std::vector< TTTrack< Ref_PixelDigi_ > > );

  /// Get the Stubs already stored away
  edm::Handle< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > > > TTStubHandle;
  edm::Handle< std::vector< TTTrack< Ref_PixelDigi_ > > > TTPatternHandle;

  iEvent.getByLabel( TTStubsInputTag, TTStubHandle );
  iEvent.getByLabel( TTPatternsInputTag, TTPatternHandle );

  /// STEP 0
  /// Prepare output
  TTTracksForOutput->clear();
  TTTracksBinForOutput->clear();

  int layer  = 0;
  int ladder = 0;
  int module = 0;

  int nbLayers = 0;

  /// Loop over Patterns
  unsigned int tkCnt = 0;
  unsigned int j     = 0;

  std::map< unsigned int , edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > > stubMap;
  

  // PLACEHOLDER
  PCATrackFitter* PCAFIT = new PCATrackFitter(nbLayers); 

  TCBuilder* TCB  = new TCBuilder(nbLayers); // Floating point
  TCBuilder* TCBb = new TCBuilder(nbLayers); // Bit-wise
  TCBb->setHardwareEmulation(true);

  /// STEP 1
  /// Loop over patterns

  //  std::cout << "Start the loop over " << TTPatternHandle->size() << " pattern(s) in order to recover the stubs" << std::endl;

  std::vector<Hit*> m_hits;
  for(unsigned int i=0;i<m_hits.size();i++) delete m_hits[i];
  m_hits.clear();

  std::vector<Track*> tracks;
  for(unsigned int i=0;i<tracks.size();i++) delete tracks[i];
  tracks.clear();

  std::vector<Track*> tracksb;
  for(unsigned int i=0;i<tracksb.size();i++) delete tracksb[i];
  tracksb.clear();

  edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >::const_iterator inputIter;
  edmNew::DetSet< TTStub< Ref_PixelDigi_ > >::const_iterator stubIter;

  std::vector< TTTrack< Ref_PixelDigi_ > >::const_iterator iterTTTrack;
  
  /// Go on only if there are Patterns from PixelDigis
  if ( TTPatternHandle->size() > 0 )
  {
    for ( iterTTTrack = TTPatternHandle->begin();
	  iterTTTrack != TTPatternHandle->end();
	  ++iterTTTrack )
    {
      edm::Ptr< TTTrack< Ref_PixelDigi_ > > tempTrackPtr( TTPatternHandle, tkCnt++ );

      //      std::cout << tkCnt << std::endl;


      j = 0;

      m_hits.clear();
      tracks.clear();
      stubMap.clear();

      /// Get everything relevant
      unsigned int seedSector = tempTrackPtr->getSector();
      nbLayers = tempTrackPtr->getWedge();

      // Get the stubs in the road

      std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_  > >, TTStub< Ref_PixelDigi_  > > > trackStubs = tempTrackPtr->getStubRefs();

      //      std::cout << "PT1" << std::endl;

      // Loop over stubs contained in the pattern to recover the info

      for(unsigned int i=0;i<trackStubs.size();i++)
      {
	++j;

	//	std::cout << "PT2" << std::endl;

	edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > tempStubRef = trackStubs.at(i);

	stubMap.insert( std::make_pair( j, tempStubRef ) );

	/// Calculate average coordinates col/row for inner/outer Cluster
	/// These are already corrected for being at the center of each pixel
	MeasurementPoint mp0 = tempStubRef->getClusterRef(0)->findAverageLocalCoordinates();
	GlobalPoint posStub  = theStackedTracker->findGlobalPosition( &(*tempStubRef) );
	
	StackedTrackerDetId detIdStub( tempStubRef->getDetId() );
	//bool isPS = theStackedTracker->isPSModule( detIdStub );
	
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
	
	// Here we rearrange the number in order to be compatible with the AM emulator
	if ( detIdStub.isBarrel() )
	{
	  layer  = detIdStub.iLayer()+4;
	  ladder = detIdStub.iPhi()-1;
	  module = detIdStub.iZ()-1;
	}
	else if ( detIdStub.isEndcap() )
	{
	  layer  = 10+detIdStub.iZ()+abs((int)(detIdStub.iSide())-2)*7;
	  ladder = detIdStub.iRing()-1;
	  module = detIdStub.iPhi()-1;
	}

	//cout << layer << " / " << ladder << " / " << module << " / " << std::endl;

	int strip  =  mp0.x();
	int tp     = -1;
	float eta  = 0;
	float phi0 = 0;
	float spt  = 0;
	float x    = posStub.x();
	float y    = posStub.y();
	float z    = posStub.z();
	float x0   = 0.;
	float y0   = 0.;
	float z0   = 0.;
	float ip   = sqrt(x0*x0+y0*y0);

        //cout << x << " / " << y << " / " << z << " / " << std::endl;
	
	Hit* h = new Hit(layer,ladder, module, segment, strip, 
			 j, tp, spt, ip, eta, phi0, x, y, z, x0, y0, z0, tempStubRef->getTriggerDisplacement()-tempStubRef->getTriggerOffset());
	m_hits.push_back(h);

      } /// End of loop over track stubs

      //      std::cout << "PT3" << std::endl;

      TCB->setSectorID(tempTrackPtr->getSector());
      TCB->fit(m_hits);
      TCBb->setSectorID(tempTrackPtr->getSector());
      TCBb->fit(m_hits);

      //      std::cout << "PT4" << std::endl;

      tracks = TCB->getTracks();
      TCB->clean();

      // Perform the PCA fit PLACEHOLDER
      PCAFIT->setSectorID(tempTrackPtr->getSector()); 
      PCAFIT->setTracks(tracks);
      PCAFIT->fit(m_hits); 
      /*
      if (PCAFIT->getTracks().size() != 0)
        tracks = PCAFIT->getTracks();
      */

      tracksb = TCBb->getTracks();
      TCBb->clean();

      //      std::cout << "PT5" << std::endl;

      // Store the tracks (no duplicate cleaning yet)
      //cout<<"Found "<<tracks.size()<<" track"<<endl;

      std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > > tempVec;

      for(unsigned int tt=0;tt<tracks.size();tt++)
      {	
	tempVec.clear();

	vector<int> stubs = tracks[tt]->getStubs();
	for(unsigned int sti=0;sti<stubs.size();sti++) tempVec.push_back( stubMap[ stubs[sti] ]);

	double pz = tracks[tt]->getCurve()/(tan(2*atan(exp(-tracks[tt]->getEta0()))));
	
	TTTrack< Ref_PixelDigi_ > tempTrack( tempVec );
	GlobalPoint POCA(0.,0.,tracks[tt]->getZ0());
	GlobalVector mom(tracks[tt]->getCurve()*cos(tracks[tt]->getPhi0()),
			 tracks[tt]->getCurve()*sin(tracks[tt]->getPhi0()),
			 pz);
	
	//std::cout << tracks[tt]->getCurve()<<" / "<<tracks[tt]->getZ0() << " / " <<tracks[tt]->getEta0()<<" / "<<tracks[tt]->getPhi0()<< std::endl;
	
	tempTrack.setSector( seedSector );
	tempTrack.setWedge( -1 );
	tempTrack.setMomentum( mom , 5);
	tempTrack.setPOCA( POCA , 5);
	//std::cout << tracks[tt]->getZ0() << " / " << POCA.z() << " / " << tempTrack.getPOCA().z() << std::endl;
	TTTracksForOutput->push_back( tempTrack );
	
	delete tracks[tt];
      }

      //      std::cout << "PT6" << std::endl;

      for(unsigned int tt=0;tt<tracksb.size();tt++)
      {	
	tempVec.clear();

	vector<int> stubs = tracksb[tt]->getStubs();
	for(unsigned int sti=0;sti<stubs.size();sti++) tempVec.push_back( stubMap[ stubs[sti] ]);

	double pz = tracksb[tt]->getCurve()/(tan(2*atan(exp(-tracksb[tt]->getEta0()))));
	
	TTTrack< Ref_PixelDigi_ > tempTrack( tempVec );
	GlobalPoint POCA(0.,0.,tracksb[tt]->getZ0());
	GlobalVector mom(tracksb[tt]->getCurve()*cos(tracksb[tt]->getPhi0()),
			 tracksb[tt]->getCurve()*sin(tracksb[tt]->getPhi0()),
			 pz);
	
	//std::cout << tracks[tt]->getCurve()<<" / "<<tracks[tt]->getZ0() << " / " <<tracks[tt]->getEta0()<<" / "<<tracks[tt]->getPhi0()<< std::endl;
	
	tempTrack.setSector( seedSector );
	tempTrack.setWedge( -1 );
	tempTrack.setMomentum( mom , 5);
	tempTrack.setPOCA( POCA , 5);
	//std::cout << tracks[tt]->getZ0() << " / " << POCA.z() << " / " << tempTrack.getPOCA().z() << std::endl;
	TTTracksBinForOutput->push_back( tempTrack );
	
	delete tracksb[tt];
      }

      //      std::cout << "PT7" << std::endl;

    } // End of loop over patterns
      
    //    std::cout << "PT8" << std::endl;

    delete(TCB);    
    delete(TCBb);  

    //PLACEHOLDER
    delete(PCAFIT);

  }

  /// Put in the event content
  iEvent.put( TTTracksForOutput, TTTrackOutputTag);
  iEvent.put( TTTracksBinForOutput, TTTrackBinaryOutputTag);
}

// DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(TrackFitTCProducer);

#endif

