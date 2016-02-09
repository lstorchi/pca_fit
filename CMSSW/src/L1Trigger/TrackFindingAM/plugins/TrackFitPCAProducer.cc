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

#include "L1Trigger/TrackFindingAM/interface/TCBuilder.h"
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
  //edm::InputTag                TTPatternsInputTag;
  edm::InputTag                TCTTrackInputTag;
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
  //TTPatternsInputTag       = iConfig.getParameter< edm::InputTag >( "TTInputPatterns" );
  TTTrackOutputTag         = iConfig.getParameter< std::string >( "TTTrackName" );       // Container of C++ TCs
  //TTTrackBinaryOutputTag   = iConfig.getParameter< std::string >( "TTTrackBinaryName" ); // Container of bit-wise TCs
  TCTTrackInputTag = iConfig.getParameter< edm::InputTag >( "TCBuilderTTracks" );

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
  std::cout << "PCA Fitter begin run !!!!" << std::endl;

}

/// End run
void TrackFitPCAProducer::endRun( const edm::Run& run, const edm::EventSetup& iSetup ) {}

/// Implement the producer
void TrackFitPCAProducer::produce( edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  /// Prepare output
  /// The temporary collection is used to store tracks
  /// before removal of duplicates
  std::auto_ptr< std::vector< TTTrack< Ref_PixelDigi_ > > > 
    TTTracksForOutput( new std::vector< TTTrack< Ref_PixelDigi_ > > );

  /// Get the Stubs already stored away
  edm::Handle< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > > > TTStubHandle;
  //edm::Handle< std::vector< TTTrack< Ref_PixelDigi_ > > > TTPatternHandle;
  edm::Handle< std::vector< TTTrack< Ref_PixelDigi_ > > > TCBuilderTTracksHandle;

  iEvent.getByLabel( TTStubsInputTag, TTStubHandle );       // Get the stubs
  //iEvent.getByLabel( TTPatternsInputTag, TTPatternHandle ); // Get the matched roads
  iEvent.getByLabel( TCTTrackInputTag, TCBuilderTTracksHandle); // Get the matched roads

  std::cout << "Inside PCA produce" << std::endl;
  /// STEP 0
  /// Prepare output
  TTTracksForOutput->clear();

  int layer  = 0;
  int ladder = 0;
  int module = 0;
  int nbLayers = 0;

  unsigned int j     = 0;
  //unsigned int tkCnt = 0;
  PCATrackFitter* pcafitter = new PCATrackFitter(nbLayers);

  std::vector<Hit*> m_hits;
  std::vector<Track*> pca_tracks, tcb_tracks;
  std::map< unsigned int , edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > > stubMap;

  //Loop over TTracks from TCBuilder 
  std::vector< TTTrack< Ref_PixelDigi_ > >::const_iterator iterTCTTrack;
  iterTCTTrack = (*TCBuilderTTracksHandle.product()).begin();
  for(;iterTCTTrack != (*TCBuilderTTracksHandle.product()).end(); ++iterTCTTrack ) 
  {
    const TTTrack< Ref_PixelDigi_ > trk = *iterTCTTrack;
    Track* tempt = new Track();
    tempt->setCurve(iterTCTTrack->getMomentum(5).perp());
    tempt->setEta0(iterTCTTrack->getMomentum(5).eta());
    tempt->setPhi0(iterTCTTrack->getMomentum(5).phi());
    tempt->setZ0(iterTCTTrack->getPOCA(5).z());
    tempt->setCharge(-(iterTCTTrack->getRInv(5)));//charge = -RInv
    std::cout << "Charge from TC track=" << tempt->getCharge() << std::endl;
    tcb_tracks.push_back(tempt);
    nbLayers = iterTCTTrack->getWedge();
    unsigned int seedSector = iterTCTTrack->getSector();
    m_hits.clear();
    tcb_tracks.clear();    
    stubMap.clear();
    j = 0;
    std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_  > >, 
      TTStub< Ref_PixelDigi_  > > > trackStubs = iterTCTTrack->getStubRefs();
    for(unsigned int i=0; i<trackStubs.size(); ++i)
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

    } /// End of loop over tcb_track stubs
    
    pcafitter->setSectorID(seedSector);
    pcafitter->setTracks(tcb_tracks, stubMap, seedSector);
    pcafitter->fit(m_hits);
    pca_tracks = pcafitter->getTracks();

    std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, 
      TTStub< Ref_PixelDigi_ > > > tempVec;
    for(unsigned int tt=0; tt<pca_tracks.size(); tt++)
    {
      tempVec.clear();
      //Stubs used for the fit 
      vector<int> stubs = pca_tracks[tt]->getStubs();
      for(unsigned int sti=0;sti<stubs.size();sti++) 
        tempVec.push_back( stubMap[ stubs[sti] ] );
      TTTrack< Ref_PixelDigi_ > tempTrack( tempVec );
      double pz = pca_tracks[tt]->getCurve()/(tan(2*atan(exp(-pca_tracks[tt]->getEta0()))));
      GlobalPoint POCA(0.,0.,pca_tracks[tt]->getZ0());
      GlobalVector mom(pca_tracks[tt]->getCurve()*cos(pca_tracks[tt]->getPhi0()),
			 pca_tracks[tt]->getCurve()*sin(pca_tracks[tt]->getPhi0()),
			 pz);
      tempTrack.setSector( seedSector );
      tempTrack.setWedge( -1 );
      tempTrack.setMomentum( mom , 5);
      tempTrack.setPOCA( POCA , 5);
  
      TTTracksForOutput->push_back( tempTrack );
	
      delete pca_tracks[tt];
    }
  }    
  iEvent.put( TTTracksForOutput, TTTrackOutputTag);
}

// DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(TrackFitPCAProducer);

#endif
