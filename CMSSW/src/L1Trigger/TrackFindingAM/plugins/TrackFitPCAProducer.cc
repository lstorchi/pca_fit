/*! \class   TrackFitPCAProducer
 *
 *  \author S Viret / L Storchi
 *  \date   2016, Mar 4
 *
 */


#ifndef TRACK_FITTER_PCA_H
#define TRACK_FITTER_PCA_H

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
  std::string                  TTTrackBinaryOutputTag;

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
  TTTrackOutputTag         = iConfig.getParameter< std::string >( "TTTrackName" );

  produces< std::vector< TTTrack< Ref_PixelDigi_ > > >( TTTrackOutputTag );
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

  /// Get the Stubs already stored away
  edm::Handle< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > > > TTStubHandle;
  edm::Handle< std::vector< TTTrack< Ref_PixelDigi_ > > > TTPatternHandle;

  iEvent.getByLabel( TTStubsInputTag, TTStubHandle );
  iEvent.getByLabel( TTPatternsInputTag, TTPatternHandle );

  /// STEP 0
  /// Prepare output
  TTTracksForOutput->clear();

  int layer  = 0;
  int ladder = 0;
  int module = 0;

  int nbLayers = 0;

  /// Loop over TCs
  unsigned int tkCnt = 0;
  unsigned int j     = 0;

  std::map< unsigned int , edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > > stubMap;
  

  PCATrackFitter* PCA  = new PCATrackFitter(nbLayers); // Floating point

  /// STEP 0
  /// Read PCAConst file
  PCA->read_float_const_filename ("../data/infnpca_const_files/barrel_tow16_pca_const.txt");
  PCA->read_float_const_filename ("../data/infnpca_const_files/barrel_tow17_pca_const.txt");
  PCA->read_float_const_filename ("../data/infnpca_const_files/barrel_tow18_pca_const.txt");
  PCA->read_float_const_filename ("../data/infnpca_const_files/barrel_tow19_pca_const.txt");
  PCA->read_float_const_filename ("../data/infnpca_const_files/barrel_tow20_pca_const.txt");
  PCA->read_float_const_filename ("../data/infnpca_const_files/barrel_tow21_pca_const.txt");
  PCA->read_float_const_filename ("../data/infnpca_const_files/barrel_tow22_pca_const.txt");
  PCA->read_float_const_filename ("../data/infnpca_const_files/barrel_tow23_pca_const.txt");
  PCA->read_float_const_filename ("../data/infnpca_const_files/barrel_tow24_pca_const.txt");
  PCA->read_float_const_filename ("../data/infnpca_const_files/barrel_tow25_pca_const.txt");
  PCA->read_float_const_filename ("../data/infnpca_const_files/barrel_tow26_pca_const.txt");
  PCA->read_float_const_filename ("../data/infnpca_const_files/barrel_tow27_pca_const.txt");
  PCA->read_float_const_filename ("../data/infnpca_const_files/barrel_tow28_pca_const.txt");
  PCA->read_float_const_filename ("../data/infnpca_const_files/barrel_tow29_pca_const.txt");
  PCA->read_float_const_filename ("../data/infnpca_const_files/barrel_tow30_pca_const.txt");
  PCA->read_float_const_filename ("../data/infnpca_const_files/barrel_tow31_pca_const.txt");

  /// STEP 1
  /// Loop over track candidates

  //  std::cout << "Start the loop over " << TTPatternHandle->size() << " pattern(s) in order to recover the stubs" << std::endl;

  std::vector<Hit*> m_hits;
  for(unsigned int i=0;i<m_hits.size();i++) delete m_hits[i];
  m_hits.clear();

  std::vector<Track*> tracks;
  for(unsigned int i=0;i<tracks.size();i++) delete tracks[i];
  tracks.clear();

  edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >::const_iterator inputIter;
  edmNew::DetSet< TTStub< Ref_PixelDigi_ > >::const_iterator stubIter;

  std::vector< TTTrack< Ref_PixelDigi_ > >::const_iterator iterTTTrack;
  
  /// Go on only if there are TCs from PixelDigis
  if ( TTPatternHandle->size() > 0 )
  {
    for ( iterTTTrack = TTPatternHandle->begin();
	  iterTTTrack != TTPatternHandle->end();
	  ++iterTTTrack )
    {
      edm::Ptr< TTTrack< Ref_PixelDigi_ > > tempTrackPtr( TTPatternHandle, tkCnt++ );

      j = 0;



      m_hits.clear();
      tracks.clear();
      stubMap.clear();

      /// Get everything relevant
      unsigned int seedSector = tempTrackPtr->getSector();

      Track* TC = new Track();
      // TODO getWedge is unsigned ???
      if (tempTrackPtr->getWedge() == 1)
        TC->setCharge(-1.0);
      else
        TC->setCharge(+1.0);
      //std::cout << "Wedge: " << tempTrackPtr->getWedge() << std::endl;
      TC->setCurve(tempTrackPtr->getMomentum(5).perp());
      TC->setEta0(tempTrackPtr->getMomentum(5).eta());
      TC->setPhi0(tempTrackPtr->getMomentum(5).phi());
      TC->setZ0(tempTrackPtr->getPOCA(5).z());


      // Get the stubs in the TC

      std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_  > >, TTStub< Ref_PixelDigi_  > > > trackStubs = tempTrackPtr->getStubRefs();

      // Loop over stubs contained in the pattern to recover the info

      for(unsigned int i=0;i<trackStubs.size();i++)
      {
	++j;

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
	
	Hit* h = new Hit(layer,ladder, module, segment, strip, 
			 j, tp, spt, ip, eta, phi0, x, y, z, x0, y0, z0, tempStubRef->getTriggerDisplacement()-tempStubRef->getTriggerOffset());
	m_hits.push_back(h);

      } /// End of loop over track stubs

      //if (tempTrackPtr->getSector()!=18) continue; // For the moment don't need to bother with the rest
      if (tempTrackPtr->getSector() < 16 || tempTrackPtr->getSector() > 31) continue; // For the moment don't need to bother with the rest

      cout<<"Dealing with TC with "<< j << " stubs in tower " << tempTrackPtr->getSector() << endl;

      PCA->cleanChi2();
      PCA->setSectorID(tempTrackPtr->getSector());
      PCA->setTrack(TC);
      PCA->fit(m_hits);
      tracks = PCA->getTracks();
      PCA->clean();

      delete TC;      

      std::vector<double> chi2v = PCA->get_chi2f();

      if (chi2v.size() == tracks.size())
      {
        // Store the tracks (no duplicate cleaning yet)
	cout<<"Found "<<tracks.size()<<" track"<<endl;
        
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
          
          tempTrack.setSector( seedSector );
          tempTrack.setWedge( tracks[tt]->getCharge() );
          tempTrack.setMomentum( mom , 5);
          tempTrack.setPOCA( POCA , 5);
          tempTrack.setChi2(chi2v[tt]);
        
          TTTracksForOutput->push_back( tempTrack );
          
          delete tracks[tt];
        }
      }
    } // End of loop over patterns

    delete(PCA);    

  }

  /// Put in the event content
  iEvent.put( TTTracksForOutput, TTTrackOutputTag);
}

// DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(TrackFitPCAProducer);

#endif

