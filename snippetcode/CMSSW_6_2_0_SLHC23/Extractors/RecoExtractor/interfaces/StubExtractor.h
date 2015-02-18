#ifndef STUBEXTRACTOR_H
#define STUBEXTRACTOR_H

/**
 * StubExtractor
 * \brief: Base class for extracting pixel info
 */


//Include RECO inf
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerDetUnit.h"

#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"

#include "DataFormats/SiPixelDetId/interface/StackedTrackerDetId.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"
#include "FWCore/Utilities/interface/InputTag.h"

//Include std C++
#include <iostream>
#include <vector>
#include <boost/shared_ptr.hpp>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "MCExtractor.h"

class StubExtractor
{

 public:

  StubExtractor(std::string CLUS_tag,std::string CLUS_name,std::string STUB_tag,std::string STUB_name, bool doTree);
  StubExtractor(TFile *a_file);
  ~StubExtractor();


  void init(const edm::EventSetup *setup);
  void writeInfo(const edm::Event *event, MCExtractor *mc); 
  void getInfo(int ievt); 

  void reset();
  void fillTree(); 
  void fillSize(int size);
  int  getSize();
  int  n_events() {return m_n_events;}

  bool isOK() {return m_OK;}

  //Getters

  int getClust1Idx(float x, float y, float z);
  int getClust2Idx(int idx1, float dist);
  int get_id(int lay,int lad,int mod,float x,float y,float z);

  int getNDigis() {return m_clus;}

  float getStub_ptGen(int i) {return sqrt(m_stub_pxGEN->at(i)*m_stub_pxGEN->at(i)
					 +m_stub_pyGEN->at(i)*m_stub_pyGEN->at(i));}
  float getStub_etaGen(int i) {return m_stub_etaGEN->at(i);}
  float getStub_strip(int i) {return m_stub_strip->at(i);}
  float getStub_x(int i) {return m_stub_x->at(i);}
  float getStub_y(int i) {return m_stub_y->at(i);}
  float getStub_z(int i) {return m_stub_z->at(i);}
  float getStub_bend(int i) {return m_stub_deltas->at(i);}

  float getStub_X0(int i) {return m_stub_X0->at(i);}
  float getStub_Y0(int i) {return m_stub_Y0->at(i);}
  float getStub_Z0(int i) {return m_stub_Z0->at(i);}
  float getStub_PHI0(int i) {return m_stub_PHI0->at(i);}



  int getStub_ID(int i) {return 1000000*m_stub_layer->at(i)+10000*m_stub_ladder->at(i)
      +100*m_stub_module->at(i)+m_stub_seg->at(i);}
  int getStub_tp(int i) {return m_stub_tp->at(i);}
  int getNStub() {return m_stub;}

 private:
  
  std::string m_CLUS_tag;
  std::string m_STUB_tag;

  std::string m_CLUS_name;
  std::string m_STUB_name;


  int n_tot_evt;
  int m_nstubs;

  bool  m_verb;
  bool  m_matchStubs;
  bool  m_posMatch;
  int   m_max_wclus;
  int   m_max_dist;
  float m_max_StubDev;
  int   m_PDG_id;
  float m_window_size;
  float m_thresh;

  TTree* m_tree;
  edm::Handle< edm::SimTrackContainer >  SimTrackHandle;
  edm::Handle< edm::SimVertexContainer > SimVtxHandle;

  edm::Handle<TrackingParticleCollection>  TPCollection ;
  edm::Handle< edm::DetSetVector<PixelDigiSimLink> > pDigiLinkColl;
  edm::DetSet<PixelDigiSimLink> pDigiLinks;

  //  edm::ESHandle<TrackerGeometry> theTrackerGeometry;
  edm::InputTag m_tag;
  bool m_OK;
  bool m_matching;
  int  m_n_events;
  double mMagneticFieldStrength ;

  /// Geometry handles etc
  edm::ESHandle< TrackerGeometry >                theTrackerGeometry;
  edm::ESHandle< StackedTrackerGeometry >         theStackedTrackerGeometry;
  const StackedTrackerGeometry*                   theStackedGeometry;
  StackedTrackerGeometry::StackContainerIterator  StackedTrackerIterator;

  edm::Handle< edmNew::DetSetVector< TTCluster< Ref_PixelDigi_ > > > PixelDigiL1TkClusterHandle;
  edm::Handle< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > > > PixelDigiL1TkStubHandle;
  edm::Handle< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > > > PixelDigiL1TkFailedStubHandle;

  edm::Handle< TTClusterAssociationMap< Ref_PixelDigi_ > > MCTruthTTClusterHandle;
  edm::Handle< TTStubAssociationMap< Ref_PixelDigi_ > >    MCTruthTTStubHandle;

  edm::Handle< std::vector< TrackingParticle > > TrackingParticleHandle;
  edm::Handle< std::vector< TrackingVertex > > TrackingVertexHandle;

  std::vector<int>                 *m_clus_used;   


  int m_clus;

  // Size of the following vectors is m_clus
  std::vector<float>  *m_clus_x;      // x pos of cluster i (in cm)
  std::vector<float>  *m_clus_y;      // y pos of cluster i (in cm)
  std::vector<float>  *m_clus_z;      // z pos of cluster i (in cm)
  std::vector<float>  *m_clus_e;      // UNUSED
  std::vector<int>    *m_clus_layer;  // layer of cluster i
  std::vector<int>    *m_clus_module; // module of cluster i
  std::vector<int>    *m_clus_ladder; // ladder/ring of cluster i
  std::vector<int>    *m_clus_seg;    // segment (column) of cluster i
  std::vector<float>  *m_clus_strip;  // strip position of cluster i (in half-strip unit) 
  std::vector<int>    *m_clus_sat;    // number of saturating strips in cluster i (digis>255 ADC)
  std::vector<int>    *m_clus_nstrips;// number of strips in cluster i
  std::vector<int>    *m_clus_matched;// number of simhits matched to the cluster
  std::vector<int>    *m_clus_PS;     // number of segments in the corresponding module
  std::vector<int>    *m_clus_nrows;  // number of strips in the corresponding module
  std::vector<int>    *m_clus_pdgID;  // list of tracking particles inducing cluster i
  std::vector<float>  *m_clus_ptGEN;  // list of simhits inducing cluster i
  std::vector<int>    *m_clus_pid;    // process id inducing cluster i (see MCExtractor.h)
  std::vector<int>    *m_clus_rank_PIX;   // Number of pixel clusters in the chip (no size limit)  
  std::vector<int>    *m_clus_rank_STR;   // Number of strip clusters in the chip (no size limit)  
  std::vector<int>    *m_clus_rank_PIX_4; // Number of pixel clusters in the chip (4 strips limit)    
  std::vector<int>    *m_clus_rank_STR_4; // Number of strip clusters in the chip (4 strips limit)    
  std::vector<int>    *m_clus_rank_PIX_8; // Number of pixel clusters in the chip (8 strips limit)      
  std::vector<int>    *m_clus_rank_STR_8; // Number of strip clusters in the chip (8 strips limit)    
  std::vector<int>    *m_clus_tp;  
  std::vector<std::vector<int> >    *m_clus_pix;   // Pixel digis linked to the cluster 
  std::vector<std::vector<int> >    *m_clus_mult;  // Overall chip cluster occupancies for the event

  int m_stub;

  // Size of the following vectors is m_stub
  std::vector<float>  *m_stub_pt;     // pt calculated for stub i (in GeV/c)
  std::vector<float>  *m_stub_ptMC;   // pt calculated using MC coordinates for stub i (in GeV/c)
  std::vector<float>  *m_stub_pxGEN;  // px generated of stub i (in GeV/c)
  std::vector<float>  *m_stub_pyGEN;  // py generated of stub i (in GeV/c)
  std::vector<float>  *m_stub_etaGEN; // eta generated of stub i (in GeV/c)
  std::vector<float>  *m_stub_X0;     // x origin of particle inducing stub i (in cm)
  std::vector<float>  *m_stub_Y0;     // y origin of particle inducing stub i (in cm)
  std::vector<float>  *m_stub_Z0;     // z origin of particle inducing stub i (in cm)
  std::vector<float>  *m_stub_PHI0;   // phi origin of particle inducing stub i (in rad)
  std::vector<int>    *m_stub_layer;  // layer of stub i
  std::vector<int>    *m_stub_module; // module of stub i
  std::vector<int>    *m_stub_ladder; // ladder/ring of stub i
  std::vector<int>    *m_stub_seg;    // segment of stub i
  std::vector<int>    *m_stub_chip;   // chip of stub i
  std::vector<float>  *m_stub_strip;  // strip of stub i (innermost module value)
  std::vector<float>  *m_stub_x;      // x pos of stub i (in cm)
  std::vector<float>  *m_stub_y;      // y pos of stub i (in cm)
  std::vector<float>  *m_stub_z;      // z pos of stub i (in cm)
  std::vector<int>    *m_stub_clust1; // bottom cluster index of stub i 
  std::vector<int>    *m_stub_clust2; // top cluster index of stub i 
  std::vector<int>    *m_stub_cw1;    // bottom cluster width (in strips) 
  std::vector<int>    *m_stub_cw2;    // top cluster width (in strips) 
  std::vector<float>  *m_stub_deltas; // corrected stub width (in half-strips)
  std::vector<float>  *m_stub_cor;    // stub wirdth applied correction (in half-strips)
  std::vector<int>    *m_stub_tp;     // index of the TP inducing the stub in the MC tree
  std::vector<int>    *m_stub_pdg;    // PDG code of the particle inducing the stub
  std::vector<int>    *m_stub_pid;    // process id inducing cluster i (see MCExtractor.h)
  std::vector<int>    *m_stub_rank;   // Rank of the stub at the FE output (in bend order, then strip)
};

#endif 
