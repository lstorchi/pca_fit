#ifndef _PCAFITTER_H_
#define _PCAFITTER_H_

#include "TrackFitter.h"

#include "FWCore/Framework/interface/EDProducer.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"

#include "L1Trigger/TrackFindingAM/interface/CMSPatternLayer.h"
#include "L1Trigger/TrackFindingAM/interface/PatternFinder.h"
#include "L1Trigger/TrackFindingAM/interface/SectorTree.h"
#include "L1Trigger/TrackFindingAM/interface/Hit.h"


#include <math.h>

#include <iomanip>
#include <map>
#include <set>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>

/**
   \brief Implementation of PCA fitter
**/
class PCATrackFitter:public TrackFitter{

 private:

  friend class boost::serialization::access;
  
  template<class Archive> void save(Archive & ar, const unsigned int version) const
  {
    ar << boost::serialization::base_object<TrackFitter>(*this);
    ar << sec_phi;
  }
  
  template<class Archive> void load(Archive & ar, const unsigned int version)
  {
    ar >> boost::serialization::base_object<TrackFitter>(*this);
    ar >> sec_phi;
  }
  
  BOOST_SERIALIZATION_SPLIT_MEMBER()

  std::vector<Track*> tracks_;
  std::map< unsigned int , edm::Ref< edmNew::DetSetVector<
            TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > > stubmap_;
  unsigned int seedsector_;

 public:

  /**
     \brief Default constructor   
  **/
  PCATrackFitter();
  /**
     \brief Constructor   
     \param nb Layers number
  **/  
  PCATrackFitter(int nb);
  ~PCATrackFitter();

  void initialize();
  void fit(vector<Hit*> hits);
  void fit();
  void setTracks(const std::vector<Track*> &, const 
      std::map< unsigned int , edm::Ref< edmNew::DetSetVector<
         TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > > &,
         const unsigned int &);

  void mergePatterns();
  void mergeTracks();

  TrackFitter* clone();
};
#endif
