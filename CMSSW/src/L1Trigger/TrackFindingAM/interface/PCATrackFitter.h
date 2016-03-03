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
   \brief PCA fitter
**/
class PCATrackFitter:public TrackFitter
{

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
  std::string cfname_;
  bool useinteger_ ;

  /* we will merge these using TypeIs... structs */
  void fit_double(vector<Hit*> hits);
  void fit_integer(vector<Hit*> hits);

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
  void setTracks(const std::vector<Track*> &); 

  void mergePatterns();
  void mergeTracks();

  void set_const_filename(const std::string & in)
  {
    cfname_ = in;
  }
  const std::string & get_const_filename () const
  {
    return cfname_;
  }

  void set_useinteger (const bool & in)
  {
    useinteger_ = in;
  }
  bool get_useinteger () const
  {
    return useinteger_;
  }

  TrackFitter* clone();
};
#endif
