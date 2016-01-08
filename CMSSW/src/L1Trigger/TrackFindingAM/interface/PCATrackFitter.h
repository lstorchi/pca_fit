#ifndef _PCATRACKFITTER_H_
#define _PCATRACKFITTER_H_

#include "TrackFitter.h"
#include "Hit.h"

#include <vector>
#include <set>
#include <map>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>

/**
   \brief Implementation of a PCA track fitter
**/
class PCATrackFitter:public TrackFitter{

 private:

  unsigned int verboseLevel_; 
  unsigned int event_counter_;
  unsigned int road_id_;

  void initialize();

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

  void mergePatterns();
  void mergeTracks();
  void fit();
  void fit(std::vector<Hit*> hits);
  TrackFitter* clone();

  void setEventCounter(unsigned int event_counter)
  {
    event_counter_ = event_counter;
  };
  
  void setRoadID(unsigned int road_id)
  {
    road_id_ = road_id;
  };
  
  void setVerboseLevel(unsigned int verboseLevel)
  {
    verboseLevel_ = verboseLevel;
  }; 

};
#endif
