#ifndef _PCAFITTER_H_
#define _PCAFITTER_H_

#include "TrackFitter.h"
#include <math.h>

#include <iomanip>
#include <map>
#include <set>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>


#include "../interface/pcaconst.hpp"

/**
   \brief PCA fitter hfile : Loriano Storchi 2016
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

  Track *track_;

  std::string cfname_;
  bool useinteger_ ;

  std::vector<pca::matrixpcaconst<double> > pcacontvct_float_;
  std::vector<pca::matrixpcaconst<int32_t> > pcacontvct_integer_;

  /* we will merge these using TypeIs... structs */
  void fit_float(vector<Hit*> hits);
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

  void fit();
  void initialize();
  void fit(vector<Hit*> hits);

  void setTrack(Track* intc); 

  void mergePatterns();
  void mergeTracks();

  /* In case a second file is read values will be appended,
   * initialize will remove all the values */
  void read_float_const_filename (const std::string & in);
  void read_integegr_const_filename (const std::string & in);

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
