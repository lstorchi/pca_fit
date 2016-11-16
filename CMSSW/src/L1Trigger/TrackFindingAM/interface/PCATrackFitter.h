#ifndef _PCAFITTER_H_
#define _PCAFITTER_H_

#include "TrackFitter.h"
#include <math.h>

#include <iomanip>
#include <map>
#include <set>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>


#include "../interface/pcaconst.h"

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

  std::string cfname_ ;
  bool useboundaries_ ;
  bool useinteger_ ;
  bool storecpt_;

  std::vector<pca::matrixpcaconst<double> > pcacontvct_float_;
  std::vector<pca::matrixpcaconst<long long int> > pcacontvct_integer_;

  std::vector<double> chi2v_;

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
  void fit(vector<Hit*> hits, int pattern_id=-1);

  void setTrack(Track* intc); 

  const std::vector<double> & get_chi2() const
  {
    return chi2v_;
  }

  void cleanChi2();
  void mergePatterns();
  void mergeTracks();

  /* In case a second file is read values will be appended,
   * initialize will remove all the values */
  void read_float_const_filename (const std::string & in);
  void read_integer_const_filename (const std::string & in);

  void set_useinteger (const bool & in)
  {
    useinteger_ = in;
  }
  bool get_useinteger () const
  {
    return useinteger_;
  }

  void set_useboundaries (const bool & in)
  {
    useboundaries_ = in;
  }
  bool get_useboundaries () const
  {
    return useboundaries_;
  }

  void set_store_chargept (const bool & in)
  {
    storecpt_ = in;
  }
  bool get_store_chargept () const
  {
    return storecpt_;
  }

  TrackFitter* clone();
};
#endif
