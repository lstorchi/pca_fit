#ifndef _PCAFITTER_H_
#define _PCAFITTER_H_

#include "TrackFitter.h"
#include <math.h>

#include <iomanip>
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


  double fenetre_b_z[3][6];
  double fenetre_b_phi[3][6];

  double fenetre_e_r[3][15];
  double fenetre_e_phi[3][15];

  std::vector< std::vector<int> > idxs_lay;
  std::vector< std::vector<int> > list_cand6;
  std::vector<float> list_scores6;
  std::vector< std::vector<int> > list_cand5;
  std::vector<float> list_scores5;
  
  unsigned int m_nLayer;

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
  void mergePatterns();
  void mergeTracks();
  void fit();
  void fit(vector<Hit*> hits);

  TrackFitter* clone();
};
#endif
