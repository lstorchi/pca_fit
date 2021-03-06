#ifndef _PCAFFUNCTYPE_INC_ 
#define _PCAFFUNCTYPE_INC_

#include <armadillo>

#include <string>
#include <vector>
#include <map>
#include <set>

#include "pcaconst.hpp"

// to be moved
#define PCA_COTTHETAIDX 0
#define PCA_Z0IDX 1

#define PCA_ONEOVERPTIDX 0
#define PCA_PHIIDX 1

#define ISBARREL 0
#define ISHYBRID 1
#define ISENDCAP 2

namespace pca
{
  bool is_a_valid_layers_seq(const std::string &, 
      int, const int, const bool);

  bool validate_barrel_sequence_5 (
      const std::string & sequence);
  
  bool check_to_read (bool, bool, int);
  
  bool check_charge_sign (int, std::set<int> &);
 
  std::string get_paramname_from_id (int);

  int numofline (const char * fname);

  void tokenize (const std::string & str, std::vector<std::string> & tokens,
      const std::string & delimiters = " ");

  void write_to_file (const char * fname, 
      const arma::mat & vec, int idx);
  
  void write_armmat (const char *, arma::mat & cmtx);
  
  void write_armvct (const char *, arma::rowvec & q);
  
  void read_armmat (const char *, arma::mat & cmtx);
  
  void read_armvct (const char *, arma::rowvec & q);

  int get_num_of_ent(const char *);

  void read_armmat_ibw (const char * fname, arma::mat & cmtx);

  void read_armvct_ibw (const char * fname, arma::rowvec & q);
  
  bool import_pca_const (
    std::vector<pca::matrixpcaconst<int32_t> > & vct,
    pca::matrixpcaconst<int32_t> & cmtx_rz, 
    pca::matrixpcaconst<int32_t> & qvec_rz, 
    pca::matrixpcaconst<int32_t> & amtx_rz, 
    pca::matrixpcaconst<int32_t> & kvec_rz, 
    pca::matrixpcaconst<int32_t> & cmtx_rphi, 
    pca::matrixpcaconst<int32_t> & qvec_rphi, 
    pca::matrixpcaconst<int32_t> & amtx_rphi, 
    pca::matrixpcaconst<int32_t> & kvec_rphi, 
    double eta, double pt, 
    int chargesignin,
    const std::string & layersid,
    const std::string & pslayersid,
    int towerid);

  bool import_pca_const (
    std::vector<pca::matrixpcaconst<double> > & vct,
    pca::matrixpcaconst<double> & cmtx_rz, 
    pca::matrixpcaconst<double> & qvec_rz, 
    pca::matrixpcaconst<double> & amtx_rz, 
    pca::matrixpcaconst<double> & kvec_rz, 
    pca::matrixpcaconst<double> & cmtx_rphi, 
    pca::matrixpcaconst<double> & qvec_rphi, 
    pca::matrixpcaconst<double> & amtx_rphi, 
    pca::matrixpcaconst<double> & kvec_rphi, 
    double eta, double pt, 
    int chargesignin,
    const std::string & layersid,
    const std::string & pslayersid,
    int towerid);
 
#if 0
  bool reading_from_file_split (const pca::pcafitter &, const char * filename, 
          arma::mat & paramin, arma::mat & coordin,
          bool useonlyeven, bool useonlyodd, 
          bool rzread, bool rphiread, double etamin,
          double etamax, double ptmin, double ptmax, 
          bool chargeoverpt, int chargesign,
          bool excludesmodule, bool usealsod0,
          bool usex0y0, int singleparam, 
          double phimin, double phimax,
          double z0min, double z0max,
          double d0min, double d0max, 
          bool usealsox0, 
          bool verbose, 
          arma::vec & ptvals,
          bool checklayersid,
          int maxnumoflayers);
#endif

  void global_to_relative (arma::mat & coordin, 
      double coord1min = std::numeric_limits<double>::infinity(), 
      double coord2min = std::numeric_limits<double>::infinity() );

  bool file_exists(const std::string& filename);

  double delta_phi(double phi1, double phi2);

  double cot (double);

  double atanall (double);

  template<typename T> bool armamat_to_pcamat (const arma::mat & in, 
      pca::matrixpcaconst<T> & out )
  {
    out.reset(in.n_rows, in.n_cols);

    for (unsigned int i=0; i<in.n_rows; ++i)
      for (unsigned int j=0; j<in.n_cols; ++j)
        out.element(i, j) = in(i, j);

    return true;
  } 

  template<typename T> bool pcamat_to_armamat (const pca::matrixpcaconst<T> & in, 
      arma::mat & out) 
  {
    out.set_size(in.n_rows(), in.n_cols());

    for (unsigned int i=0; i<in.n_rows(); ++i)
      for (unsigned int j=0; j<in.n_cols(); ++j)
        out(i, j) = in(i, j);

    return true;
  } 

  template<typename T> bool pcamat_to_armarowvec (const pca::matrixpcaconst<T> & in, 
      arma::rowvec & out) 
  {
    unsigned int dim = in.n_rows();
    if (dim == 1)
      dim= in.n_cols();

    out.set_size(dim);

    unsigned int idx = 0;
    for (unsigned int i=0; i<in.n_rows(); ++i)
    {
      for (unsigned int j=0; j<in.n_cols(); ++j)
      {
        out(idx) = in(i, j);
        ++idx;
      }
    }

    return true;
  } 

  void dump_element (const matrixpcaconst<double> & in, 
    std::ostream & out);
};
  
#endif
