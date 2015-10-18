#ifndef _PCAFFUNCTYPE_INC_ 
#define _PCAFFUNCTYPE_INC_

#include <armadillo>

#include <string>
#include <vector>
#include <map>

// to be moved
#define ONEOVERPTIDX 0
#define PHIIDX 1
#define D0IDX 2
#define COTTHETAIDX 3
#define Z0IDX 4

#define SPLIT_COTTHETAIDX 0
#define SPLIT_Z0IDX 1

#define SPLIT_ONEOVERPTIDX 0
#define SPLIT_PHIIDX 1

#define SPLIT_X0IDX 0
#define SPLIT_Y0IDX 1

#define SPLIT_D0IDX 2
#define SPLIT_X0IDX_NS 2

namespace pca
{
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

  void global_to_relative (arma::mat & coordin);

  bool file_exists(const std::string& filename);

  double delta_phi(double phi1, double phi2);

  double cot (double);

  double atanall (double);

};
  
#endif
