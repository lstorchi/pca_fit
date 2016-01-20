#ifndef _PCAFFUNCTYPE_INC_ 
#define _PCAFFUNCTYPE_INC_

#include <armadillo>

#include <string>
#include <vector>
#include <map>
#include <set>

// to be moved
#define PCA_COTTHETAIDX 0
#define PCA_Z0IDX 1

#define PCA_ONEOVERPTIDX 0
#define PCA_PHIIDX 1

namespace pca
{
  bool is_a_valid_layers_seq(const std::string &, 
      int, const bool, const bool);

  bool validate_barrel_sequence_5 (bool excludesmodule, 
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

};
  
#endif
