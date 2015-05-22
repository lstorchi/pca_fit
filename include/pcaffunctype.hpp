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
#define COTTETHAIDX 3
#define Z0IDX 4

#define SPLIT_COTTETHAIDX 0
#define SPLIT_Z0IDX 1

#define SPLIT_ONEOVERPTIDX 0
#define SPLIT_PHIIDX 1

#define SPLIT_X0IDX 0
#define SPLIT_Y0IDX 1

#define SPLIT_D0IDX 2

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

  void reading_from_file (const char * filename, 
          arma::mat & paramin, arma::mat & coordin,
          arma::mat & layer, arma::mat & ladder, 
          arma::mat & module, 
          std::map<std::string, int> & subsectors, 
          std::map<std::string, int> & subladders,
          std::vector<std::string> & subsectorslist,
          std::vector<std::string> & subladderslist,
          int num_of_ent, bool usesegid, 
          bool useonlyeven, bool useonlyodd);

  bool reading_from_file_split (const pca::pcafitter &, const char * filename, 
          arma::mat & paramin, arma::mat & coordin,
          int num_of_ent, bool useonlyeven, bool useonlyodd, 
          bool rzread, bool rphiread, double etamin,
          double etamax, double ptmin, double ptmax, 
          bool chargeoverpt, int chargesign,
          bool excludesmodule, bool usealsod0,
          bool usex0y0, int singleparam, 
          double phimin, double phimax,
          double z0min, double z0max,
          double d0min, double d0max);

  bool file_exists(const std::string& filename);

  double delta_phi(double phi1, double phi2);

  double cot (double);

  double atanall (double);

};
  
#endif
