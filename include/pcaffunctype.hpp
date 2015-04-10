#ifndef _PCAFFUNCTYPE_INC_ 
#define _PCAFFUNCTYPE_INC_

#include <armadillo>
#include <vector>
#include <map>

namespace pca
{
  int numofline (const char * fname);

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
};
  
#endif