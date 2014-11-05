#ifndef _PCAFITTER_PRIVATE_INC_ 
#define _PCAFITTER_PRIVATE_INC_

#define ENTDIM 8
#define COORDIM (ENTDIM-2)
#define DIMPERCOORD 3

#define PARAMDIM 5
#define PTIDX 0
#define PHIIDX 1
#define D0IDX 2
#define ETAIDX 3
#define Z0IDX 4

#include <armadillo>
#include <vector>
#include <map>

namespace pcafitter
{
  std::string get_version_string();

  std::string paramidxtostring (int i);

  int numofline (const char * fname);
  
  void writetofile (const char * fname, 
      const arma::mat & vec, int idx);

  void writearmmat (const char *, arma::mat & cmtx);

  void writearmvct (const char *, arma::rowvec & q);

  void readarmmat (const char *, arma::mat & cmtx);

  void readarmvct (const char *, arma::rowvec & q);

  void computeparameters (const arma::mat & cmtx, 
    const arma::rowvec & q, 
    const arma::mat & coordslt, 
    double * ptcmp, double * phicmp, 
    double * etacmp, double * z0cmp, 
    double * d0cmp);
  
  void readingfromfile (const char * filename, 
      arma::mat & paramin, arma::mat & coordin,
      arma::mat & layer, arma::mat & ladder, 
      arma::mat & module, 
      std::map<std::string, int> & subsectors, 
      std::map<std::string, int> & subladders,
      std::vector<std::string> & subsectorslist,
      std::vector<std::string> & subladderslist,
      int num_of_ent);
  
  void select_bigger_sub (
      const std::map<std::string, int> & sublist, 
      bool verbose, int & maxnumber, std::string & slctsubsec );

  void extract_sub (const std::vector<std::string> & sublist, 
      const std::string & slctsubsec, 
      const arma::mat & paramin,
      const arma::mat & coordin, 
      arma::mat & param,
      arma::mat & coord);

  void compute_pca_constants (
      const arma::mat & param, 
      const arma::mat & coord,
      arma::mat & cmtx, 
      arma::rowvec & q);
}

#endif
