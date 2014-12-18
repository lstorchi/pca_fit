#ifndef _PCAFITTER_INC_ 
#define _PCAFITTER_INC_

#define ENTDIM 8
#define NUMOFLAYER (ENTDIM-2)

#define PARAMDIM 5

#define PTIDX 0
#define PHIIDX 1
#define D0IDX 2
#define TETHAIDX 3
#define Z0IDX 4

#include <armadillo>
#include <vector>
#include <map>

// lstorchi: this should be moved into include as soon as the utility code should be
//           a class... This is not a proper class, just a way to collect some functions 

class pcafitter 
{
  public:
    pcafitter();
    ~pcafitter();

    static std::string get_version_string();
    
    static std::string paramidx_to_string (int i);

    // remove it
    static int numofline (const char * fname);

    int get_coordim () const;
    void set_coordim (int);
    
    void write_to_file (const char * fname, 
        const arma::mat & vec, int idx);
    
    void write_armmat (const char *, arma::mat & cmtx);
    
    void write_armvct (const char *, arma::rowvec & q);
    
    void read_armmat (const char *, arma::mat & cmtx);
    
    void read_armvct (const char *, arma::rowvec & q);
    
    void compute_parameters (const arma::mat & cmtx, 
      const arma::rowvec & q, 
      const arma::mat & coordslt, 
      double * oneoverptcmp, double * phicmp, 
      double * etacmp, double * z0cmp, 
      double * d0cmp);
    
    void reading_from_file (const char * filename, 
        arma::mat & paramin, arma::mat & coordin,
        arma::mat & layer, arma::mat & ladder, 
        arma::mat & module, 
        std::map<std::string, int> & subsectors, 
        std::map<std::string, int> & subladders,
        std::vector<std::string> & subsectorslist,
        std::vector<std::string> & subladderslist,
        int num_of_ent, bool usesegid);
    
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

  private:
    int coordim_;
};

#endif
