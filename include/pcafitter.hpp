#ifndef _PCAFITTER_INC_ 
#define _PCAFITTER_INC_

#include <armadillo>
#include <vector>
#include <map>

// lstorchi: this should be moved into include as soon as the utility code should be
//           a class... This is not a proper class, just a way to collect some functions 
//           as soon as have time please rewrite it and move main_computation here ... 

namespace pca
{
  class pcafitter 
  {
    public:
      pcafitter();
      ~pcafitter();
  
      static std::string get_version_string();
      
      bool set_paramidx (int i, const std::string &);
      std::string paramidx_to_string (int i);
  
      int get_coordim () const;
      void set_coordim (int);
  
      int get_paramdim () const;
      void set_paramdim (int);
         
      bool compute_parameters (
          const arma::mat & cmtx, 
          const arma::rowvec & q, 
          const arma::mat & amtx,
          const arma::mat & vmtx,
          const arma::rowvec & kvct, 
          const arma::mat & coord, 
#ifdef INTBITEWISE         
          int16_t ** paraptr,
#else
          double ** paraptr,
#endif
          int paramdim,
          arma::rowvec & chi2values);
           
      void select_bigger_sub (
          const std::map<std::string, int> & sublist, 
          bool verbose, int & maxnumber, std::string & slctsubsec );
      
      void extract_sub (const std::vector<std::string> & sublist, 
          const std::string & slctsubsec, 
          const arma::mat & paramin,
          const arma::mat & coordin, 
          arma::mat & param,
          arma::mat & coord);

      
      bool compute_pca_constants (
          const arma::mat & param, 
          const arma::mat & coord,
          arma::mat & cmtx, 
          arma::rowvec & q,
          arma::mat & vmtx, 
          arma::mat & amtx, 
          arma::rowvec & kivec,
          int verbositylevel);

      const std::string & get_errmsg () const;
      int get_errnum() const;
  
    private:
  
      int errnum_;
      int coordim_;
      int paramdim_;

      std::string errmsg_;

      std::map<int, std::string> paramname_;

      void set_errmsg (int , const std::string &);
      void reset_error ();
  };
};

#endif
