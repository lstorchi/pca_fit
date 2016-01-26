#ifndef _PCAROOTFILEREADER_INC_ 
#define _PCAROOTFILEREADER_INC_

#include <armadillo>
#include <vector>
#include <map>

// lstorchi: this should be moved into include as soon as the utility code should be
//           a class... This is not a proper class, just a way to collect some functions 
//           as soon as have time please rewrite it and move main_computation here ... 

namespace pca
{
  class rootfilereader
  {
    private:
      // lstorchi: also with INTBITEWISE I would suggest to use 
      //           double here and convert afterwards 
      typedef struct 
      {
        int dim;
        std::vector<double> x, y, z;
        std::vector<int> layer, ladder, module, segid;

        double x0, y0, z0, d0, pt, phi, eta;
        int pdg;
      } track_str;

    public:
      rootfilereader();
      ~rootfilereader();
  
      void reset();

      std::string get_filename();
      void set_filename(const std::string &);
      
      const std::string & get_errmsg () const;
      int get_errnum() const;

      void set_rzplane (bool);
      bool get_rzplane () const;
      void set_rphiplane (bool);
      bool get_rphiplane () const;

      void set_savecheckfiles (bool in);
      bool get_savecheckfiles () const;

      void set_chargeoverpt (bool);
      bool get_chargeoverpt () const;

      void set_checklayersids (bool);
      bool get_checklayersids () const;

      void set_excludesmodule (bool);
      bool get_excludesmodule () const;

      void set_performlinearinterpolation (bool);
      bool get_performlinearinterpolation () const;

      void set_printoutstdinfo (bool);
      bool get_printoutstdinfo () const;

      void set_verbose (bool);
      bool get_verbose () const;

      void set_specificseq (const char *);
      const std::string & get_specificseq () const;

      void set_etalimits (double & min, double & max);
      void set_philimits (double & min, double & max);
      void set_ptlimits (double & min, double & max);
      void set_d0limits (double & min, double & max);
      void set_z0limits (double & min, double & max);
      void set_maxnumoflayers (int val);
      void set_chargesign (int val);
      void set_maxnumoftracks (unsigned int);

      void get_etalimits (double & min, double & max) const;
      void get_philimits (double & min, double & max) const;
      void get_ptlimits (double & min, double & max) const;
      void get_d0limits (double & min, double & max) const;
      void get_z0limits (double & min, double & max) const;
      int get_maxnumoflayers () const;
      int get_chargesign () const;
      unsigned int get_maxnumoftracks () const;

      bool reading_from_root_file (const pca::pcafitter & fitter, 
          arma::mat & paramin, arma::mat & coordin, 
          arma::vec & ptvalsout);

    private:
  
      std::vector<track_str> tracks_vct_;

      int errnum_;
      std::string errmsg_, filename_;

      bool rzplane_, rphiplane_, chargeoverpt_, 
          excludesmodule_, verbose_, checklayersids_, 
          savecheckfiles_, printoutstdinfo_, 
          performlinearinterpolation_;
      double etamin_, etamax_, phimin_, phimax_,
            z0min_, z0max_, d0min_, d0max_, 
            ptmin_, ptmax_;
      int maxnumoflayers_, chargesign_;
      unsigned int maxnumoftracks_;
      std::string specificseq_;

      bool extract_data (const pca::pcafitter & fitter, 
          arma::mat & paramin, arma::mat & coordin, 
          arma::vec & ptvalsout);

      bool check_if_withinranges (const int & charge, 
          const double & eta, const double & phi, 
          const double & z0, const double & d0, 
          const double & pt, const std::string &) const;

      bool linearinterpolation ();

      void set_errmsg (int , const std::string &);
      void reset_error ();
  };
};

#endif
