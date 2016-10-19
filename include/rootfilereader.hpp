#ifndef _PCAROOTFILEREADER_INC_ 
#define _PCAROOTFILEREADER_INC_

#include <armadillo>
#include <vector>
#include <map>
#include <set>

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
        std::vector<int16_t> i_x, i_y, i_z;
        std::vector<int> layer, ladder, module, segid;

        double x0, y0, z0, d0, pt, phi, eta;
        int pdg;

        std::string layersids;
      } track_str;

      typedef struct 
      {
        int dim;
        std::vector<double> r, phii, z;
        std::vector<int32_t> i_r, i_phii, i_z;
        std::vector<int> layer;

        double x0, y0, z0, d0, pt, phi, eta;
        int pdg;

        std::string layersids;
      } track_rphiz_str;

    public:
      rootfilereader();
      ~rootfilereader();
  
      void reset();

      void set_region_type (int);
      int get_region_type () const;

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

      void set_fkfiveoutofsix (bool, int);
      bool get_fkfiveoutofsix (int &) const;

      void set_performlinearinterpolation (bool);
      bool get_performlinearinterpolation () const;

      void set_printoutstdinfo (bool);
      bool get_printoutstdinfo () const;

      void set_verbose (bool);
      bool get_verbose () const;

      void set_multiple_charge_pe (bool in)
      {
        multiple_charge_pe_ = in;
      };
      bool get_multiple_charge_pe () const
      {
        return multiple_charge_pe_;
      };

      void set_towid (int);
      int get_towid () const;

      void set_specificseq (const char *);
      const std::string & get_specificseq () const;

      const std::string & get_actualseq () const;
      const std::set<std::string> & get_actualseq_set () const;
      std::string get_actualseq_set_string () const;
      
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

      void set_use3layers(std::set<int> &);
      void reset_use3layers ();

      bool reading_from_root_file (const pca::pcafitter & fitter, 
          arma::mat & paramin, arma::mat & coordin, 
          arma::vec & ptvalsout, arma::vec & etavalout, 
          std::vector<std::string> &, std::vector<std::string> &);
      
      bool info_from_root_file (unsigned int & numevent, 
          double & xmin, double & xmax, double & ymin, double & ymax,
          double & zmin, double & zmax, double & etamin, double & etamax,
          double & ptmin, double & ptmax, double & phimin, double & phimax, 
          double & x0min, double & x0max, double & y0min, double & y0max,
          double & z0min, double & z0max);

      // add back to fitted phi values
      double get_rotation_angle () const
      {
        return sec_phi_;
      }

      // need for const generation
      void apply_rotation_to_phi (bool in)
      {
        applyrotationtophi_ = in;
      }
      bool get_apply_rotation_to_phi () const
      {
        return applyrotationtophi_;
      }

      // need for const generation
      void apply_rotation_to_xy (bool in)
      {
        applyrotationtoxy_ = in;
      }
      bool get_apply_rotation_to_xy () const
      {
        return applyrotationtoxy_;
      }

    private:
  
      std::vector<track_str> tracks_vct_;

      int tow_;
      double sec_phi_;
      bool applyrotationtophi_;
      bool applyrotationtoxy_;

      int errnum_;
      std::string errmsg_, filename_;

      std::string layersid_;
      std::set<std::string> layersid_set_;

      bool use3layers_;
      std::set<int> tlayers_;

      bool rzplane_, rphiplane_, chargeoverpt_, 
          excludesmodule_, verbose_, checklayersids_, 
          savecheckfiles_, printoutstdinfo_, 
          performlinearinterpolation_,
          fkfiveoutofsix_;
      double etamin_, etamax_, phimin_, phimax_,
            z0min_, z0max_, d0min_, d0max_, 
            ptmin_, ptmax_;
      int maxnumoflayers_, chargesign_, layeridtorm_, 
          regiontype_;
      unsigned int maxnumoftracks_;
      std::string specificseq_;

      bool multiple_charge_pe_;

      bool extract_data (const pca::pcafitter & fitter, 
          arma::mat & paramin, arma::mat & coordin, 
          arma::vec & ptvalsout, arma::vec & etavalout, 
          std::vector<std::string> &, std::vector<std::string> &);

      bool check_if_withinranges (const int & charge, 
          const double & eta, const double & phi, 
          const double & z0, const double & d0, 
          const double & pt, const std::string &) const;

      bool linearinterpolation ();

      bool from_root_file (int & i, unsigned int & countlayerswithdupid,
          std::vector<float> & stubx, 
          std::vector<float> & stuby, std::vector<float> & stubz, 
          std::vector<float> & pt, std::vector<float> & x0, 
          std::vector<float> & y0, std::vector<float> & z0, 
          std::vector<float> & eta, std::vector<float> & phi, 
          std::vector<float> & pdg, std::vector<int> & moduleid, 
          std::set<std::string> & layersids_set, 
          std::set<int> & layeridlist, std::ofstream & ss, 
          std::ofstream & ptfile, 
          std::ofstream & phifile, std::ofstream & d0file, 
          std::ofstream & etafile, std::ofstream & z0file, 
          std::ofstream & sstrack);

      bool remove_layer ();
      bool remove_last_layer ();

      bool linearinterpolationrphiz (std::vector<track_rphiz_str> &);

      bool convertorphiz (std::vector<track_rphiz_str> &);

      void set_errmsg (int , const std::string &);
      void reset_error ();
  };
};

#endif
