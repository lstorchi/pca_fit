#include <algorithm>
#include <iostream>
#include <iterator>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h" 
#include "TBasket.h"

// can be included in nay case if -std=c++11
#ifdef INTBITEWISE
#include "stdint.h"
#endif 

// Loriano: let's try Armadillo quick code 
#include <armadillo>
#include <cassert>

#include <set>

#include <pcafitter.hpp>
#include <pcaffunctype.hpp>
#include <rootfilereader.hpp>

#include <sys/stat.h>

using namespace pca;

namespace 
{
  bool check_charge (const int inval, const int chargesign)
  {
    if (chargesign == 0)
      return true;
    else 
    {
      if ((chargesign > 0) && (inval > 0))
        return true;

      if ((chargesign < 0) && (inval < 0))
        return true;
    }

    return false;
  }

  bool check_sequence (const std::string & layersid, 
      const std::string & specificseq)
  {
    if ((specificseq == "") || (specificseq == layersid))
      return true;

    return false;
  }
}

rootfilereader::rootfilereader () 
{
  reset();
}

rootfilereader::~rootfilereader ()
{
}

void rootfilereader::reset()
{
  rzplane_ = false;
  rphiplane_ = false; 
  chargeoverpt_ = true;
  excludesmodule_ = false; 
  verbose_ = false; 
  checklayersids_ = false;
  savecheckfiles_ = true;
  printoutstdinfo_ = true;

  etamin_ = -INFINITY; 
  etamax_ = INFINITY; 
  phimin_ = -INFINITY;
  phimax_ = INFINITY; 
  ptmin_ = -INFINITY; 
  ptmax_ = INFINITY;
  z0min_ = -INFINITY; 
  z0max_ = INFINITY;
  d0min_ = -INFINITY; 
  d0max_ = INFINITY; 
  maxnumoflayers_ = 6;
  chargesign_ = 0;
  maxnumoftracks_ = INFINITY;
  specificseq_ = "";
  performlinearinterpolation_ = false;

  reset_error();
  filename_ = "";
}

void rootfilereader::set_printoutstdinfo (bool in)
{
  printoutstdinfo_ =  in;
}

bool rootfilereader::get_printoutstdinfo () const
{
  return printoutstdinfo_;
}

void rootfilereader::set_filename (const std::string & in)
{
  filename_ = in;
}

void rootfilereader::set_savecheckfiles (bool in)
{
  savecheckfiles_ = in;
}

bool rootfilereader::get_savecheckfiles () const
{
  return savecheckfiles_;
}

void rootfilereader::get_etalimits (double & min, double & max) const
{
  min = etamin_;
  max = etamax_;
}

void rootfilereader::get_philimits (double & min, double & max) const
{
  min = phimin_;
  max = phimax_;
}

void rootfilereader::get_ptlimits (double & min, double & max) const
{
  min = ptmin_;
  max = ptmax_;
}

void rootfilereader::get_d0limits (double & min, double & max) const
{
  min = d0min_;
  max = d0max_;
}

void rootfilereader::get_z0limits (double & min, double & max) const
{
  min = z0min_;
  max = z0max_;
}

int rootfilereader::get_maxnumoflayers () const
{
  return maxnumoflayers_;
}

void rootfilereader::set_chargeoverpt (bool in)
{
  chargeoverpt_ = in;
}

bool rootfilereader::get_chargeoverpt () const
{
  return chargeoverpt_;
}

int rootfilereader::get_chargesign () const
{
  return chargesign_;
}

void rootfilereader::set_etalimits (double & min, double & max)
{
  etamin_ = min;
  etamax_ = max;
}

void rootfilereader::set_philimits (double & min, double & max)
{
  phimin_ = min;
  phimax_ = max;
}

void rootfilereader::set_ptlimits (double & min, double & max)
{
  ptmin_ = min;
  ptmax_ = max;
}

void rootfilereader::set_d0limits (double & min, double & max)
{
  d0min_ = min;
  d0max_ = max;
}

void rootfilereader::set_z0limits (double & min, double & max)
{
  z0min_ = min;
  z0max_ = max;
}

void rootfilereader::set_maxnumoflayers (int val)
{
  maxnumoflayers_ = val;
}

void rootfilereader::set_chargesign (int val)
{
  bool setval = true;
  set_chargeoverpt (setval);

  chargesign_ = val;
}

void rootfilereader::set_rzplane (bool in)
{
  rzplane_ = in;
}

bool rootfilereader::get_rzplane () const
{
  return rzplane_;
}

void rootfilereader::set_rphiplane (bool in)
{
  rphiplane_ = in;
}

bool rootfilereader::get_rphiplane () const
{
  return rphiplane_;
}

void rootfilereader::set_checklayersids (bool in)
{
  checklayersids_ = in;
}

bool rootfilereader::get_checklayersids () const
{
  return checklayersids_;
}

void rootfilereader::set_excludesmodule (bool in)
{
  excludesmodule_ = in;
}

bool rootfilereader::get_excludesmodule () const
{
  return excludesmodule_;
}

void rootfilereader::set_verbose (bool in)
{
  verbose_ = in;
}

bool rootfilereader::get_verbose () const
{
  return verbose_;
}

const std::string & rootfilereader::get_errmsg () const
{
  return errmsg_;
}

int rootfilereader::get_errnum() const
{
  return errnum_;
} 

unsigned int rootfilereader::get_maxnumoftracks() const
{
  return maxnumoftracks_;
} 

void rootfilereader::set_maxnumoftracks(unsigned int in)
{
  maxnumoftracks_ = in;
} 

void rootfilereader::set_specificseq (const char * in)
{
  specificseq_ = in;
}

const std::string & rootfilereader::get_specificseq () const
{
  return specificseq_;
}

void rootfilereader::set_performlinearinterpolation (bool in)
{
  performlinearinterpolation_ = in;
}

bool rootfilereader::get_performlinearinterpolation () const
{
  return performlinearinterpolation_;
}

bool rootfilereader::reading_from_root_file (
    const pca::pcafitter & fitter, arma::mat & paramin, 
    arma::mat & coordin, arma::vec & ptvalsout)
{
  TFile* inputFile = TFile::Open(filename_.c_str());

  std::ofstream ss, ssext, ssext1, ssext2, ptfile, 
    phifile, d0file, etafile, z0file, sstrack;

  if (rzplane_ && rphiplane_) 
  {
    set_errmsg (1, "Cannot use together rz and rphi plane");
    return false;
  }

  if (!rzplane_ && !rphiplane_) 
  {
    set_errmsg (1, "Select rz or rphi plane");
    return false;
  }

  if (savecheckfiles_)
  {
    ss.open("bankstub.txt");
    ssext.open("bankstub_notequal.txt");
    ssext1.open("bankstub_less_layers.txt");
    ssext2.open("bankstub_more_layers.txt");

    ptfile.open("pt_bankstubs.txt");
    phifile.open("phi_bankstubs.txt");
    d0file.open("d0_bankstubs.txt");
    etafile.open("eta_bankstubs.txt");
    z0file.open("z0_bankstubs.txt");

    sstrack.open("bankstub_filtered.txt");
  }

  TChain* TT = (TChain*) inputFile->Get("BankStubs");

  std::vector<int> moduleid, * p_moduleid; 
  p_moduleid = &moduleid;

  TT->SetBranchAddress("STUB_modid", &p_moduleid); // QA come determino layerid e altro ? 
                                                   //    devo caricare la geometria ?
  std::vector<float> stubx, * p_stubx, stuby, * p_stuby, stubz, * p_stubz,
    pt, * p_pt, x0, * p_x0, y0, * p_y0, z0, * p_z0, eta, * p_eta,
    phi, * p_phi;
  std::vector<float> pdg, * p_pdg;

  p_stubx = &stubx;
  p_stuby = &stuby;
  p_stubz = &stubz;
  p_pt = &pt;
  p_x0 = &x0;
  p_y0 = &y0;
  p_z0 = &z0;
  p_eta = &eta;
  p_phi = &phi;
  p_pdg = &pdg;

  TT->SetBranchAddress("STUB_x", &p_stubx);
  TT->SetBranchAddress("STUB_y", &p_stuby);
  TT->SetBranchAddress("STUB_z", &p_stubz);

  TT->SetBranchAddress("STUB_ptGEN", &p_pt);
  TT->SetBranchAddress("STUB_X0", &p_x0);
  TT->SetBranchAddress("STUB_Y0", &p_y0);
  TT->SetBranchAddress("STUB_Z0", &p_z0);
  TT->SetBranchAddress("STUB_etaGEN", &p_eta);
  TT->SetBranchAddress("STUB_PHI0", &p_phi);
  TT->SetBranchAddress("STUB_pdg", &p_pdg);

  unsigned int countevt = 0;
  Int_t nevent = TT->GetEntries(); 
  if (savecheckfiles_)
    ss << "We got " << nevent << " events in BankStubs" << std::endl; 

  if (printoutstdinfo_)
    std::cout << "Total num of events: " << nevent << std::endl;

  std::set<int> layeridlist;
  unsigned int countlayerswithdupid = 0;

  for (Int_t i=0; i<nevent; ++i) 
  { 
     TT->GetEntry(i);
     
     assert (moduleid.size() == stubx.size());
     assert (moduleid.size() == stuby.size());
     assert (moduleid.size() == stubz.size());
     assert (moduleid.size() == pt.size());
     assert (moduleid.size() == x0.size());
     assert (moduleid.size() == y0.size());
     assert (moduleid.size() == z0.size());
     assert (moduleid.size() == eta.size());
     assert (moduleid.size() == phi.size());
     assert (moduleid.size() == pdg.size());

     bool allAreEqual = ((std::find_if(z0.begin() + 1, z0.end(), 
        std::bind1st(std::not_equal_to<int>(), z0.front())) == z0.end()) &&
                        (std::find_if(x0.begin() + 1, x0.end(), 
        std::bind1st(std::not_equal_to<int>(), x0.front())) == x0.end()) &&
                        (std::find_if(y0.begin() + 1, y0.end(), 
        std::bind1st(std::not_equal_to<int>(), y0.front())) == y0.end()) &&
                        (std::find_if(pt.begin() + 1, pt.end(), 
        std::bind1st(std::not_equal_to<int>(), pt.front())) == pt.end()) &&
                        (std::find_if(eta.begin() + 1, eta.end(), 
        std::bind1st(std::not_equal_to<int>(), eta.front())) == eta.end()) &&
                        (std::find_if(phi.begin() + 1, phi.end(), 
        std::bind1st(std::not_equal_to<int>(), phi.front())) == phi.end()) &&
                        (std::find_if(pdg.begin() + 1, pdg.end(),
        std::bind1st(std::not_equal_to<int>(), pdg.front())) == pdg.end()));

     if ((moduleid.size() == (unsigned int) maxnumoflayers_)  && 
         allAreEqual) // QA nel caso dei BankStubs questo check e' utile ?
     {
       rootfilereader::track_str single_track;

       double d0val;
       //d0val = (y0[0]-(tan(phi[0])*x0[0]))*cos(phi[0]);
       d0val = y0[0]*cos(phi[0])-x0[0]*sin(phi[0]);
       //double d0val = x0[0];

       if (savecheckfiles_)
       {
         ptfile << pt[0] << std::endl;
         phifile << phi[0] << std::endl;
         d0file << d0val << std::endl;
         etafile << eta[0] << std::endl;
         z0file << z0[0] << std::endl;
         
         ss << i+1 << " " << moduleid.size() << std::endl;
       }

       single_track.dim = (int)moduleid.size();

       std::ostringstream osss;
       std::set<int> layeridset;
 
       int j = 0;
       for (; j<(int)moduleid.size(); ++j)
       {
         if (savecheckfiles_)
           ss << stubx[j] << " " << stuby[j] << " " <<
             stubz[j] << " ";
         
         single_track.x.push_back(stubx[j]);
         single_track.y.push_back(stuby[j]);
         single_track.z.push_back(stubz[j]);
         
         int value = moduleid[j];
         int layer = value/1000000;
         value = value-layer*1000000;
         int ladder = value/10000;
         value = value-ladder*10000;
         int module = value/100;
         value = value-module*100;
         int segid = value; // QA is just this ? from the source code seems so, I need to / by 10 ?

         osss << layer;

         if (savecheckfiles_)
           ss << layer << " " << ladder << " " << 
             module << " " << segid << " " << pdg[j] << std::endl;
         
         single_track.layer.push_back(layer);
         single_track.ladder.push_back(ladder);
         single_track.module.push_back(module);
         single_track.segid.push_back(segid);

         layeridset.insert(layer);
         layeridlist.insert(layer);
       }
       --j;

       if (savecheckfiles_)
         ss << pt[j]<< " "  <<
           phi[j] << " " << d0val << " " 
           << eta[j] << " " << z0[j] << " " <<
           x0[j] << " " << y0[j] << std::endl;

       single_track.pt = pt[j];
       single_track.pdg = pdg[j];
       single_track.phi = phi[j];
       single_track.d0 = d0val;
       single_track.eta = eta[j];
       single_track.x0 = x0[j];
       single_track.y0 = y0[j];
       single_track.z0 = z0[j];
       single_track.layersids = osss.str();
 
       if (layeridset.size() != (unsigned int)single_track.dim)
         ++countlayerswithdupid;

       if (check_if_withinranges (pdg[j], 
             eta[j], phi[j], d0val, z0[j], 
             pt[j], osss.str()))
       {
         tracks_vct_.push_back(single_track);

         if (savecheckfiles_)
         {
           sstrack << tracks_vct_.size() << " events " << std::endl;
           sstrack << i+1 << " " << moduleid.size() << std::endl;

           for (int i=0; i<(int)moduleid.size(); ++i)
           {
             sstrack << single_track.x[i] << " " 
                     << single_track.y[i] << " " 
                     << single_track.z[i] << " "
                     << single_track.layer[i] << " " 
                     << single_track.ladder[i] << " " 
                     << single_track.module[i] << " "
                     << single_track.segid[i] << " "
                     << single_track.pdg << std::endl; 
           }

           sstrack << single_track.pt << " "  
                   << single_track.phi << " " 
                   << single_track.d0 << " " 
                   << single_track.eta << " " 
                   << single_track.z0 << " " 
                   << single_track.x0 << " " 
                   << single_track.y0 << std::endl;
         }
       }

       countevt++;
     }
     else if ((moduleid.size() > (unsigned int) maxnumoflayers_) && allAreEqual)
     {
       if (savecheckfiles_)
       {
         double d0val;
         d0val = y0[0]*cos(phi[0])-x0[0]*sin(phi[0]);

         ptfile << pt[0] << std::endl;
         phifile << phi[0] << std::endl;
         d0file << d0val << std::endl;
         etafile << eta[0] << std::endl;
         z0file << z0[0] << std::endl;

         ssext2 << i+1 << " " << moduleid.size() << std::endl;

         int j = 0;
         for (; j<(int)moduleid.size(); ++j)
         {
           ssext2 << stubx[j] << " " << stuby[j] << " " <<
              stubz[j] << " ";
         
           int value = moduleid[j];
           int layer = value/1000000;
           value = value-layer*1000000;
           int ladder = value/10000;
           value = value-ladder*10000;
           int module = value/100;
           value = value-module*100;
           int segid = value; // QA is just this ? from the source code seems so, I need to / by 10 ?
           
           ssext2 << layer << " " << ladder << " " << 
             module << " " << segid << " " << pdg[j] << std::endl;
         }
         --j;
         
         ssext2 << pt[j]<< " "  <<
           phi[j] << " " << d0val << " " 
           << eta[j] << " " << z0[j] << " " <<
           x0[j] << " " << y0[j] << std::endl;
       }

       countevt++;
     }
     else if ((moduleid.size() < (unsigned int) maxnumoflayers_) && allAreEqual)
     {
       if (savecheckfiles_)
       {
         double d0val;
         d0val = y0[0]*cos(phi[0])-x0[0]*sin(phi[0]);

         ptfile << pt[0] << std::endl;
         phifile << phi[0] << std::endl;
         d0file << d0val << std::endl;
         etafile << eta[0] << std::endl;
         z0file << z0[0] << std::endl;

         ssext1 << i+1 << " " << moduleid.size() << std::endl;

         int j = 0;
         for (; j<(int)moduleid.size(); ++j)
         {
           ssext1 << stubx[j] << " " << stuby[j] << " " <<
              stubz[j] << " ";
         
           int value = moduleid[j];
           int layer = value/1000000;
           value = value-layer*1000000;
           int ladder = value/10000;
           value = value-ladder*10000;
           int module = value/100;
           value = value-module*100;
           int segid = value; // QA is just this ? from the source code seems so, I need to / by 10 ?
           
           ssext1 << layer << " " << ladder << " " << 
             module << " " << segid << " " << pdg[j] << std::endl;
         }
         --j;
         
         ssext1 << pt[j]<< " "  <<
           phi[j] << " " << d0val << " " 
           << eta[j] << " " << z0[j] << " " <<
           x0[j] << " " << y0[j] << std::endl;
       }

       countevt++;
     }
     else
     {
       if (savecheckfiles_)
       {
         double d0val;
         d0val = y0[0]*cos(phi[0])-x0[0]*sin(phi[0]);

         ptfile << pt[0] << std::endl;
         phifile << phi[0] << std::endl;
         d0file << d0val << std::endl;
         etafile << eta[0] << std::endl;
         z0file << z0[0] << std::endl;

         ssext << i+1 << " " << moduleid.size() << std::endl;

         int j = 0;
         for (; j<(int)moduleid.size(); ++j)
         {
          if (savecheckfiles_)
            ssext << stubx[j] << " " << stuby[j] << " " <<
               stubz[j] << " ";

          int value = moduleid[j];
          int layer = value/1000000;
          value = value-layer*1000000;
          int ladder = value/10000;
          value = value-ladder*10000;
          int module = value/100;
          value = value-module*100;
          int segid = value; // QA is just this ? from the source code seems so, I need to / by 10 ?
          
          ssext << layer << " " << ladder << " " << 
            module << " " << segid << " " << pdg[j] << std::endl;
         }
         --j;

         ssext << pt[j]<< " "  <<
           phi[j] << " " << d0val << " " 
           << eta[j] << " " << z0[j] << " " <<
           x0[j] << " " << y0[j] << std::endl;
       }

       countevt++;
     }

     if (countevt >= maxnumoftracks_)
       break;
  }

  inputFile->Close();

  if (savecheckfiles_)
  {
    ptfile.close();
    phifile.close();
    d0file.close();
    etafile.close();
    z0file.close();
    ss.close();
    ssext.close();
    ssext1.close();
    ssext2.close();
    sstrack.close();
  }

  if (printoutstdinfo_)
  {
    std::cout << "Layers IDs list: " << std::endl;
    std::set<int>::iterator lids = layeridlist.begin();
    for (; lids != layeridlist.end(); ++lids)
      std::cout << " " << (*lids) << std::endl;
    std::cout << std::endl;

    std::cout << "Event with DupIds: " << countlayerswithdupid << std::endl;
  }

  return rootfilereader::extract_data (fitter, 
    paramin, coordin, ptvalsout);
}

///////////////////////////////////////////////////////////////////////////////
//                                PRIVATE                                    //
///////////////////////////////////////////////////////////////////////////////

bool rootfilereader::check_if_withinranges (const int & charge, 
    const double & eta, const double & phi, const double & z0, 
    const double & d0, const double & pt, 
    const std::string & layersid) const
{
  bool isbarrel = true;

  if (is_a_valid_layers_seq(layersid, maxnumoflayers_, 
        isbarrel, checklayersids_))
    if (check_sequence (layersid, specificseq_))
      if (check_charge (charge, chargesign_))
        if ((eta <= etamax_) && (eta >= etamin_))
          if ((pt <= ptmax_) && (pt >= ptmin_))
            if ((phi <= phimax_) && (phi >= phimin_))
              if ((d0 <= d0max_) && (d0 >= d0min_))
                if ((z0 <= z0max_) && (z0 >= z0min_))
                  return true;

  return false;
}

void rootfilereader::set_errmsg (int num, const std::string & msg)
{ 
  errnum_ = num;
  errmsg_ = msg;
}

void rootfilereader::reset_error ()
{
  errnum_ = 0;
  errmsg_ = "";
}

bool rootfilereader::extract_data (const pca::pcafitter & fitter, 
    arma::mat & paramin, arma::mat & coordin, arma::vec & ptvalsout)
{
  if (printoutstdinfo_)
    std::cout << "Extracted  " << tracks_vct_.size() << " tracks " << std::endl;

  coordin.resize(tracks_vct_.size(), fitter.get_coordim());
  paramin.resize(tracks_vct_.size(), fitter.get_paramdim());
  ptvalsout.resize(tracks_vct_.size());

  int excludesmodval = 0;

  if (maxnumoflayers_ == 5)
    excludesmodval = 1;
  else if (maxnumoflayers_ ==  6)
    excludesmodval = 2;

  if (performlinearinterpolation_)
  {
    if (!linearinterpolation ())
      return false;

    excludesmodval = 2;
  }

  /* leave the code as it was */
  int counter = 0;
  std::vector<track_str>::const_iterator track = tracks_vct_.begin();
  for (; track != tracks_vct_.end(); ++track)
  {
    for (int j = 0; j < track->dim; ++j)
    {
      if (excludesmodule_)
        if (j > excludesmodval)
          continue;
      
      double ri = sqrt(pow(track->x[j], 2.0) + 
          pow (track->y[j], 2.0));

      if (rzplane_)
      {
        coordin(counter, j*2) = track->z[j];
        coordin(counter, j*2+1) = ri;
      }
      else if (rphiplane_)
      {
        double phii = acos(track->x[j]/ri);
       
        coordin(counter, j*2) = phii;
        coordin(counter, j*2+1) = ri;
      }
    }

    ptvalsout(counter) = track->pt;

    if (rzplane_)
    {
      paramin(counter, PCA_Z0IDX) = track->z0;
      // lstorchi: I use this to diretcly convert input parameters into
      //     better parameters for the fitting 
      // eta = -ln[tan(theta / 2)]
      // theta = 2 * arctan (e^(-eta))
      // cotan (theta) = cotan (2 * arctan (e^(-eta)))
      paramin(counter, PCA_COTTHETAIDX) =  
        cot(2.0 * atan (exp (-1.0e0 * track->eta)));
      //double theta = atan(1.0 /  paramread(counter, PCA_COTTHETAIDX));
      //std::cout << etaread << " " << theta * (180/M_PI) << std::endl;
      //just to visualize pseudorapidity 
      counter++;
    }
    else if (rphiplane_)
    {
      paramin(counter, PCA_PHIIDX) = track->phi;
      // use 1/pt
      if (chargeoverpt_)
      {
        if (chargesign_ < 0)
          paramin(counter, PCA_ONEOVERPTIDX) = -1.0e0 / track->pt;
        else
          paramin(counter, PCA_ONEOVERPTIDX) = 1.0e0 / track->pt;
      }
      else
        paramin(counter, PCA_ONEOVERPTIDX) = 1.0e0 / track->pt;

      ++counter;
    }

    if (verbose_ && printoutstdinfo_)
    {
      std::cout << "ETA : " << track->eta << std::endl;
      std::cout << "PT  : " << track->pt << std::endl;
      std::cout << "PHI : " << track->phi << std::endl;
      std::cout << "D0  : " << track->d0 << std::endl;
      std::cout << "Z0  : " << track->z0 << std::endl;
    }
  }

  return true;
}

bool rootfilereader::linearinterpolation ()
{
  if (maxnumoflayers_ == 5)
  {
    std::vector<track_str>::iterator track = tracks_vct_.begin();
    for (; track != tracks_vct_.end(); ++track)
    {
      if (track->layersids == "678910")
      {
        set_errmsg (1, "TODO not yet implemented");
        return false;
      }
      else if (track->layersids == "578910")
      {
        double doverv = 0.4596;

        double v1 = track->x[1] - track->x[0]; 
        double v2 = track->y[1] - track->y[0]; 
        double v3 = track->z[1] - track->z[0]; 

        double pd1 = track->x[0] + doverv * v1;
        double pd2 = track->y[0] + doverv * v2;
        double pd3 = track->z[0] + doverv * v3;

        std::vector<double>::iterator it = track->x.begin();
        ++it;
        track->x.insert(it, pd1);

        it = track->y.begin();
        ++it;
        track->y.insert(it, pd2);

        it = track->z.begin();
        ++it;
        track->z.insert(it, pd3);

        std::vector<int>::iterator iit = track->layer.begin();
        ++iit;
        track->layer.insert(iit, 6);

        iit = track->ladder.begin();
        ++iit;
        track->ladder.insert(iit, -1);

        iit = track->module.begin();
        ++iit;
        track->module.insert(iit, -1);
 
        iit = track->segid.begin();
        ++iit;
        track->segid.insert(iit, -1);

        track->layersids == "5678910";
      }
      else if (track->layersids == "568910")
      {
        set_errmsg (1, "TODO not yet implemented");
        return false;
      }
      else if (track->layersids == "567910")
      {
        set_errmsg (1, "TODO not yet implemented");
        return false;
      }
      else if (track->layersids ==  "567810")
      {
        set_errmsg (1, "TODO not yet implemented");
        return false;
      }
      else if (track->layersids == "56789")
      {
        set_errmsg (1, "TODO not yet implemented");
        return false;
      }
    }
  }
  else 
  {
    set_errmsg (1, "Can work only using 5 layers out of six");
    return false;
  }

  return true;
}
