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

// can be included in any case if -std=c++11
#include "stdint.h"

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

  bool is_avalid_layerid (int regid, bool excludesmodule, int layerid)
  {
    if (regid == ISBARREL)
    {
      if (excludesmodule)
      {
        if ((layerid >= 5) && (layerid <= 7))
          return true;
      }
      else
      {
        if ((layerid >= 5) && (layerid <= 10))
          return true;
      }

      return false;
    }
    else
    {
      return false;
    }
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
  verbose_ = false; 
  regiontype_ = ISBARREL;
  rphiplane_ = false; 
  chargeoverpt_ = true;
  excludesmodule_ = false; 
  checklayersids_ = false;
  savecheckfiles_ = true;
  printoutstdinfo_ = true;
  fkfiveoutofsix_ = false;

  useintbitewise_ = false;
  use3layers_ = false;
  tlayers_.clear();

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
  maxnumoftracks_ = (unsigned int)  INFINITY;
  specificseq_ = "";
  performlinearinterpolation_ = false;
  layeridtorm_ = -1;

  layersid_ = "";

  tracks_vct_.clear();

  reset_error();
  filename_ = "";

  tow_ = -1;
  sec_phi_ = 0.0;
  applyrotationtophi_ = false;
  applyrotationtoxy_ = false;
}

void rootfilereader::set_region_type (int in)
{
  regiontype_ = in;
}

int rootfilereader::get_region_type () const
{
  return regiontype_;
}

void rootfilereader::set_useintbitewise (bool in)
{
  useintbitewise_ = in;
}

bool rootfilereader::get_useintbitewise () const
{
  return useintbitewise_;
}

void rootfilereader::set_printoutstdinfo (bool in)
{
  printoutstdinfo_ =  in;
}

bool rootfilereader::get_printoutstdinfo () const
{
  return printoutstdinfo_;
}

void rootfilereader::set_fkfiveoutofsix (bool in, int ini)
{
  fkfiveoutofsix_ = in;
  layeridtorm_ = ini;
}

bool rootfilereader::get_fkfiveoutofsix (int & out) const
{
  out = layeridtorm_;
  return fkfiveoutofsix_;
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

void rootfilereader::set_towid (int towid)
{
  tow_ = towid;

  double sec_phi = (tow_ % 8) * M_PI / 4.0 - 0.4;
  sec_phi_ = sec_phi;
}

int rootfilereader::get_towid () const
{
  return tow_;
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

const std::string & rootfilereader::get_actualseq () const
{
  return layersid_;
}

void rootfilereader::set_performlinearinterpolation (bool in)
{
  performlinearinterpolation_ = in;
}

bool rootfilereader::get_performlinearinterpolation () const
{
  return performlinearinterpolation_;
}

bool rootfilereader::info_from_root_file (unsigned int & numevent, 
    double & xmin, double & xmax, double & ymin, double & ymax,
    double & zmin, double & zmax, double & etamin, double & etamax,
    double & ptmin, double & ptmax, double & phimin, double & phimax, 
    double & x0min, double & x0max, double & y0min, double & y0max,
    double & z0min, double & z0max)
{
  TFile* inputFile = TFile::Open(filename_.c_str());

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
  numevent = (unsigned int) nevent;

  if (printoutstdinfo_)
    std::cout << "Total num of events: " << nevent << std::endl;

  std::set<int> layeridlist;

  std::set<std::string> layersids_set;
  bool thefirst = true;

  for (Int_t i=0; i<nevent; ++i) 
  { 
     TT->GetEntry(i);
     
     if ((moduleid.size() == stubx.size()) &&
         (moduleid.size() == stuby.size()) && 
         (moduleid.size() == stubz.size()) && 
         (moduleid.size() == pt.size()) && 
         (moduleid.size() == x0.size()) && 
         (moduleid.size() == y0.size()) && 
         (moduleid.size() == z0.size()) && 
         (moduleid.size() == eta.size()) && 
         (moduleid.size() == phi.size()) && 
         (moduleid.size() == pdg.size()))
     {
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
       
       if (allAreEqual) // QA nel caso dei BankStubs questo check e' utile ?
       {
         int j = 0;
         for (; j<(int)moduleid.size(); ++j)
         {
           if (thefirst && (j == 0))
           {
             xmin = xmax = stubx[j];
             ymin = ymax = stuby[j];
             zmin = zmax = stubz[j];
           }
           else
           {
             if (xmin > stubx[j])
               xmin = stubx[j];
             if (ymin > stuby[j])
               ymin = stuby[j];
             if (zmin > stubz[j])
               zmin = stubz[j];
             
             if (xmax < stubx[j])
               xmax = stubx[j];
             if (ymax < stuby[j])
               ymax = stuby[j];
             if (zmax < stubz[j])
               zmax = stubz[j];
           }
       
           /*
           int value = moduleid[j];
           int layer = value/1000000;
           value = value-layer*1000000;
           int ladder = value/10000;
           value = value-ladder*10000;
           int module = value/100;
           value = value-module*100;
           int segid = value; 
           */
       
         }
       
         j = 0;
       
         if (thefirst)
         {
           ptmin = ptmax = pt[j];
           etamin = etamax = eta[j];
           phimin = phimax = phi[j];
       
           x0min = x0max = x0[j];
           y0min = y0max = y0[j];
           z0min = z0max = z0[j];

           thefirst = false;
         }
         else
         {
           if (ptmin > pt[j])
             ptmin = pt[j];
           if (etamin > eta[j])
             etamin = eta[j];
           if (phimin > phi[j])
             phimin = phi[j];
           if (x0min > x0[j])
             x0min = x0[j];
           if (y0min > y0[j])
             y0min = y0[j];
           if (z0min > z0[j])
             z0min = z0[j];
       
           if (ptmax < pt[j])
             ptmax = pt[j];
           if (etamax < eta[j])
             etamax = eta[j];
           if (phimax < phi[j])
             phimax = phi[j];
           if (x0max < x0[j])
             x0max = x0[j];
           if (y0max < y0[j])
             y0max = y0[j];
           if (z0max < z0[j])
             z0max = z0[j];
         }
       }
     }
 
     countevt++;

     stubx.clear();
     stuby.clear();
     stubz.clear();
     pt.clear();
     x0.clear();
     y0.clear();
     z0.clear();
     eta.clear();
     phi.clear();
     pdg.clear();

     if (countevt >= maxnumoftracks_)
       break;
  }

  inputFile->Close();

  return true;
}
 

bool rootfilereader::reading_from_root_file (
    const pca::pcafitter & fitter, arma::mat & paramin, 
    arma::mat & coordin, arma::vec & ptvalsout, 
    arma::vec & etavalout)
{
  TFile* inputFile = TFile::Open(filename_.c_str());

  std::ofstream ss, ssext, ptfile, phifile, d0file, 
    etafile, z0file, sstrack;

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

  if (fkfiveoutofsix_)
  {
    /*
    if (use3layers_)
    {
      set_errmsg (11, "Not yet implemented");
      return false;
    }
    */

    if (!is_avalid_layerid (regiontype_, excludesmodule_, layeridtorm_) )
    {
      set_errmsg (11, "Invalid layer to remove");
      return false;
    }
  }

  if (savecheckfiles_)
  {
    ss.open("bankstub.txt");
    ssext.open("bankstub_notequal.txt");

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

  std::set<std::string> layersids_set;

  for (Int_t i=0; i<nevent; ++i) 
  { 
     TT->GetEntry(i);
     
     if ((moduleid.size() == stubx.size())
        && (moduleid.size() == stuby.size())
        && (moduleid.size() == stubz.size())
        && (moduleid.size() == pt.size())
        && (moduleid.size() == x0.size())
        && (moduleid.size() == y0.size())
        && (moduleid.size() == z0.size())
        && (moduleid.size() == eta.size())
        && (moduleid.size() == phi.size())
        && (moduleid.size() == pdg.size()))
     {
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
       
       //if ((moduleid.size() == (unsigned int) maxnumoflayers_)  &&
       // cannot perfomr this check in hybrid and maybe endcap
       
       if (allAreEqual) // QA nel caso dei BankStubs questo check e' utile ?
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
       
           single_track.i_x.push_back((int16_t) 10 * stubx[j]);
           single_track.i_y.push_back((int16_t) 10 * stuby[j]);
           single_track.i_z.push_back((int16_t) 10 * stubz[j]);
       
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
         {
           ++countlayerswithdupid;
           // track with duplicated layid are removed 
         }
         else
         {
           if (regiontype_ == ISBARREL)
           {
             if (moduleid.size() == (unsigned int) maxnumoflayers_)
             {
               // do not copy duplicated 
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
                             << single_track.pdg << " " 
                             << single_track.i_x[i] << " " 
                             << single_track.i_y[i] << " " 
                             << single_track.i_z[i] << " "
                             << std::endl; 
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
             }
           }
           else if (regiontype_ == ISHYBRID)
           {
             if (moduleid.size() >= (unsigned int) maxnumoflayers_)
             {
               layersids_set.insert(single_track.layersids);
               if (check_if_withinranges (pdg[j], 
                     eta[j], phi[j], d0val, z0[j], 
                     pt[j], osss.str()))
               {
                 if (moduleid.size() == (unsigned int) maxnumoflayers_)
                 {
                   // do not copy duplicated 
                   tracks_vct_.push_back(single_track);
                 }
                 else if (moduleid.size() == (unsigned int) (maxnumoflayers_ + 1))
                 {
                   // need to remove last layer
                   single_track.x.pop_back();
                   single_track.y.pop_back();
                   single_track.z.pop_back();
                 
                   single_track.i_x.pop_back();
                   single_track.i_y.pop_back();
                   single_track.i_z.pop_back();
                 
                   single_track.layer.pop_back();
                   single_track.ladder.pop_back();
                   single_track.module.pop_back();
                   single_track.segid.pop_back();
       
                   single_track.dim = 6;
       
                   std::ostringstream ossl;
                   for (int i=0; i<(int)single_track.layer.size(); ++i)
                     ossl << single_track.layer[i];
       
                   single_track.layersids = ossl.str();
                 
                   tracks_vct_.push_back(single_track);
                 }
                 else if (moduleid.size() == (unsigned int) (maxnumoflayers_ + 2))
                 {
                   // need to remove last 2 layers
                   single_track.x.pop_back();
                   single_track.y.pop_back();
                   single_track.z.pop_back();
                 
                   single_track.i_x.pop_back();
                   single_track.i_y.pop_back();
                   single_track.i_z.pop_back();
                 
                   single_track.layer.pop_back();
                   single_track.ladder.pop_back();
                   single_track.module.pop_back();
                   single_track.segid.pop_back();
       
                   single_track.x.pop_back();
                   single_track.y.pop_back();
                   single_track.z.pop_back();
                 
                   single_track.i_x.pop_back();
                   single_track.i_y.pop_back();
                   single_track.i_z.pop_back();
                 
                   single_track.layer.pop_back();
                   single_track.ladder.pop_back();
                   single_track.module.pop_back();
                   single_track.segid.pop_back();
                 
                   single_track.dim = 6;
       
                   std::ostringstream ossl;
                   for (int i=0; i<(int)single_track.layer.size(); ++i)
                     ossl << single_track.layer[i];
       
                   single_track.layersids = ossl.str();
       
                   /*
                   std::set<int>::iterator iit = layeridset.begin();
                   for (; iit != layeridset.end(); ++iit)
                   {
                     std::cerr << *iit << " : ";
                   }
                   std::cerr << single_track.eta << std::endl;
                   std::cerr << " " << single_track.pt << std::endl;
                   std::cerr << std::endl;
                   */
       
                   tracks_vct_.push_back(single_track);
                 }
                 else 
                 {
                   std::set<int>::iterator iit = layeridset.begin();
                   for (; iit != layeridset.end(); ++iit)
                   {
                     std::cerr << *iit << " : ";
                   }
                   std::cerr << std::endl;
       
                   set_errmsg (1, "HYBRID maybe > 8 layers ");
                   return false;
                 }
                 
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
                             << single_track.pdg << " " 
                             << single_track.i_x[i] << " " 
                             << single_track.i_y[i] << " " 
                             << single_track.i_z[i] << " "
                             << std::endl; 
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
             }
             else
             {
               // there could be a 5 layers track ?
       
             }
           }
           else if (regiontype_ == ISENDCAP)
           {
             if (moduleid.size() >= (unsigned int) maxnumoflayers_)
             {
               layersids_set.insert(single_track.layersids);
             }
           }
         }
       }
     }

     countevt++;

     stubx.clear();
     stuby.clear();
     stubz.clear();
     pt.clear();
     x0.clear();
     y0.clear();
     z0.clear();
     eta.clear();
     phi.clear();
     pdg.clear();

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

    std::cout << "Sequences: " << std::endl;
    std::set<std::string>::iterator lidi = layersids_set.begin();
    for (; lidi != layersids_set.end(); ++lidi)
      std::cout << *lidi << std::endl;
  }

  if (regiontype_ == ISENDCAP)
  {
    set_errmsg (1, "ENDCAP not yet implemented ");
    return false;
  }
 

  return rootfilereader::extract_data (fitter, 
    paramin, coordin, ptvalsout, etavalout);
}

///////////////////////////////////////////////////////////////////////////////
//                                PRIVATE                                    //
///////////////////////////////////////////////////////////////////////////////

bool rootfilereader::convertorphiz (std::vector<track_rphiz_str> & 
    rphiz_tracks)
{
  if (tow_ == -1)
  {
    set_errmsg (-15, "Need to set towid");
    return false;
  }

  if (useintbitewise_) 
  {
    set_errmsg (-16, "Integer phi rotation not yet implemented");
    return false;
  }

  double ci = cos(-sec_phi_);
  double si = sin(-sec_phi_);

  std::vector<track_str>::const_iterator track = tracks_vct_.begin();
  for (; track != tracks_vct_.end(); ++track)
  {
    track_rphiz_str single_track;

    single_track.dim = track->dim;
    single_track.layer = track->layer;
    single_track.x0 = track->x0;
    single_track.y0 = track->y0;
    single_track.z0 = track->z0;
    single_track.d0 = track->d0;
    single_track.pt = track->pt;
    single_track.phi = track->phi;
    single_track.eta = track->eta;
    single_track.pdg = track->pdg;
    single_track.layersids = track->layersids;

    for (int j = 0; j < track->dim; ++j)
    {
      double x = track->x[j];
      double y = track->y[j];
 
      if (applyrotationtoxy_ )
      {
        x = track->x[j] * ci - track->y[j] * si;
        y = track->x[j] * si + track->y[j] * ci;
      }

      double ri = sqrt(pow(x, 2.0) + pow (y, 2.0));
      //double phii = acos(track->x[j]/ri);
      double phii = atan2(y, x);

      int32_t i_ri = sqrt(pow(track->i_x[j], 2.0) + 
          pow (track->i_y[j], 2.0));
      //int32_t i_phii = 50000*acos((double) track->i_x[j]/
      //    (double) i_ri);
      int32_t i_phii = 50000*atan2((double)track->i_y[j],
          (double)track->i_x[j]);

      single_track.z.push_back(track->z[j]);
      single_track.r.push_back(ri);
      single_track.phii.push_back(phii);

      single_track.i_z.push_back(track->i_z[j]);
      single_track.i_r.push_back(i_ri);
      single_track.i_phii.push_back(i_phii);
    }

    rphiz_tracks.push_back(single_track);
  }

  return true;
}

void rootfilereader::set_use3layers(std::set<int> & in)
{
  use3layers_ = true;
  tlayers_ = in;
}

void rootfilereader::reset_use3layers ()
{
  use3layers_ = false;
  tlayers_.clear();
}


bool rootfilereader::check_if_withinranges (const int & charge, 
    const double & eta, const double & phi, const double & z0, 
    const double & d0, const double & pt, 
    const std::string & layersid) const
{
  if (is_a_valid_layers_seq(layersid, maxnumoflayers_, 
        regiontype_, checklayersids_))
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
    arma::mat & paramin, arma::mat & coordin, arma::vec & ptvalsout, 
    arma::vec & etavalout)
{
  if (printoutstdinfo_)
    std::cout << "Extracted  " << tracks_vct_.size() << " tracks " << std::endl;

  if (tracks_vct_.size() == 0)
    return false;

  if (fkfiveoutofsix_)
  {
    maxnumoflayers_ = 5;

    if (!remove_layer ())
      return false;
  }

  if (excludesmodule_)
  {
    if (fitter.get_coordim() != (maxnumoflayers_ - 3) * 2)
    {
      set_errmsg (1, "Wrong coord dim");
      return false;
    }
  }
  else if (use3layers_)
  {
    if (fitter.get_coordim() != 3 * 2)
    {
      set_errmsg (1, "Wrong coord dim");
      return false;
    }
  }
  else
  {
    if (fitter.get_coordim() != maxnumoflayers_ * 2)
    {
      set_errmsg (1, "Wrong coord dim");
      return false;
    }
  }

  //std::cout << fitter.get_coordim() << std::endl;
  
  coordin.resize(tracks_vct_.size(), fitter.get_coordim());
  paramin.resize(tracks_vct_.size(), fitter.get_paramdim());
  ptvalsout.resize(tracks_vct_.size());
  etavalout.resize(tracks_vct_.size());

  int excludesmodval = 0;

  if (maxnumoflayers_ == 5)
    excludesmodval = 1;
  else if (maxnumoflayers_ ==  6)
    excludesmodval = 2;

#if 0 

  if (performlinearinterpolation_)
  {
    if (!linearinterpolation ())
      return false;

    excludesmodval = 2;
  }

  std::vector<track_rphiz_str> rphiz_tracks;
  if (!convertorphiz (rphiz_tracks))
    return false;

#else

  std::vector<track_rphiz_str> rphiz_tracks;
  if (!convertorphiz (rphiz_tracks))
    return false;

  if (performlinearinterpolation_)
  {
    if (!linearinterpolationrphiz (rphiz_tracks))
      return false;

    excludesmodval = 2;
  }

#endif

  /* leave the code as it was */
  int counter = 0;
  std::vector<track_rphiz_str>::const_iterator track = rphiz_tracks.begin();
  for (; track != rphiz_tracks.end(); ++track)
  {
    int jidx = 0;
    std::ostringstream osss;
    std::string actuallayersids = "";
    for (int j = 0; j < track->dim; ++j)
    {
      if (excludesmodule_)
        if (j > excludesmodval)
          continue;

      if (use3layers_)
        if (tlayers_.find(track->layer[j]) == tlayers_.end())
          continue;

      if (useintbitewise_)
      {
        if (rzplane_)
        {
          coordin(counter, jidx*2) = track->i_z[j];
          coordin(counter, jidx*2+1) = track->i_r[j];
        }
        else if (rphiplane_)
        {
          coordin(counter, jidx*2) = track->i_phii[j];
          coordin(counter, jidx*2+1) = track->i_r[j];
        }
      }
      else
      { 
        if (rzplane_)
        {
          coordin(counter, jidx*2) = track->z[j];
          coordin(counter, jidx*2+1) = track->r[j];
        }
        else if (rphiplane_)
        {
          coordin(counter, jidx*2) = track->phii[j];
          coordin(counter, jidx*2+1) = track->r[j];
        }
      }

      ++jidx;

      osss << track->layer[j] << ":";
    }

    actuallayersids = osss.str();
    actuallayersids.erase(actuallayersids.end()-1);

    // this should be removed if we decide to use a single 
    // set for each possible combination 
    if (checklayersids_ && (regiontype_ == ISBARREL))
    {
      if (layersid_ == "")
      {
        layersid_ = actuallayersids;
      }
      else
      {
        if (actuallayersids != layersid_)
        {
          set_errmsg(-1, "different layers id sequence " + 
              actuallayersids + " vs " + layersid_);
          return false;
        }
      }
    }

    ptvalsout(counter) = track->pt;
    etavalout(counter) = track->eta;

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
      double phiorig = track->phi;
      if (applyrotationtophi_ )
      {
        phiorig += sec_phi_;
        phiorig = fmod(phiorig + M_PI, 2 * M_PI) - M_PI;
      }

      paramin(counter, PCA_PHIIDX) = phiorig;

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

  // restore old value
  if (fkfiveoutofsix_)
    maxnumoflayers_ = 6;

  return true;
}

bool rootfilereader::remove_layer()
{
  std::vector<track_str>::iterator track = tracks_vct_.begin();
  for (; track != tracks_vct_.end(); ++track)
  {
    std::ostringstream osss;
    for (int j = 0; j < track->dim; ++j)
      if (track->layer[j] != layeridtorm_)
        osss << track->layer[j];
  
    track->layersids = osss.str();

    for (int j = 0; j < track->dim; ++j)
    {
      if (track->layer[j] == layeridtorm_)
      {
        track->x.erase(track->x.begin()+j);
        track->y.erase(track->y.begin()+j);
        track->z.erase(track->z.begin()+j);

        track->i_x.erase(track->i_x.begin()+j);
        track->i_y.erase(track->i_y.begin()+j);
        track->i_z.erase(track->i_z.begin()+j);

        track->layer.erase(track->layer.begin()+j);
        track->ladder.erase(track->ladder.begin()+j);
        track->module.erase(track->module.begin()+j);
        track->segid.erase(track->segid.begin()+j);
        track->dim--;

        break;
      }
    }
  }

  return true;
}

bool rootfilereader::linearinterpolationrphiz (
    std::vector<track_rphiz_str> & rphiz_tracks)
{

  if (!useintbitewise_)
  {
    if (maxnumoflayers_ == 5)
    {
      std::vector<track_rphiz_str>::iterator track = rphiz_tracks.begin();
      for (; track != rphiz_tracks.end(); ++track)
      {
        if (rphiplane_)
        {
          if (track->layersids == "678910")
          {
            set_errmsg (1, "TODO not yet implemented");
            return false;
          }
          else if (track->layersids == "578910")
          {
            // simplest approach we will need to store a single scalar 
            // and maybe a second one to remove the bias 
            double doverv = 0.4596;
            
            double v1 = track->r[1] - track->r[0]; 
            double v2 = track->phii[1] - track->phii[0]; 
            
            double pd1 = track->r[0] + doverv * v1;
            double pd2 = track->phii[0] + doverv * v2;
            
            std::vector<double>::iterator it = track->r.begin();
            ++it;
            track->r.insert(it, pd1);
            
            it = track->phii.begin();
            ++it;
            track->phii.insert(it, pd2);
            
            std::vector<int>::iterator iit = track->layer.begin();
            ++iit;
            track->layer.insert(iit, 6);
            
            track->layersids == "5678910";
            track->dim = 6;
    
            return true;
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
        else if (rzplane_)
        {
          if (track->layersids == "678910")
          {
            set_errmsg (1, "TODO not yet implemented");
            return false;
          }
          else if (track->layersids == "578910")
          {
            // simplest approach we will need to store a single scalar 
            // and maybe a second one to remove the bias 
            double doverv = 0.4596;
            
            double v1 = track->r[1] - track->r[0]; 
            double v2 = track->z[1] - track->z[0]; 
            
            double pd1 = track->r[0] + doverv * v1;
            double pd2 = track->z[0] + doverv * v2;
            
            std::vector<double>::iterator it = track->r.begin();
            ++it;
            track->r.insert(it, pd1);
            
            it = track->z.begin();
            ++it;
            track->z.insert(it, pd2);
            
            std::vector<int>::iterator iit = track->layer.begin();
            ++iit;
            track->layer.insert(iit, 6);
            
            track->layersids == "5678910";
            track->dim = 6;
    
            return true;
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
    }
    else 
    {
      set_errmsg (1, "Can work only using 5 layers out of six");
      return false;
    }
    
    return true;
  }
  else
  {
    set_errmsg (1, "INTDITEWISE not yet implemented");
    return false;
  }
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
        // simplest approach we will need to store a single scalar 
        // and maybe a second one to remove the bias 
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
        track->dim = 6;
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
