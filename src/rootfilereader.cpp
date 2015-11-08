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

rootfilereader::rootfilereader () 
{
  rzplane_ = false;
  rphiplane_ = false; 
  chargeoverpt_ = false;
  excludesmodule_ = false; 
  verbose_ = false; 
  checklayersids_ = false;
  useodd_ = false;
  useeven_ = false;

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
  maxnumoflayers_ = INFINITY;
  chargesign_ = 0;
  maxnumoftracks_ = INFINITY;

  reset_error();
  filename_ = "";
}

rootfilereader::~rootfilereader ()
{
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

void rootfilereader::set_maxnumoflayers (int & val)
{
  maxnumoflayers_ = val;
}

void rootfilereader::set_chargesign (int & val)
{
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

void rootfilereader::set_useodd (bool in)
{
  useodd_ = in;
}

bool rootfilereader::get_useodd () const
{
  return useodd_;
}

void rootfilereader::set_useeven (bool in)
{
  useeven_ = in;
}

bool rootfilereader::get_useeven () const
{
  return useeven_;
}

void rootfilereader::set_chargeoverpt (bool in)
{
  chargeoverpt_ = in;
}

bool rootfilereader::get_chargeoverpt () const
{
  return chargeoverpt_;
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

bool rootfilereader::reading_from_root_file (
    const pca::pcafitter & fitter, arma::mat & paramin, 
    arma::mat & coordin, arma::vec & ptvalsout)
{
  TFile* inputFile = TFile::Open(filename_.c_str());

  std::ofstream ss("bankstub.txt");
  std::ofstream ssext("bankstub_notequal.txt");
  std::ofstream ssext1("bankstub_lesst6layers.txt");

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
  ss << "We got " << nevent << " events in BankStubs" << std::endl; 

  std::ofstream ptfile("pt_BankStubs.txt");
  std::ofstream phifile("phi_BankStubs.txt");
  std::ofstream d0file("d0_BankStubs.txt");
  std::ofstream etafile("eta_BankStubs.txt");
  std::ofstream z0file("z0_BankStubs.txt");

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
        std::bind1st(std::not_equal_to<int>(), phi.front())) == phi.end()));


     if ((moduleid.size() == 6)  && allAreEqual) // QA nel caso dei BankStubs questo check e' utile ?
     {
       double d0val;
       //d0val = (y0[0]-(tan(phi[0])*x0[0]))*cos(phi[0]);
       d0val = y0[0]*cos(phi[0])-x0[0]*sin(phi[0]);
       //double d0val = x0[0];

       ptfile << pt[0] << std::endl;
       phifile << phi[0] << std::endl;
       d0file << d0val << std::endl;
       etafile << eta[0] << std::endl;
       z0file << z0[0] << std::endl;

       ss << i+1 << " " << moduleid.size() << std::endl;

       int j = 0;
       for (; j<(int)moduleid.size(); ++j)
       {
#ifdef INTBITEWISE         
        //Can we provide these scale factors from outside
        int16_t stubX = stubx[j]*10;
        int16_t stubY = stuby[j]*10;
        int16_t stubZ = stubz[j]*10;

        ss << stubX << " " << stubY << " " <<
           stubZ << " ";
#else
        ss << stubx[j] << " " << stuby[j] << " " <<
          stubz[j] << " ";
#endif

        int value = moduleid[j];
        int layer = value/1000000;
        value = value-layer*1000000;
        int ladder = value/10000;
        value = value-ladder*10000;
        int module = value/100;
        value = value-module*100;
        int segid = value; // QA is just this ? from the source code seems so, I need to / by 10 ?

        ss << layer << " " << ladder << " " << 
          module << " " << segid << " " << pdg[j] << std::endl;
       }
       --j;

       ss << pt[j]<< " "  <<
         phi[j] << " " << d0val << " " 
         << eta[j] << " " << z0[j] << " " <<
         x0[j] << " " << y0[j] << std::endl;

       countevt++;
     }
     else if ((moduleid.size() > 6) && allAreEqual)
     {
       double d0val;
       //d0val = (y0[0]-(tan(phi[0])*x0[0]))*cos(phi[0]);
       d0val = y0[0]*cos(phi[0])-x0[0]*sin(phi[0]);
       //double d0val = x0[0];

       ptfile << pt[0] << std::endl;
       phifile << phi[0] << std::endl;
       d0file << d0val << std::endl;
       etafile << eta[0] << std::endl;
       z0file << z0[0] << std::endl;

       ss << i+1 << " " << moduleid.size() << std::endl;

       int j = 0;
       for (; j<(int)moduleid.size(); ++j)
       {
#ifdef INTBITEWISE
        int16_t stubX = stubx[j]*10;
        int16_t stubY = stuby[j]*10;
        int16_t stubZ = stubz[j]*10;

        ss << stubX << " " << stubY << " " <<
           stubZ << " ";
#else
        ss << stubx[j] << " " << stuby[j] << " " <<
           stubz[j] << " ";
#endif
        int value = moduleid[j];
        int layer = value/1000000;
        value = value-layer*1000000;
        int ladder = value/10000;
        value = value-ladder*10000;
        int module = value/100;
        value = value-module*100;
        int segid = value; // QA is just this ? from the source code seems so, I need to / by 10 ?

        ss << layer << " " << ladder << " " << 
          module << " " << segid << " " << pdg[j] << std::endl;
       }
       --j;

       ss << pt[j]<< " "  <<
         phi[j] << " " << d0val << " " 
         << eta[j] << " " << z0[j] << " " <<
         x0[j] << " " << y0[j] << std::endl;

       countevt++;
     }
     else if ((moduleid.size() < 6) && allAreEqual)
     {
       double d0val;
       //d0val = (y0[0]-(tan(phi[0])*x0[0]))*cos(phi[0]);
       d0val = y0[0]*cos(phi[0])-x0[0]*sin(phi[0]);
       //double d0val = x0[0];

       ptfile << pt[0] << std::endl;
       phifile << phi[0] << std::endl;
       d0file << d0val << std::endl;
       etafile << eta[0] << std::endl;
       z0file << z0[0] << std::endl;

       ssext1 << i+1 << " " << moduleid.size() << std::endl;

       int j = 0;
       for (; j<(int)moduleid.size(); ++j)
       {
#ifdef INTBITEWISE
        int16_t stubX = stubx[j]*10;
        int16_t stubY = stuby[j]*10;
        int16_t stubZ = stubz[j]*10;

        ssext1 << stubX << " " << stubY << " " <<
           stubZ << " ";
#else
        ssext1 << stubx[j] << " " << stuby[j] << " " <<
           stubz[j] << " ";
#endif
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

       countevt++;
     }
     else
     {
       double d0val;
       //d0val = (y0[0]-(tan(phi[0])*x0[0]))*cos(phi[0]);
       d0val = y0[0]*cos(phi[0])-x0[0]*sin(phi[0]);
       //double d0val = x0[0];

       ptfile << pt[0] << std::endl;
       phifile << phi[0] << std::endl;
       d0file << d0val << std::endl;
       etafile << eta[0] << std::endl;
       z0file << z0[0] << std::endl;

       ssext << i+1 << " " << moduleid.size() << std::endl;

       int j = 0;
       for (; j<(int)moduleid.size(); ++j)
       {
#ifdef INTBITEWISE
        int16_t stubX = stubx[j]*10;
        int16_t stubY = stuby[j]*10;
        int16_t stubZ = stubz[j]*10;

        ssext << stubX << " " << stubY << " " <<
           stubZ << " ";
#else
        ssext << stubx[j] << " " << stuby[j] << " " <<
           stubz[j] << " ";
#endif
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

       countevt++;
     }

     if (countevt >= maxnumoftracks_)
       break;
  }

  ptfile.close();
  phifile.close();
  d0file.close();
  etafile.close();
  z0file.close();

  inputFile->Close();

  ss.close();
  ssext.close();
  ssext1.close();

  return true;
}

/*****************************************************************************
 *                                PRIVATE                                    *
 *****************************************************************************/

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

/*
bool pca::reading_from_root_file (const pca::pcafitter & fitter, 
     const char * filename, 
     arma::mat & paramin, arma::mat & coordin, 
     bool useonlyeven, bool useonlyodd,
     bool rzplane, bool rphiplane, 
     double etamin, double etamax, 
     double ptmin, double ptmax,
     bool chargeoverpt, int chargesign, 
     bool excludesmodule, bool usealsod0, 
     bool usex0y0, int singleparam,
     double phimin, double phimax,
     double z0min, double z0max,
     double d0min, double d0max, 
     bool usealsox0,
     bool verbose,
     arma::vec & ptvalsout,
     bool checklayersids,
     int maxnumoflayers)
{
  std::vector<std::string> layersids;

  int extdim = 9;
  std::string line;

  arma::mat paramread;
  arma::mat coordread;

  arma::vec etavals;
  arma::vec phivals;
  arma::vec ptvals;
  arma::vec d0vals;
  arma::vec z0vals;

  // non performante but easy to go 
  int num_of_ent = pca::get_num_of_ent(filename);
  std::cout << "Num of events: " << num_of_ent << std::endl;

  coordread.set_size(num_of_ent, fitter.get_coordim());
  paramread.set_size(num_of_ent, fitter.get_paramdim());
  etavals.set_size(num_of_ent);
  phivals.set_size(num_of_ent);
  ptvals.set_size(num_of_ent);
  d0vals.set_size(num_of_ent);
  z0vals.set_size(num_of_ent);

  if (useonlyeven && useonlyodd)
  {
    useonlyeven = false;
    useonlyodd = false;
  }

  std::ifstream mytfp;
  mytfp.open (filename, std::ios::in);

  std::set<int> layeridlist;
  unsigned int countlayerswithdupid = 0;

  int counter = 0;
  std::getline (mytfp, line);
  for (int i = 0; i < num_of_ent; ++i)
  {
    int fake1, realsize;
    mytfp >> fake1 >> realsize ;

    std::ostringstream osss;
    std::set<int> pidset, layeridset;
    
    for (int j = 0; j < realsize; ++j)
    {
      int a, b, c, segid, pid;
#ifdef INTBITEWISE
      int16_t x, y, z;
#else
      double x, y, z;
#endif

      if (chargeoverpt)
      {
        mytfp >> x >> 
                 y >> 
                 z >> 
                 a >> b >> c >> segid >> pid; // segid I am reading because can be used as local ccordinate ?
                                              // in case of l1tkstubs is the tp value here 
                                              // pid is particle id (charge)
      }
      else
      {
        mytfp >> x >> 
                 y >> 
                 z >> 
                 a >> b >> c >> segid; // segid I am reading because can be used as local ccordinate ?
                                              // in case of l1tkstubs is the tp value here 
                                              // pid is particle id (charge)
      }

      osss << a;
      layeridset.insert(a);
      layeridlist.insert(a);

      if (j < maxnumoflayers)
      {
        if (excludesmodule)
          if (j > 2)
            continue;
        
        pidset.insert(pid);
        
        if (check_charge_sign(chargesign, pidset))
        {
          if (check_to_read (useonlyeven,useonlyodd,i))
          {
#ifdef INTBITEWISE
            int16_t ri = sqrt(pow(x, 2.0) + pow (y, 2.0));
#else
            double ri = sqrt(pow(x, 2.0) + pow (y, 2.0));
#endif 

            if (rzplane)
            {
              coordread(counter, j*2) = z;
              coordread(counter, j*2+1) = ri;
        
              // fake to use XY or similar plane
              //coordread(counter, j*2) = x;
              //coordread(counter, j*2+1) = y;
        
            }
            else if (rphiplane)
            {
//recast x and r for arccos calculation. Though arccos operation not permitted in integer representation.
//Check X-Y view instead
#ifdef INTBITEWISE
              int16_t phii = 10*acos((double) x/(double) ri);
#else
              double phii = acos(x/ri);
#endif
          
              coordread(counter, j*2) = phii;
              coordread(counter, j*2+1) = ri;

              // quick switch to XY view
              //coordread(counter, j*2) = x;
              //coordread(counter, j*2+1) = y;
 
            }
          }
        }
      }
    }

    if (layeridset.size() != (unsigned int)realsize)
    {
      ++countlayerswithdupid;
      //std::cerr << "Layer id duplicated " << std::endl;
    }

    //Need to change for Integer Representation , but for the time 
    //being its alright since bankstub.txt has int16_t
    double ptread, phiread, d0read, etaread, z0read,
           x0read, y0read;

    mytfp >> ptread >> 
             phiread >> 
             d0read >> 
             etaread >> 
             z0read >>
             x0read >>
             y0read;

    if (check_to_read (useonlyeven,useonlyodd,i))
    {
      layersids.push_back(osss.str());
      etavals(counter) = etaread;
      phivals(counter) = phiread;
      ptvals(counter) = ptread;
      z0vals(counter) = z0read;
      d0vals(counter) = d0read;

      if ((singleparam >= 1) && (singleparam <= 7))
      {
        if (check_charge_sign(chargesign, pidset))
        {
          switch (singleparam)
          {
            case (1):
              paramread(counter, 0) =  cot(2.0 * atan (exp (-1.0e0 * etaread)));
              break;
            case (2):
              if (chargeoverpt)
              {
                if (pidset.size() != 1)
                {
                  std::cerr << "pid values differ" << std::endl;
                  return false;
                }
             
                if (*(pidset.begin()) < 0)
                  paramread(counter, 0) = -1.0e0 / ptread;
                else
                  paramread(counter, 0) = 1.0e0 / ptread;
              }
              else
                paramread(counter, 0) = 1.0e0 / ptread;
 
              break;
            case (3):
              paramread(counter, 0) = z0read;
              break;
            case (4):
              paramread(counter, 0) = phiread;
              break;
            case (5):
              paramread(counter, 0) = x0read;
              break;
            case (6):
              paramread(counter, 0) = y0read;
              break;
            case (7):
              paramread(counter, 0) = d0read;
              break;
            default:
              paramread(counter, 0) = 0.0e0;
              std::cerr << "wrong paramindex value" << std::endl;
              break;
          }
        }
      }
      else
      {
        if (rzplane)
        {
          if (usex0y0)
          {
            paramread(counter, SPLIT_X0IDX) = x0read;
            paramread(counter, SPLIT_Y0IDX) = y0read;
          }
          else
          {
            paramread(counter, SPLIT_Z0IDX) = z0read;
            // lstorchi: I use this to diretcly convert input parameters into
            //     better parameters for the fitting 
            // eta = -ln[tan(theta / 2)]
            // theta = 2 * arctan (e^(-eta))
            // cotan (theta) = cotan (2 * arctan (e^(-eta)))
            paramread(counter, SPLIT_COTTHETAIDX) =  cot(2.0 * atan (exp (-1.0e0 * etaread)));
            //double theta = atan(1.0 /  paramread(counter, SPLIT_COTTHETAIDX));
            //std::cout << etaread << " " << theta * (180/M_PI) << std::endl;
            //just to visualize pseudorapidity 
          }
        
          if (usealsod0)
            paramread(counter, SPLIT_D0IDX) = d0read;
          else if (usealsox0)
            paramread(counter, SPLIT_X0IDX_NS) = x0read;
        }
        else if (rphiplane)
        {
          if (check_charge_sign(chargesign, pidset))
          {
            if (usex0y0)
            {
              paramread(counter, SPLIT_X0IDX) = x0read;
              paramread(counter, SPLIT_Y0IDX) = y0read;
            }
            else
            {
              paramread(counter, SPLIT_PHIIDX) = phiread;
              // use 1/pt
              if (chargeoverpt)
              {
                if (pidset.size() != 1)
                {
                  std::cerr << "pid values differ" << std::endl;
                  return false;
                }
             
                if (*(pidset.begin()) < 0)
                  paramread(counter, SPLIT_ONEOVERPTIDX) = -1.0e0 / ptread;
                else
                  paramread(counter, SPLIT_ONEOVERPTIDX) = 1.0e0 / ptread;
              }
              else
                paramread(counter, SPLIT_ONEOVERPTIDX) = 1.0e0 / ptread;
            }
        
            if (usealsod0)
              paramread(counter, SPLIT_D0IDX) = d0read;
            else if (usealsox0)
              paramread(counter, SPLIT_X0IDX_NS) = x0read;
          }
        }
      }

      if (check_charge_sign(chargesign, pidset))
        ++counter;
    }
  }

  assert (coordread.n_rows == paramread.n_rows);
  assert (etavals.n_rows == paramread.n_rows);
  assert (phivals.n_rows == paramread.n_rows);
  assert (ptvals.n_rows == paramread.n_rows);
  assert (z0vals.n_rows == paramread.n_rows);
  assert (d0vals.n_rows == paramread.n_rows);
  assert (layersids.size() == paramread.n_rows);

  extdim = 0;
  for (int i=0; i<(int)etavals.n_rows; ++i)
    if (is_a_valid_layers_seq(layersids[i], checklayersids))
      if ((etavals(i) <= etamax) && (etavals(i) >= etamin))
        if ((ptvals(i) <= ptmax) && (ptvals(i) >= ptmin))
          if ((phivals(i) <= phimax) && (phivals(i) >= phimin))
            if ((d0vals(i) <= d0max) && (d0vals(i) >= d0min))
              if ((z0vals(i) <= z0max) && (z0vals(i) >= z0min))
                ++extdim;

  coordin.resize(extdim, fitter.get_coordim());
  paramin.resize(extdim, fitter.get_paramdim());
  ptvalsout.resize(extdim);

  //std::cout << "extdim: " << extdim << std::endl;

  counter = 0;
  for (int i=0; i<(int)etavals.n_rows; ++i)
  {
    if (is_a_valid_layers_seq(layersids[i], checklayersids))
    {
      if ((etavals(i) <= etamax) && (etavals(i) >= etamin))
      {
        if ((ptvals(i) <= ptmax) && (ptvals(i) >= ptmin))
        {
          if ((phivals(i) <= phimax) && (phivals(i) >= phimin))
          {
            if ((d0vals(i) <= d0max) && (d0vals(i) >= d0min))
            {
              if ((z0vals(i) <= z0max) && (z0vals(i) >= z0min))
              {
                if (verbose)
                {
                  std::cout << "ETA : " << etavals(i) << std::endl;
                  std::cout << "PT  : " << ptvals(i) << std::endl;
                  std::cout << "PHI : " << phivals(i) << std::endl;
                  std::cout << "D0  : " << d0vals(i) << std::endl;
                  std::cout << "Z0  : " << z0vals(i) << std::endl;
                }
                
                for (int j=0; j<(int)paramread.n_cols; ++j)
                {
                  paramin(counter, j) = paramread(i, j);
                  //std::cout << "Track: " << counter+1 << std::endl;
                  //std::cout <<  paramin(counter, j) << " from " << 
                  //  paramread(i, j) << std::endl;
                }
                
                for (int j=0; j<(int)coordread.n_cols; ++j)
                  coordin(counter, j) = coordread(i, j);
                
                ptvalsout(counter) = ptvals(i);
                
                ++counter;
              }
            }
          }
        }
      }
    }
  }

  //std::cout << "counter: " << counter << std::endl;
  
  std::cout << "Layers IDs list: " << std::endl;
  std::set<int>::iterator lids = layeridlist.begin();
  for (; lids != layeridlist.end(); ++lids)
    std::cout << " " << (*lids) << std::endl;
  std::cout << std::endl;

  std::cout << "Event with DupIds: " << countlayerswithdupid << std::endl;

  mytfp.close();

  return true;
}
*/

