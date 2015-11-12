#include <iostream>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>

#ifdef INTBITEWISE
#include "stdint.h"
#endif 

// Loriano: let's try Armadillo quick code 
#include <armadillo>
#include <cassert>

#include <set>

#include <pcafitter.hpp>
#include <pcaffunctype.hpp>

#include <sys/stat.h>

//#define DEBUG 1

#define STDDIM_PRV_VLS 1

namespace
{
  bool is_a_valid_layers_seq(const std::string & in, bool tocheck)
  {
    static const char * valid_layers_seq[STDDIM_PRV_VLS];
    valid_layers_seq[0] = "5678910";

    if (tocheck)
    {
      for (int i=0; i<STDDIM_PRV_VLS; i++)
        if (in == valid_layers_seq[i])
          return true;

      return false;
    }
  
    return true;
  }

  bool check_to_read (bool useonlyeven, bool useonlyodd, int i)
  {
    if (useonlyeven || useonlyodd)
    {
      if (useonlyeven)
      {
        if (!((i+1) % 2))
          return true;
      }

      if (useonlyodd)
      {
        if ((i+1) % 2)
          return true;
      }
    }
    else 
      return true;
        
    return false;
  }

  bool check_charge_sign (int chargesign, std::set<int> & pidset)
  {
    if (pidset.size() == 1)
    {
      if (chargesign == 0)
        return true;

      // make it clear 
      if ((chargesign == -1) && (*pidset.begin() < 0))
        return true;

      if ((chargesign == 1) && (*pidset.begin() > 0))
        return true;
    }

    return false;
  }
}


using namespace pca;

void pca::tokenize (const std::string & str, std::vector<std::string> & tokens,
        const std::string & delimiters)
{
  // Skip delimiters at beginning.
  std::string::size_type lastPos = 
  str.find_first_not_of(delimiters, 0);
  std::string::size_type pos = 
  str.find_first_of(delimiters, lastPos);
                 
  while (std::string::npos != pos || 
         std::string::npos != lastPos)
  {
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    lastPos = str.find_first_not_of(delimiters, pos);
    pos = str.find_first_of(delimiters, lastPos);
  }
}

double pca::cot (double x)
{
  return tan(M_PI_2 - x);
}

bool pca::file_exists(const std::string& filename)
{
  struct stat buf;
  if (stat(filename.c_str(), &buf) != -1)
    return true;
              
  return false;
}

double pca::delta_phi(double phi1, double phi2) // http://cmslxr.fnal.gov/source/DataFormats/Math/interface/deltaPhi.h
{ 
  double result = phi1 - phi2;
  
  while (result > M_PI) result -= 2*M_PI;
  while (result <= -M_PI) result += 2*M_PI;
  
  return result;
}


void pca::read_armmat (const char * fname, arma::mat & cmtx)
{
  int n, m;
  std::ifstream myfilec(fname, std::ios::binary | std::ios::in);
  myfilec.read((char *)&n, sizeof(int));
  myfilec.read((char *)&m, sizeof(int));
  cmtx.set_size(n, m);
  for (int i=0; i<(int)cmtx.n_rows; i++)
  {
    for (int j=0; j<(int)cmtx.n_cols; j++)
    {
      double v;
      myfilec.read((char *)&v, sizeof(v));

#ifdef INTBITEWISE      
      cmtx(i, j) = (int16_t) v;
#else
      cmtx(i, j) = v;
#endif
    }
  }

  myfilec.close();
}

void pca::read_armvct (const char * fname, arma::rowvec & q)
{
  int n;
  std::ifstream myfileq(fname, std::ios::binary);
  myfileq.read((char*)&(n), sizeof(n));
  q.set_size(n);
  for (int i=0; i<(int)q.n_cols; i++)
  {
    double v;
    myfileq.read((char*)&v, sizeof(v));

#ifdef INTBITEWISE
    q(i) = (int16_t) v;
#else
    q(i) = v;
#endif
  }

  myfileq.close();
}

void pca::write_armmat (const char * fname, arma::mat & cmtx)
{
  std::ofstream myfilec(fname, std::ios::binary);
  myfilec.write((const char*)&(cmtx.n_rows), sizeof(cmtx.n_rows));
  myfilec.write((const char*)&(cmtx.n_cols), sizeof(cmtx.n_cols));
  for (int i=0; i<(int)cmtx.n_rows; i++)
    for (int j=0; j<(int)cmtx.n_cols; j++)
      myfilec.write((const char*)&(cmtx(i, j)), sizeof(cmtx(i, j)));
  myfilec.close();
}

void pca::write_armvct (const char * fname, arma::rowvec & q)
{
  std::ofstream myfileq(fname, std::ios::binary);
  myfileq.write((const char*)&(q.n_cols), sizeof(q.n_cols));
  for (int i=0; i<(int)q.n_cols; i++)
    myfileq.write((const char*)&(q(i)), sizeof(q(i)));
  myfileq.close();
}

void pca::write_to_file (const char * fname, 
    const arma::mat & vec, int idx) 
{
  std::ofstream myfile(fname);
    
  for (int i=0; i<(int)vec.n_rows; i++)
    myfile << vec(i, idx) << std::endl;
    
  myfile.close();
}

int pca::numofline (const char * fname) 
{ 
  int number_of_lines = 0;
  std::string line;
  std::ifstream myfile(fname);
  
  while (std::getline(myfile, line))
    ++number_of_lines;

  myfile.close();
                          
  return number_of_lines;
}

int pca::get_num_of_ent (const char * fname)
{
  int num = 0;

  std::string line;
  std::ifstream myfile(fname);
  
  // read first fake line
  std::getline(myfile, line);
  while (std::getline(myfile, line))
  {
    std::istringstream splitstr(line);
    std::vector<std::string> tokens;

    std::copy(std::istream_iterator<std::string>(splitstr),
         std::istream_iterator<std::string>(),
         back_inserter(tokens));

    if (tokens.size() == 2)
    {
      num++;

      int dim = atoi(tokens[1].c_str());
      for (int i=0; i<dim+1; i++)
      {
        // read all lines within events 
        std::getline(myfile, line);
      }
    }
    else
    {
      //should add an error event handler 
      std::cerr << "Error in file format" << std::endl;
      return 0;
    }
  }

  myfile.close();

  return num;
}

void pca::global_to_relative (arma::mat & coordin,
    double coord1min, double coord2min)
{
  if (coordin.n_cols%2 == 0)
  {
    double min1coord = std::numeric_limits<double>::infinity();
    double min2coord = std::numeric_limits<double>::infinity();

    if ((coord1min == std::numeric_limits<double>::infinity()) &&
        (coord2min == std::numeric_limits<double>::infinity()))
    {
      for (int i=0; i<(int)coordin.n_cols; i+=2)
      {
        double minval = coordin.col(i).min();
        if (minval < min1coord)
          min1coord = minval;
        minval = coordin.col(i+1).min();
        if (minval < min2coord)
          min2coord = minval;
      }
    }
    else
    {
      min1coord = coord1min;
      min2coord = coord2min;
    }

    std::cout << "Min values: " << min1coord << " and " << 
      min2coord << std::endl;

    for (int i=0; i<(int)coordin.n_rows; ++i)
    {
      for (int j=0; j<(int)coordin.n_cols; j+=2)
      {
        coordin(i, j) -= min1coord;
        coordin(i, j+1) -= min2coord;
      }
    }
  }

}

bool pca::reading_from_file_split (const pca::pcafitter & fitter, 
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
      double x, y, z;

      mytfp >> x >> 
               y >> 
               z >> 
               a >> b >> c >> segid >> pid; 

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
          double ri = sqrt(pow(x, 2.0) + pow (y, 2.0));

          if (rzplane)
          {
            coordread(counter, j*2) = z;
            coordread(counter, j*2+1) = ri;
        
          }
          else if (rphiplane)
          {
            double phii = acos(x/ri);
          
            coordread(counter, j*2) = phii;
            coordread(counter, j*2+1) = ri;
          }
        }
      }
    }

    if (layeridset.size() != (unsigned int)realsize)
    {
      ++countlayerswithdupid;
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

    if (check_charge_sign(chargesign, pidset))
    {
      layersids.push_back(osss.str());
      etavals(counter) = etaread;
      phivals(counter) = phiread;
      ptvals(counter) = ptread;
      z0vals(counter) = z0read;
      d0vals(counter) = d0read;

      if (rzplane)
      {
        paramread(counter, SPLIT_Z0IDX) = z0read;
        paramread(counter, SPLIT_COTTHETAIDX) =  cot(2.0 * atan (exp (-1.0e0 * etaread)));
      }
      else if (rphiplane)
      {
        paramread(counter, SPLIT_PHIIDX) = phiread;
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
                  paramin(counter, j) = paramread(i, j);
                
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


  std::cout << "Layers IDs list: " << std::endl;
  std::set<int>::iterator lids = layeridlist.begin();
  for (; lids != layeridlist.end(); ++lids)
    std::cout << " " << (*lids) << std::endl;
  std::cout << std::endl;

  std::cout << "Event with DupIds: " << countlayerswithdupid << std::endl;

  mytfp.close();

  return true;
}

std::string pca::get_paramname_from_id (int id)
{
  switch (id)
  {
    case (1):
      return "eta";
    case (2):
      return "pt";
    case (3):
      return "z0";
    case (4):
      return "phi";
    case (5):
      return "x0";
    case (6):
      return "y0";
    case (7):
      return "d0";
    default:
      return "";
  }

  return "";
}
