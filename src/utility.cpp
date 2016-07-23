#include <iostream>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>

#include "stdint.h"

// Loriano: let's try Armadillo quick code 
#include <armadillo>
#include <cassert>

#include <set>

#include <pcafitter.hpp>
#include <pcaffunctype.hpp>

#include <sys/stat.h>

//#define DEBUG 1

#define STDDIM_BARREL_PRV_VLS_6 1
#define STDDIM_BARREL_PRV_VLS_5 6
#define STDDIM_BARREL_PRV_VLS_5_NS 3

#define STDDIM_HYBRID_PRV_VLS_6 6
#define STDDIM_HYBRID_PRV_VLS_5 35


namespace 
{
  static const char * 
    valid_barrel_layers_seq_6[STDDIM_BARREL_PRV_VLS_6] = {"5678910"};

  static const char * 
    valid_barrel_layers_seq_5[STDDIM_BARREL_PRV_VLS_5] = {"678910",
                                                          "578910",
                                                          "568910",
                                                          "567910",
                                                          "567810",
                                                          "56789" };

  static const char * 
    valid_hybrid_layers_seq_6[STDDIM_HYBRID_PRV_VLS_6] = {"5678910",
                                                          "5678918",
                                                          "56781819",
                                                          "567181920",
                                                          "5618192021",
                                                          "51819202122"};

  static const char * 
    valid_hybrid_layers_seq_5[STDDIM_HYBRID_PRV_VLS_5] = {"678910",
                                                          "578910",
                                                          "568910",
                                                          "567910",
                                                          "567810",
                                                          "56789", 
                                                          "678918",
                                                          "578918",
                                                          "568918",
                                                          "567918",
                                                          "567818",
                                                          "6781819",
                                                          "5781819",
                                                          "5681819",
                                                          "5671819",
                                                          "567819",
                                                          "567818",
                                                          "67181920",
                                                          "57181920",
                                                          "56181920",
                                                          "5671920",
                                                          "5671820",
                                                          "5671819",
                                                          "618192021",
                                                          "518192021",
                                                          "56192021",
                                                          "56182021",
                                                          "56181921",
                                                          "56181920",
                                                          "1819202122",
                                                          "519202122",
                                                          "518202122",
                                                          "518192122",
                                                          "518192022",
                                                          "518192021"};
}

using namespace pca;

bool pca::validate_barrel_sequence_5 (const std::string & sequence)
{
  for (int i=0; i<STDDIM_BARREL_PRV_VLS_5; i++)
    if (sequence.compare(valid_barrel_layers_seq_5[i]) == 0)
      return true;

  return false;
}

bool pca::is_a_valid_layers_seq(const std::string & in, int maxnumoflayers, 
    const int regiontype, const bool tocheck)
{
  if (tocheck)
  {
    if (regiontype == ISBARREL)
    {
      if (maxnumoflayers == 6)
      {
        for (int i=0; i<STDDIM_BARREL_PRV_VLS_6; i++)
          if (in.compare(valid_barrel_layers_seq_6[i]) == 0)
            return true;
      }
      else if (maxnumoflayers == 5)
      {
        for (int i=0; i<STDDIM_BARREL_PRV_VLS_5; i++)
          if (in.compare(valid_barrel_layers_seq_5[i]) == 0)
            return true;
      }
    }
    else if (regiontype == ISHYBRID)
    {
      if (maxnumoflayers == 6)
      {
        for (int i=0; i<STDDIM_HYBRID_PRV_VLS_6; i++)
          if (in.compare(valid_hybrid_layers_seq_6[i]) == 0)
            return true;
      }
      else if (maxnumoflayers == 5)
      {
        for (int i=0; i<STDDIM_HYBRID_PRV_VLS_5; i++)
          if (in.compare(valid_hybrid_layers_seq_5[i]) == 0)
            return true;
      }
    }

    return false;
  }

  return true;
}

bool pca::check_to_read (bool useonlyeven, bool useonlyodd, int i)
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

bool pca::check_charge_sign (int chargesign, std::set<int> & pidset)
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

void pca::read_armmat_ibw (const char * fname, arma::mat & cmtx)
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
      cmtx(i, j) = (int32_t) v;
    }
  }

  myfilec.close();
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
      cmtx(i, j) = v;
    }
  }

  myfilec.close();
}

void pca::read_armvct_ibw (const char * fname, arma::rowvec & q)
{
  int n;
  std::ifstream myfileq(fname, std::ios::binary);
  myfileq.read((char*)&(n), sizeof(n));
  q.set_size(n);
  for (int i=0; i<(int)q.n_cols; i++)
  {
    double v;
    myfileq.read((char*)&v, sizeof(v));
    q(i) = (int32_t) v;
  }

  myfileq.close();
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
    q(i) = v;
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

#if 0 

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

#ifdef DEBUG
  double x0min = -1.0e0 * INFINITY, x0max = +1.0e0 * INFINITY;
  double y0min = -1.0e0 * INFINITY, y0max = +1.0e0 * INFINITY;

  /*
  x0min = -0.05;
  x0max =  0.05;
  y0min = -0.05;
  y0max =  0.05;
  */

  arma::mat xyzvals;
#endif

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

#ifdef DEBUG
  /* to set x0 and y0 range */
  xyzvals.set_size(num_of_ent, 3*maxnumoflayers);
  arma::vec x0vals;
  arma::vec y0vals;
  x0vals.set_size(num_of_ent);
  y0vals.set_size(num_of_ent);
#endif

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

#if 0
    if (realsize > maxnumoflayers)
    {
      std::cerr << "WARNING: More than 6 layers" << std::endl;
      std::cerr << "         Not checking missing stub TODO" << std::endl;
    }

    if (realsize < maxnumoflayers)
    {
      std::cerr << "WARNING: less than ",maxnumoflayers," layers" << std::endl;
      std::cerr << "         wont use it" << std::endl;
      continue;
    }
#endif

    std::ostringstream osss;
    std::set<int> pidset, layeridset;
    
    for (int j = 0; j < realsize; ++j)
    {
      int a, b, c, segid, pid;
      double x, y, z;

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
	    int32_t ri = sqrt(pow(x, 2.0) + pow (y, 2.0));
	    z  = (int32_t) z;
#else
	    double ri = sqrt(pow(x, 2.0) + pow (y, 2.0));
#endif
	    
#ifdef DEBUG
            xyzvals(counter, (3*j)+0) = x;
            xyzvals(counter, (3*j)+1) = y;
            xyzvals(counter, (3*j)+2) = z;
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

#ifdef INTBITEWISE
	      int32_t phii = 50000*acos((double) x/(double) ri);
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

#ifdef DEBUG
      x0vals(counter) = x0read;
      y0vals(counter) = y0read;
#endif

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
            paramread(counter, PCA_X0IDX) = x0read;
            paramread(counter, PCA_Y0IDX) = y0read;
          }
          else
          {
            paramread(counter, PCA_Z0IDX) = z0read;
            // lstorchi: I use this to diretcly convert input parameters into
            //     better parameters for the fitting 
            // eta = -ln[tan(theta / 2)]
            // theta = 2 * arctan (e^(-eta))
            // cotan (theta) = cotan (2 * arctan (e^(-eta)))
            paramread(counter, PCA_COTTHETAIDX) =  cot(2.0 * atan (exp (-1.0e0 * etaread)));
            //double theta = atan(1.0 /  paramread(counter, PCA_COTTHETAIDX));
            //std::cout << etaread << " " << theta * (180/M_PI) << std::endl;
            //just to visualize pseudorapidity 
          }
        
          if (usealsod0)
            paramread(counter, PCA_D0IDX) = d0read;
          else if (usealsox0)
            paramread(counter, PCA_X0IDX_NS) = x0read;
        }
        else if (rphiplane)
        {
          if (check_charge_sign(chargesign, pidset))
          {
            if (usex0y0)
            {
              paramread(counter, PCA_X0IDX) = x0read;
              paramread(counter, PCA_Y0IDX) = y0read;
            }
            else
            {
              paramread(counter, PCA_PHIIDX) = phiread;
              // use 1/pt
              if (chargeoverpt)
              {
                if (pidset.size() != 1)
                {
                  std::cerr << "pid values differ" << std::endl;
                  return false;
                }
             
                if (*(pidset.begin()) < 0)
                  paramread(counter, PCA_ONEOVERPTIDX) = -1.0e0 / ptread;
                else
                  paramread(counter, PCA_ONEOVERPTIDX) = 1.0e0 / ptread;
              }
              else
                paramread(counter, PCA_ONEOVERPTIDX) = 1.0e0 / ptread;
            }
        
            if (usealsod0)
              paramread(counter, PCA_D0IDX) = d0read;
            else if (usealsox0)
              paramread(counter, PCA_X0IDX_NS) = x0read;
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

#ifdef DEBUG  
  assert (x0vals.n_rows == paramread.n_rows);
  assert (y0vals.n_rows == paramread.n_rows);
#endif

  extdim = 0;
  for (int i=0; i<(int)etavals.n_rows; ++i)
    if (is_a_valid_layers_seq(layersids[i], checklayersids))
      if ((etavals(i) <= etamax) && (etavals(i) >= etamin))
        if ((ptvals(i) <= ptmax) && (ptvals(i) >= ptmin))
          if ((phivals(i) <= phimax) && (phivals(i) >= phimin))
            if ((d0vals(i) <= d0max) && (d0vals(i) >= d0min))
              if ((z0vals(i) <= z0max) && (z0vals(i) >= z0min))
#ifdef DEBUG              
                if ((x0vals(i) <= x0max) && (x0vals(i) >= x0min))
                  if ((y0vals(i) <= y0max) && (y0vals(i) >= y0min))
#endif
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
#ifdef DEBUG              
                if ((x0vals(i) <= x0max) && (x0vals(i) >= x0min))
                {
                  if ((y0vals(i) <= y0max) && (y0vals(i) >= y0min))
                  {
#endif
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
                   
#ifdef DEBUG                    
                    //for (int j=0; j<maxnumoflayers; ++j)
                    for (int j=0; j<1; ++j)
                      std::cerr << xyzvals(i,(3*j)+0) << " " << xyzvals(i,(3*j)+1) << " " 
                        << xyzvals(i,(3*j)+2) << std::endl;
#endif             
                   
                    
                    ++counter;
#ifdef DEBUG                  
                  }
                }
#endif
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
#endif

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

/* quick and very dirty */
void pca::dump_element (const pca::matrixpcaconst<double> & in, 
    std::ostream & out)
{
  out.precision(6);
  for (unsigned int i = 0; i<in.n_rows(); ++i)
  {
    for (unsigned int j = 0; j<in.n_cols(); ++j)
      out << std::scientific << in(i, j) << " ";
    out << std::endl;
  }
}
