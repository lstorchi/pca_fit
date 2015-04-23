#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>

// Loriano: let's try Armadillo quick code 
#include <armadillo>
#include <cassert>

#include <set>

#include <pcafitter.hpp>
#include <pcaffunctype.hpp>

#include <sys/stat.h>

namespace
{
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
  myfilec.read((char *)&n, sizeof(n));
  myfilec.read((char *)&m, sizeof(m));
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

void pca::reading_from_file (const char * filename,
    arma::mat & paramin, arma::mat & coordin, 
    arma::mat & layer, arma::mat & ladder, 
    arma::mat & module, 
    std::map<std::string, int> & subsectors, 
    std::map<std::string, int> & subladders,
    std::vector<std::string> & subsectorslist,
    std::vector<std::string> & subladderslist,
    int num_of_ent, bool usesegid,
    bool useonlyeven, bool useonlyodd)
{

  std::string line;
  std::ifstream mytfp;
  mytfp.open (filename, std::ios::in);

  if (useonlyeven && useonlyodd)
  {
    useonlyeven = false;
    useonlyodd = false;
  }

  int counter = 0;
  std::getline (mytfp, line);
  for (int i = 0; i < num_of_ent; ++i)
  {
    int fake1, fake2;
    // valori aggiunti solo di controllo 
    mytfp >> fake1 >> fake2 ;
#ifdef DEBUG    
    std::cout << fake1 << " " << fake2 << std::endl;
#endif
    std::ostringstream osss, ossl;
    osss << std::setfill('0');
    ossl << std::setfill('0');
    
    for (int j = 0; j < NUMOFLAYER; ++j)
    {
      int a, b, c, segid, pid;
      double x, y, z;

      mytfp >> x >> 
               y >> 
               z >> 
               a >> b >> c >> segid >> pid; // segid I am reading because can be used as local ccordinate ?
                                            // in case of l1tkstubs is the tp value here 
 
      if (check_to_read (useonlyeven,useonlyodd,i))
      {
        if (usesegid)
          coordin(counter, j) = (double) segid;
        else
        {
          coordin(counter, j*3) = x;
          coordin(counter, j*3+1) = y;
          coordin(counter, j*3+2) = z;
        }
 
        layer(counter, j) = a;
        ladder(counter, j) = b;
        module(counter, j) = c;
        
        osss << std::setw(2) << layer(counter, j);
        osss << std::setw(2) << ladder(counter, j);
        if (j != NUMOFLAYER-1)
          osss<<"-";
        
        ossl << std::setw(2) << layer(counter, j);
        ossl << std::setw(2) << ladder(counter, j);
        ossl << std::setw(2) << module(counter, j);
        if (j != NUMOFLAYER-1)
          ossl<<"-";
      }
    }
      
      
    if (check_to_read (useonlyeven,useonlyodd,i))
    {
      
      subsectorslist.push_back(osss.str());
      subladderslist.push_back(ossl.str());
      
      std::map<std::string, int>::iterator its = subsectors.find(osss.str());
      if (its == subsectors.end())
        subsectors[osss.str()] = 1;
      else 
        subsectors[osss.str()] += 1;
      
      std::map<std::string, int>::iterator itl = subladders.find(ossl.str());
      if (itl == subladders.end())
        subladders[ossl.str()] = 1;
      else 
        subladders[ossl.str()] += 1;
    }
    
    double ptread, phiread, d0read, etaread, z0read;

    mytfp >> ptread >> 
             phiread >> 
             d0read >> 
             etaread >> 
             z0read;
    
    if (check_to_read (useonlyeven,useonlyodd,i))
    {
      paramin(counter, PHIIDX) = phiread;
      paramin(counter, D0IDX) = d0read;
      paramin(counter, Z0IDX) = z0read;
      // lstorchi: I use this to diretcly convert input parameters into
      //     better parameters for the fitting 
      // eta = -ln[tan(tetha / 2)]
      // tetha = 2 * arctan (e^(-eta))
      // cotan (tetha) = cotan (2 * arctan (e^(-eta)))
      paramin(counter, COTTETHAIDX) =  cot(2.0 * atan (exp (-1.0e0 * etaread)));
      // use 1/pt 
      paramin(counter, ONEOVERPTIDX) = 1.0e0 / ptread;

      ++counter;
    }
  }

  mytfp.close();
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

bool pca::reading_from_file_split (const pca::pcafitter & fitter, 
     const char * filename, 
     arma::mat & paramin, arma::mat & coordin, 
     int num_of_ent, bool useonlyeven, bool useonlyodd,
     bool rzplane, bool rphiplane, 
     double etamin, double etamax, 
     bool chargeoverpt, int chargesign)
{
  int extdim = 9;
  std::string line;
  std::ifstream mytfp;
  mytfp.open (filename, std::ios::in);

  arma::mat paramread;
  arma::mat coordread;
  arma::vec etavals;

  coordread.set_size(num_of_ent, fitter.get_coordim());
  paramread.set_size(num_of_ent, fitter.get_paramdim());
  etavals.set_size(num_of_ent);

  if (useonlyeven && useonlyodd)
  {
    useonlyeven = false;
    useonlyodd = false;
  }

  int counter = 0;
  std::getline (mytfp, line);
  for (int i = 0; i < num_of_ent; ++i)
  {
    int fake1, fake2;
    // valori aggiunti solo di controllo 
    mytfp >> fake1 >> fake2 ;
#ifdef DEBUG    
    std::cout << fake1 << " " << fake2 << std::endl;
#endif
    std::ostringstream osss, ossl;
    osss << std::setfill('0');
    ossl << std::setfill('0');

    std::set<int> pidset;
    
    for (int j = 0; j < NUMOFLAYER; ++j)
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


      pidset.insert(pid);

      if (check_charge_sign(chargesign, pidset))
      {
        if (check_to_read (useonlyeven,useonlyodd,i))
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
      
    double ptread, phiread, d0read, etaread, z0read;

    mytfp >> ptread >> 
             phiread >> 
             d0read >> 
             etaread >> 
             z0read;

    if (check_to_read (useonlyeven,useonlyodd,i))
    {
      etavals(counter) = etaread;

      if (rzplane)
      {
        paramread(counter, SPLIT_Z0IDX) = z0read;
        // lstorchi: I use this to diretcly convert input parameters into
        //     better parameters for the fitting 
        // eta = -ln[tan(tetha / 2)]
        // tetha = 2 * arctan (e^(-eta))
        // cotan (tetha) = cotan (2 * arctan (e^(-eta)))
        paramread(counter, SPLIT_COTTETHAIDX) =  cot(2.0 * atan (exp (-1.0e0 * etaread)));
        //double tetha = atan(1.0 /  paramread(counter, SPLIT_COTTETHAIDX));
        //std::cout << etaread << " " << tetha * (180/M_PI) << std::endl;
        //just to visualize pseudorapidity 
      }
      else if (rphiplane)
      {
        if (check_charge_sign(chargesign, pidset))
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
      }

      if (check_charge_sign(chargesign, pidset))
        ++counter;
    }
  }

  assert (coordread.n_rows == paramread.n_rows);
  assert (etavals.n_rows == paramread.n_rows);

  extdim = 0;
  for (int i=0; i<(int)etavals.n_rows; ++i)
    if ((etavals(i) <= etamax) && (etavals(i) >= etamin))
      ++extdim;

  coordin.resize(extdim, fitter.get_coordim());
  paramin.resize(extdim, fitter.get_paramdim());

  //std::cout << "extdim: " << extdim << std::endl;

  counter = 0;
  for (int i=0; i<(int)etavals.n_rows; ++i)
  {
    if ((etavals(i) <= etamax) && (etavals(i) >= etamin))
    {
      for (int j=0; j<(int)paramread.n_cols; ++j)
      {
        paramin(counter, j) = paramread(i, j);
        //std::cout << "Track: " << counter+1 << std::endl;
        //std::cout <<  paramin(counter, j) << " from " << 
        //  paramread(i, j) << std::endl;
      }


      for (int j=0; j<(int)coordread.n_cols; ++j)
        coordin(counter, j) = coordread(i, j);

      ++counter;
    }
  }

  //std::cout << "counter: " << counter << std::endl;

  mytfp.close();

  return true;
}

