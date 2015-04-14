#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>

// Loriano: let's try Armadillo quick code 
#include <armadillo>
#include <cassert>

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
}


using namespace pca;

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
      int a, b, c, segid;
      double x, y, z;

      mytfp >> x >> 
               y >> 
               z >> 
               a >> b >> c >> segid; // segid I am reading because can be used as local ccordinate ?
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
    
    double valread, valr, phiread, d0read, z0read;

    mytfp >> valr >> 
             phiread >> 
             d0read >> 
             valread >> 
             z0read;
    
    if (check_to_read (useonlyeven,useonlyodd,i))
    {
      paramin(counter, PHIIDX) = phiread;
      paramin(counter, D0IDX) = d0read;
      paramin(counter, Z0IDX) = z0read;
      // lstorchi: I use this to diretcly convert input parameters into
      //     better parameters for the fitting 
      // cot (tetha/2) = 1 / e^(-eta)
      paramin(counter, TETHAIDX) = 1.0e0 / exp (-1.0e0 * valread);
      // use 1/pt 
      paramin(counter, PTIDX) = 1.0e0 / valr;

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

void pca::reading_from_file_split (const char * filename, 
     arma::mat & paramin, arma::mat & coordin, 
     int num_of_ent, bool useonlyeven, bool useonlyodd,
     bool rzplane, bool rphiplane)
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
      int a, b, c, segid;
      double x, y, z;

      mytfp >> x >> 
               y >> 
               z >> 
               a >> b >> c >> segid; // segid I am reading because can be used as local ccordinate ?
                                       // in case of l1tkstubs is the tp value here 
 
      if (check_to_read (useonlyeven,useonlyodd,i))
      {
        double ri = sqrt(pow(x, 2.0) + pow (y, 2.0));

        if (rzplane)
        {
          coordin(counter, j*2) = z;
          coordin(counter, j*2+1) = ri;
        }
        else if (rphiplane)
        {
          double phii = asin(y/ri);

          coordin(counter, j*2+1) = ri;
          coordin(counter, j*2) = phii;
        }
      }
    }
      
      
    double valread, valr, phiread, d0read, z0read;

    mytfp >> valr >> 
             phiread >> 
             d0read >> 
             valread >> 
             z0read;

    if (check_to_read (useonlyeven,useonlyodd,i))
    {
      if (rzplane)
      {
        paramin(counter, SPLIT_Z0IDX) = z0read;
        // lstorchi: I use this to diretcly convert input parameters into
        //     better parameters for the fitting 
        // cot (tetha/2) = 1 / e^(-eta)
        paramin(counter, SPLIT_COTTETHAIDX) = 1.0e0 / exp (-1.0e0 * valread);
      }
      else if (rphiplane)
      {
        paramin(counter, SPLIT_PHIIDX) = phiread;
        paramin(counter, SPLIT_PTIDX) = 1.0e0 / valr;
      }

      ++counter;
    }
  }

  mytfp.close();

}

