#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>

// Loriano: let's try Armadillo quick code 
#include <armadillo>
#include <cassert>

#include <getopt.h>
#include <unistd.h>
#include <alloca.h>

#include <version.h>
#include <pcafitter.hpp>

// lstorchi: here all the basic routines, in principles can be used 
//           to start building a proper class.

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

///////////////////////////////////////////////////////////////////////////////
//  PUBLIC
///////////////////////////////////////////////////////////////////////////////

pcafitter::pcafitter()
{
  coordim_ = 6*3;
  paramdim_ = 5;
}

pcafitter::~pcafitter()
{
}

int pcafitter::get_paramdim () const
{
  return paramdim_;
}

void pcafitter::set_paramdim (int i)
{
  paramdim_ = i;
}

int pcafitter::get_coordim () const
{
  return coordim_;
}

void pcafitter::set_coordim (int i)
{
  coordim_ = i;
}

std::string pcafitter::get_version_string( )
{
  std::ostringstream version;
  version << PCAFITTER_MAJOR_VER << "." << PCAFITTER_MINOR_VER << "." << 
    PCAFITTER_PATCH_VER ;

  return version.str();
}

std::string pcafitter::paramidx_to_string (int i)
{
  switch (i)
  {
    case PTIDX:
      return "oneoverpt";
      break;
    case PHIIDX:
      return "phi";
      break;
    case TETHAIDX:
      return "cot(tetha/2)";
      break;
    case Z0IDX:
      return "z0";
      break;
    case D0IDX:
      return "d0";
      break;
    default:
      return "";
      break;
  }
}

void pcafitter::compute_parameters (const arma::mat & cmtx, 
    const arma::rowvec & q, 
    const arma::mat & coord, 
    double * oneoverptcmp, double * phicmp, 
    double * etacmp, double * z0cmp, 
    double * d0cmp)
{
  for (int i=0; i<(int)coord.n_rows; ++i)
  {
    // 1 / pt 
    oneoverptcmp[i] = q(PTIDX);
    for (int k=0; k<coordim_; ++k)
      oneoverptcmp[i] += cmtx(PTIDX,k)*coord(i,k);

    // phi
    phicmp[i] = q(PHIIDX);
    for (int k=0; k<coordim_; ++k)
      phicmp[i] += cmtx(PHIIDX,k)*coord(i,k);

    // eta
    etacmp[i] = q(TETHAIDX);
    for (int k=0; k<coordim_; ++k)
      etacmp[i] += cmtx(TETHAIDX,k)*coord(i,k);

    // z0
    z0cmp[i] = q(Z0IDX);
    for (int k=0; k<coordim_; ++k)
      z0cmp[i] += cmtx(Z0IDX,k)*coord(i,k);

    // d0
    d0cmp[i] = q(D0IDX);
    for (int k=0; k<coordim_; ++k)
      d0cmp[i] += cmtx(D0IDX,k)*coord(i,k);
  }
}


void pcafitter::read_armmat (const char * fname, arma::mat & cmtx)
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

void pcafitter::read_armvct (const char * fname, arma::rowvec & q)
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

void pcafitter::write_armmat (const char * fname, arma::mat & cmtx)
{
  std::ofstream myfilec(fname, std::ios::binary);
  myfilec.write((const char*)&(cmtx.n_rows), sizeof(cmtx.n_rows));
  myfilec.write((const char*)&(cmtx.n_cols), sizeof(cmtx.n_cols));
  for (int i=0; i<(int)cmtx.n_rows; i++)
    for (int j=0; j<(int)cmtx.n_cols; j++)
      myfilec.write((const char*)&(cmtx(i, j)), sizeof(cmtx(i, j)));
  myfilec.close();
}

void pcafitter::write_armvct (const char * fname, arma::rowvec & q)
{
  std::ofstream myfileq(fname, std::ios::binary);
  myfileq.write((const char*)&(q.n_cols), sizeof(q.n_cols));
  for (int i=0; i<(int)q.n_cols; i++)
    myfileq.write((const char*)&(q(i)), sizeof(q(i)));
  myfileq.close();
}

void pcafitter::write_to_file (const char * fname, 
    const arma::mat & vec, int idx) 
{
  std::ofstream myfile(fname);
    
  for (int i=0; i<(int)vec.n_rows; i++)
    myfile << vec(i, idx) << std::endl;
    
  myfile.close();
}

void pcafitter::reading_from_file (const char * filename,
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

void pcafitter::select_bigger_sub (
    const std::map<std::string, int> & sublist, 
    bool verbose, int & maxnumber, std::string & slctsubsec ) 
{
  maxnumber = -1;
  std::map<std::string, int>::const_iterator it = sublist.begin();
  for (; it != sublist.end(); ++it)
  {
    if (it->second > maxnumber)
    {
      maxnumber = it->second;
      slctsubsec = it->first;
    }
  
    if (verbose)
      std::cout << "Subsector " << it->first << " has " << 
          it->second << " tracks " << std::endl;
  } 
}

void pcafitter::extract_sub (const std::vector<std::string> & sublist, 
    const std::string & slctsubsec, 
    const arma::mat & paramin,
    const arma::mat & coordin, 
    arma::mat & param,
    arma::mat & coord)
{
  int dim = 0;

  std::vector<std::string>::const_iterator it = sublist.begin();
  for (; it != sublist.end(); ++it)
    if (*it == slctsubsec)
      dim++;

  coord.set_size(dim,coordim_);
  param.set_size(dim,paramdim_);

  int k = 0;
  it = sublist.begin();
  for (int i=0; it != sublist.end(); ++it, ++i)
  {
    if (*it == slctsubsec)
    {
      for (int j = 0; j < coordim_; ++j)
        coord(k,j) = coordin(i,j);
      for (int j = 0; j < paramdim_; ++j)
        param(k,j) = paramin (i,j);
      ++k;
    }
  }
}

void pcafitter::compute_pca_constants (
    const arma::mat & param, 
    const arma::mat & coord,
    arma::mat & cmtx, 
    arma::rowvec & q)
{
  arma::mat v = arma::zeros<arma::mat>(coordim_,coordim_);
  v = arma::cov(coord);

#ifdef DEBUG
  /* correlation matrix ricorda su dati standardizzati coincide con la matrice 
   *   di covarianza : 
   *   z = x -<x> / sigma */
  arma::mat corr = arma::zeros<arma::mat>(coordim_,coordim_);
  
  for (int i=0; i<(coordim_); ++i)
    for (int j=0; j<(coordim_); ++j)
      corr(i,j) = v(i,j) / sqrt(v(i,i)*v(j,j));
  
  std::cout << "Correlation matrix: " << std::endl;
  std::cout << corr;
#endif

  arma::mat vi = arma::zeros<arma::mat>(coordim_,coordim_);
  vi = arma::inv(v); 
  //vi = v.i();

#ifdef DEBUG
  std::cout << "inverse by cov matrix: " << std::endl;
  std::cout << v * vi ;
#endif

  // and so on ...
  arma::mat hcap = arma::zeros<arma::mat>(coordim_,paramdim_);
  arma::rowvec paramm = arma::zeros<arma::rowvec>(paramdim_);
  arma::rowvec coordm = arma::zeros<arma::rowvec>(coordim_);

  //hcap = arma::cov()

  /*
  for (int i=0; i<(int)coord.n_cols; ++i)
  {
    for (int j=0; j<(int)coord.n_rows; ++j)
      coordm(i) += coord(j,i);
    coordm(i) /= (double) coord.n_rows;
  } 
  std::cout << coordm << std::endl;
  */ 

  coordm = mean(coord, 0);

  /*
  for (int i=0; i<(int)coord.n_rows; ++i)
  {
    paramm(PTIDX) += param(i, PTIDX);
    paramm(PHIIDX) += param(i, PHIIDX);
    paramm(D0IDX) += param(i, D0IDX);
    paramm(TETHAIDX) += param(i, TETHAIDX);
    paramm(Z0IDX) += param(i, Z0IDX);
  }

  paramm(PTIDX) /= (double) coord.n_rows;
  paramm(PHIIDX) /= (double) coord.n_rows;
  paramm(D0IDX) /= (double) coord.n_rows;
  paramm(TETHAIDTETHAIDX(double) coord.n_rows;
  paramm(Z0IDX) /= (double) coord.n_rows;

  std::cout << paramm(PTIDX) << " " << mean(param, 0) << std::endl;
  */

  paramm = mean(param, 0);

  hcap = arma::cov(coord, param);

  //std::cout << hcap << std::endl;

  for (int i=0; i<paramdim_; ++i)
    for (int l=0; l<coordim_; ++l)
      for (int m=0; m<coordim_; ++m)
        cmtx(i,l) += vi(l,m) * hcap (m,i);

#ifdef DEBUG
  std::cout << "C matrix: " << std::endl;
  std::cout << cmtx;
#endif 

  for (int i=0; i<paramdim_; ++i)
  {
    q(i) = paramm(i);
    for (int l=0; l<coordim_; ++l)
      q(i) -= cmtx(i,l)*coordm(l);
  }

#ifdef DEBUG
  std::cout << "Q vector: " << std::endl;
  for (int i=0; i<paramdim_; ++i)
    std::cout << q(i) << std::endl;
#endif

}

int pcafitter::numofline (const char * fname) 
{ 
  int number_of_lines = 0;
  std::string line;
  std::ifstream myfile(fname);
  
  while (std::getline(myfile, line))
    ++number_of_lines;

  myfile.close();
                          
  return number_of_lines;
}

///////////////////////////////////////////////////////////////////////////////
//  PRIVATE
///////////////////////////////////////////////////////////////////////////////

