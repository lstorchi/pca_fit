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
#include <pcaffunctype.hpp>

// lstorchi: here all the basic routines, in principles can be used 
//           to start building a proper class.

///////////////////////////////////////////////////////////////////////////////
//  PUBLIC
///////////////////////////////////////////////////////////////////////////////

using namespace pca;

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
  paramname_.clear();
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

bool pcafitter::set_paramidx (int i, const std::string & val)
{ 
  reset_error();

  if (i < paramdim_)
  {
    paramname_[i] = val;

    return true;
  }

 set_errmsg (1, "invalid index");

 return false; 
}

std::string pcafitter::paramidx_to_string (int i)
{
  std::map<int, std::string>::iterator it;
  if ((it = paramname_.find(i)) != paramname_.end())
    return it->second;

  return "";
}

bool pcafitter::compute_parameters (const arma::mat & cmtx, 
    const arma::rowvec & q, 
    const arma::mat & coord, 
    double ** paraptr, int paramdim)
{
  reset_error();

  if (paramdim != paramdim_)
  {
    set_errmsg (2, "invalid number of parameters");
    return false;
  }

  for (int j=0; j<paramdim; ++j)
  {
    double *ptr = paraptr[j];
    for (int i=0; i<(int)coord.n_rows; ++i)
    {
      ptr[i] = q(j);
      for (int k=0; k<coordim_; ++k)
        ptr[i] += cmtx(j,k)*coord(i,k);
    }
  }

  return true;
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

const std::string & pcafitter::get_errmsg () const
{
  return errmsg_;
}

int pcafitter::get_errnum() const
{
  return errnum_;
} 

///////////////////////////////////////////////////////////////////////////////
//  PRIVATE
///////////////////////////////////////////////////////////////////////////////

 
void pcafitter::set_errmsg (int num, const std::string & msg)
{
  errnum_ = num;
  errmsg_ = msg;
}

void pcafitter::reset_error ()
{
  errnum_ = 0;
  errmsg_ = "";
}
