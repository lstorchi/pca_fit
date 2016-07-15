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
  useintbitewise_ = false;
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

bool pcafitter::compute_parameters (
    const arma::mat & cmtx, 
    const arma::rowvec & q, 
    const arma::mat & amtx,
    const arma::rowvec & kvct, 
    const arma::mat & coord, 
    int32_t ** paraptr,
    int paramdim, 
    arma::rowvec & chi2values)
{
  reset_error();

  if (!useintbitewise_)
  {
    set_errmsg (21, "this method can be used only with inbitewise");
    return false;
  }

  if (paramdim != paramdim_)
  {
    set_errmsg (2, "invalid number of parameters");
    return false;
  }

  for (int j=0; j<paramdim; ++j)
  {
    int32_t *ptr = paraptr[j];
    for (int i=0; i<(int)coord.n_rows; ++i)
    {
      ptr[i] = q(j);
      for (int k=0; k<coordim_; ++k)
      {
	ptr[i] += (int32_t) (double (cmtx(j,k)*coord(i,k)) /50.0);  
        //Probably no need to convert to double and then to int32
	//ptr[i] += ((cmtx(j,k)*coord(i,k))/1000.0);
      }
    }
  }
    
  std::cout << "coord.n_rows : " << coord.n_rows << std::endl;
  std::cout << "coord.n_cols : " << coord.n_cols << std::endl;

  for (int b=0; b<(int)coord.n_rows; ++b) // loop over tracks 
  {
    int32_t chi2 = 0.0;

    for (int i=0; i<coordim_-paramdim_; ++i)
    {
      int32_t val = 0;
      
      /* sum over j */
      for (int j=0; j<coordim_; ++j)
        val += amtx(i,j) * coord(b,j);

      /* wrong sign */
      val -= kvct(i);

      chi2 += val*val;
    }

    chi2values(b) = chi2;
  }

  return true;
}

bool pcafitter::compute_parameters (
    const arma::mat & cmtx, 
    const arma::rowvec & q, 
    const arma::mat & amtx,
    const arma::rowvec & kvct, 
    const arma::mat & coord, 
    double ** paraptr,
    int paramdim,
    arma::rowvec & chi2values,
    arma::rowvec & chi2values1,
    const arma::mat & vinv,
    const arma::rowvec & coordm)
{
  reset_error();

  if (useintbitewise_)
  {
    set_errmsg (21, "invalid number of parameters");
    return false;
  }

  if (paramdim != paramdim_)
  {
    set_errmsg (2, "use this method only with float");
    return false;
  }

  for (int j=0; j<paramdim; ++j)
  {
    double *ptr = paraptr[j];
    for (int i=0; i<(int)coord.n_rows; ++i)
    {
      ptr[i] = q(j);
      for (int k=0; k<coordim_; ++k)
	ptr[i] += (cmtx(j,k)*coord(i,k));
    }
  }
    
  /*
  for (int k=0; k<coordim_; ++k)
  {
    std::cout << coord[k].mean() << " " 
       << coord[k].stddev() << std::endl;
    std::cout << coordm(k) << std::endl;
  }
  */

  if ((vinv.n_rows != 0) && (vinv.n_cols != 0))
  {
    for (int k=0; k<(int)coord.n_rows; ++k)
    {
      for (int i=0; i<coordim_; ++i)
      {
        for (int j=0; j<coordim_; ++j)
        {
          chi2values1(k) += (coord(k,i) - coordm(i)) * 
            vinv(i, j) * (coord(k,j) - coordm(j));
        }
      }
  
      //std::cout << "chi2 using eq 10 pg 112 " << k << " ==> " 
      //   << chi2values1(k)  << std::endl;
    }
  }

  std::cout << "coord.n_rows : " << coord.n_rows << std::endl;
  std::cout << "coord.n_cols : " << coord.n_cols << std::endl;

  for (int b=0; b<(int)coord.n_rows; ++b) // loop over tracks 
  {
    double chi2 = 0.0;

    for (int i=0; i<coordim_-paramdim_; ++i)
    {
      double val = 0.0;
      
      /* sum over j */
      for (int j=0; j<coordim_; ++j)
        val += amtx(i,j) * coord(b,j);

      /* wrong sign */
      val -= kvct(i);

      //std::cout << "  val:" << val << std::endl;

      chi2 += val*val;
    }

    chi2values(b) = chi2;

    //std::cout << "chi2 using eq 11 pg 112 " << b << " ==> " 
    //  << chi2  << std::endl;
  }

  return true;
}

bool pcafitter::compute_parameters (
    std::vector<pca::matrixpcaconst<double> > & allconst,
    const arma::mat & coord, 
    const arma::mat & paramslt,
    const std::string & layersid,
    const std::string & pslayersid,
    int towerid,
    double ** paraptr,
    int paramdim, bool rphiplane,
    arma::rowvec & chi2values)
{

  reset_error();

  if (useintbitewise_)
  {
    set_errmsg (21, "not yet implemented");
    return false;
  }

  if (paramdim != paramdim_)
  {
    set_errmsg (2, "use this method only with float");
    return false;
  }

  std::cout << "coord.n_rows : " << coord.n_rows << std::endl;
  std::cout << "coord.n_cols : " << coord.n_cols << std::endl;

  for (int b=0; b<(int)coord.n_rows; ++b) // loop over tracks 
  {
    arma::mat cmtx;
    arma::rowvec q; 
    arma::mat amtx;
    arma::rowvec kvct; 
  
    pca::matrixpcaconst<double> cmtx_c(0, 0), 
      qvec_c(0, 0), amtx_c(0, 0), kvec_c(0, 0);
 
    double theta = atan(1.0e0 / paramslt(b, PCA_COTTHETAIDX));
    double etaorig = 0.0e0;
    double tantheta2 = tan (theta/2.0e0);
    if (tantheta2 < 0.0)
      etaorig = 1.0e0 * log (-1.0e0 * tantheta2);
    else
      etaorig = -1.0e0 * log (tantheta2);

    double qoverptorig = paramslt(b, PCA_ONEOVERPTIDX);
    int chargesignin = ((qoverptorig < 0.0) ? -1 : 1);
    double pt = ((double)chargesignin) / qoverptorig;

    if (rphiplane)
    {
      if (!import_pca_const (allconst, cmtx_c, qvec_c, 
            amtx_c, kvec_c, etaorig, pt, 
            chargesignin, layersid,  
            towerid, pca::FLOATPT, pca::RPHI))
        return false;
    }
    else 
    {
      chargesignin = 0;

      if (!import_pca_const (allconst, cmtx_c, qvec_c, 
            amtx_c, kvec_c, etaorig, pt, 
            chargesignin, pslayersid, 
            towerid, pca::FLOATPT, pca::RZ))
        return false;
    }

    pcamat_to_armamat (cmtx_c, cmtx);
    pcamat_to_armamat (amtx_c, amtx);
    pcamat_to_armarowvec (qvec_c, q);
    pcamat_to_armarowvec (kvec_c, kvct);

    //std::cout << "Here" << std::endl;

    for (int i=0; i<(int)paramdim; ++i)
    {
      double *ptr = paraptr[i];

      ptr[b] = q(i);
      for (int k=0; k<coordim_; ++k)
        ptr[b] += cmtx(i,k)*coord(b,k);
    }

    double chi2 = 0.0;

    for (int i=0; i<coordim_-paramdim_; ++i)
    {
      double val = 0.0;
      
      /* sum over j */
      for (int j=0; j<coordim_; ++j)
        val += amtx(i,j) * coord(b,j);

      /* wrong sign */
      val -= kvct(i);

      chi2 += val*val;
    }

    chi2values(b) = chi2;
  }

  return true;
}


bool pcafitter::compute_parameters_cgpca (
    std::vector<pca::matrixpcaconst<double> > & cgconst,
    std::vector<pca::matrixpcaconst<double> > & allconst,
    const arma::mat & coord, 
    const arma::mat & paramslt,
    const std::string & layersid,
    const std::string & pslayersid,
    int towerid,
    double ** paraptr,
    int paramdim, bool rphiplane,
    arma::rowvec & chi2values)
{

  reset_error();

  if (useintbitewise_)
  {
    set_errmsg (21, "not yet implemented");
    return false;
  }

  if (paramdim != paramdim_)
  {
    set_errmsg (2, "use this method only with float");
    return false;
  }

  std::cout << "coord.n_rows : " << coord.n_rows << std::endl;
  std::cout << "coord.n_cols : " << coord.n_cols << std::endl;

  for (int b=0; b<(int)coord.n_rows; ++b) // loop over tracks 
  {
    arma::mat cmtx;
    arma::rowvec q; 
    arma::mat amtx;
    arma::rowvec kvct; 

    if (coord.n_cols != 12)
    {
      set_errmsg (21, "not yet implemented for 5 oof 6");
      return false;
    }
  
    double theta = atan(1.0e0 / paramslt(b, PCA_COTTHETAIDX));
    double etaorig = 0.0e0;
    double tantheta2 = tan (theta/2.0e0);
    if (tantheta2 < 0.0)
      etaorig = 1.0e0 * log (-1.0e0 * tantheta2);
    else
      etaorig = -1.0e0 * log (tantheta2);

    double qoverptorig = paramslt(b, PCA_ONEOVERPTIDX);
    int chargesignin = ((qoverptorig < 0.0) ? -1 : 1);
    double pt = ((double)chargesignin) / qoverptorig;

    pca::matrixpcaconst<double> cmtx_c(0, 0), 
      qvec_c(0, 0), amtx_c(0, 0), kvec_c(0, 0);

    if (rphiplane)
    {
      /* CG PCA */
      pca::matrixpcaconst<double> cgcmtx_c(0, 0), 
        cgqvec_c(0, 0), cgamtx_c(0, 0), cgkvec_c(0, 0);

      //std::string layids = "5:8:10";
      std::string layids = "5:6:7:8:9:10";

      if (!import_pca_const (cgconst, cgcmtx_c, cgqvec_c, 
            cgamtx_c, cgkvec_c, etaorig, pt, 
            chargesignin, layids,  
            towerid, pca::FLOATPT, pca::RPHI))
        return false;

      arma::mat cgcmtx;
      arma::rowvec cgq; 
 
      pcamat_to_armamat (cgcmtx_c, cgcmtx);
      pcamat_to_armarowvec (cgqvec_c, cgq);

      //std::cout << "cgcmtx.n_rows : " << cgcmtx.n_rows << std::endl;
      //std::cout << "cgcmtx.n_cols : " << cgcmtx.n_cols << std::endl;

      /* 3 layers
      double cgqoverpt = 0.0; 
      cgqoverpt = cgq(0);
      // lay 5
      cgqoverpt += cgcmtx(0,0)*coord(b,0);
      cgqoverpt += cgcmtx(0,1)*coord(b,1);
      // lay 8
      cgqoverpt += cgcmtx(0,2)*coord(b,6);
      cgqoverpt += cgcmtx(0,3)*coord(b,7);
      // lay 10
      cgqoverpt += cgcmtx(0,4)*coord(b,10);
      cgqoverpt += cgcmtx(0,5)*coord(b,11);
      // 5     6     7     8    9      10 
      //0 1 / 2 3 / 4 5 / 6 7 / 8 9 / 10 11 
      */

      /* 5 layers */
      double cgqoverpt = 0.0; 
      cgqoverpt = cgq(0);
      // lay 5
      cgqoverpt += cgcmtx(0,0)*coord(b,0);
      cgqoverpt += cgcmtx(0,1)*coord(b,1);
      cgqoverpt += cgcmtx(0,2)*coord(b,2);
      cgqoverpt += cgcmtx(0,3)*coord(b,3);
      cgqoverpt += cgcmtx(0,4)*coord(b,4);
      cgqoverpt += cgcmtx(0,5)*coord(b,5);
      cgqoverpt += cgcmtx(0,6)*coord(b,6);
      cgqoverpt += cgcmtx(0,7)*coord(b,7);
      cgqoverpt += cgcmtx(0,8)*coord(b,8);
      cgqoverpt += cgcmtx(0,9)*coord(b,9);
      cgqoverpt += cgcmtx(0,10)*coord(b,10);
      cgqoverpt += cgcmtx(0,11)*coord(b,11);

      //std::cout << ((double)chargesignin)/cgqoverpt << " vs " << pt << std::endl; 

      double estpt = ((double)chargesignin)/cgqoverpt;

      if (estpt < 3.0)
      {
        std::cout << "CG PCA low 3.o GeV " << estpt << std::endl;
        estpt = 3.0 ;
      }

      if (estpt > 200.0)
      {
        std::cout << "CG PCA gt 200.o GeV " << estpt << std::endl;
        estpt = 200.0 ;
      }

      if (!import_pca_const (allconst, cmtx_c, qvec_c, 
            amtx_c, kvec_c, etaorig, estpt, 
            chargesignin, layersid,  
            towerid, pca::FLOATPT, pca::RPHI))
      {
        std::cerr << chargesignin << " " <<  
          layersid << " " << towerid << " " << 
          estpt << " " << etaorig << std::endl;
        return false;
      }
    }
    else 
    {
      chargesignin = 0;

      if (!import_pca_const (allconst, cmtx_c, qvec_c, 
            amtx_c, kvec_c, etaorig, pt, 
            chargesignin, pslayersid, 
            towerid, pca::FLOATPT, pca::RZ))
        return false;
    }

    pcamat_to_armamat (cmtx_c, cmtx);
    pcamat_to_armamat (amtx_c, amtx);
    pcamat_to_armarowvec (qvec_c, q);
    pcamat_to_armarowvec (kvec_c, kvct);

    //std::cout << "Here" << std::endl;

    for (int i=0; i<(int)paramdim; ++i)
    {
      double *ptr = paraptr[i];

      ptr[b] = q(i);
      for (int k=0; k<coordim_; ++k)
        ptr[b] += cmtx(i,k)*coord(b,k);
    }

    double chi2 = 0.0;

    for (int i=0; i<coordim_-paramdim_; ++i)
    {
      double val = 0.0;
      
      /* sum over j */
      for (int j=0; j<coordim_; ++j)
        val += amtx(i,j) * coord(b,j);

      /* wrong sign */
      val -= kvct(i);

      chi2 += val*val;
    }

    chi2values(b) = chi2;
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

bool pcafitter::compute_pca_constants (
    const arma::mat & param, 
    const arma::mat & coord,
    arma::mat & cmtx, 
    arma::rowvec & q, 
    arma::mat & vinv,
    arma::mat & amtx,
    arma::rowvec & kivec,
    arma::rowvec & coordm,
    int verbositylevel)
{

  if (this->get_paramdim() <= 0)
    return false; 
  
  if (this->get_coordim() <= 0)
    return false;

  // ordered 
  arma::vec eigval;
  // by row or by column ?
  arma::mat eigvec;

  if (verbositylevel == 1)
    std::cout << "Compute covariance mtx" << std::endl;

  arma::mat v = arma::zeros<arma::mat>(coordim_, coordim_);
  arma::vec eigvaltmp = arma::zeros<arma::vec>(coordim_);

  eigvec = arma::zeros<arma::mat>(coordim_, coordim_);
  eigval = arma::zeros<arma::vec>(coordim_);

  v = arma::cov(coord);
  //v = arma::cor(coord);

  if (verbositylevel == 1)
    std::cout << "Eigensystem" << std::endl;

  arma::eig_sym(eigvaltmp, eigvec, v);
  //arma::mat coeff;
  //arma::mat score;
  //arma::princomp(coeff, score, v);

  /* stored column by column  
  for (int i=0; i<(int)eigvaltmp.n_rows; ++i)
    std::cout << i << " ==> " << 
      eigvaltmp(i) << std::endl;

  for (int i=0; i<(int)eigvec.n_rows; ++i)
  {
    std::cout << i << " eigvector                     ==> ";
    for (int j=0; j<(int)eigvec.n_rows; ++j)
      std::cout << eigvec(j,i) << " ";
    std::cout << std::endl;

    std::cout << i << " eigvector * eigvalue          ==> ";
    for (int j=0; j<(int)eigvec.n_rows; ++j)
      std::cout << eigvaltmp(i) * eigvec(j,i) << " ";
    std::cout << std::endl;

    arma::vec resvec = arma::zeros<arma::vec>(this->get_coordim());
    std::cout << i << " eigvector * mat               ==> ";
    for (int j=0; j<(int)eigvec.n_rows; ++j)
    {
      double lam = 0.0;
      for (int k=0; k<(int)eigvec.n_rows; ++k)
        lam += v(j, k) * eigvec(k,i);
      std::cout << lam << " ";
      resvec(i) = lam;
    }
    std::cout << std::endl;
  }

  arma::mat deltas = eigvec.t() * v * eigvec;
  deltas.print(std::cout);
  */

  /* compute A matrix */
  amtx.resize(this->get_coordim()-this->get_paramdim(), 
      this->get_coordim());
  for (int i=0; i<this->get_coordim()-this->get_paramdim(); ++i)
    for (int j=0; j<this->get_coordim(); ++j)
      amtx(i, j) = eigvec(j, i) / sqrt(eigvaltmp(i));

  /*
  coord.print (std::cout);
  for (int i=0; i<(int)coord.n_cols; ++i)
    std::cout << arma::mean(coord.col(i)) << " ";
  std::cout << std::endl;
  */

  /* compute kivec */
  kivec.resize(this->get_coordim()-this->get_paramdim());
  for (int i=0; i<this->get_coordim()-this->get_paramdim(); ++i)
  {
    kivec(i) = 0.0e0;
    for (int k=0; k<this->get_coordim(); ++k)
      kivec(i) += amtx(i, k) * arma::mean(coord.col(k));
  }
 
  for (int i=0; i<this->get_coordim(); ++i)
    eigval(i) = eigvaltmp(this->get_coordim()-i-1);

  double totval = 0.0e0;
  for (int i=0; i<this->get_coordim(); ++i)
    totval += eigval(i);

  if (verbositylevel == 1)
    std::cout << "Eigenvalues: " << std::endl;

  double totvar = 0.0e0; 
  for (int i=0; i<this->get_coordim(); ++i)
  {
    if (i < this->get_paramdim())
      totvar += 100.0e0*(eigval(i)/totval);

    if (verbositylevel == 1)
      std::cout << i+1 << " ==> " << 100.0e0*(eigval(i)/totval) 
              << "% value: " << eigval(i) <<  std::endl;
  }

  if (verbositylevel == 1)
    std::cout << this->get_paramdim() << " eigenvalues: " << totvar << std::endl;

  v = arma::zeros<arma::mat>(coordim_,coordim_);
  v = arma::cov(coord);

#ifdef DEBUG
  /* correlation matrix ricorda su dati standardizzati coincide con la matrice 
   *   di covarianza : 
   *   z = x - <x> / sigma */
  arma::mat corr = arma::zeros<arma::mat>(coordim_,coordim_);
  
  for (int i=0; i<(coordim_); ++i)
    for (int j=0; j<(coordim_); ++j)
      corr(i,j) = v(i,j) / sqrt(v(i,i)*v(j,j));
  
  std::cout << "Correlation matrix: " << std::endl;
  std::cout << corr;
#endif

  //arma::mat vi = arma::zeros<arma::mat>(coordim_,coordim_);
  vinv = arma::inv(v); 
  //vi = v.i();

#ifdef DEBUG
  std::cout << "inverse by cov matrix: " << std::endl;
  std::cout << v * vinv;
#endif

  // and so on ...
  arma::mat hcap = arma::zeros<arma::mat>(coordim_,paramdim_);
  arma::rowvec paramm = arma::zeros<arma::rowvec>(paramdim_);
  coordm.resize(coordim_);
  coordm = arma::zeros<arma::rowvec>(coordim_);

  coordm = mean(coord, 0);
  paramm = mean(param, 0);

  /*
  for (int i=0; i<(int)coord.n_cols; ++i)
  {
    for (int j=0; j<(int)coord.n_rows; ++j)
      coordm(i) += coord(j,i);
    coordm(i) /= (double) coord.n_rows;
  } 
  std::cout << coordm << std::endl;

  for (int i=0; i<(int)coord.n_rows; ++i)
  {
    paramm(PTIDX) += param(i, PTIDX);
    paramm(PHIIDX) += param(i, PHIIDX);
    paramm(D0IDX) += param(i, D0IDX);
    paramm(THETAIDX) += param(i, THETAIDX);
    paramm(Z0IDX) += param(i, Z0IDX);
  }

  paramm(PTIDX) /= (double) coord.n_rows;
  paramm(PHIIDX) /= (double) coord.n_rows;
  paramm(D0IDX) /= (double) coord.n_rows;
  paramm(THETAIDTHETAIDX(double) coord.n_rows;
  paramm(Z0IDX) /= (double) coord.n_rows;

  std::cout << paramm(PTIDX) << " " << mean(param, 0) << std::endl;
  */

  hcap = arma::cov(coord, param);

  //std::cout << hcap << std::endl;

  for (int i=0; i<paramdim_; ++i)
    for (int l=0; l<coordim_; ++l)
      for (int m=0; m<coordim_; ++m)
        cmtx(i,l) += vinv(l,m) * hcap (m,i);

  for (int i=0; i<paramdim_; ++i)
  {
    q(i) = paramm(i);
    for (int l=0; l<coordim_; ++l)
      q(i) -= cmtx(i,l)*coordm(l);
  }

  if (useintbitewise_)
  {
    //cmtx and q constants are still in float point. Multiply them to bring in integer range.
    //They are casted as int32_t while reading from pca::read_armmat and pca::read_armvct in fit pca step.
    
    for (int i=0; i<paramdim_; ++i)
    {
      for (int l=0; l<coordim_; ++l)
      {
        //R-Z Factors
        //if (i == 0) cmtx(i,l) *= 15000000;
        //else if (i == 1) cmtx(i,l) *= 1000000;
        //if (i == 0) cmtx(i,l) *= 30000;
        //else if (i == 1) cmtx(i,l) *= 2000;
        //R-Phi Fcators
        if (i == 0) 
          cmtx(i,l) *= 60000;  //if (i == 0) cmtx(i,l) *= 60000;
        else if (i == 1) 
          cmtx(i,l) *= 8000; //else if (i == 1) cmtx(i,l) *= 8000;
        
        if (l%2 == 1) 
          cmtx(i, l) *= 50; //Further scale up the r corresponding constants
      }
    }
    
    for (int i=0; i<paramdim_; ++i)
    {
      //R-Z Factors
      //if (i == 0) q(i) *=15000000;
      //else if (i == 1) q(i) *= 1000000;
      //R-Phi Factors
      if (i == 0)      
        q(i) = (q(i)*60000000 - 14000000);  //if (i == 0)      q(i) = (q(i)*60000000 - 14000000);//- for Pt+ & + for Pt-
      else if (i == 1) 
        q(i) = (q(i)*8000000 - 16000000);
    }
    
    // TODO the same for the chi2 related matrix and kvec
  }

  if (useintbitewise_)
  {
    arma::imat icmtx = arma::conv_to<arma::imat>::from(cmtx);
    std::cout << "C matrix: " << std::endl;
    std::cout << icmtx;
  
    arma::ivec iq = arma::conv_to<arma::ivec>::from(q);
    std::cout << "Q vector: " << std::endl;
    for (int i=0; i<paramdim_; ++i)
      std::cout << iq(i) << std::endl;
  }
  else
  {
    std::cout << "C matrix: " << std::endl;
    std::cout << cmtx;
  
    std::cout << "Q vector: " << std::endl;
    for (int i=0; i<paramdim_; ++i)
      std::cout << q(i) << std::endl;
  }

  return true;
}

const std::string & pcafitter::get_errmsg () const
{
  return errmsg_;
}

int pcafitter::get_errnum() const
{
  return errnum_;
} 

void pcafitter::set_useintbitewise (bool in)
{
  useintbitewise_ = in;
}

bool pcafitter::get_useintbitewise () const
{
  return useintbitewise_;
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
