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

#include <pcafitter_private.hpp>

void pcafitter::readarmmat (const char * fname, arma::mat & cmtx)
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

void pcafitter::readarmvct (const char * fname, arma::rowvec & q)
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

void pcafitter::writearmmat (const char * fname, arma::mat & cmtx)
{
  std::ofstream myfilec(fname, std::ios::binary);
  myfilec.write((const char*)&(cmtx.n_rows), sizeof(cmtx.n_rows));
  myfilec.write((const char*)&(cmtx.n_cols), sizeof(cmtx.n_cols));
  for (int i=0; i<(int)cmtx.n_rows; i++)
    for (int j=0; j<(int)cmtx.n_cols; j++)
      myfilec.write((const char*)&(cmtx(i, j)), sizeof(cmtx(i, j)));
  myfilec.close();
}

void pcafitter::writearmvct (const char * fname, arma::rowvec & q)
{
  std::ofstream myfileq(fname, std::ios::binary);
  myfileq.write((const char*)&(q.n_cols), sizeof(q.n_cols));
  for (int i=0; i<(int)q.n_cols; i++)
    myfileq.write((const char*)&(q(i)), sizeof(q(i)));
  myfileq.close();
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

void pcafitter::writetofile (const char * fname, 
    const arma::rowvec & vec) 
{
  std::ofstream myfile(fname);
    
  for (int i=0; i<(int)vec.n_cols; i++)
    myfile << i << " " << vec(i) << std::endl;
    
  myfile.close();
}

void pcafitter::readingfromfile (const char * filename, 
    arma::rowvec & ptin, arma::rowvec & phiin, 
    arma::rowvec & d0in, arma::rowvec & etain, 
    arma::rowvec & z0in, arma::mat & coordin, 
    arma::mat & layer, arma::mat & ladder, 
    arma::mat & module, 
    std::map<std::string, int> & subsectors, 
    std::map<std::string, int> & subladders,
    std::vector<std::string> & subsectorslist,
    std::vector<std::string> & subladderslist,
    int num_of_ent)
{

  std::string line;
  std::ifstream mytfp;
  mytfp.open (filename, std::ios::in);

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
    
    for (int j = 0; j < COORDIM; ++j)
    {
      int a, b, c;
      mytfp >> coordin(i, j*3) >> 
               coordin(i, j*3+1) >> 
               coordin(i, j*3+2) >> 
               a >> b >> c; 
    
      layer(i, j) = a;
      ladder(i, j) = b;
      module(i, j) = c;
      
      osss << std::setw(2) << layer(i, j);
      osss << std::setw(2) << ladder(i, j);
      if (j != COORDIM-1)
        osss<<"-";

      ossl << std::setw(2) << layer(i, j);
      ossl << std::setw(2) << ladder(i, j);
      ossl << std::setw(2) << module(i, j);
      if (j != COORDIM-1)
        ossl<<"-";

    }

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
    
    mytfp >> ptin(i) >> 
             phiin(i) >> 
             d0in(i) >> 
             etain(i) >> 
             z0in(i);
  }

  mytfp.close();
}

void pcafitter::select_bigger_sub (
    const std::map<std::string, int> & sublist, 
    bool verbose, int & maxnumber, std::string & slctsubsec ) 
{
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
    const arma::rowvec & ptin, 
    const arma::rowvec & phiin, 
    const arma::rowvec & d0in, 
    const arma::rowvec & etain, 
    const arma::rowvec & z0in, 
    const arma::mat & coordin, 
    arma::rowvec & pt, 
    arma::rowvec & phi,
    arma::rowvec & d0,
    arma::rowvec & eta,
    arma::rowvec & z0,
    arma::mat & coord)
{
  int dim = 0;

  std::vector<std::string>::const_iterator it = sublist.begin();
  for (; it != sublist.end(); ++it)
    if (*it == slctsubsec)
      dim++;

  coord.set_size(dim,3*COORDIM);
  pt.set_size(dim);
  phi.set_size(dim);
  d0.set_size(dim);
  eta.set_size(dim);
  z0.set_size(dim);

  int k = 0;
  it = sublist.begin();
  for (int i=0; it != sublist.end(); ++it, ++i)
  {
    if (*it == slctsubsec)
    {
      for (int j = 0; j < 3*COORDIM; ++j)
        coord(k,j) = coordin(i,j);

      pt(k) = ptin(i);
      phi(k) = phiin(i);
      d0(k) = d0in(i);
      eta(k) = etain(i);
      z0(k) = z0in(i);

      ++k;
    }
  }
}

void pcafitter::compute_pca_constants (
    const arma::rowvec & pt, 
    const arma::rowvec & phi,
    const arma::rowvec & d0,
    const arma::rowvec & eta,
    const arma::rowvec & z0,
    const arma::mat & coord,
    arma::mat & cmtx, 
    arma::rowvec & q)
{
  arma::mat v = arma::zeros<arma::mat>(3*COORDIM,3*COORDIM);
  v = arma::cov(coord);

#ifdef DEBUG
  /* correlation matrix ricorda su dati standardizzati coincide con la matrice 
   *   di covarianza : 
   *   z = x -<x> / sigma */
  arma::mat corr = arma::zeros<arma::mat>(3*COORDIM,3*COORDIM);
  
  for (int i=0; i<(3*COORDIM); ++i)
    for (int j=0; j<(3*COORDIM); ++j)
      corr(i,j) = v(i,j) / sqrt(v(i,i)*v(j,j));
  
  std::cout << "Correlation matrix: " << std::endl;
  std::cout << corr;
#endif

  arma::mat vi = arma::zeros<arma::mat>(3*COORDIM,3*COORDIM);
  vi = v.i(); 

#ifdef DEBUG
  std::cout << "inverse by cov matrix: " << std::endl;
  std::cout << v * vi ;
#endif
  
  // and so on ...
  arma::mat hcap = arma::zeros<arma::mat>(3*COORDIM,PARAMDIM);
  arma::rowvec paramm = arma::zeros<arma::rowvec>(PARAMDIM);
  arma::rowvec coordm = arma::zeros<arma::rowvec>(3*COORDIM);
  double sum = 1.0e0;

  for (int l=0; l<(int)coord.n_rows; ++l) 
  {
    sum += 1.0e0;
    for (int i=0; i<(3*COORDIM); ++i)
      coordm(i) += (coord(l,i)-coordm(i))/sum;
    
    paramm(0) += (pt(l)-paramm(0))/sum;
    paramm(1) += (phi(l)-paramm(0))/sum;
    paramm(2) += (d0(l)-paramm(0))/sum;
    paramm(3) += (eta(l)-paramm(0))/sum;
    paramm(4) += (z0(l)-paramm(0))/sum;
    
    for (int i=0; i<(3*COORDIM); ++i)
    {
      hcap(i,0) += ((coord(l,i) - coordm(i))*
                    (pt(l) - paramm(0))-
                    (sum-1.0e0)*hcap(i,0)/sum)/(sum-1.0e0);

      hcap(i,1) += ((coord(l,i) - coordm(i))*
                    (phi(l) - paramm(1))-
                    (sum-1.0e0)*hcap(i,1)/sum)/(sum-1.0e0);

      hcap(i,2) += ((coord(l,i) - coordm(i))*
                    (d0(l) - paramm(2))-
                    (sum-1.0e0)*hcap(i,2)/sum)/(sum-1.0e0);

      hcap(i,3) += ((coord(l,i) - coordm(i))*
                    (eta(l) - paramm(3))-
                    (sum-1.0e0)*hcap(i,3)/sum)/(sum-1.0e0);

      hcap(i,4) += ((coord(l,i) - coordm(i))*
                    (z0(l) - paramm(4))-
                    (sum-1.0e0)*hcap(i,4)/sum)/(sum-1.0e0);

    }
  }

  for (int i=0; i<PARAMDIM; ++i)
    for (int l=0; l<(3*COORDIM); ++l)
      for (int m=0; m<(3*COORDIM); ++m)
        cmtx(i,l) += vi(l,m) * hcap (m,i);

#ifdef DEBUG
  std::cout << "C matrix: " << std::endl;
  std::cout << cmtx;
#endif 

  for (int i=0; i<PARAMDIM; ++i)
  {
    q(i) = paramm(i);
    for (int l=0; l<(3*COORDIM); ++l)
      q(i) -= cmtx(i,l)*coordm[l];
  }

#ifdef DEBUG
  std::cout << "Q vector: " << std::endl;
  for (int i=0; i<PARAMDIM; ++i)
    std::cout << q(i) << std::endl;
#endif

}
