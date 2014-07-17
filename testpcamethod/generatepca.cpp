#include <iostream>
#include <fstream>
#include <string>

// Loriano: let's try Armadillo quick code 
#include <armadillo>

#define ENTDIM 8
#define COORDIM (ENTDIM-2)
#define PARAMDIM 5

namespace 
{
  int numofline (const char * fname) 
  { 
    int number_of_lines = 0;
    std::string line;
    std::ifstream myfile(fname);
    
    while (std::getline(myfile, line))
      ++number_of_lines;

    myfile.close();
                            
    return number_of_lines;
  }
}

int main (int argc, char ** argv)
{
  if (argc != 2)
  {
    std::cerr << "usage: " << argv[0] << " coordinatesfile " << std::endl;
    return 1;
  }

  int num_of_line = numofline(argv[1]);
  std::cout << "file has " << num_of_line << " line " << std::endl;
  int num_of_ent = (num_of_line-1)/ENTDIM;
  std::cout << "  " << num_of_ent << " entries " << std::endl;

  // non perfomante ma easy to go
  arma::mat paraminp = arma::zeros<arma::mat>(num_of_ent,PARAMDIM);
  arma::mat coordinp = arma::zeros<arma::mat>(num_of_ent,3*COORDIM);

  arma::mat tracker;
  tracker.set_size(num_of_ent,3*COORDIM);

  // leggere file coordinate tracce simulate plus parametri
  std::string line;
  std::ifstream mytfp;
  mytfp.open (argv[1], std::ios::in);

  std::getline (mytfp, line);
  //std::cout << line << std::endl;
  
  for (int i = 0; i < num_of_ent; ++i)
  {
    int fake1, fake2;
    mytfp >> fake1 >> fake2 ;
#ifdef DEBUG    
    std::cout << fake1 << " " << fake2 << std::endl;
#endif
    for (int j = 0; j < COORDIM; ++j)
    {
      int a, b, c;
      mytfp >> coordinp(i, j*3) >> 
               coordinp(i, j*3+1) >> 
               coordinp(i, j*3+2) >> 
               a >> b >> c; 

               tracker(i, j*3) = a;
               tracker(i, j*3+1) = b;
               tracker(i, j*3+2) = c;
    }
    mytfp >> paraminp(i,0) >> 
             paraminp(i,1) >> 
             paraminp(i,2) >> 
             paraminp(i,3) >> 
             paraminp(i,4);
  }

  mytfp.close();

#ifdef DEBUG
  for (int i = 0; i < num_of_ent; ++i)
  {
    for (int j = 0; j < COORDIM; ++j)
    {
      std::cout << coordinp(i, j*3) << " " <<
                   coordinp(i, j*3+1) << " " <<
                   coordinp(i, j*3+2) << std::endl;
    }
    std::cout << paraminp(i,0) << " " <<
                 paraminp(i,1) << " " <<
                 paraminp(i,2) << " " <<
                 paraminp(i,3) << " " <<
                 paraminp(i,4) << std::endl;
  }
#endif

  arma::mat param;
  arma::mat coord;

  // to be used to select inly a ladder .... 
  coord = coordinp;
  param = paraminp;
  
  // projection 
  arma::mat score;
  // ordered 
  arma::vec eigval;
  // by row or by column ?
  arma::mat eigvec;

  arma::princomp(eigvec, score, eigval, coord);

  //std::cout << score.n_rows << " " << score.n_cols << std::endl;

  //for (int i=0; i<score.n_rows; ++i)
  //  std::cout << score(i,0)  << " " 
  //    << score(i,1) << " " << score(i,2) << std::endl;

  double totval = 0.0e0;
  for (int i=0; i<(3*COORDIM); ++i)
    totval += eigval(i);

  std::cout << "Eigenvalues: " << std::endl;
  double totvar = 0.0e0; 
  for (int i=0; i<(3*COORDIM); ++i)
  {
    if (i < PARAMDIM)
      totvar += 100.0e0*(eigval(i)/totval);

    std::cout << i+1 << " ==> " << 100.0e0*(eigval(i)/totval) 
      << "% value: " << eigval(i) <<  std::endl;
  }
  std::cout << "PARAMDIM eigenvalues: " << totvar << std::endl;

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
  arma::mat paramm = arma::zeros<arma::mat>(PARAMDIM);
  arma::mat coordm = arma::zeros<arma::mat>(3*COORDIM);
  double sum = 1.0e0;
  
  for (int l=0; l<num_of_ent; ++l) 
  {
    sum += 1.0e0;
    for (int i=0; i<(3*COORDIM); ++i)
      coordm(i) += (coord(l,i)-coordm(i))/sum;

    for (int i=0; i<PARAMDIM; ++i)
      paramm(i) += (param(l,i)-paramm(i))/sum;

    for (int i=0; i<(3*COORDIM); ++i)
    {
      for (int j=0; j<PARAMDIM; ++j)
      {
        hcap(i,j) += ((coord(l,i) - coordm(i))*
                      (param(l,j) - paramm(j))-
                      (sum-1.0e0)*hcap(i,j)/sum)/(sum-1.0e0);
      }
    }
  }

  arma::mat cmtx = arma::zeros<arma::mat>(PARAMDIM,3*COORDIM);

  for (int i=0; i<PARAMDIM; ++i)
    for (int l=0; l<(3*COORDIM); ++l)
      for (int m=0; m<(3*COORDIM); ++m)
        cmtx(i,l) += vi(l,m) * hcap (m,i);

#ifdef DEBUG
  std::cout << "C matrix: " << std::endl;
  std::cout << cmtx;
#endif 

  arma::mat q = arma::zeros<arma::mat>(PARAMDIM);

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

  //test back
  arma::running_stat<double> chi2stats;
  arma::running_stat<double> pc[PARAMDIM];
  for (int l=0; l<num_of_ent; ++l)
  {

    std::cout << "Track: " << l+1 << std::endl;

    for (int i=0; i<PARAMDIM; ++i)
    {
      double p = q(i);
      for (int k=0; k<(3*COORDIM); ++k)
        p += cmtx(i,k)*coord(l,k);
      
      std::cout << "   computed "  << p << " real " << param(l,i) << std::endl;
    
      pc[i](fabs(p - param(l,i))/(fabs(p + param(l,i))/2.0));
    }

    /* chi**2 */
    double chi2 = 0.0e0;
    for (int i=0; i<(3*COORDIM); ++i)
      for (int j=0; j<(3*COORDIM); ++j)
        chi2 += (coord(l,i) - coordm(i)) * 
                vi(i, j) * (coord(l,j) - coordm(j)); 
    std::cout << "   chi2: " << chi2 << std::endl;
    chi2stats(chi2);
  }
  std::cout << "chi2 mean   = " << chi2stats.mean() << std::endl;
  std::cout << "chi2 stdev  = " << chi2stats.stddev()  << std::endl;
  std::cout << "chi2 min    = " << chi2stats.min() << std::endl;
  std::cout << "chi2 max    = " << chi2stats.max() << std::endl;

  for (int i=0; i<PARAMDIM; ++i)
  {
    std::cout << pc[i].mean() << " " << pc[i].stddev() << std::endl;

    arma::running_stat<double> stats;

    for (int l=0; l<num_of_ent; ++l) 
      stats(param(l,i));

    std::cout << "  mean   = " << stats.mean() << std::endl;
    std::cout << "  stdev  = " << stats.stddev()  << std::endl;
    std::cout << "  min    = " << stats.min() << std::endl;
    std::cout << "  max    = " << stats.max() << std::endl;
  }

  return 0;
}
