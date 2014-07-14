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
  double ** param_mtx = new double *[num_of_ent];
  double ** coord_mtx = new double *[num_of_ent];
  for (int i = 0; i < num_of_ent; ++i)
  {
    coord_mtx[i] = new double[3*COORDIM]; 
    param_mtx[i] = new double[PARAMDIM];
  }

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
      mytfp >> coord_mtx[i][j*3] >> 
               coord_mtx[i][j*3+1] >> 
               coord_mtx[i][j*3+2] >> 
               a >> b >> c; 
    }
    mytfp >> param_mtx[i][0] >> 
             param_mtx[i][1] >> 
             param_mtx[i][2] >> 
             param_mtx[i][3] >> 
             param_mtx[i][4];
  }

  mytfp.close();

#ifdef DEBUG
  for (int i = 0; i < num_of_ent; ++i)
  {
    for (int j = 0; j < COORDIM; ++j)
    {
      std::cout << coord_mtx[i][j*3] << " " <<
                   coord_mtx[i][j*3+1] << " " <<
                   coord_mtx[i][j*3+2] << std::endl;
    }
    std::cout << param_mtx[i][0] << " " <<
                 param_mtx[i][1] << " " <<
                 param_mtx[i][2] << " " <<
                 param_mtx[i][3] << " " <<
                 param_mtx[i][4] << std::endl;
  }
#endif

  arma::mat v = arma::zeros<arma::mat>(num_of_ent,3*COORDIM);
  for (int i = 0; i < num_of_ent; ++i)
    for (int j = 0; j < 3*COORDIM; ++j)
      v(i,j) = coord_mtx[i][j];

  arma::mat coeff;
  arma::mat score;
  arma::vec latent;
  arma::vec tsquared;

  arma::princomp(coeff, score, latent, tsquared, v);

  std::cout << "Eigenvalue: " << std::endl;
  std::cout << latent;

  double sum = 1.0e0;
  arma::mat coordm = arma::zeros<arma::mat>(3*COORDIM);
  arma::mat hca = arma::zeros<arma::mat>(3*COORDIM,3*COORDIM);

  for (int l=0; l<num_of_ent; ++l) 
  {
    sum += 1.0e0;
    for (int i=0; i<(3*COORDIM); ++i)
      coordm(i) += (coord_mtx[l][i]-coordm(i))/sum;

    for (int i=0; i<(3*COORDIM); ++i)
    {
      for (int j=0; j<(3*COORDIM); ++j)
      {
        hca(i,j) += ((coord_mtx[l][i] - coordm(i))*
                     (coord_mtx[l][j] - coordm(j))-
                     (sum-1.0e0)*hca(i,j)/sum)/(sum-1.0e0);
      }
    }
  }

  /* correlation matrix ricorda su dati standardizzati coincide con la matrice 
 *   di covarianza : 
 *   z = x -<x> / sigma */
  arma::mat corr = arma::zeros<arma::mat>(3*COORDIM,3*COORDIM);

  for (int i=0; i<(3*COORDIM); ++i)
    for (int j=0; j<(3*COORDIM); ++j)
      corr(i,j) = hca(i,j) / sqrt(hca(i,i)*hca(j,j));

  std::cout << "Correlation matrix: " << std::endl;
  std::cout << corr;

  for (int i=0; i<(3*COORDIM); ++i)
    for (int j=i+1; j<(3*COORDIM); ++j)
      if (hca(i,j) != hca(j,i))
        std::cout << i << " " << j << " " << 
          hca(i,j) << " ERROR" << std::endl;;

  arma::vec eigval = arma::zeros<arma::mat>(3*COORDIM);
  arma::mat eigvec = arma::zeros<arma::mat>(3*COORDIM,3*COORDIM);

  arma::eig_sym(eigval, eigvec, hca);

  double totval = 0.0e0;
  for (int i=0; i<(3*COORDIM); ++i)
    totval += eigval(i);

  std::cout << "Eigenvalues: " << std::endl;
  int j = 1;
  double totvar = 0.0e0; 
  for (int i=(3*COORDIM-1); i>=0; --i)
  {
    if (j <= PARAMDIM)
      totvar += 100.0e0*(eigval(i)/totval);
    ++j;

    std::cout << i+1 << " ==> " << 100.0e0*(eigval(i)/totval) 
      << " ==> " << eigval(i) <<  std::endl;

  }
  std::cout << "PARAMDIM eigenvalues: " << totvar << std::endl;

  arma::mat hcai = arma::zeros<arma::mat>(3*COORDIM,3*COORDIM);
  hcai = hca.i(); 

//#ifdef DEBUG
  //std::cout << hca * hcai ;
//#endif
  
  // and so on ...
  arma::mat paramm = arma::zeros<arma::mat>(PARAMDIM);
  arma::mat hcap = arma::zeros<arma::mat>(3*COORDIM,PARAMDIM);

  coordm.fill(0.0e0);
  sum = 1.0e0;
  
  for (int l=0; l<num_of_ent; ++l) 
  {
    sum += 1.0e0;
    for (int i=0; i<(3*COORDIM); ++i)
      coordm(i) += (coord_mtx[l][i]-coordm(i))/sum;

    for (int i=0; i<PARAMDIM; ++i)
      paramm(i) += (param_mtx[l][i]-paramm(i))/sum;

    for (int i=0; i<(3*COORDIM); ++i)
    {
      for (int j=0; j<PARAMDIM; ++j)
      {
        hcap(i,j) += ((coord_mtx[l][i] - coordm(i))*
                      (param_mtx[l][j] - paramm(j))-
                      (sum-1.0e0)*hcap(i,j)/sum)/(sum-1.0e0);
      }
    }
  }

  /* correlation matrix
  double pstdev[PARAMDIM];
  for (int i=0; i<PARAMDIM; ++i)
  {
    arma::running_stat<double> stats;

    for (int l=0; l<num_of_ent; ++l) 
      stats(param_mtx[l][i]);

#ifdef DEBUG
    std::cout << "mean   = " << stats.mean() << std::endl;
    std::cout << "         "  << paramm(i) << std::endl;
    std::cout << "stdev  = " << stats.stddev()  << std::endl;
#endif

    pstdev[i] = stats.stddev();
  }

  for (int i=0; i<(3*COORDIM); ++i)
    for (int j=0; j<PARAMDIM; ++j)
      hcap(i,j) = hcap(i,j) / (cstdev[i]*pstdev[j]);
  */

  arma::mat cmtx = arma::zeros<arma::mat>(PARAMDIM,3*COORDIM);

  for (int i=0; i<PARAMDIM; ++i)
    for (int l=0; l<(3*COORDIM); ++l)
      for (int m=0; m<(3*COORDIM); ++m)
        cmtx(i,l) += hcai(l,m) * hcap (m,i);

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

  arma::mat a = arma::zeros<arma::mat>((3*COORDIM),(3*COORDIM));
  for (int i=0; i<(3*COORDIM); ++i)
    for (int j=0; j<(3*COORDIM); ++j)
      a(i,j) = eigvec(i,j)/sqrt(eigval(i));

  arma::mat k = arma::zeros<arma::mat>(3*COORDIM);
  for (int i=0; i<(3*COORDIM); ++i)
    for (int j=0; j<(3*COORDIM); ++j)
      k(i) += a(i,j)*coordm(j);


  //test back
  arma::running_stat<double> chi2stats;
  arma::running_stat<double> pc[PARAMDIM];
  for (int l=0; l<num_of_ent; ++l)
  {

    for (int i=0; i<PARAMDIM; ++i)
    {
      double p = q(i);
      for (int k=0; k<(3*COORDIM); ++k)
        p += cmtx(i,k)*coord_mtx[l][k];
    
      pc[i](fabs(p - param_mtx[l][i])/(fabs(p + param_mtx[l][i])/2.0));
    }

    /* chi**2 */
    double chi2 = 0.0e0;
    for (int i=0; i<(3*COORDIM)-PARAMDIM; ++i)
    {
      double v = k(i);
    
      for (int j=0; j<(3*COORDIM); ++j)
        v += a(i,j) * coord_mtx[l][j];
      
      chi2 += (v*v);
    } 
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
      stats(param_mtx[l][i]);

    std::cout << "  mean   = " << stats.mean() << std::endl;
    std::cout << "  stdev  = " << stats.stddev()  << std::endl;
    std::cout << "  min    = " << stats.min() << std::endl;
    std::cout << "  max    = " << stats.max() << std::endl;
  }

  for (int i = 0; i < num_of_ent; ++i)
  {
    delete(coord_mtx[i]);
    delete(param_mtx[i]);
  }
  delete(coord_mtx);
  delete(param_mtx);

  return 0;
}
