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
 

  double sum = 1.0e0;
  double coordm[3*COORDIM];
  double hc[3*COORDIM][3*COORDIM];

  std::fill_n (coordm, (3*COORDIM), 0.0e0);
  std::fill_n (&(hc[0][0]), (3*COORDIM)*(3*COORDIM), 0.0e0);

  for (int l=0; l<num_of_ent; ++l) 
  {
    sum += 1.0e0;
    for (int i=0; i<(3*COORDIM); ++i)
      coordm[i] += (coord_mtx[l][i]-coordm[i])/sum;

    for (int i=0; i<(3*COORDIM); ++i)
    {
      for (int j=0; j<(3*COORDIM); ++j)
      {
        hc[i][j] += ((coord_mtx[l][i] - coordm[i])*
                     (coord_mtx[l][j] - coordm[j])-
                     (sum-1.0e0)*hc[i][j]/sum)/(sum-1.0e0);
      }
    }
  }

  double cstdev[3*COORDIM];

  for (int i=0; i<(3*COORDIM); ++i)
  {
    arma::running_stat<double> stats;

    for (int l=0; l<num_of_ent; ++l) 
      stats(coord_mtx[l][i]);

    std::cout << "mean = " << stats.mean() << std::endl;
    std::cout << "       "  << coordm[i] << std::endl;

    std::cout << "var  = " << stats.stddev()  << std::endl;

    cstdev[i] = stats.stddev();
  }

  for (int l=0; l<num_of_ent; ++l) 
  {
    for (int i=0; i<(3*COORDIM); ++i)
      coordm[i] += (coord_mtx[l][i]-coordm[i])/sum;

    for (int i=0; i<(3*COORDIM); ++i)
    {
      for (int j=0; j<(3*COORDIM); ++j)
      {
        hc[i][j] += ((coord_mtx[l][i] - coordm[i])*
                     (coord_mtx[l][j] - coordm[j])-
                     (sum-1.0e0)*hc[i][j]/sum)/(sum-1.0e0);
      }
    }
  }

/*
  for (int i=0; i<(3*COORDIM); ++i)
  {
    for (int j=0; j<(3*COORDIM); ++j)
    {
      hc[i][j] = hc[i][j]/(cstdev[i]*cstdev[j]);
    }
  }
*/

  for (int i=0; i<(3*COORDIM); ++i)
    for (int j=i+1; j<(3*COORDIM); ++j)
      if (hc[i][j] != hc[j][i])
        std::cout << i << " " << j << " " << 
          hc[i][j] << " ERROR" << std::endl;;

#ifdef DEBUG
  double coordm_d[3*COORDIM];
  std::fill_n (coordm_d, (3*COORDIM), 0.0e0);
  for (int i=0; i<(3*COORDIM); ++i)
  {
    for (int l=0; l<num_of_ent; ++l)
      coordm_d[i] += coord_mtx[l][i];
    coordm_d[i] = coordm_d[i] / (double)num_of_ent;
    //std::cout << coordm_d[i] << " " << coordm[i] << std::endl;
  }

  double hc_d[3*COORDIM][3*COORDIM];
  std::fill_n (&(hc_d[0][0]), (3*COORDIM)*(3*COORDIM), 0.0e0);
  for (int i=0; i<(3*COORDIM); ++i)
  {
    for (int j=0; j<(3*COORDIM); ++j)
    {
      for (int l=0; l<num_of_ent; ++l)
      {
        hc_d[i][j] += (coord_mtx[l][i] - coordm[i])*
                      (coord_mtx[l][j] - coordm[j]);
      }
      hc_d[i][j] = hc_d[i][j] / (double)num_of_ent;
      //std::cout << hc_d[i][j] << " " << hc[i][j] << std::endl;
    }
  }

  arma::mat hca_d = arma::zeros<arma::mat>(3*COORDIM,3*COORDIM);
  for (int i=0; i<(3*COORDIM); ++i)
    for (int j=0; j<(3*COORDIM); ++j)
      hca_d(i,j) = hc_d[i][j];

  arma::vec eigval_d;
  arma::mat eigvec_d;

  arma::eig_sym(eigval_d, eigvec_d, hca_d);

  for (int i=(3*COORDIM-1); i>=0; --i)
  {
    std::cout << i+1 << " ==> " << eigval_d(i) << std::endl;
  }
  std::cout << "====" << std::endl;

#endif

  arma::mat hca = arma::zeros<arma::mat>(3*COORDIM,3*COORDIM);
  for (int i=0; i<(3*COORDIM); ++i)
    for (int j=0; j<(3*COORDIM); ++j)
      hca(i,j) = hc[i][j];

  arma::vec eigval;
  arma::mat eigvec;

  arma::eig_sym(eigval, eigvec, hca);

#ifdef DEBUG

  double xnew[3*COORDIM];

  std::fill_n (xnew, (3*COORDIM), 0.0e0);

  for (int l=0; l<num_of_ent; ++l) 
  {
    for (int i=0; i<(3*COORDIM); ++i)
    {
      xnew[i] = 0.0e0;
      for (int j=0; j<(3*COORDIM); ++j)
        xnew[i] += eigvec(i, j) * (coord_mtx[l][i] - coordm[i]);
    }
    
    for (int i=0; i<(3*COORDIM); ++i)
    {
      std::cout << i << " --> " << xnew[i] << std::endl;
    }


    for (int i=0; i<PARAMDIM; ++i)
      std::cout << "      " <<  param_mtx[i][i] << std::endl;

  }

#endif

  double totval = 0.0e0;
  for (int i=0; i<(3*COORDIM); ++i)
    totval += eigval(i);

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

  //std::cout << hca * hcai ;

  //exit(1);
  
  // and so on ...
  double paramm[PARAMDIM];
  arma::mat hcap = arma::zeros<arma::mat>(3*COORDIM,PARAMDIM);

  std::fill_n(coordm, (3*COORDIM), 0.0e0 );
  std::fill_n(paramm, PARAMDIM, 0.0e0 );
  sum = 1.0e0;
  
  for (int l=0; l<num_of_ent; ++l) 
  {
    sum += 1.0e0;
    for (int i=0; i<(3*COORDIM); ++i)
      coordm[i] += (coord_mtx[l][i]-coordm[i])/sum;

    for (int i=0; i<PARAMDIM; ++i)
      paramm[i] += (param_mtx[l][i]-paramm[i])/sum;

    for (int i=0; i<(3*COORDIM); ++i)
    {
      for (int j=0; j<PARAMDIM; ++j)
      {
        hcap(i,j) += ((coord_mtx[l][i] - coordm[i])*
                      (param_mtx[l][j] - paramm[j])-
                      (sum-1.0e0)*hcap(i,j)/sum)/(sum-1.0e0);
      }
    }
  }

  double pstdev[PARAMDIM];

  for (int i=0; i<PARAMDIM; ++i)
  {
    arma::running_stat<double> stats;

    for (int l=0; l<num_of_ent; ++l) 
      stats(param_mtx[l][i]);

    std::cout << "mean = " << stats.mean() << std::endl;
    std::cout << "       "  << paramm[i] << std::endl;

    std::cout << "var  = " << stats.stddev()  << std::endl;

    pstdev[i] = stats.stddev();
  }



  /*
  for (int i=0; i<(3*COORDIM); ++i)
    for (int j=0; j<PARAMDIM; ++j)
      hcap(i,j) = hcap(i,j) / (cstdev[i]*pstdev[j]);
  */

  arma::mat cmtx = arma::zeros<arma::mat>(PARAMDIM,3*COORDIM);

  for (int i=0; i<PARAMDIM; ++i)
    for (int l=0; l<(3*COORDIM); ++l)
      for (int m=0; m<(3*COORDIM); ++m)
        cmtx(i,l) += hcai(l,m) * hcap (m,i);
 
  std::cout << "C matrix: " << std::endl;
  std::cout << cmtx;

  double q[PARAMDIM];
  std::fill_n(q, PARAMDIM, 0.0e0 );

  //std::cout << cmtx;

  for (int i=0; i<PARAMDIM; ++i)
  {
    q[i] = paramm[i];
    for (int l=0; l<(3*COORDIM); ++l)
      q[i] -= cmtx(i,l)*coordm[l];
  }

  std::cout << "Q vector: " << std::endl;
  for (int i=0; i<PARAMDIM; ++i)
    std::cout << q[i] << std::endl;

  //test back
  for (int i=0; i<PARAMDIM; ++i)
  {
    double p = q[i];
    for (int l=0; l<(3*COORDIM); ++l)
      p += cmtx(i,l)*coord_mtx[0][l];
    
    std::cout << "P " << i+1 << " ==> " << p << std::endl;
    std::cout << "  " << param_mtx[0][i] << std::endl;
  }

  // documento Annovi
  // calcolo matrice di correlazione traccie HC 
  // diagonalizzo HC e determino A, matrice 5 autovettori principali (5 componenti pricipali)
  // A matrice rotazione che mi permette di calcolare la traslazione usando i paamtri di tracce 
  // simulate. 

  // documento ATLAS 
  // calcolare V e data V calcolare inversa V e quindi C, matrice di rotazione 
  // data C determinare il vettore di traslazione q  
  // c plus q costanti PCA 
  // write constants in a file

  for (int i = 0; i < num_of_ent; ++i)
  {
    delete(coord_mtx[i]);
    delete(param_mtx[i]);
  }
  delete(coord_mtx);
  delete(param_mtx);

  return 0;
}
