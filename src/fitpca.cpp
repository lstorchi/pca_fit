#include <iostream>

#include <armadillo>
#include <cassert>

#include <pcafitter_private.hpp>

int main (int argc, char ** argv)
{
  if (argc != 4)
  {
    std::cerr << "usage: " << argv[0] << " coordinatesfile cmtx qvec" << std::endl;
    return 1;
  }

  // leggere file coordinate tracce e file costanti PCA
  // N righe di 9 double sono le coordinate
  // matrice C e vettore q sono le costanti
  
  arma::mat cmtx;
  arma::rowvec q;

  std::cout << "Read constant to file" << std::endl;
  pcafitter::readarmmat(argv[2], cmtx);
  pcafitter::readarmvct(argv[3], q);

  std::cout << "Reading data from " << argv[1] << " file " << std::endl;
  int num_of_line = pcafitter::numofline(argv[1]);
  std::cout << "file has " << num_of_line << " line " << std::endl;
  int num_of_ent = (num_of_line-1)/ENTDIM;
  std::cout << "file has " << num_of_ent << " entries " << std::endl;

  arma::rowvec pt = arma::zeros<arma::rowvec>(num_of_ent);
  arma::rowvec phi = arma::zeros<arma::rowvec>(num_of_ent);
  arma::rowvec d0 = arma::zeros<arma::rowvec>(num_of_ent);
  arma::rowvec eta = arma::zeros<arma::rowvec>(num_of_ent);
  arma::rowvec z0 = arma::zeros<arma::rowvec>(num_of_ent);

  arma::mat layer, ladder, module, coord;
  layer.set_size(num_of_ent,COORDIM);
  ladder.set_size(num_of_ent,COORDIM);
  module.set_size(num_of_ent,COORDIM);
  coord.set_size(num_of_ent,3*COORDIM);

  std::map<std::string, int> subsectors, subladders;
  std::vector<std::string> subladderslist, subsectorslist;

  // leggere file coordinate tracce simulate plus parametri
  pcafitter::readingfromfile (argv[1], pt, phi, d0, eta, z0, coord, 
      layer, ladder, module, subsectors, subladders, 
      subsectorslist, subladderslist, num_of_ent);
 


  return 0;
}
