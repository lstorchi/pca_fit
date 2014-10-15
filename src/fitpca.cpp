#include <iostream>

int main (int argc, char ** argv)
{
  if (argc != 3)
  {
    std::cerr << "usage: " << argv[0] << " coordinatesfile constants" << std::endl;
    return 1;
  }

  // leggere file coordinate tracce e file costanti PCA
  // N righe di 9 double sono le coordinate
  // matrice C e vettore q sono le costanti

  // calcolare i parametri per ogni traccia
  

  return 0;
}
