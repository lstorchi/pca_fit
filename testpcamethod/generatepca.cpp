#include <iostream>

int main (int argc, char ** argv)
{
  if (argc != 2)
  {
    std::cerr << "usage: " << argv[0] << " coordinatesfile " << std::endl;
    return 1;
  }

  // leggere file coordinate tracce simulate plus parametri
  // N righe di 14 double
  

  // calcolare V e data V calcolare inversa V e quindi C, matrice di rotazione 
  
  // data C determinare il vettore di traslazione q  

  // c plus q costanti PCA 
  // write constants in a file


  return 0;
}
