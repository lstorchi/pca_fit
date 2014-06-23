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


  return 0;
}
