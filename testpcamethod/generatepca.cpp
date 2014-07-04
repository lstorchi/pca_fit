#include <iostream>
#include <fstream>
#include <string>

#define ENTDIM 8

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
    coord_mtx[i] = new double[3*(ENTDIM-2)]; 
    param_mtx[i] = new double[5];
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
    std::cout << fake1 << " " << fake2 << std::endl;
    for (int j = 0; j < ENTDIM-2; ++j)
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
    for (int j = 0; j < ENTDIM-2; ++j)
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
