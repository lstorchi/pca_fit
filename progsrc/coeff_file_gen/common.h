#include "pcaconst.hpp"

// questo e' il rapporto fra scaling factor di r e phi so pca::risf/pca::pisf
double rpdval_get ();
double rzdval_get ();

int check_val(double val, int width, std::string mstr); 
int get_missing_layer(pca::matrixpcaconst<double> c);
