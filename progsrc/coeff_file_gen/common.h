#include "pcaconst.hpp"

// questo e' il rapporto fra scaling factor di r e phi so pca::risf/pca::pisf
const long long int rpdval = 32;

int check_val(double val, int width, std::string mstr); 
int get_missing_layer(pca::matrixpcaconst<double> c);
