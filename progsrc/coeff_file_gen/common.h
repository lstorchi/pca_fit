#include "pcaconst.hpp"

const int mult_factor = 1e6;
const int const_mult_factor = mult_factor*1024;
const int chisq_mult_factor = 3e4;
const int chisq_const_mult_factor = chisq_mult_factor*1024;
const int const_w = 25;
const int add_const_w = 36;

int check_val(double val, int width, std::string mstr); 
int get_missing_layer(pca::matrixpcaconst<double> c);
