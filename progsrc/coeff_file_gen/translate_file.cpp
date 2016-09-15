#include <iostream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <fstream>
#include <bitset>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

#include "common.h"

const int mult_factor = 1e6;
const int const_mult_factor = mult_factor*1024;
const int chisq_mult_factor = 1e4;
const int chisq_const_mult_factor = chisq_mult_factor*1024;
const int const_w = 25;
const int add_const_w = 36;

using namespace pca;
using namespace std;

int main (int argc, char *argv[]) 
{
  if (argc != 2)
  {
    std::cerr << "usage: " << argv[0] << " pca_const_file.txt " << std::endl;
    return 1;
  }

  std::vector< matrixpcaconst<double> > all_constants;
  bool a = read_pcaconst_from_file<double>(all_constants, argv[1]);

  // write rphi param
  std::vector<double> constants_rphi_cmtx;
  for (matrixpcaconst<double> c : all_constants) 
  {
    if (c.get_plane_type() == plane_type::RPHI)
    {
      if (c.get_const_type() == const_type::CMTX) 
      {
	for (int i = 0; i < c.n_rows(); ++i)
        {
	  for (int j = 0; j < c.n_cols(); ++j)
          {
	    if ( (get_missing_layer(c) >= 0) &&
		 (get_missing_layer(c) <= (j/2))
		 )
	      constants_rphi_cmtx.push_back(((get_missing_layer(c) == (j/2))?0:c.element(i, j-2)));
	    else
	      constants_rphi_cmtx.push_back(c.element(i, j));
          }
        }
      }
    }
  }

  std::vector<double> constants_rphi_qvec;
  for (matrixpcaconst<double> c : all_constants) 
  {
    if (c.get_plane_type() == plane_type::RPHI)
    {
      if (c.get_const_type() == const_type::QVEC) 
      {
	for (int i = 0; i < c.n_cols(); ++i)
        {
	  constants_rphi_qvec.push_back(c.element(0, i));
        }
      }
    }
  }

  return 0;
}
