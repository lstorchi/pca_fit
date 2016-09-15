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

using namespace pca;
using namespace std;

int check_val(double val, int width, string mstr) 
{
  if (abs(val) > pow(2, width-1))
    cout << "overflow at " << mstr << endl;
}

int get_missing_layer(matrixpcaconst<double> c) 
{
  int missing_layer;
  missing_layer = -1;

  string layers_str;
  layers_str = c.get_layersids();
  
  if (c.get_plane_type() == plane_type::RPHI) 
  {
    for (int i = 0; i < 6; ++i)
    {
      if ( (layers_str[i*2] != ('5'+i) && (i != 5) ) ||
	   ( (layers_str[i*2] != '1') && (i == 5) ) )
	{
	  missing_layer = i;
	  break;
	}
    }
  } 
  else 
  {
    for (int i = 0; i < 3; ++i)
    {
      if ( layers_str[i*2] != ('5'+i))
      {
        missing_layer = i;
	break;
      }
    }
  }

  return missing_layer;
}

