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

double rpdval_get ()
{
  double v  = pca::pisf/pca::risf;

  return v;
}

double rzdval_get ()
{
  double v = pca::risf/pca::zisf;
   
  return v;
}

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


void get_eta_ranges (const int towerid, double & etamin, double & etamax)
{
  switch (towerid)
  {
    case 0:  etamin = -2.4; etamax = -1.2; break;
    case 1:  etamin = -2.4; etamax = -1.2; break;
    case 2:  etamin = -2.4; etamax = -1.2; break;
    case 3:  etamin = -2.4; etamax = -1.2; break;
    case 4:  etamin = -2.4; etamax = -1.2; break;
    case 5:  etamin = -2.4; etamax = -1.2; break;
    case 6:  etamin = -2.4; etamax = -1.2; break;
    case 7:  etamin = -2.4; etamax = -1.2; break;
    case 8:  etamin = -1.7; etamax = -0.4; break;
    case 9:  etamin = -1.7; etamax = -0.4; break;
    case 10: etamin = -1.7; etamax = -0.4; break;
    case 11: etamin = -1.7; etamax = -0.4; break;
    case 12: etamin = -1.7; etamax = -0.4; break;
    case 13: etamin = -1.7; etamax = -0.4; break;
    case 14: etamin = -1.7; etamax = -0.4; break;
    case 15: etamin = -1.7; etamax = -0.4; break;
    case 16: etamin = -0,6; etamax = 0.4; break;
    case 17: etamin = -0,6; etamax = 0.4; break;
    case 18: etamin = -0,6; etamax = 0.4; break;
    case 19: etamin = -0,6; etamax = 0.4; break;
    case 20: etamin = -0,6; etamax = 0.4; break;
    case 21: etamin = -0,6; etamax = 0.4; break;
    case 22: etamin = -0,6; etamax = 0.4; break;
    case 23: etamin = -0,6; etamax = 0.4; break;
    case 24: etamin = -0.4; etamax = 0.6; break;
    case 25: etamin = -0.4; etamax = 0.6; break;
    case 26: etamin = -0.4; etamax = 0.6; break;
    case 27: etamin = -0.4; etamax = 0.6; break;
    case 28: etamin = -0.4; etamax = 0.6; break;
    case 29: etamin = -0.4; etamax = 0.6; break;
    case 30: etamin = -0.4; etamax = 0.6; break;
    case 31: etamin = -0.4; etamax = 0.6; break;
    case 32: etamin = 0.4 ; etamax = 1.7; break;
    case 33: etamin = 0.4 ; etamax = 1.7; break;
    case 34: etamin = 0.4 ; etamax = 1.7; break;
    case 35: etamin = 0.4 ; etamax = 1.7; break;
    case 36: etamin = 0.4 ; etamax = 1.7; break;
    case 37: etamin = 0.4 ; etamax = 1.7; break;
    case 38: etamin = 0.4 ; etamax = 1.7; break;
    case 39: etamin = 0.4 ; etamax = 1.7; break;
    case 40: etamin = 1.2 ; etamax = 2.4; break;
    case 41: etamin = 1.2 ; etamax = 2.4; break;
    case 42: etamin = 1.2 ; etamax = 2.4; break;
    case 43: etamin = 1.2 ; etamax = 2.4; break;
    case 44: etamin = 1.2 ; etamax = 2.4; break;
    case 45: etamin = 1.2 ; etamax = 2.4; break;
    case 46: etamin = 1.2 ; etamax = 2.4; break;
    case 47: etamin = 1.2 ; etamax = 2.4; break;
  }
}
