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

int main (int argc, char *argv[]) 
{
  if (argc != 2)
  {
    std::cerr << "usage: " << argv[0] << " pca_const_file.txt " << std::endl;
    return EXIT_FAILURE;
  }

  std::string out_fname;

  std::vector< matrixpcaconst<double> > all_constants;
  bool a = read_pcaconst_from_file<double>(all_constants, argv[1]);

  if (!a)
  {
    std::cerr << "Error while reading from " << argv[1] << std::endl;
    return EXIT_FAILURE; 
  }

  out_fname = argv[1];
  out_fname.append ("_int");

  fstream fin(out_fname.c_str());
  if (fin)
  {
    std::cout << "File " << out_fname << " found, deleting..." << std::endl;
    fin.close();

    remove(out_fname.c_str());
  } 

  double rpdval = rpdval_get(); 
  double rzdval = rzdval_get();
  
  for (matrixpcaconst<double> c : all_constants) 
  {
    matrixpcaconst<long long int> ci(c.n_rows(), c.n_cols());

    ci.set_const_type(c.get_const_type());
    ci.set_towerid(c.get_towerid());
    ci.set_layersids(c.get_layersids());
    ci.set_sector_type(c.get_sector_type());
    ci.set_plane_type(c.get_plane_type());
    ci.set_chargesign(c.get_chargesign());

    double ptmin, ptmax, etamin, etamax;
    c.get_ptrange(ptmin, ptmax);
    c.get_etarange(etamin, etamax);
    ci.set_ptrange(ptmin, ptmax);
    ci.set_etarange(etamin, etamax);

    ci.set_ttype(ttype::INTEGPT);

    if (c.get_plane_type() == plane_type::RPHI)
    {
      if (c.get_const_type() == const_type::CMTX) 
      {
        int counter = 0;
	for (int i = 0; i < c.n_rows(); ++i)
        {
	  for (int j = 0; j < c.n_cols(); ++j)
          {
            long long int tmp = mult_factor*c.element(i,j);
            /* TODO: unclear ask, where counter is the global index of element */
            if ( counter%2 == 0) tmp /= rpdval;
            check_val(tmp, const_w, "rphi matrix consts");
            ci.element(i, j) = tmp;
            counter++;
          }
        }
      }
      else if (c.get_const_type() == const_type::QVEC)
      {
	for (int i = 0; i < c.n_cols(); ++i)
        {
          long long int tmp = const_mult_factor*c.element(0,i);
          check_val (tmp, add_const_w, "rphi add consts");
          ci.element(0, i) = tmp;
        }
      }
      else if (c.get_const_type() == const_type::AMTX) 
      {
        int counter = 0;
	for (int i = 0; i < c.n_rows(); ++i)
        {
	  for (int j = 0; j < c.n_cols(); ++j)
          {
            long long int tmp = chisq_mult_factor*c.element(i,j);
            /* TODO: unclear ask, where counter is the global index of element */
            if ( counter%2 == 0) tmp /= rpdval;
            check_val(tmp, const_w, "rphi matrix consts");
            ci.element(i, j) = tmp;
            counter++;
          }
        }
      }
      else if (c.get_const_type() == const_type::KVEC)
      {
	for (int i = 0; i < c.n_cols(); ++i)
        {
          long long int tmp = chisq_const_mult_factor*c.element(0,i);
          check_val (tmp, add_const_w, "rphi add consts");
          ci.element(0, i) = tmp;
        }
      }
      else
      {
        std::cerr << "Unknown Const Type" << std::endl;
        return EXIT_FAILURE;
      }
    }
    else if (c.get_plane_type() == plane_type::RZ)
    {
      if (c.get_const_type() == const_type::CMTX) 
      {
        int counter = 0;
	for (int i = 0; i < c.n_rows(); ++i)
        {
	  for (int j = 0; j < c.n_cols(); ++j)
          {
            long long int tmp = mult_factor*c.element(i,j);
            /* TODO: unclear ask, where counter is the global index of element */
            if ( counter%2 == 0) tmp *= rzdval;
            check_val(tmp, const_w, "rz matrix consts");
            ci.element(i, j) = tmp;
            counter++;
          }
        }
      }
      else if (c.get_const_type() == const_type::QVEC)
      {
	for (int i = 0; i < c.n_cols(); ++i)
        {
          long long int tmp = const_mult_factor*c.element(0,i);
          check_val (tmp, add_const_w, "rz add consts");
          ci.element(0, i) = tmp;
        }
      }
      else if (c.get_const_type() == const_type::AMTX) 
      {
        int counter = 0;
	for (int i = 0; i < c.n_rows(); ++i)
        {
	  for (int j = 0; j < c.n_cols(); ++j)
          {
            long long int tmp = chisq_mult_factor*c.element(i,j);
            /* TODO: unclear ask, where counter is the global index of element */
            if ( counter%2 == 0) tmp *= rzdval;
            check_val(tmp, const_w, "rz matrix consts");
            ci.element(i, j) = tmp;
            counter++;
          }
        }
      }
      else if (c.get_const_type() == const_type::KVEC)
      {
	for (int i = 0; i < c.n_cols(); ++i)
        {
          long long int tmp = chisq_const_mult_factor*c.element(0,i);
          check_val (tmp, add_const_w, "rz add consts");
          ci.element(0, i) = tmp;
        }
      }
      else
      {
        std::cerr << "Unknown Const Type" << std::endl;
        return EXIT_FAILURE;
      }
    }
    else 
    {
      std::cerr << "Unknown Plane Type" << std::endl;
      return EXIT_FAILURE;
    }

    write_pcaconst_to_file (ci, out_fname.c_str());
  }

  return EXIT_SUCCESS;
}
