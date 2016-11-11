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
  if (argc != 3)
  {
    std::cerr << "usage: " << argv[0] << " pca_const_file_1.txt pca_const_file_2.txt" << std::endl;
    std::cerr << "       the second file is added to the first " << std::endl;
    return EXIT_FAILURE;
  }

  std::string out_fname;

  std::vector< matrixpcaconst<double> > all_constants_1;
  bool a = read_pcaconst_from_file<double>(all_constants_1, argv[1]);
  std::vector< matrixpcaconst<double> > all_constants_2;
  bool b = read_pcaconst_from_file<double>(all_constants_2, argv[2]);

  if ((!a) || (!b))
  {
    std::cerr << "Error while reading from " << argv[1] << " or " << 
      argv[2] << std::endl;
    return EXIT_FAILURE; 
  }

  out_fname = argv[1];
  out_fname.append ("_allranges");

  fstream fin(out_fname.c_str());
  if (fin)
  {
    std::cout << "File " << out_fname << " found, deleting..." << std::endl;
    fin.close();

    remove(out_fname.c_str());
  } 

  double rpdval = rpdval_get(); 
  double rzdval = rzdval_get();
  
  std::vector< matrixpcaconst<double> > cout;

  for (matrixpcaconst<double> c : all_constants_1) 
  {
    matrixpcaconst<double> ci(c.n_rows(), c.n_cols());

    ci.set_const_type(c.get_const_type());
    ci.set_towerid(c.get_towerid());
    ci.set_layersids(c.get_layersids().c_str());
    ci.set_sector_type(c.get_sector_type());
    ci.set_plane_type(c.get_plane_type());
    ci.set_chargesign(c.get_chargesign());
    ci.set_ttype(c.get_ttype());

    double ptmin, ptmax, etamin, etamax;
    c.get_ptrange(ptmin, ptmax);
    c.get_etarange(etamin, etamax);
    ci.set_ptrange(ptmin, ptmax);
    ci.set_etarange(etamin, etamax);

    bool allzero = false;

    if (c.get_plane_type() == plane_type::RPHI)
    {
      if (c.get_const_type() == const_type::CMTX) 
      {
        int czeros = 0;
	for (int i = 0; i < c.n_rows(); ++i)
        {
	  for (int j = 0; j < c.n_cols(); ++j)
          {
            ci.element(i, j) = c.element(i,j);
            if (ci.element(i, j) == 0.0)
              ++czeros;
          }
        }

        if (czeros == (c.n_rows() * c.n_cols()))
          allzero = true;
      }
      else if (c.get_const_type() == const_type::QVEC)
      {
        int czeros = 0;
	for (int i = 0; i < c.n_cols(); ++i)
        {
          ci.element(0, i) = c.element(0,i);
          if (ci.element(0, i) == 0.0)
            ++czeros;
        }

        if (czeros == c.n_cols())
          allzero = true;
      }
      else if (c.get_const_type() == const_type::AMTX) 
      {
        int czeros = 0;
	for (int i = 0; i < c.n_rows(); ++i)
        {
	  for (int j = 0; j < c.n_cols(); ++j)
          {
            ci.element(i, j) = c.element(i,j);
            if (ci.element(i, j) == 0.0)
              ++czeros;
          }
        }

        if (czeros == (c.n_rows() * c.n_cols()))
          allzero = true;
      }
      else if (c.get_const_type() == const_type::KVEC)
      {
        int czeros = 0;
	for (int i = 0; i < c.n_cols(); ++i)
        {
          ci.element(0, i) = c.element(0,i);
          if (ci.element(0, i) == 0.0)
            ++czeros; 
        }

        if (czeros == c.n_cols())
          allzero = true;
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
        int czeros = 0;
	for (int i = 0; i < c.n_rows(); ++i)
        {
	  for (int j = 0; j < c.n_cols(); ++j)
          {
            ci.element(i, j) = c.element(i,j);
            if (ci.element(i, j) == 0.0)
              ++czeros;
          }
        }

        if (czeros == c.n_cols())
          allzero = true;
      }
      else if (c.get_const_type() == const_type::QVEC)
      {
        int czeros = 0;
	for (int i = 0; i < c.n_cols(); ++i)
        {
          ci.element(0, i) = c.element(0,i);
          if (ci.element(0, i) == 0.0)
            ++czeros;
        }

        if (czeros == c.n_cols())
          allzero = true;
      }
      else if (c.get_const_type() == const_type::AMTX) 
      {
        int czeros = 0;
	for (int i = 0; i < c.n_rows(); ++i)
        {
	  for (int j = 0; j < c.n_cols(); ++j)
          {
            ci.element(i, j) = c.element(i,j);
            if (ci.element(i, j) == 0.0)
              ++czeros;
          }
        }

        if (czeros == c.n_cols())
          allzero = true;
      }
      else if (c.get_const_type() == const_type::KVEC)
      {
        int czeros = 0;
	for (int i = 0; i < c.n_cols(); ++i)
        {
          ci.element(0, i) = c.element(0,i);
          if (ci.element(0, i) == 0.0)
            ++czeros;
        }

        if (czeros == c.n_cols())
          allzero = true;
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

    if (allzero)
      std::cerr << "all value are zero " << 
       pca::matrixpcaconst<double>::const_type_to_string(
           ci.get_const_type()) << std::endl;

    cout.push_back(ci);
  }

  for (matrixpcaconst<double> c : all_constants_2) 
  {
    bool jump = false;
    for (matrixpcaconst<double> ctcmp : cout) 
    {
      if (ctcmp == c)
      {
        jump = true;
        break;
      }
    }

    if (! jump)
    {
      matrixpcaconst<double> ci(c.n_rows(), c.n_cols());

      ci.set_const_type(c.get_const_type());
      ci.set_towerid(c.get_towerid());
      ci.set_layersids(c.get_layersids().c_str());
      ci.set_sector_type(c.get_sector_type());
      ci.set_plane_type(c.get_plane_type());
      ci.set_chargesign(c.get_chargesign());
      ci.set_ttype(c.get_ttype());
     
      double ptmin, ptmax, etamin, etamax;
      c.get_ptrange(ptmin, ptmax);
      c.get_etarange(etamin, etamax);
      ci.set_ptrange(ptmin, ptmax);
      ci.set_etarange(etamin, etamax);

      if (c.get_plane_type() == plane_type::RPHI)
      {
        if (c.get_const_type() == const_type::CMTX) 
        {
          for (int i = 0; i < c.n_rows(); ++i)
            for (int j = 0; j < c.n_cols(); ++j)
              ci.element(i, j) = c.element(i,j);
        }
        else if (c.get_const_type() == const_type::QVEC)
        {
          for (int i = 0; i < c.n_cols(); ++i)
            ci.element(0, i) = c.element(0,i);
        }
        else if (c.get_const_type() == const_type::AMTX) 
        {
          for (int i = 0; i < c.n_rows(); ++i)
            for (int j = 0; j < c.n_cols(); ++j)
              ci.element(i, j) = c.element(i,j);
        }
        else if (c.get_const_type() == const_type::KVEC)
        {
          for (int i = 0; i < c.n_cols(); ++i)
            ci.element(0, i) = c.element(0,i);
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
          for (int i = 0; i < c.n_rows(); ++i)
            for (int j = 0; j < c.n_cols(); ++j)
              ci.element(i, j) = c.element(i,j);
        }
        else if (c.get_const_type() == const_type::QVEC)
        {
          for (int i = 0; i < c.n_cols(); ++i)
            ci.element(0, i) = c.element(0,i);
        }
        else if (c.get_const_type() == const_type::AMTX) 
        {
          int counter = 0;
          for (int i = 0; i < c.n_rows(); ++i)
            for (int j = 0; j < c.n_cols(); ++j)
              ci.element(i, j) = c.element(i,j);
        }
        else if (c.get_const_type() == const_type::KVEC)
        {
          for (int i = 0; i < c.n_cols(); ++i)
            ci.element(0, i) = c.element(0,i);
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

      cout.push_back(ci);
    }
  }

  quick_write_pcaconst_to_file (cout, out_fname.c_str());

  return EXIT_SUCCESS;
}
