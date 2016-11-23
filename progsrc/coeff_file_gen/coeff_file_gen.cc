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
    return 1;
  }

  std::vector< matrixpcaconst<double> > all_constants;
  bool a = read_pcaconst_from_file<double>(all_constants, argv[1]);
  if (!a)
  {
    cout << "Error loading pca constants, exiting\n";
    return 1;
  }

  double rpdval = rpdval_get(); 
  double rzdval = rzdval_get();

  int maxN;

  // load rphi param
  std::vector<double> constants_rphi_cmtx;
  for (matrixpcaconst<double> c : all_constants) 
  {
    if (c.get_plane_type() == plane_type::RPHI)
    {
      if (c.get_const_type() == const_type::CMTX) 
      {
        if (get_missing_layer(c) == -1)
          maxN = 12;
	else
	  maxN = 10;
	for (int i = 0; i < 2; ++i)
        {
	  for (int j = 0; j < maxN; ++j)
            constants_rphi_cmtx.push_back(c.element(i, j));
          if (maxN != 12)
          {
	    constants_rphi_cmtx.push_back(0);
	    constants_rphi_cmtx.push_back(0);
	  }
        }
      }
    }
  }

  std::vector<double> constants_rphi_qvec;
  for (matrixpcaconst<double> c : all_constants) 
    if (c.get_plane_type() == plane_type::RPHI)
      if (c.get_const_type() == const_type::QVEC) 
	for (int i = 0; i < 2; ++i)
	  constants_rphi_qvec.push_back(c.element(0, i));

  cout << "Matrix sizes:" << endl;
  cout << "RPhi CMTX:\t" << constants_rphi_cmtx.size() << "\tQVEC:\t"  << constants_rphi_qvec.size() << endl;

  // write rphi param
  ofstream ofs_rphi_param("rphi_patt_coeff.coe", ofstream::out);
  ofs_rphi_param << "memory_initialization_radix=2;\nmemory_initialization_vector=\n";

  ofstream ofs_rphi_param_nb("rphi_patt_coeff", ofstream::out);

  int base_addr = 0;
  long long int tmp;
  for (int ptbins = 0; ptbins < 14; ++ptbins)
  {
    for (int mlcombs = 0; mlcombs < 7; ++mlcombs) 
    {
      // In the coefficients file, the constants are arranged like so:
      // Pt bin: 0, missing layer: -1, charge:  1
      // Pt bin: 0, missing layer: -1, charge: -1
      // Pt bin: 0, missing layer:  0, charge:  1
      // Pt bin: 0, missing layer:  0, charge: -1
      // Pt bin: 0, missing layer:  1, charge:  1
      // Pt bin: 0, missing layer:  1, charge: -1
      //                      ...
      // Pt bin: 7, missing layer:  5, charge: -1
      // thus the base_addr addressing scheme
      // On the memory address space, coefficients are arranged in the same way
      base_addr = (ptbins & 0xfe)*7 + mlcombs*2 + (ptbins & 0x1);

      for (int i = 2-1; i>=0; --i)
      {
	tmp = const_mult_factor*constants_rphi_qvec[base_addr*2 + i];
        ofs_rphi_param_nb << tmp << " ";
	ofs_rphi_param << format_and_check_val<add_const_w>(tmp, "rphi add consts");
      }

      for (int i = 2*12-1; i>=0; --i) 
      {
	tmp = mult_factor*constants_rphi_cmtx[base_addr*2*12 + i];
        ofs_rphi_param_nb << tmp << std::endl;

	if (i%2==0)
	  tmp /= rpdval;
	ofs_rphi_param << format_and_check_val<const_w>(tmp, "rphi matrix consts");
      }
      ofs_rphi_param << ',' << endl;
      ofs_rphi_param_nb << ',' << endl;
    }
  }
  
  ofs_rphi_param.close();
  ofs_rphi_param_nb.close();

  // load rphi chisq
  std::vector<double> constants_rphi_amtx;
  for (matrixpcaconst<double> c : all_constants) 
  {
    if (c.get_plane_type() == plane_type::RPHI)
    {
      if (c.get_const_type() == const_type::AMTX) 
      {
	//cout << c.n_rows() << ' ' << c.n_cols() << endl;
	for (int i = 0; i < 10; ++i) 
        {
	  if ( (get_missing_layer(c) == -1) ||
	       ( (get_missing_layer(c) >= 0) &&
		 (i < 8)
		 )
	       ) 
          {
            if (get_missing_layer(c) == -1)
	      maxN = 12;
	    else
	      maxN = 10;
	    for (int j = 0; j < maxN; ++j)
              constants_rphi_amtx.push_back(c.element(i, j));
            if (maxN != 12) {
	      constants_rphi_amtx.push_back(0);
	      constants_rphi_amtx.push_back(0);
	    }
	  } 
          else 
          {
	    for (int j = 0; j < 12; ++j) 
	      constants_rphi_amtx.push_back(0);
	  }
	}
      }
    }
  }

  std::vector<double> constants_rphi_kvec;
  for (matrixpcaconst<double> c : all_constants) 
  {
    if (c.get_plane_type() == plane_type::RPHI)
    {
      if (c.get_const_type() == const_type::KVEC) 
      {
	//cout << c.n_rows() << ' ' << c.n_cols() << endl;
	for (int i = 0; i < 10; ++i)
        {
	  if ( (get_missing_layer(c) == -1) ||
	       ( (get_missing_layer(c) >= 0) &&
		 (i < 8)
		 )
	       ) 
	    constants_rphi_kvec.push_back(c.element(0, i));
	  else
	    constants_rphi_kvec.push_back(0);
        }
      }
    }
  }

  cout << "RPhi AMTX:\t" << constants_rphi_amtx.size() << "\tKVEC:\t" << constants_rphi_kvec.size() << endl;

  // write rphi chisq
  ofstream ofs_rphi_chisq("rphi_patt_coeff_chisq.coe", ofstream::out);
  ofs_rphi_chisq << "memory_initialization_radix=2;\nmemory_initialization_vector=\n";

  for (int ptbins = 0; ptbins < 14; ++ptbins)
  {
    for (int mlcombs = 0; mlcombs < 7; ++mlcombs) 
    {
      base_addr = (ptbins & 0xfe)*7 + mlcombs*2 + (ptbins & 0x1);

      for (int i = 10-1; i>=0; --i) 
      {
	tmp = chisq_const_mult_factor*constants_rphi_kvec[base_addr*10 + i];
	ofs_rphi_chisq << format_and_check_val<add_const_w>(tmp, "rphi chisq add consts");
      }

      for (int i = 10*12-1; i>=0; --i) 
      {
	tmp = chisq_mult_factor*constants_rphi_amtx[base_addr*10*12 + i];
	if (i%2==0)
	  tmp /= rpdval;
	ofs_rphi_chisq << format_and_check_val<const_w>(tmp, "rphi chisq matrix consts");
      }
      ofs_rphi_chisq << ',' << endl;
    }
  }
  
  ofs_rphi_chisq.close();

  // load rz param
  std::vector<double> constants_rz_cmtx;
  for (matrixpcaconst<double> c : all_constants) 
  {
    if (c.get_plane_type() == plane_type::RZ)
    {
      if (c.get_const_type() == const_type::CMTX)
      {
        if (get_missing_layer(c) == -1)
	  maxN = 6;
	else
	  maxN = 4;
	for (int i = 0; i < 2; ++i)
        {
	  for (int j = 0; j < maxN; ++j)
	    constants_rz_cmtx.push_back(c.element(i, j));
	  if (maxN != 6)
          {
	      constants_rz_cmtx.push_back(0);
	      constants_rz_cmtx.push_back(0);
	  }
	}
      }
    }
  }

  std::vector<double> constants_rz_qvec;
  for (matrixpcaconst<double> c : all_constants) 
    if (c.get_plane_type() == plane_type::RZ)
      if (c.get_const_type() == const_type::QVEC) 
	for (int i = 0; i < 2; ++i)
	  constants_rz_qvec.push_back(c.element(0, i));

  cout << "RZ CMTX:\t" << constants_rz_cmtx.size() << "\tQVEC:\t" << constants_rz_qvec.size() << endl;

  // write rz param
  ofstream ofs_rz_param("rz_patt_coeff.coe", ofstream::out);
  ofs_rz_param << "memory_initialization_radix=2;\nmemory_initialization_vector=\n";

  for (int etabins = 0; etabins < 20; ++etabins)
  {
    for (int mlcombs = 0; mlcombs < 4; ++mlcombs) 
    {
      // In the coefficients file, the constants are arranged like so:
      // eta bin: 0, missing layer:  0
      // eta bin: 0, missing layer:  1
      // eta bin: 0, missing layer:  2
      // eta bin: 0, missing layer: -1
      // eta bin: 1, missing layer:  0
      // eta bin: 1, missing layer:  1
      // eta bin: 1, missing layer:  2
      // eta bin: 1, missing layer: -1
      //             ...
      // In the memory we want to write coefficients like this:
      // eta bin: 0, missing layer: -1
      // eta bin: 0, missing layer:  0
      // eta bin: 0, missing layer:  1
      // eta bin: 0, missing layer:  2
      base_addr = etabins*4 + ((mlcombs+3)%4);

      for (int i = 2-1; i>=0; --i)
      {
	tmp = const_mult_factor*constants_rz_qvec[base_addr*2 + i];
	ofs_rz_param << format_and_check_val<add_const_w>(tmp, "rz add consts");
      }

      for (int i = 2*6-1; i>=0; --i) 
      {
	tmp = mult_factor*constants_rz_cmtx[base_addr*2*6 + i];
	if (i%2==0)
	  tmp *= rzdval;
	ofs_rz_param << format_and_check_val<const_w>(tmp, "rz matrix consts");
      }
      ofs_rz_param << ',' << endl;
    }
  }
  
  ofs_rz_param.close();

  // load rz chisq
  std::vector<double> constants_rz_amtx;
  for (matrixpcaconst<double> c : all_constants) 
  {
    if (c.get_plane_type() == plane_type::RZ)
    {
      if (c.get_const_type() == const_type::AMTX) 
      {
	//cout << c.n_rows() << ' ' << c.n_cols() << endl;
	for (int i = 0; i < 4; ++i) 
        {
	  if ( (get_missing_layer(c) == -1) ||
	       ( (get_missing_layer(c) >= 0) &&
		 (i < 2)
		 )
	       ) 
          {
            if (get_missing_layer(c) == -1)
	      maxN = 6;
	    else
	      maxN = 4;
	    for (int j = 0; j < maxN; ++j)
	      constants_rz_amtx.push_back(c.element(i, j));
	    if (maxN != 6)
            {
	      constants_rz_amtx.push_back(0);
	      constants_rz_amtx.push_back(0);
	    }
	  } 
          else 
          {
	    for (int j = 0; j < 6; ++j) 
	      constants_rz_amtx.push_back(0);
	  }
	}
      }
    }
  }

  std::vector<double> constants_rz_kvec;
  for (matrixpcaconst<double> c : all_constants) 
  {
    if (c.get_plane_type() == plane_type::RZ)
    {
      if (c.get_const_type() == const_type::KVEC) 
      {
	//cout << c.n_rows() << ' ' << c.n_cols() << endl;
	for (int i = 0; i < 4; ++i)
        {
	  if ( (get_missing_layer(c) == -1) ||
	       ( (get_missing_layer(c) >= 0) &&
		 (i < 2)
		 )
	       ) 
	    constants_rz_kvec.push_back(c.element(0, i));
          else
	    constants_rz_kvec.push_back(0);
        }
      }
    }
  }

  cout << "RZ AMTX:\t" << constants_rz_amtx.size() << "\tKVEC:\t" << constants_rz_kvec.size() << endl;

  // write rz chisq
  ofstream ofs_rz_chisq("rz_patt_coeff_chisq.coe", ofstream::out);
  ofs_rz_chisq << "memory_initialization_radix=2;\nmemory_initialization_vector=\n";

  for (int etabins = 0; etabins < 20; ++etabins)
  {
    for (int mlcombs = 0; mlcombs < 4; ++mlcombs) 
    {
      base_addr = etabins*4 + ((mlcombs+3)%4);

      for (int i = 4-1; i>=0; --i) 
      {
	tmp = chisq_const_mult_factor*constants_rz_kvec[base_addr*4 + i];
	ofs_rz_chisq << format_and_check_val<add_const_w>(tmp, "rz chisq add consts");
      }

      for (int i = 4*6-1; i>=0; --i) 
      {
	tmp = chisq_mult_factor*constants_rz_amtx[base_addr*4*6 + i];
	if (i%2==0)
	  tmp *= rzdval;
	ofs_rz_chisq << format_and_check_val<const_w>(tmp, "rz chisq matrix consts");
      }
      ofs_rz_chisq << ',' << endl;
    }
  }
  
  ofs_rz_chisq.close();

  return 0;
}
