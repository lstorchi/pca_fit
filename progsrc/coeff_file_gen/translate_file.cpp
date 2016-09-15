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
	for (int i = 0; i < 2; ++i)
        {
	  for (int j = 0; j < 12; ++j)
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
    if (c.get_plane_type() == plane_type::RPHI)
      if (c.get_const_type() == const_type::QVEC) 
	for (int i = 0; i < 2; ++i)
	  constants_rphi_qvec.push_back(c.element(0, i));

  cout << constants_rphi_cmtx.size() << ' ' << constants_rphi_qvec.size() << endl;
  
  ofstream ofs_rphi_param("rphi_patt_coeff.coe", ofstream::out);
  ofs_rphi_param << "memory_initialization_radix=2;\nmemory_initialization_vector=\n";

  ofstream ofs_rphi_param_nb("rphi_patt_coeff", ofstream::out);

  int base_addr = 0;
  long long int tmp;
  for (int ptbins = 0; ptbins < 14; ++ptbins)
  {
    for (int mlcombs = 0; mlcombs < 7; ++mlcombs) 
    {
      base_addr = (ptbins & 0xfe)*7 + mlcombs*2 + (ptbins & 0x1);
      tmp = const_mult_factor*constants_rphi_qvec[base_addr*2+1];
      check_val(tmp, add_const_w, "rphi add consts");
      ofs_rphi_param_nb << tmp << " ";
      bitset<add_const_w> tmp_bin1(tmp);
      ofs_rphi_param << tmp_bin1;

      tmp = const_mult_factor*constants_rphi_qvec[base_addr*2];
      check_val(tmp, add_const_w, "rphi add consts");
      ofs_rphi_param_nb << tmp << " ";
      bitset<add_const_w> tmp_bin2(tmp);
      ofs_rphi_param << tmp_bin2;

      for (int i = 2*12-1; i>=0; --i) 
      {
	tmp = mult_factor*constants_rphi_cmtx[base_addr*2*12 + i];
        ofs_rphi_param_nb << tmp << std::endl;

	if (i%2==0)
	  tmp /= 64;
	check_val(tmp, const_w, "rphi matrix consts");
	bitset<const_w> tmp_bin3(tmp);
	ofs_rphi_param << tmp_bin3;
      }
      ofs_rphi_param << ',' << endl;
      ofs_rphi_param_nb << ',' << endl;
    }
  }
  
  ofs_rphi_param.close();
  ofs_rphi_param_nb.close();

  // write rphi chisq
  //int prevsize;
  std::vector<double> constants_rphi_amtx;
  for (matrixpcaconst<double> c : all_constants) 
  {
    if (c.get_plane_type() == plane_type::RPHI)
    {
      if (c.get_const_type() == const_type::AMTX) 
      {
	//cout << c.n_rows() << ' ' << c.n_cols() << endl;
	//prevsize = constants_rphi_amtx.size();
	for (int i = 0; i < 10; ++i) 
        {
	  if ( (get_missing_layer(c) == -1) ||
	       ( (get_missing_layer(c) >= 0) &&
		 (i < 8)
		 )
	       ) 
          {
	    for (int j = 0; j < 12; ++j) 
            {
	      if ( (get_missing_layer(c) >= 0) &&
		   (get_missing_layer(c) <= (j/2))
		   )
		constants_rphi_amtx.push_back(((get_missing_layer(c) == (j/2))?0:c.element(i, j-2)));
	      else
		constants_rphi_amtx.push_back(c.element(i, j));
	    }
	  } 
          else 
          {
	    for (int j = 0; j < 12; ++j) 
	      constants_rphi_amtx.push_back(0);
	  }
	}
	//cout << constants_rphi_amtx.size() - prevsize << endl;
      }
    }
  }

  //int prevsize;
  std::vector<double> constants_rphi_kvec;
  for (matrixpcaconst<double> c : all_constants) 
  {
    if (c.get_plane_type() == plane_type::RPHI)
    {
      if (c.get_const_type() == const_type::KVEC) 
      {
	//prevsize = constants_rphi_kvec.size();
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
	//cout << constants_rphi_kvec.size() - prevsize << endl;
      }
    }
  }

  cout << constants_rphi_amtx.size() << ' ' << constants_rphi_kvec.size() << endl;
  
  ofstream ofs_rphi_chisq("rphi_patt_coeff_chisq.coe", ofstream::out);
  ofs_rphi_chisq << "memory_initialization_radix=2;\nmemory_initialization_vector=\n";

  for (int ptbins = 0; ptbins < 14; ++ptbins)
  {
    for (int mlcombs = 0; mlcombs < 7; ++mlcombs) 
    {
      base_addr = (ptbins & 0xfe)*7 + mlcombs*2 + (ptbins & 0x1);

      for (int i = 10-1; i>=0; --i) 
      {
	tmp = chisq_const_mult_factor*constants_rphi_kvec[base_addr*10+i];
	check_val(tmp, add_const_w, "rphi chisq add consts");
	bitset<add_const_w> tmp_bin1(tmp);
	ofs_rphi_chisq << tmp_bin1;
      }

      for (int i = 10*12-1; i>=0; --i) 
      {
	tmp = chisq_mult_factor*constants_rphi_amtx[base_addr*10*12 + i];
	if (i%2==0)
	  tmp /= 64;
	check_val(tmp, const_w, "rphi chisq matrix consts");
	bitset<const_w> tmp_bin3(tmp);
	ofs_rphi_chisq << tmp_bin3;
      }
      ofs_rphi_chisq << ',' << endl;
    }
  }
  
  ofs_rphi_chisq.close();

  // write rz param
  std::vector<double> constants_rz_cmtx;
  for (matrixpcaconst<double> c : all_constants) 
  {
    if (c.get_plane_type() == plane_type::RZ)
    {
      if (c.get_const_type() == const_type::CMTX) 
      {
	for (int i = 0; i < 2; ++i)
        {
	  for (int j = 0; j < 6; ++j)
          {
	    if ( (get_missing_layer(c) >= 0) &&
		 (get_missing_layer(c) <= (j/2))
		 )
	      constants_rz_cmtx.push_back(((get_missing_layer(c) == (j/2))?0:c.element(i, j-2)));
	    else
	      constants_rz_cmtx.push_back(c.element(i, j));
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

  cout << constants_rz_cmtx.size() << ' ' << constants_rz_qvec.size() << endl;
  
  ofstream ofs_rz_param("rz_patt_coeff.coe", ofstream::out);
  ofs_rz_param << "memory_initialization_radix=2;\nmemory_initialization_vector=\n";

  for (int etabins = 0; etabins < 20; ++etabins)
  {
    for (int mlcombs = 0; mlcombs < 4; ++mlcombs) 
    {
      base_addr = etabins*4 + ((mlcombs+3)%4);
      tmp = const_mult_factor*constants_rz_qvec[base_addr*2+1];
      check_val(tmp, add_const_w, "rz add consts");
      bitset<add_const_w> tmp_bin1(tmp);
      ofs_rz_param << tmp_bin1;

      tmp = const_mult_factor*constants_rz_qvec[base_addr*2];
      check_val(tmp, add_const_w, "rz add consts");
      bitset<add_const_w> tmp_bin2(tmp);
      ofs_rz_param << tmp_bin2;

      for (int i = 2*6-1; i>=0; --i) 
      {
	tmp = mult_factor*constants_rz_cmtx[base_addr*2*6 + i];
	if (i%2==0)
	  tmp *= 4;
	check_val(tmp, const_w, "rz matrix consts");
	bitset<const_w> tmp_bin3(tmp);
	ofs_rz_param << tmp_bin3;
      }
      ofs_rz_param << ',' << endl;
    }
  }
  
  ofs_rz_param.close();

  // write rz chisq
  //int prevsize;
  std::vector<double> constants_rz_amtx;
  for (matrixpcaconst<double> c : all_constants) 
  {
    if (c.get_plane_type() == plane_type::RZ)
    {
      if (c.get_const_type() == const_type::AMTX) 
      {
	//cout << c.n_rows() << ' ' << c.n_cols() << endl;
	//prevsize = constants_rz_amtx.size();
	for (int i = 0; i < 4; ++i) 
        {
	  if ( (get_missing_layer(c) == -1) ||
	       ( (get_missing_layer(c) >= 0) &&
		 (i < 2)
		 )
	       ) 
          {
	    for (int j = 0; j < 6; ++j) 
            {
	      if ( (get_missing_layer(c) >= 0) &&
		   (get_missing_layer(c) <= (j/2))
		   )
		constants_rz_amtx.push_back(((get_missing_layer(c) == (j/2))?0:c.element(i, j-2)));
	      else
		constants_rz_amtx.push_back(c.element(i, j));
	    }
	  } 
          else 
          {
	    for (int j = 0; j < 6; ++j) 
	      constants_rz_amtx.push_back(0);
	  }
	}
	//cout << constants_rz_amtx.size() - prevsize << endl;
      }
    }
  }

  //int prevsize;
  std::vector<double> constants_rz_kvec;
  for (matrixpcaconst<double> c : all_constants) 
  {
    if (c.get_plane_type() == plane_type::RZ)
    {
      if (c.get_const_type() == const_type::KVEC) 
      {
	//prevsize = constants_rz_kvec.size();
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
	//cout << constants_rz_kvec.size() - prevsize << endl;
      }
    }
  }

  cout << constants_rz_amtx.size() << ' ' << constants_rz_kvec.size() << endl;
  
  ofstream ofs_rz_chisq("rz_patt_coeff_chisq.coe", ofstream::out);
  ofs_rz_chisq << "memory_initialization_radix=2;\nmemory_initialization_vector=\n";

  for (int etabins = 0; etabins < 20; ++etabins)
  {
    for (int mlcombs = 0; mlcombs < 4; ++mlcombs) 
    {
      base_addr = etabins*4 + ((mlcombs+3)%4);

      for (int i = 4-1; i>=0; --i) 
      {
	tmp = chisq_const_mult_factor*constants_rz_kvec[base_addr*4+i];
	check_val(tmp, add_const_w, "rz chisq add consts");
	bitset<add_const_w> tmp_bin1(tmp);
	ofs_rz_chisq << tmp_bin1;
      }

      for (int i = 4*6-1; i>=0; --i) 
      {
	tmp = chisq_mult_factor*constants_rz_amtx[base_addr*4*6 + i];
	if (i%2==0)
	  tmp *= 4;
	check_val(tmp, const_w, "rz chisq matrix consts");
	bitset<const_w> tmp_bin3(tmp);
	ofs_rz_chisq << tmp_bin3;
      }
      ofs_rz_chisq << ',' << endl;
    }
  }
  
  ofs_rz_chisq.close();

  return 0;
}
