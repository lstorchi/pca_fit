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

#include "pcaconst.hpp"

using namespace pca;
using namespace std;

template <typename T>
void dump_element (const pca::matrixpcaconst<T> & in, std::ostream & out)
{
  out.precision(6);
  for (unsigned int i = 0; i<in.n_rows(); ++i)
  {
    for (unsigned int j = 0; j<in.n_cols(); ++j)
      out << std::scientific << in(i, j) << " ";
    out << std::endl;
  }
}

bool check_val(double val, int width, bool sign = true) 
{
  if (sign)
  {
    if (abs(val) > pow(2, width-1))
      std::cerr << "value: " << (long long int) val << " width: " << width << std::endl;

    return (abs(val) > pow(2, width-1));
  }

  if (abs(val) > pow(2, width))
    std::cerr << "value: " << (long long int) val << " width: " << width << std::endl;

  return (abs(val) > pow(2, width));
}

bool read_integer_const_filename (const std::string & in, 
    std::vector<pca::matrixpcaconst<long long int> > & pcacontvct_integer)
{
  std::cout << "Reading " << in << std::endl;
  if (!pca::read_pcaconst_from_file (pcacontvct_integer, in.c_str()))
  {
    std::cerr << "Error while reading constant from " << in << " read only " <<                   
      pcacontvct_integer.size() << std::endl;
    return false;
  }
        
  std::vector<pca::matrixpcaconst<long long int> >::const_iterator it = 
    pcacontvct_integer.begin();
  for (; it != pcacontvct_integer.end(); ++it)            
  {
    if (it->get_ttype() != pca::INTEGPT)
    {        
      std::cerr << "Wrong PCAconst type for " << it->get_towerid() <<      
        " input file " <<  in << std::endl;      
      return false;                                        
    }                    
  }
            
  std::cout << "Done" << std::endl;

  return true;
}


int main (int argc, char *argv[]) 
{
  if (argc != 2)
  {
    std::cerr << "usage: " << argv[0] << " pca_const_file.txt " << std::endl;
    return EXIT_FAILURE;
  }

  std::string in_fname = argv[1];
  std::vector<pca::matrixpcaconst<long long int> > pcacontvct_integer;
  
  if (read_integer_const_filename (in_fname, pcacontvct_integer))
  {
    pca::matrixpcaconst<long long int> cmtx_rz(0, 0);
    pca::matrixpcaconst<long long int> qvec_rz(0, 0); 
    pca::matrixpcaconst<long long int> amtx_rz(0, 0); 
    pca::matrixpcaconst<long long int> kvec_rz(0, 0); 
    pca::matrixpcaconst<long long int> cmtx_rphi(0, 0); 
    pca::matrixpcaconst<long long int> qvec_rphi(0, 0); 
    pca::matrixpcaconst<long long int> amtx_rphi(0, 0); 
    pca::matrixpcaconst<long long int> kvec_rphi(0, 0); 

    double eta_est, pt_est, charge;
    int tow;

    pca::matrixpcaconst<long long int> zrv(1, 12), phirv(1, 12);
    std::string layersid, pslayersid;

    bool rz5oof6, rphi5oof6;

    /* TEST FW 
    |---------+----------------|
    | value   | scaling factor |
    |---------+----------------|
    | r       |           1024 | OK
    | phi     |          65536 | OK
    | z       |            256 | OK
    | pt      |            256 | ==> ?
    | eta     |           1024 | ==> ?
    | results |       1e5*1024 | OK
    |---------+----------------|

    ----------------- Track #1 --------------------

    |----------+------+-----|
    | Pattern# |   pt | eta |
    |----------+------+-----|
    | 0x05C5A2 | 1172 | 171 |
    |----------+------+-----|
    */

    eta_est = 0.1669921875;
    pt_est = 4.578125;
    charge = +1;

    tow = 25;
 
    /*
    |--------+--------+-------+------|
    | Layer# |      r |   phi |    z |
    |--------+--------+-------+------|
    |      1 |  34971 | 43112 | -534 | missing layer is 5
    |      2 |  51069 | 44438 |  141 |
    |      8 |  69286 | 45914 |  751 |
    |      9 |  91956 | 47752 | 2038 |
    |     10 | 108987 | 49128 | 3267 |
    |--------+--------+-------+------|
    */

    rz5oof6 = true;
    rphi5oof6 = true;

    layersid = "6:7:8:9:10"; 
    pslayersid = "6:7";

    zrv(0, 0) = -534; zrv(0, 1) = 43112;
    zrv(0, 2) =  141; zrv(0, 3) = 44438;
    //zrv(0, 4) = 733; zrv(0, 5) = 53194;
    
    phirv(0, 0) = 43112; phirv(0, 1) =  34971;
    phirv(0, 2) = 44438; phirv(0, 3) =  51069;
    phirv(0, 4) = 45914; phirv(0, 5) =  69286;
    phirv(0, 6) = 47752; phirv(0, 7) =  91956;
    phirv(0, 8) = 49128; phirv(0, 9) = 108987;
    //phirv(0, 10) = 35082; phirv(0, 11) = 111634;

    double sec_phi = (tow%8) * M_PI / 4.0 - 0.4;

    if (pca::import_pca_const (pcacontvct_integer, 
                          cmtx_rz, 
                          qvec_rz, 
                          amtx_rz, 
                          kvec_rz, 
                          cmtx_rphi, 
                          qvec_rphi, 
                          amtx_rphi, 
                          kvec_rphi, 
                          eta_est, 
                          pt_est, 
                          charge,
                          layersid, 
                          pslayersid, 
                          tow, 
                          pca::INTEGPT))
    {
      std::cout << "CMTX RZ: " << std::endl;
      dump_element(cmtx_rz, std::cout);
      
      std::cout << "QVEC RZ: " << std::endl;
      dump_element(qvec_rz, std::cout);
      
      std::cout << "CMTX RPHI: " << std::endl;
      dump_element(cmtx_rphi, std::cout);
      
      std::cout << "QVEC RPHI: " << std::endl;
      dump_element(qvec_rphi, std::cout);

      long long int cottheta = 0;
      long long int z0 = 0;
      
      cottheta = qvec_rz(0,0);
      z0 = qvec_rz(0,1);
      if (check_val((double) cottheta, pca::add_const_w))
        std::cerr << "Overflow in cottheta add_const_w " << std::endl;
      if (check_val((double) z0, pca::add_const_w))
        std::cerr << "Overflow in z0 add_const_w" << std::endl;
      for (int i=0; i<(int)cmtx_rz.n_cols(); ++i)
      {
        cottheta += cmtx_rz(0, i) * zrv(0, i);
        z0 += cmtx_rz(1, i) * zrv(0, i);
        if (check_val((double) z0, pca::add_const_w))
          std::cerr << "Overflow in z0 add_const_w" << std::endl;
        if (check_val((double) cottheta, pca::add_const_w))
          std::cerr << "Overflow in cottheta add_const_w " << std::endl;
      }
      if (check_val((double) z0, pca::result_w))
        std::cerr << "Overflow in z0 result_w" << std::endl;
      if (check_val((double) cottheta, pca::result_w))
        std::cerr << "Overflow in cottheta result_w" << std::endl;

      std::cout << "int cotetha:    " << cottheta << std::endl;
      std::cout << "int z0:         " << z0 << std::endl;

      double d_cottheta = (double) ((double)cottheta / (double) pca::const_mult_factor);
      double d_z0 = (double) ((double)z0 / (double)pca::const_mult_factor);

      double eta = 0.0e0;
      double theta = atan(1.0e0 / d_cottheta); 
      double tantheta2 = tan (theta/2.0e0); 
      if (tantheta2 < 0.0)
        eta = 1.0e0 * log (-1.0e0 * tantheta2);
      else
        eta = -1.0e0 * log (tantheta2);
      
      long long int coverpt = 0; // pt
      long long int phi = 0;
      
      coverpt = qvec_rphi(0,0);
      phi = qvec_rphi(0,1);
      if (check_val((double) coverpt, pca::add_const_w))
        std::cerr << "Overflow in coverpt add_const_w " << std::endl;
      if (check_val((double) phi, pca::add_const_w))
        std::cerr << "Overflow in phi add_const_w" << std::endl;
      for (int i=0; i<(int)cmtx_rphi.n_cols(); ++i)
      {
        coverpt += cmtx_rphi(0, i) * phirv(0, i);
        phi += cmtx_rphi(1, i) * phirv(0, i);
        if (check_val((double) coverpt, pca::add_const_w))
          std::cerr << "Overflow in coverpt add_const_w " << std::endl;
        if (check_val((double) phi, pca::add_const_w))
          std::cerr << "Overflow in phi add_const_w" << std::endl;
      }
      if (check_val((double) coverpt, pca::result_w))
        std::cerr << "Overflow in coverpt result_w" << std::endl;
      if (check_val((double) phi, pca::result_w))
        std::cerr << "Overflow in phi result_w" << std::endl;

      std::cout << "int coverpt:    " << coverpt << std::endl;
      std::cout << "int phi:        " << phi << std::endl;

      double d_pt = (double) charge / ((double) coverpt / (double) pca::const_mult_factor);
      double d_phi = (double) ((double) phi / (double) pca::const_mult_factor);

      if ((tow == 19) || (tow == 20) ||
          (tow == 27) || (tow == 28) ||
          (tow == 11) || (tow == 12) ||
          (tow == 35) || (tow == 36))
        d_phi -= sec_phi;

      int coordim = 6, paramdim = 2;
      if (rz5oof6)
        coordim = 4;

      long long int chi2rz = 0.0;
      for (int i=0; i<coordim-paramdim; ++i)
      {
        long long int val = 0.0;
                                    
        for (int j=0; j<coordim; ++j)
        {
          val += amtx_rz(i,j) * zrv(0, j);
          if (check_val((double) val, pca::add_const_w))
            std::cerr << "Overflow in val chi2rz 1 add_const_w " << std::endl;
        }

        val -= kvec_rz(0, i);
        if (check_val((double) val, pca::add_const_w))
          std::cerr << "Overflow in val chi2rz 2 add_const_w " << std::endl;
        
        chi2rz += val*val;
        if (check_val((double) chi2rz, pca::result_w, false))
          std::cerr << "Overflow in chi2rz 3 result_w " << std::endl;
      }

      coordim = 12, paramdim = 2;
      if (rphi5oof6)
        coordim = 10;

      long long int chi2rphi = 0.0;
      for (int i=0; i<coordim-paramdim; ++i)
      {
        long long int val = 0.0;
                                    
        for (int j=0; j<coordim; ++j)
        {
          val += amtx_rphi(i,j) * phirv(0, j);
          if (check_val((double) val, pca::add_const_w))
            std::cerr << "Overflow in val chi2rphi 1 add_const_w " << std::endl;
        }

        val -= kvec_rphi(0, i);
        if (check_val((double) val, pca::add_const_w))
          std::cerr << "Overflow in val chi2rphi 2 add_const_w " << std::endl;
        
        chi2rphi += val*val;
        if (check_val((double) chi2rphi, pca::result_w, false))
          std::cerr << "Overflow in val chi2rphi 3 result_w " << std::endl;
      }

      if (check_val((double) chi2rz, pca::result_w, false))
        std::cerr << "Overflow in chi2rz result_w" << std::endl;
      if (check_val((double) chi2rphi, pca::result_w, false))
        std::cerr << "Overflow in chi2rphi result_w" << std::endl;

      double d_chi2rz = (double) chi2rz / (double) pow(pca::chisq_const_mult_factor, 2);
      double d_chi2rphi = (double) chi2rphi / (double) pow(pca::chisq_const_mult_factor, 2);

      std::cerr << "d_chi2rz: " << d_chi2rz << " d_chi2rphi: " << d_chi2rphi  << std::endl;

      std::cout << "pt:         " << d_pt << std::endl;
      std::cout << "phi:        " << d_phi << std::endl; 
      std::cout << "eta:        " << eta << std::endl;
      std::cout << "z0:         " << d_z0 << std::endl;

      if (rz5oof6)
        std::cout << "chi2rz:     " << d_chi2rz/2.0 << std::endl;
      else
        std::cout << "chi2rz:     " << d_chi2rz/4.0 << std::endl;

      if (rphi5oof6)
        std::cout << "chi2rphi:   " << d_chi2rphi/8.0 << std::endl;
      else
        std::cout << "chi2rphi:   " << d_chi2rphi/10.0 << std::endl;
    }
  }

  return EXIT_SUCCESS;
}
