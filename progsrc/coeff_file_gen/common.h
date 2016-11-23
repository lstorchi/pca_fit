#include "pcaconst.hpp"
#include <bitset>
#include <typeinfo>

// questo e' il rapporto fra scaling factor di r e phi so pca::risf/pca::pisf
double rpdval_get ();
double rzdval_get ();

template<int width>
std::bitset<width> format_and_check_val(double val, std::string mstr)
{
  std::bitset<width> tmp(val);
  if (abs(val) > pow(2, width-1))
    std::cout << "overflow at " << mstr << std::endl;
  return tmp;
}

int get_missing_layer(pca::matrixpcaconst<double> c);
void get_eta_ranges (const int towerid, double & etamin, double & etamax);

template<typename T>
bool quick_write_pcaconst_to_file (const std::vector<pca::matrixpcaconst<T> > & vct, 
    const char * filename)
{
  std::ofstream outf;
  outf.open(filename);
  
  outf << vct.size() << std::endl;
  for (unsigned int i = 0; i != vct.size(); ++i)
  {
    double min, max;
    outf << pca::matrixpcaconst<T>::const_type_to_string(vct[i].get_const_type()) 
      << std::endl;
    outf << vct[i].get_towerid() << std::endl;
    if (vct[i].get_layersids() == "")
      outf << "unspecified" << std::endl;
    else
      outf << vct[i].get_layersids() << std::endl;
    outf << pca::matrixpcaconst<T>::sector_type_to_string(vct[i].get_sector_type()) 
      << std::endl;
    outf << pca::matrixpcaconst<T>::plane_type_to_string(vct[i].get_plane_type()) 
      << std::endl;
    outf << pca::matrixpcaconst<T>::ttype_to_string(vct[i].get_ttype()) 
      << std::endl;
    vct[i].get_ptrange(min, max);
    outf.precision(10);
    outf << std::scientific << min << " " << max << std::endl;
    vct[i].get_etarange(min, max);
    outf << std::scientific << min << " " << max << std::endl;
    outf << vct[i].get_chargesign () << std::endl;
    outf << vct[i].n_rows() << std::endl;
    outf << vct[i].n_cols() << std::endl;
  
    outf.precision(10);
    if ((typeid(T) == typeid(double)) || 
        (typeid(T) == typeid(float)))
    {
      for (unsigned int j = 0; j<vct[i].n_rows(); ++j)
        for (unsigned int k = 0; k<vct[i].n_cols(); ++k)
          outf << std::scientific << vct[i](j, k) << std::endl;
    }
    else
    {
      for (unsigned int j = 0; j<vct[i].n_rows(); ++j)
        for (unsigned int k = 0; k<vct[i].n_cols(); ++k)
          outf << vct[i](j, k) << std::endl;
    }
 
  }
  
  outf.close();

  return true;
}

