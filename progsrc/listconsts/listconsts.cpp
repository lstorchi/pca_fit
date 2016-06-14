#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>
#include <set>

// Loriano: let's try Armadillo quick code 
#include <armadillo>
#include <cassert>

#include <getopt.h>
#include <unistd.h>
#include <alloca.h>

#include <pcafitter.hpp>
#include <pcaffunctype.hpp>
#include <rootfilereader.hpp>

#include "TROOT.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"

#ifdef INTBITEWISEFIT
#include "stdint.h"
#endif

// lstorchi: basi code to fit tracks, using the PCA constants generated 
//           by the related generatepca

bool read_pca_const (const std::string & cfname)
{
  std::vector<pca::matrixpcaconst<double> > vct;
  if (read_pcaconst_from_file (vct, cfname.c_str()))
  {
    int totaldim = 0, muplus = 0, muneg = 0;
    std::vector<pca::matrixpcaconst<double> >::const_iterator it = 
      vct.begin();
    std::cout << "plane_type const_type towerid sector_type " << 
      "ttype layerseq ptmin ptmax etamin etamax chargesign dim" << std::endl;
    for (; it != vct.end(); ++it)
    {
      double ptmin, ptmax, etamin, etamax;
      std::string layerseq;
      int chargesign;

      it->get_ptrange(ptmin, ptmax);
      it->get_etarange(etamin, etamax);
      chargesign = it->get_chargesign();
      layerseq = it->get_layersids();

      std::cout 
        << pca::matrixpcaconst<double>::plane_type_to_string(it->get_plane_type()) << " " 
        << pca::matrixpcaconst<double>::const_type_to_string(it->get_const_type()) << " "
        << it->get_towerid() << " "
        << pca::matrixpcaconst<double>::sector_type_to_string(it->get_sector_type()) << " "
        << pca::matrixpcaconst<double>::ttype_to_string(it->get_ttype()) << " ";
      std::cout << layerseq << " " << ptmin << " " <<  ptmax << " " 
        << etamin << " "<< etamax << " " << chargesign << " ";
      if (pca::matrixpcaconst<double>::plane_type_to_string(it->get_plane_type()) == "RPHI")
      {
        if (chargesign > 0)
          muplus++;
        else if (chargesign < 0)
          muneg++;
        else 
        {
          std::cerr << "error charge " << std::endl;
          exit(1);
        }
      }
      int dim = it->n_rows() * it->n_cols();
      std::cout << dim << std::endl;
      totaldim += dim;
    }

    std::cout << "Total values: " << totaldim << std::endl;
    std::cout << "mu+: " << muplus << " mu-: " << muneg << std::endl;
    if (muplus != muneg)
      std::cerr << "Error  in mu+ mu- number " << std::endl;

    return true;
  }

  return false;
}

void usage (char * name)
{
  std::cerr << "usage: " << name << " constfile " << std::endl;
  std::cerr << std::endl;
  std::cerr << " -h, --help                       : display this help and exit" << std::endl;

  exit(1);
}


# ifndef __CINT__
int main (int argc, char ** argv)
{
  gROOT->ProcessLine("#include <vector>");

  pca::pcafitter fitter;

  std::string cfname = "pca_const.txt";

  while (1)
  {
    int c, option_index;
    static struct option long_options[] = {
      {"help", 0, NULL, 'h'},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "h", 
        long_options, &option_index);

    if (c == -1)
      break;

    switch (c)
    {
      case 'h':
        usage (argv[0]);
        break;
      default:
        usage (argv[0]);
        break;
    } 
  }

  if (optind >= argc) 
    usage (argv[0]);

  char * filename = (char *) alloca (strlen(argv[optind]) + 1);
  strcpy (filename, argv[optind]);

  if (pca::file_exists(filename))
  {
    std::cout << "Reading " << filename << std::endl;

    if (!read_pca_const (filename))
    {
      std::cerr << "Error in reading constants from file" << std::endl;
      return EXIT_FAILURE;
    }
  }
  else
  {
    std::cerr << filename << " does not exist" << std::endl;
    return 1;
  }

  return EXIT_SUCCESS;
}
#endif
