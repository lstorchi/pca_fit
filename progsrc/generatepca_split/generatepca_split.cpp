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
#include <sys/stat.h>

#include <getopt.h>
#include <unistd.h>
#include <alloca.h>

#include <pcafitter.hpp>
#include <pcaffunctype.hpp>
#include <rootfilereader.hpp>

#include "TROOT.h"

#define MINDIMLINIT 25

// lstorchi: basic quick code to generate PCA constants

namespace
{
  bool file_exists(const std::string& filename)
  {
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
      return true;
                
    return false;
  }
}

void usage (char * name)
{
  std::cerr << "usage: " << name << " [options] rootcoordinatesfile " << std::endl;
  std::cerr << std::endl;
  std::cerr << " -h, --help                      : display this help and exit" << std::endl;
  std::cerr << " -v, --version                   : print version and exit" << std::endl;
  std::cerr << " -V, --verbose                   : verbose mode on" << std::endl;
  std::cerr << " -l, --correlation               : compute and print correlation" << std::endl;
  std::cerr << " -p, --dump-allcoords            : dump all stub coordinates to a file" << std::endl;
  std::cerr << std::endl;
  std::cerr << " -z, --rz-plane                  : use rz plane view (fit eta and z0)" << std::endl;
  std::cerr << " -r, --rphi-plane                : use r-phi plane view (fit pt and phi)" << std::endl;
  std::cerr << " -a, --relative                  : use relative coordinates (compute min values)" << std::endl;
  std::cerr << " -b, --relative-values=[v1;v2]   : use relative coordinates (using v1 (phi or z) and v2 (r) as min)" 
    << std::endl;
  std::cerr << std::endl;
  std::cerr << " -k, --check-layersids           : check exact layers sequence (is_a_valid_layers_seq for seq list)" 
    << std::endl;
  std::cerr << " -f, --five-hits=[\"sequence\"]    : build constants for 5 / 6, specify the sequence, " << std::endl;
  std::cerr << "                                   if -x is used you should specify three layers sequence"  << std::endl;
  std::cerr << " -g, --charge-sign=[+/-]         : use only + particle or - paricle (again both planes) " << std::endl;
  std::cerr << " -t, --eta-range=\"etamin;etamax\" : specify the eta range to use " << std::endl;
  std::cerr << " -n, --pt-range=\"ptmin;ptmax\"    : specify the pt range to use " << std::endl;
  std::cerr << " -m, --phi-range=\"phimin;phimax\" : specify the phi range to use " << std::endl;
  std::cerr << " -o, --z0-range=\"z0min;z0max\"    : specify the z0 range to use " << std::endl;
  std::cerr << " -u, --d0-range=\"d0min;d0max\"    : specify the d0 range to use " << std::endl;
  std::cerr << std::endl;
  std::cerr << " -x, --exclude-s-module          : exclude S-module (last three layer) so 6 coordinates " << 
    "inseatd of 12 " << std::endl;                                  

  exit(1);
}

void perform_main_computation (const arma::mat & coord, 
    const arma::mat & param, 
    const std::string & cfname, 
    const std::string & qfname, 
    const std::string & afname,
    const std::string & vfname, 
    const std::string & kfname, 
    const std::string & cmfname ,
    pca::pcafitter & fitter, bool verbose)
{
  std::cout << fitter.get_paramdim() << " X " << fitter.get_coordim() << std::endl;

  arma::mat cmtx = arma::zeros<arma::mat>(fitter.get_paramdim(),
      fitter.get_coordim());
  arma::rowvec q = arma::zeros<arma::rowvec>(fitter.get_paramdim());
  arma::mat vmtx = arma::zeros<arma::mat>(fitter.get_coordim(),
      fitter.get_coordim());
  arma::mat amtx = arma::zeros<arma::mat>(
      fitter.get_coordim()-fitter.get_paramdim(),
      fitter.get_coordim());

  int verbositylevel = 1;
  if (verbose) 
    verbositylevel = 2;

  arma::rowvec kivec, coordmvec;
  std::cout << "Compute PCA constants " << std::endl;
  if (!fitter.compute_pca_constants (param,
         coord, cmtx, q, vmtx, amtx, kivec, coordmvec, 
         verbositylevel))
  {
    std::cerr << "compute_pca_constants error" << std::endl;
    return;
  }

  std::cout << "Write constant to file" << std::endl;
  pca::write_armmat(cfname.c_str(), cmtx);
  pca::write_armvct(qfname.c_str(), q);
  pca::write_armmat(afname.c_str(), amtx);
  pca::write_armmat(vfname.c_str(), vmtx);
  pca::write_armvct(kfname.c_str(), kivec);
  pca::write_armvct(cmfname.c_str(), coordmvec);
}

# ifndef __CINT__
int main (int argc, char ** argv)
{
  gROOT->ProcessLine("#include <vector>");

  pca::pcafitter fitter; 

  bool rzplane = false;
  bool rphiplane = false;
  bool correlation = false;
  bool savecheckfiles = false;
  bool checklayersids = false;
  bool printallcoords = false;
  bool userelativecoord = false;

  int chargesign = 0;
  int numoflayers =  6;

  double etamin = -1.0e0 * INFINITY, etamax = +1.0e0 * INFINITY;
  double ptmin = -1.0e0 * INFINITY, ptmax = +1.0e0 * INFINITY;
  double phimin = -1.0e0 * INFINITY, phimax = +1.0e0 * INFINITY;
  double z0min = -1.0e0 * INFINITY, z0max = +1.0e0 * INFINITY;
  double d0min = -1.0e0 * INFINITY, d0max = +1.0e0 * INFINITY;
  double coord1min = std::numeric_limits<double>::infinity();
  double coord2min = std::numeric_limits<double>::infinity();

  std::vector<std::string> tokens;
  std::string sequence;

  bool excludesmodule = false;
  bool verbose = false;

  while (1)
  {
    int c, option_index;
    static struct option long_options[] = {
      {"help", 0, NULL, 'h'},
      {"version", 0, NULL, 'v'},
      {"verbose", 0, NULL, 'V'},
      {"correlation", 0, NULL, 'l'},
      {"dump-allcoords", 0, NULL, 'p'},
      {"charge-sign", 1, NULL, 'g'},
      {"rz-plane", 0, NULL, 'z'},
      {"rphi-plane", 0, NULL, 'r'},
      {"exclude-s-module", 0, NULL, 'x'},
      {"pt-range", 1, NULL, 'n'},
      {"eta-range", 1, NULL, 't'},
      {"phi-range", 1, NULL, 'm'},
      {"z0-range", 1, NULL, 'o'},
      {"d0-range", 1, NULL, 'u'},
      {"check-layersids", 1, NULL, 'k'},
      {"relative", 0, NULL, 'a'},
      {"five-hits", 1, NULL, 'f'},
      {"relative-values", 1, NULL, 'b'},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "aVlkxhvpzrb:g:t:n:m:o:u:f:", long_options, &option_index);

    if (c == -1)
      break;

    switch (c)
    {
      case 'f':
        numoflayers = 5;
        sequence = optarg;
        break;
      case 'b':
        userelativecoord = true;
        tokens.clear();
        pca::tokenize (optarg, tokens, ";");
        if (tokens.size() != 2)
          usage (argv[0]);

        coord1min = atof(tokens[0].c_str());
        coord2min = atof(tokens[1].c_str());

        break;
 
      case 'a':
        userelativecoord = true;
        break;
      case 'k':
        checklayersids = true;
        break;
      case 'l':
        correlation = true;
        break;
      case 'm':
        tokens.clear();
        pca::tokenize (optarg, tokens, ";");
        if (tokens.size() != 2)
          usage (argv[0]);

        phimin = atof(tokens[0].c_str());
        phimax = atof(tokens[1].c_str());

        break;
      case 'o':
        tokens.clear();
        pca::tokenize (optarg, tokens, ";");
        if (tokens.size() != 2)
          usage (argv[0]);

        z0min = atof(tokens[0].c_str());
        z0max = atof(tokens[1].c_str());

        break;
      case 'u':
        tokens.clear();
        pca::tokenize (optarg, tokens, ";");
        if (tokens.size() != 2)
          usage (argv[0]);

        d0min = atof(tokens[0].c_str());
        d0max = atof(tokens[1].c_str());

        break;
      case 'V':
        verbose = true;
        break;
      case 'x':
        excludesmodule = true;
        break;
      case 'n':
        tokens.clear();
        pca::tokenize (optarg, tokens, ";");
        if (tokens.size() != 2)
          usage (argv[0]);

        ptmin = atof(tokens[0].c_str());
        ptmax = atof(tokens[1].c_str());

        break;
      case 't':
        tokens.clear();
        pca::tokenize (optarg, tokens, ";");
        if (tokens.size() != 2)
          usage (argv[0]);

        etamin = atof(tokens[0].c_str());
        etamax = atof(tokens[1].c_str());

        break;
      case 'g':
        if (strlen(optarg) > 1)
          usage (argv[0]);
        
        if (*optarg == '-')
          chargesign = -1;
        else if (*optarg == '+')
          chargesign = +1;
        else
          usage (argv[0]);

        break;
      case 'z':
        rzplane = true;
        break;
      case 'r':
        rphiplane = true;
        break;
      case 'p':
        printallcoords = true;
        break;
      case 'h':
        usage (argv[0]);
        break;
      case 'v':
        std::cout << "Version: " << pca::pcafitter::get_version_string() << std::endl;
        exit(1);
        break;
      default:
        usage (argv[0]);
        break;
    } 
  }

  if (numoflayers == 5)
  {
    if (!pca::validate_barrel_sequence_5 (excludesmodule, 
          sequence))
    {
      std::cerr << "Wrong sequence" << std::endl;
      return EXIT_FAILURE;
    }
  }

  TODO 

  if (optind >= argc) 
    usage (argv[0]);

  if ((rzplane && rphiplane) ||
      (!rzplane && !rphiplane))
  {
    std::cerr << "r-phi or r-z plane ?" << std::endl;
    usage (argv[0]);
  }

  if (excludesmodule)
    fitter.set_coordim (2*3);
  else
    fitter.set_coordim (2*6);

  fitter.set_paramdim(2);

  if (rzplane)
  {
    if (!fitter.set_paramidx(PCA_COTTHETAIDX, "cot(theta)"))
    {
      std::cerr << fitter.get_errmsg() << std::endl;
      return EXIT_FAILURE;
    }
    if (!fitter.set_paramidx(PCA_Z0IDX, "z0"))
    {
      std::cerr << fitter.get_errmsg() << std::endl;
      return EXIT_FAILURE;
    }
  }
  else if (rphiplane)
  {
    if (!fitter.set_paramidx(PCA_PHIIDX, "phi"))
    {
      std::cerr << fitter.get_errmsg() << std::endl;
      return EXIT_FAILURE;
    }
    
    if (!fitter.set_paramidx(PCA_ONEOVERPTIDX, "q/pt"))
    {
      std::cerr << fitter.get_errmsg() << std::endl;
      return EXIT_FAILURE;
    }
  }

  char * filename = (char *) alloca (strlen(argv[optind]) + 1);
  strcpy (filename, argv[optind]);

  // leggere file coordinate tracce simulate plus parametri
  if (!file_exists(filename))
  {
    std::cerr << "Inout file does not exist" << std::endl;
    return EXIT_FAILURE;
  }
                  
  arma::mat coordin, paramin;
  arma::vec ptvals;

  std::cout << "Reading data from " << filename << " file " << std::endl;

  pca::rootfilereader rootrdr;

  rootrdr.set_filename(filename);

  int mxnumoflayers = 6;
  rootrdr.set_maxnumoflayers(mxnumoflayers);
 
  rootrdr.set_rzplane(rzplane);
  rootrdr.set_rphiplane(rphiplane);
  rootrdr.set_etalimits(etamin, etamax);
  rootrdr.set_ptlimits(ptmin, ptmax);
  rootrdr.set_chargesign(chargesign);
  rootrdr.set_excludesmodule(excludesmodule);
  rootrdr.set_philimits(phimin, phimax);
  rootrdr.set_z0limits(z0min, z0max);
  rootrdr.set_d0limits(d0min, d0max);
  rootrdr.set_verbose(verbose);
  rootrdr.set_checklayersids(checklayersids);

  rootrdr.set_savecheckfiles(savecheckfiles);

  if (!rootrdr.reading_from_root_file (fitter, paramin, coordin, 
        ptvals))
  {
    std::cerr << rootrdr.get_errmsg() << std::endl;
    return EXIT_FAILURE;
  }

  if (userelativecoord)
    pca::global_to_relative(coordin, coord1min, coord2min);

  if ((coordin.n_rows == 0) || (paramin.n_rows == 0))
  {
    std::cout << "No tracks" << std::endl;
    return EXIT_FAILURE;
  }

  /* Try correlation */ 
  if (coordin.n_rows  != paramin.n_rows)
  {
    std::cerr << "num of rows should be the same" << std::endl;
    return EXIT_FAILURE;
  }

  if (correlation)
  {
    for (int i=0; i<(int)paramin.n_cols; ++i)
    {
      double avgval = 0.0;
      std::cout << "Corralation param " << i << " coord ";
      for (int j=0; j<(int)coordin.n_cols; ++j)
      {
        arma::vec x, y;
        x.set_size(coordin.n_rows);
        y.set_size(coordin.n_rows);
    
        for (int k=0; k<(int)coordin.n_rows; ++k)
        {
          x(k) = paramin(k,i);
          y(k) = coordin(k,j); 
        }
    
        double corrval;
        arma::mat corrmat = arma::cor(x,y);
        corrval = corrmat(0,0);
        avgval += corrval;
        std::cout << corrval << " "; 
    
      }
    
      std::cout << "(" << avgval/coordin.n_cols << ")" << std::endl;
    }

    for (int i=0; i<(int)paramin.n_cols; ++i)
    {
      double avgval = 0.0;
      std::cout << "Corralation param " << i << " param ";
      for (int j=0; j<(int)paramin.n_cols; ++j)
      {
        if (j != i)
        {
          arma::vec x, y;
          x.set_size(paramin.n_rows);
          y.set_size(paramin.n_rows);
          
          for (int k=0; k<(int)paramin.n_rows; ++k)
          {
            x(k) = paramin(k,i);
            y(k) = paramin(k,j); 
          }
          
          double corrval;
          arma::mat corrmat = arma::cor(x,y);
          corrval = corrmat(0,0);
          avgval += corrval;
          std::cout << corrval << " "; 
        }
      }
    
      std::cout << "(" << avgval/paramin.n_cols << ")" << std::endl;
    }
  }

  std::cout << "Using " << paramin.n_rows << " tracks" << std::endl;
  std::cout << "Writing parameters to files" << std::endl;

  std::ostringstream cfname, qfname, afname, vfname, kfname, coordmfname; 

  if (rzplane)
  {
    if (savecheckfiles)
    {
      pca::write_to_file("cottheta.txt", paramin, PCA_COTTHETAIDX);
      pca::write_to_file("z0.txt", paramin, PCA_Z0IDX);
    }

    cfname << "c.rz.bin";
    qfname << "q.rz.bin";
    afname << "a.rz.bin";
    vfname << "v.rz.bin";
    kfname << "k.rz.bin";
    coordmfname << "cm.rz.bin";
  }
  else if (rphiplane)
  {
    if (savecheckfiles)
    {
      pca::write_to_file("phi.txt", paramin, PCA_PHIIDX);
      pca::write_to_file("oneoverpt.txt", paramin, PCA_ONEOVERPTIDX);
    }

    cfname << "c.rphi.bin";
    qfname << "q.rphi.bin";
    afname << "a.rphi.bin";
    vfname << "v.rphi.bin";
    kfname << "k.rphi.bin";
    coordmfname << "cm.rphi.bin";
  }

  if (printallcoords)
  {
    std::cout << "Printout coordinates " << std::endl;
    std::ofstream myfilect("allcoords.txt");
    for (int i=0; i<(int)coordin.n_rows; ++i)
      for (int j=0; j<fitter.get_coordim(); j=j+2)
        myfilect << coordin(i, j) << " " << 
                    coordin(i, j+1) << std::endl;
    myfilect.close();
  }

  perform_main_computation (coordin, paramin,
      cfname.str(), qfname.str(), afname.str() ,
      vfname.str(), kfname.str(), coordmfname.str(),
      fitter, verbose);

  return EXIT_SUCCESS;
}
#endif
