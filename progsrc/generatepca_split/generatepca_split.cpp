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
  std::cerr << "usage: " << name << " [options] coordinatesfile " << std::endl;
  std::cerr << std::endl;
  std::cerr << " -h, --help                      : display this help and exit" << std::endl;
  std::cerr << " -v, --version                   : print version and exit" << std::endl;
  std::cerr << " -j, --jump-tracks               : generate the constants using only even tracks" << std::endl;
  std::cerr << " -p, --dump-allcoords            : dump all stub coordinates to a file" << std::endl;
  std::cerr << " -z, --rz-plane                  : use rz plane view" << std::endl;
  std::cerr << " -r, --rphi-plane                : use r-phi plane view" << std::endl;
  std::cerr << " -e, --not-use-charge            : do not read charge from coordinatesfile, by default " << std::endl;
  std::cerr << "                                   we will  use it if rphi-plane has been selected" << std::endl; 
  std::cerr << " -g, --charge-sign=[+/-]         : use only + particle or - paricle " << std::endl;
  std::cerr << " -t, --eta-range=\"etamin;etamax\" : specify the eta range to use " << std::endl;
  std::cerr << " -n, --pt-range=\"ptmin;ptmax\"    : specify the pt range to use " << std::endl;
  std::cerr << " -x, --exclude-s-module          : exclude S-module (last three layer) so 6 coordinates inseatd of 12 " 
    << std::endl;                                  
  std::cerr << " -d, --use-d0                    : use also d0 param in r-phi plane " << std::endl;

  exit(1);
}

void perform_main_computation (const arma::mat & coord, const arma::mat & param, 
    const std::string & cfname, const std::string & qfname,
    pca::pcafitter & fitter)
{
  // ordered 
  arma::vec eigval;
  // by row or by column ?
  arma::mat eigvec;

  std::cout << "Compute correlation mtx" << std::endl;
  arma::mat coordm = arma::zeros<arma::mat>(fitter.get_coordim());
  arma::mat hca = arma::zeros<arma::mat>(fitter.get_coordim(),
      fitter.get_coordim());
  arma::vec eigvaltmp = arma::zeros<arma::vec>(fitter.get_coordim());

  eigvec = arma::zeros<arma::mat>(fitter.get_coordim(),
      fitter.get_coordim());
  eigval = arma::zeros<arma::vec>(fitter.get_coordim());

  hca = arma::cov(coord);

  std::cout << "Eigensystem" << std::endl;
  arma::eig_sym(eigvaltmp, eigvec, hca);
 
  for (int i=0; i<fitter.get_coordim(); ++i)
    eigval(i) = eigvaltmp(fitter.get_coordim()-i-1);

  double totval = 0.0e0;
  for (int i=0; i<fitter.get_coordim(); ++i)
    totval += eigval(i);

  std::cout << "Eigenvalues: " << std::endl;
  double totvar = 0.0e0; 
  for (int i=0; i<fitter.get_coordim(); ++i)
  {
    if (i < fitter.get_paramdim())
      totvar += 100.0e0*(eigval(i)/totval);

    std::cout << i+1 << " ==> " << 100.0e0*(eigval(i)/totval) 
              << "% value: " << eigval(i) <<  std::endl;
  }
  std::cout << "PARAMDIM eigenvalues: " << totvar << std::endl;

  std::cout << fitter.get_paramdim() << " X " << fitter.get_coordim() << std::endl;

  arma::mat cmtx = arma::zeros<arma::mat>(fitter.get_paramdim(),
      fitter.get_coordim());
  arma::rowvec q = arma::zeros<arma::rowvec>(fitter.get_paramdim());

  std::cout << "Compute PCA constants " << std::endl;
  fitter.compute_pca_constants (param,
      coord, cmtx, q);

  std::cout << "Write constant to file" << std::endl;
  pca::write_armmat(cfname.c_str(), cmtx);
  pca::write_armvct(qfname.c_str(), q);
}

int main (int argc, char ** argv)
{
  pca::pcafitter fitter; 

  bool useonlyeven = false;
  bool printallcoords = false;
  bool rzplane = false;
  bool rphiplane = false;
  bool usecharge = true;
  bool usealsod0 = false;

  int chargesign = 0;

  double etamin = -1.0e0 * INFINITY, etamax = +1.0e0 * INFINITY;
  double ptmin = -1.0e0 * INFINITY, ptmax = +1.0e0 * INFINITY;

  std::vector<std::string> tokens;

  bool excludesmodule = false;

  while (1)
  {
    int c, option_index;
    static struct option long_options[] = {
      {"help", 0, NULL, 'h'},
      {"version", 0, NULL, 'v'},
      {"jump-tracks", 0, NULL, 'j'},
      {"dump-allcoords", 0, NULL, 'p'},
      {"rz-plane", 0, NULL, 'z'},
      {"rphi-plane", 0, NULL, 'r'},
      {"not-use-charge", 0, NULL, 'e'},
      {"charge-sign", 1, NULL, 'g'},
      {"eta-range", 1, NULL, 't'},
      {"exclude-s-module", 0, NULL, 'x'},
      {"use-d0", 0, NULL, 'd'},
      {"pt-range", 1, NULL, 'n'},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "dxehvjpzrg:t:n:", long_options, &option_index);

    if (c == -1)
      break;

    switch (c)
    {
      case 'd':
        usealsod0 = true;
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
      case 'j':
        useonlyeven = true;
        break;
      case 'h':
        usage (argv[0]);
        break;
      case 'v':
        std::cout << "Version: " << pca::pcafitter::get_version_string() << std::endl;
        exit(1);
        break;
      case 'e':
        usecharge = false;
        break;
      default:
        usage (argv[0]);
        break;
    } 
  }

  if (optind >= argc) 
    usage (argv[0]);

  if ((rzplane && rphiplane) ||
      (!rzplane && !rphiplane))
  {
    std::cerr << "r-phi or r-z plane ?" << std::endl;
    usage (argv[0]);
  }

  // R-z
  if (excludesmodule && rzplane)
    fitter.set_coordim (2*3);
  else
    fitter.set_coordim (2*6);

  if (usealsod0 && rphiplane)
    fitter.set_paramdim(3);
  else
    fitter.set_paramdim(2);

  if (rzplane)
  {
    if (!fitter.set_paramidx(SPLIT_COTTETHAIDX, "cot(tetha)"))
    {
      std::cerr << fitter.get_errmsg() << std::endl;
      return EXIT_FAILURE;
    }
    if (!fitter.set_paramidx(SPLIT_Z0IDX, "z0"))
    {
      std::cerr << fitter.get_errmsg() << std::endl;
      return EXIT_FAILURE;
    }
  }
  else if (rphiplane)
  {
    if (!fitter.set_paramidx(SPLIT_PHIIDX, "phi"))
    {
      std::cerr << fitter.get_errmsg() << std::endl;
      return EXIT_FAILURE;
    }

    if (usealsod0)
    {
      if (!fitter.set_paramidx(SPLIT_D0IDX, "d0"))
      {
        std::cerr << fitter.get_errmsg() << std::endl;
        return EXIT_FAILURE;
      }
    }

    if (usecharge)
    {
      if (!fitter.set_paramidx(SPLIT_ONEOVERPTIDX, "q/pt"))
      {
        std::cerr << fitter.get_errmsg() << std::endl;
        return EXIT_FAILURE;
      }
    }
    else
    {
      if (!fitter.set_paramidx(SPLIT_ONEOVERPTIDX, "1/pt"))
      {
        std::cerr << fitter.get_errmsg() << std::endl;
        return EXIT_FAILURE;
      }
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
                  
  int num_of_line = pca::numofline(filename);
  std::cout << "file has " << num_of_line << " line " << std::endl;
  int num_of_ent_read = (num_of_line-1)/ENTDIM;

  int num_of_ent = num_of_ent_read;

  if (useonlyeven)
  {
    if (num_of_ent_read % 2)
      num_of_ent = (num_of_ent_read-1)/2;
    else
      num_of_ent = num_of_ent_read/2;
  }

  std::cout << "file has " << num_of_ent << " entries " << std::endl;

  // non perfomante ma easy to go
  arma::mat coordin, paramin;
  coordin.set_size(num_of_ent, fitter.get_coordim());
  paramin.set_size(num_of_ent, fitter.get_paramdim());

  // leggere file coordinate tracce simulate plus parametri
  std::cout << "Reading data from " << filename << " file " << std::endl;

  if (!pca::reading_from_file_split (fitter, filename, paramin, coordin, 
         num_of_ent_read, useonlyeven, false, rzplane, rphiplane, 
         etamin, etamax, ptmin, ptmax, usecharge, chargesign, excludesmodule, 
         usealsod0))
    return EXIT_FAILURE;

  std::cout << "Using " << paramin.n_rows << " tracks" << std::endl;

  std::cout << "Writing parameters to files" << std::endl;

  std::ostringstream cfname, qfname; 

  if (rzplane)
  {
    pca::write_to_file("cottetha.txt", paramin, SPLIT_COTTETHAIDX);
    pca::write_to_file("z0.txt", paramin, SPLIT_Z0IDX);
    cfname << "c.rz.bin";
    qfname << "q.rz.bin";

  }
  else if (rphiplane)
  {
    pca::write_to_file("phi.txt", paramin, SPLIT_PHIIDX);
    pca::write_to_file("oneoverpt.txt", paramin, SPLIT_ONEOVERPTIDX);
    cfname << "c.rphi.bin";
    qfname << "q.rphi.bin";
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
      cfname.str(), qfname.str(), fitter);

  return EXIT_SUCCESS;
}
