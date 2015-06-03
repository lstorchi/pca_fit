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
  std::cerr << " -V, --verbose                   : verbose mode on" << std::endl;
  std::cerr << " -l, --correlation               : compute and print correlation" << std::endl;
  std::cerr << " -p, --dump-allcoords            : dump all stub coordinates to a file" << std::endl;
  std::cerr << " -e, --not-use-charge            : do not read charge from coordinatesfile " << std::endl;
  std::cerr << std::endl;
  std::cerr << " -z, --rz-plane                  : use rz plane view (fit eta and z0)" << std::endl;
  std::cerr << " -r, --rphi-plane                : use r-phi plane view (fit pt and phi)" << std::endl;
  std::cerr << std::endl;
  std::cerr << " -j, --jump-tracks               : generate the constants using only even tracks" << std::endl;
  std::cerr << " -g, --charge-sign=[+/-]         : use only + particle or - paricle (again both planes) " << std::endl;
  std::cerr << " -t, --eta-range=\"etamin;etamax\" : specify the eta range to use " << std::endl;
  std::cerr << " -n, --pt-range=\"ptmin;ptmax\"    : specify the pt range to use " << std::endl;
  std::cerr << " -m, --phi-range=\"phimin;phimax\" : specify the phi range to use " << std::endl;
  std::cerr << " -o, --z0-range=\"z0min;z0max\"    : specify the z0 range to use " << std::endl;
  std::cerr << " -u, --d0-range=\"d0min;d0max\"    : specify the d0 range to use " << std::endl;
  std::cerr << std::endl;
  std::cerr << " -x, --exclude-s-module          : exclude S-module (last three layer) so 6 coordinates " << 
    "inseatd of 12 " << std::endl;                                  
  std::cerr << " -d, --use-d0                    : use also d0 param in both planes " << std::endl;
  std::cerr << " -X, --use-x0                    : use also x0 param in both planes " << std::endl;
  std::cerr << " -f, --fit-x0y0                  : use and fit x0 and y0 param in both planes instead of " << std::endl;
  std::cerr << "                                   eta, pt, z0, phi " << std::endl;
  std::cerr << " -s, --fit-single-param=[num]    : use and fit X param in both planes  " << std::endl;
  std::cerr << "                                   1=eta, 2=pt, 3=z0, 4=phi, 5=x0, 6=y0, 7=d0 " << std::endl;


  exit(1);
}

void perform_main_computation (const arma::mat & coord, const arma::mat & param, 
    const std::string & cfname, const std::string & qfname,
    pca::pcafitter & fitter, bool verbose)
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
      coord, cmtx, q, verbose);

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
  bool usex0y0 = false;
  bool usesingleparam = false;
  bool usealsox0 = false;
  bool correlation = false;

  int chargesign = 0;

  int singleparam=-1;

  double etamin = -1.0e0 * INFINITY, etamax = +1.0e0 * INFINITY;
  double ptmin = -1.0e0 * INFINITY, ptmax = +1.0e0 * INFINITY;
  double phimin = -1.0e0 * INFINITY, phimax = +1.0e0 * INFINITY;
  double z0min = -1.0e0 * INFINITY, z0max = +1.0e0 * INFINITY;
  double d0min = -1.0e0 * INFINITY, d0max = +1.0e0 * INFINITY;

  std::vector<std::string> tokens;

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
      {"jump-tracks", 0, NULL, 'j'},
      {"dump-allcoords", 0, NULL, 'p'},
      {"not-use-charge", 0, NULL, 'e'},
      {"charge-sign", 1, NULL, 'g'},
      {"rz-plane", 0, NULL, 'z'},
      {"rphi-plane", 0, NULL, 'r'},
      {"exclude-s-module", 0, NULL, 'x'},
      {"use-d0", 0, NULL, 'd'},
      {"use-x0", 0, NULL, 'X'},
      {"fit-x0y0", 0, NULL, 'f'}, 
      {"fit-single-param", 1, NULL, 's'}, 
      {"pt-range", 1, NULL, 'n'},
      {"eta-range", 1, NULL, 't'},
      {"phi-range", 1, NULL, 'm'},
      {"z0-range", 1, NULL, 'o'},
      {"d0-range", 1, NULL, 'u'},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "VXlfdxehvjpzrg:t:n:s:m:o:u:", long_options, &option_index);

    if (c == -1)
      break;

    switch (c)
    {
      case 'l':
        correlation = true;
        break;
      case 'X':
        usealsox0 = true;
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
      case 's':
        singleparam=atoi(optarg);
        if (!((singleparam >= 1) && (singleparam <= 7)))
          usage(argv[0]);

        usesingleparam = true;

        break;
      case 'f':
        usex0y0 = true;
        break;
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

  if ((usealsox0) && (usealsod0))
  {
    std::cerr << "use also x0 and d0 cannot be used together" << std::endl;
    usage (argv[0]);
  }

  if (usealsox0 && usex0y0)
  {
    std::cerr << "use also x0 and x0y0 cannot be used together" << std::endl;
    usage (argv[0]);
  }


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

  if (usesingleparam)
    fitter.set_paramdim(1);
  else
  {
    if ((usealsod0) || (usealsox0))
      fitter.set_paramdim(3);
    else
      fitter.set_paramdim(2);
  }

  if (rzplane)
  {
    if (usesingleparam)
    {
      if (!fitter.set_paramidx(0, pca::get_paramname_from_id(singleparam).c_str()))
      {
        std::cerr << fitter.get_errmsg() << std::endl;
        return EXIT_FAILURE;
      }
    }
    else
    {
      if (usex0y0)
      {
        if (!fitter.set_paramidx(SPLIT_X0IDX, "x0"))
        {
          std::cerr << fitter.get_errmsg() << std::endl;
          return EXIT_FAILURE;
        }
        if (!fitter.set_paramidx(SPLIT_Y0IDX, "y0"))
        {
          std::cerr << fitter.get_errmsg() << std::endl;
          return EXIT_FAILURE;
        }
      }
      else
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
      
      if (usealsod0)
      {
        if (!fitter.set_paramidx(SPLIT_D0IDX, "d0"))
        {
          std::cerr << fitter.get_errmsg() << std::endl;
          return EXIT_FAILURE;
        }
      }
      else if (usealsox0)
      {
        if (!fitter.set_paramidx(SPLIT_X0IDX_NS, "x0"))
        {
          std::cerr << fitter.get_errmsg() << std::endl;
          return EXIT_FAILURE;
        }
      }
    }
  }
  else if (rphiplane)
  {
    if (usesingleparam)
    {
      if (!fitter.set_paramidx(0, pca::get_paramname_from_id(singleparam).c_str()))
      {
        std::cerr << fitter.get_errmsg() << std::endl;
        return EXIT_FAILURE;
      }
    }
    else
    {
      if (usealsod0)
      {
        if (!fitter.set_paramidx(SPLIT_D0IDX, "d0"))
        {
          std::cerr << fitter.get_errmsg() << std::endl;
          return EXIT_FAILURE;
        }
      }
      else if (usealsox0)
      {
        if (!fitter.set_paramidx(SPLIT_X0IDX_NS, "x0"))
        {
          std::cerr << fitter.get_errmsg() << std::endl;
          return EXIT_FAILURE;
        }
      }
      
      if (usex0y0)
      {
        if (!fitter.set_paramidx(SPLIT_X0IDX, "x0"))
        {
          std::cerr << fitter.get_errmsg() << std::endl;
          return EXIT_FAILURE;
        }
        if (!fitter.set_paramidx(SPLIT_Y0IDX, "y0"))
        {
          std::cerr << fitter.get_errmsg() << std::endl;
          return EXIT_FAILURE;
        }
      }
      else
      {
        if (!fitter.set_paramidx(SPLIT_PHIIDX, "phi"))
        {
          std::cerr << fitter.get_errmsg() << std::endl;
          return EXIT_FAILURE;
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
         usealsod0, usex0y0, singleparam, phimin, phimax, z0min, z0max,
         d0min, d0max, usealsox0, verbose))
    return EXIT_FAILURE;

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
  }

  std::cout << "Using " << paramin.n_rows << " tracks" << std::endl;

  std::cout << "Writing parameters to files" << std::endl;

  std::ostringstream cfname, qfname; 

  if (rzplane)
  {
    if (usesingleparam)
    {
      std::string name = pca::get_paramname_from_id(singleparam)+".txt"; 
      pca::write_to_file(name.c_str(), paramin, 0);
    }
    else
    {
      if (usex0y0)
      {
        pca::write_to_file("x0.txt", paramin, SPLIT_X0IDX);
        pca::write_to_file("y0.txt", paramin, SPLIT_Y0IDX);
      }
      else
      {
        pca::write_to_file("cottetha.txt", paramin, SPLIT_COTTETHAIDX);
        pca::write_to_file("z0.txt", paramin, SPLIT_Z0IDX);
      }

      if (usealsod0)
        pca::write_to_file("d0.txt", paramin, SPLIT_D0IDX);
      else if (usealsox0)
        pca::write_to_file("x0.txt", paramin, SPLIT_X0IDX_NS);
    }

    cfname << "c.rz.bin";
    qfname << "q.rz.bin";
  }
  else if (rphiplane)
  {
    if (usesingleparam)
    {
      std::string name = pca::get_paramname_from_id(singleparam)+".txt"; 
      pca::write_to_file(name.c_str(), paramin, 0);
    }
    else
    {
      if (usex0y0)
      {
        pca::write_to_file("x0.txt", paramin, SPLIT_X0IDX);
        pca::write_to_file("y0.txt", paramin, SPLIT_Y0IDX);
      }
      else
      {
        pca::write_to_file("phi.txt", paramin, SPLIT_PHIIDX);
        pca::write_to_file("oneoverpt.txt", paramin, SPLIT_ONEOVERPTIDX);
      }

      if (usealsod0)
        pca::write_to_file("d0.txt", paramin, SPLIT_D0IDX);
      else if (usealsox0)
        pca::write_to_file("x0.txt", paramin, SPLIT_X0IDX_NS);
 
    }

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
      cfname.str(), qfname.str(), fitter,
      verbose);

  return EXIT_SUCCESS;
}
