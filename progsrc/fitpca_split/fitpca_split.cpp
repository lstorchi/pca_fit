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

// lstorchi : can be included in any case 
#ifdef INTBITEWISE
// lstorchi: should we use cstdint and -std=c++11 ? 
#include "stdint.h"
#endif

// lstorchi: basi code to fit tracks, using the PCA constants generated 
//           by the related generatepca


bool build_and_compare (arma::mat & paramslt, arma::mat & coordslt, 
     arma::mat & cmtx, arma::rowvec & q, arma::mat & amtx, arma::mat & vmtx, 
     arma::rowvec & k, bool verbose, pca::pcafitter & fitter, 
     bool rzplane, bool rphiplane, arma::vec & ptvals)
{

  double ** ptrs;
  ptrs = new double* [fitter.get_paramdim()];

  double * cothetacmp = NULL, * z0cmp = NULL, * qoverptcmp = NULL, 
    * phicmp = NULL;
 
  if (rzplane)
  {
    cothetacmp = new double [(int)coordslt.n_rows];
    z0cmp = new double [(int)coordslt.n_rows];
    ptrs[PCA_COTTHETAIDX] = cothetacmp;
    ptrs[PCA_Z0IDX] = z0cmp;
  }
  else if (rphiplane)
  {
    qoverptcmp = new double [(int)coordslt.n_rows];
    phicmp = new double [(int)coordslt.n_rows];
    ptrs[PCA_ONEOVERPTIDX] = qoverptcmp;
    ptrs[PCA_PHIIDX] = phicmp;
  }

  arma::rowvec chi2values;
  chi2values.resize(coordslt.n_rows);

  if (!fitter.compute_parameters (cmtx, q, amtx, vmtx, k,
        coordslt, ptrs, fitter.get_paramdim(), 
        chi2values))
  {
    std::cerr << fitter.get_errmsg() << std::endl;
    return false;
  }

  std::ofstream myfilechi2("chi2results.txt");

  myfilechi2 << "chi2_value" << std::endl;
  for (int i=0; i<(int) chi2values.n_cols; ++i)
    myfilechi2 << chi2values(i) << std::endl;
  myfilechi2.close();

  delete [] ptrs; 

  std::ostringstream fname;
  fname << "results.txt";

  arma::running_stat<double> pcrelative[fitter.get_paramdim()];
  arma::running_stat<double> pcabsolute[fitter.get_paramdim()];
  std::ofstream myfile(fname.str().c_str());

  if (rzplane)
  {
    myfile << "pt eta_orig eta_fitt diff z0_orig z0_fitt diff" << std::endl; 
    
    for (int i=0; i<(int)coordslt.n_rows; ++i)
    {
      double thetacmp = atan(1.0e0 / cothetacmp[i]) ; 
      double etacmps = 0.0e0, tantheta2;
      tantheta2 = tan (thetacmp/2.0e0); 
      if (tantheta2 < 0.0)
        etacmps = 1.0e0 * log (-1.0e0 * tantheta2);
      else
        etacmps = -1.0e0 * log (tantheta2);

      double theta = atan(1.0e0 / paramslt(i, PCA_COTTHETAIDX));
      double etaorig = 0.0e0;
      tantheta2 = tan (theta/2.0e0);
      if (tantheta2 < 0.0)
        etaorig = 1.0e0 * log (-1.0e0 * tantheta2);
      else
        etaorig = -1.0e0 * log (tantheta2);
      
      pcrelative[PCA_COTTHETAIDX]((etacmps - etaorig)/
          etaorig);
      pcrelative[PCA_Z0IDX]((z0cmp[i] - paramslt(i, PCA_Z0IDX))/
          paramslt(i, PCA_Z0IDX));
      
      pcabsolute[PCA_COTTHETAIDX](etacmps - etaorig);
      pcabsolute[PCA_Z0IDX](z0cmp[i] - paramslt(i, PCA_Z0IDX));
        
      myfile << ptvals(i) << " " <<
        etaorig << " " << etacmps << " " <<
        (etacmps - etaorig) << " " <<
        paramslt(i, PCA_Z0IDX) << " " << z0cmp[i] << " " <<
        (z0cmp[i] - paramslt(i, PCA_Z0IDX)) << std::endl;
      
      if (verbose)
      {
        std::cout << "For track : " << i+1 << std::endl;
        std::cout << " cotheta      fitt " << cothetacmp[i] << std::endl;
        std::cout << " cotheta      orig " << paramslt(i, PCA_COTTHETAIDX) << std::endl;
        std::cout << " theta rad    fitt " << thetacmp << std::endl;
        std::cout << " theta rad    orig " << theta << std::endl;
        std::cout << " theta deg    fitt " << thetacmp*(180.0e0/M_PI) << std::endl;
        std::cout << " theta deg    orig " << theta*(180.0e0/M_PI) << std::endl;
        std::cout << " eta          fitt " << etacmps << std::endl;
        std::cout << " eta          orig " << etaorig << std::endl;
        std::cout << " z0           fitt " << z0cmp[i] << std::endl;
        std::cout << " z0           orig " << paramslt(i, PCA_Z0IDX) << std::endl;
      }
    }
  }
  else if (rphiplane)
  {
    std::ofstream myfile(fname.str().c_str());
      myfile << "pt q/pt_orig q/pt_fitt diff phi_orig phi_fitt diff" << std::endl; 
        
    for (int i=0; i<(int)coordslt.n_rows; ++i)
    {
      double qoverptorig = paramslt(i, PCA_ONEOVERPTIDX);
      double qoverptcmp = qverptcmp[i];

      pcrelative[PCA_PHIIDX]((phicmp[i] - paramslt(i, PCA_PHIIDX))/
          paramslt(i, PCA_PHIIDX));
      pcrelative[PCA_ONEOVERPTIDX]((qoverptcmp - qoverptorig)/
          qoverptorig);
    
      pcabsolute[PCA_PHIIDX](phicmp[i] - paramslt(i, PCA_PHIIDX));
      pcabsolute[PCA_ONEOVERPTIDX](qoverptcmp - qoverptorig);
    
      myfile << ptvals(i) << " " <<
        qoverptorig << " " << qoverptcmp << " " <<
        (qoverptcmp - qoverptorig) <<  " " <<
        paramslt(i, PCA_PHIIDX) << " " << phicmp[i] << " " <<
        (phicmp[i] - paramslt(i, PCA_PHIIDX)) << std::endl;
    
      if (verbose)
      {
        std::cout << "For track : " << i+1 << std::endl;
        std::cout << " q/pt         fitt " << qoverptcmp << std::endl;
        std::cout << " q/pt         orig " << qoverptorig << std::endl;
        std::cout << " phi          fitt " << phicmp[i] << std::endl;
        std::cout << " phi          orig " << paramslt(i, PCA_PHIIDX) << std::endl;
      }
    }
  }

  myfile.close();

  for (int i=0; i<fitter.get_paramdim(); ++i)
  {
     std::cout << "For " << fitter.paramidx_to_string(i) << " error " << 
       pcabsolute[i].mean() << " " << pcabsolute[i].stddev() << std::endl;

     std::cout << "For " << fitter.paramidx_to_string(i) << " error " << 
       100.0 * pcrelative[i].mean() << " % " << 100.0 * pcrelative[i].stddev() << 
       " % " << std::endl;
  }

  if (rzplane)
  {
    delete [] cothetacmp;
    delete [] z0cmp;
  }
  else if (rphiplane)
  {
    delete [] qoverptcmp;
    delete [] phicmp;
  }

  return true;
} 

void usage (char * name)
{
  std::cerr << "usage: " << name << " [options] coordinatesfile " << std::endl;
  std::cerr << std::endl;
  std::cerr << " -h, --help                      : display this help and exit" << std::endl;
  std::cerr << " -V, --verbose                   : verbose option on" << std::endl;
  std::cerr << " -v, --version                   : print version and exit" << std::endl;
  std::cerr << " -c, --cmtx=[fillename]          : CMTX filename [default is c.[rz/rphi].bin]" << std::endl;
  std::cerr << " -q, --qvct=[fillename]          : QVCT filename [default is q.[rz/rphi].bin]" << std::endl;
  std::cerr << " -c, --amtx=[fillename]          : AMTX filename [default is a.[rz/rphi].bin]" << std::endl;
  std::cerr << " -y, --kvct=[fillename]          : KVCT filename [default is k.[rz/rphi].bin]" << std::endl;
  std::cerr << " -q, --vmtx=[fillename]          : VMTX filename [default is v.[rz/rphi].bin]" << std::endl;
  std::cerr << std::endl;
  std::cerr << " -z, --rz-plane                  : use rz plane view (fit eta and z0)" << std::endl;
  std::cerr << " -r, --rphi-plane                : use r-phi plane view (fit ot and phi)" << std::endl;
  std::cerr << " -a, --relative                  : use relative coordinates (compute min values)" << std::endl;
  std::cerr << " -b, --relative-values=[v1;v2]   : use relative coordinates (using v1 (phi or z) and v2 (r) as min)" 
    << std::endl;
  std::cerr << std::endl;
  std::cerr << " -k, --check-layersids           : check exact layers sequence (is_a_valid_layers_seq for seq list)" 
    << std::endl;
  std::cerr << " -g, --charge-sign=[+/-]         : use only + particle or - paricle (again both planes)" << std::endl;
  std::cerr << " -t, --eta-range=\"etamin;etamax\" : specify the eta range to use " << std::endl;
  std::cerr << " -n, --pt-range=\"ptmin;ptmax\"    : specify the pt range to use " << std::endl;
  std::cerr << " -m, --phi-range=\"phimin;phimax\" : specify the phi range to use " << std::endl;
  std::cerr << " -o, --z0-range=\"z0min;z0max\"    : specify the z0 range to use " << std::endl;
  std::cerr << " -u, --d0-range=\"d0min;d0max\"    : specify the d0 range to use " << std::endl;
  std::cerr << std::endl;
  std::cerr << " -x, --exclude-s-module          : exclude S-module (last three layer) so 6 " << 
    "coordinates inseatd of 12 " << std::endl;

  exit(1);
}


# ifndef __CINT__
int main (int argc, char ** argv)
{
  gROOT->ProcessLine("#include <vector>");

  pca::pcafitter fitter;

  std::string qfname = "";
  std::string cfname = "";
  std::string afname = "";
  std::string vfname = "";
  std::string kfname = "";
  std::string subsec = "";
  std::string sublad = "";
  bool verbose = false;
  bool rzplane = false;
  bool rphiplane = false;
  bool checklayersids = false;
  bool savecheckfiles = false;
  bool userelativecoord = false;

  double etamin = -1.0e0 * INFINITY, etamax = +1.0e0 * INFINITY;
  double ptmin = -1.0e0 * INFINITY, ptmax = +1.0e0 * INFINITY;
  double phimin = -1.0e0 * INFINITY, phimax = +1.0e0 * INFINITY;
  double z0min = -1.0e0 * INFINITY, z0max = +1.0e0 * INFINITY;
  double d0min = -1.0e0 * INFINITY, d0max = +1.0e0 * INFINITY;
  double coord1min = std::numeric_limits<double>::infinity();
  double coord2min = std::numeric_limits<double>::infinity();

  int chargesign = 0;

  std::vector<std::string> tokens;

  bool excludesmodule = false;

  while (1)
  {
    int c, option_index;
    static struct option long_options[] = {
      {"help", 0, NULL, 'h'},
      {"cmtx", 1, NULL, 'c'},
      {"amtx", 1, NULL, 'A'},
      {"vmtx", 1, NULL, 'B'},
      {"qvct", 1, NULL, 'q'},
      {"kvct", 1, NULL, 'y'},
      {"verbose", 0, NULL, 'V'},
      {"version", 0, NULL, 'v'},
      {"jump-tracks", 0, NULL, 'j'},
      {"rz-plane", 0, NULL, 'z'},
      {"rphi-plane", 0, NULL, 'r'},
      {"charge-sign", 1, NULL, 'g'},
      {"exclude-s-module", 0, NULL, 'x'},
      {"pt-range", 1, NULL, 'n'},
      {"eta-range", 1, NULL, 't'},
      {"phi-range", 1, NULL, 'm'},
      {"z0-range", 1, NULL, 'o'},
      {"d0-range", 1, NULL, 'u'},
      {"check-layersids", 1, NULL, 'k'},
      {"relative", 0, NULL, 'a'},
      {"relative-values", 1, NULL, 'b'},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "kxzrhaVy:b:A:B:t:g:c:q:n:s:m:o:u", 
        long_options, &option_index);

    if (c == -1)
      break;

    switch (c)
    {
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
      case 'V':
        verbose = true;
        break;
      case 'v':
        std::cout << "Version: " << pca::pcafitter::get_version_string() << std::endl;
        exit(1);
        break;
      case 'h':
        usage (argv[0]);
        break;
      case'c':
        cfname = optarg;
        break;
      case 'q':
        qfname = optarg;
        break;
      case'A':
        afname = optarg;
        break;
      case 'B':
        vfname = optarg;
        break;
      case 'y':
        kfname = optarg;
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

  // leggere file coordinate tracce e file costanti PCA
  // N righe di 9 double sono le coordinate
  // matrice C e vettore q sono le costanti
  
  arma::mat cmtx, amtx, vmtx;
  arma::rowvec q, k;

  // leggere file coordinate tracce simulate plus parametri
  if (!pca::file_exists(filename))
  {
    std::cerr << "Inout file does not exist" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Reading data from " << filename << " file " << std::endl;

  arma::mat coord, param;
  arma::vec ptvals;

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

  rootrdr.set_savecheckfiles(false);

  if (!rootrdr.reading_from_root_file (fitter, param, coord, 
        ptvals))
  {
    std::cerr << rootrdr.get_errmsg() << std::endl;
    return EXIT_FAILURE;
  }

  if (userelativecoord)
    pca::global_to_relative(coord, coord1min, coord2min);

  std::cout << "Using " << param.n_rows << " tracks" << std::endl;

  if (rzplane)
  {
    if (savecheckfiles)
    {
      pca::write_to_file("tofit_cottheta.txt", param, PCA_COTTHETAIDX);
      pca::write_to_file("tofit_z0.txt", param, PCA_Z0IDX);
    }

    if (cfname == "")
      cfname = "c.rz.bin";

    if (qfname == "")
      qfname = "q.rz.bin";

    if (afname == "")
      afname = "a.rz.bin";

    if (vfname == "")
      vfname = "v.rz.bin";

    if (kfname == "")
      kfname = "k.rz.bin";
  }
  else if (rphiplane)
  {
    if (savecheckfiles)
    {
      pca::write_to_file("tofit_phi.txt", param, PCA_PHIIDX);
      pca::write_to_file("tofit_oneoverpt.txt", param, PCA_ONEOVERPTIDX);
    }

    if (cfname == "")
      cfname = "c.rphi.bin";

    if (qfname == "")
      qfname = "q.rphi.bin";

    if (afname == "")
      afname = "a.rphi.bin";

    if (vfname == "")
      vfname = "v.rphi.bin";

    if (kfname == "")
      kfname = "k.rphi.bin";
  }

  if (pca::file_exists(afname.c_str()))
  {
    std::cout << "Reading " << afname << std::endl;
    pca::read_armmat(afname.c_str(), amtx);
  }
  else
  {
    std::cerr << afname << " does not exist" << std::endl;
    return 1;
  }

  if (pca::file_exists(cfname.c_str()))
  {
    std::cout << "Reading " << cfname << std::endl;
    pca::read_armmat(cfname.c_str(), cmtx);
  }
  else
  {
    std::cerr << cfname << " does not exist" << std::endl;
    return 1;
  }

  if (pca::file_exists(vfname.c_str()))
  {
    std::cout << "Reading " << vfname << std::endl;
    pca::read_armmat(vfname.c_str(), vmtx);
  }
  else
  {
    std::cerr << vfname << " does not exist" << std::endl;
    return 1;
  }

  if (pca::file_exists(kfname.c_str()))
  {
    std::cout << "Reading " << kfname << std::endl;
    pca::read_armvct(kfname.c_str(), k);
  }
  else
  {
    std::cerr << kfname << " does not exist" << std::endl;
    return 1;
  }

  if (pca::file_exists(qfname.c_str()))
  {
    std::cout << "Reading " << qfname << std::endl;
    pca::read_armvct(qfname.c_str(), q);
  }
  else
  {
    std::cerr << qfname << " does not exist" << std::endl;
    return 1;
  }

  if (!build_and_compare (param, coord, cmtx, q, amtx, vmtx, k,
        verbose, fitter, rzplane, rphiplane, ptvals))
    return EXIT_FAILURE;

#ifdef INTBITEWISE
  //To cross check whether constants have been read correctly in int16_t mode
  std::cout << "Constants Used: C matrix: " << std::endl;
  std::cout << cmtx;
  std::cout << "Constants Used: q matrix: " << std::endl;
  std::cout << q;
#endif

  return EXIT_SUCCESS;
}
#endif
