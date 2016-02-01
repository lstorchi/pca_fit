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

bool build_and_compare (arma::mat & paramslt, arma::mat & coordslt, 
     arma::mat & cmtx, arma::rowvec & q, arma::mat & amtx, arma::mat & vmtx, 
     arma::rowvec & k, arma::rowvec & cm, bool verbose, pca::pcafitter & fitter, 
     bool rzplane, bool rphiplane, arma::vec & ptvals)
{
  int nbins = 100;

#ifdef INTBITEWISEFIT
  int32_t ** ptrs;
  ptrs = new int32_t* [fitter.get_paramdim()];

  int32_t * cothetacmp = NULL, * z0cmp = NULL, * qoverptcmp = NULL,
    * phicmp = NULL;
#else
  double ** ptrs;
  ptrs = new double* [fitter.get_paramdim()];

  double * cothetacmp = NULL, * z0cmp = NULL, * qoverptcmp = NULL,
    * phicmp = NULL;
#endif
  
  if (rzplane)
  {
#ifdef INTBITEWISEFIT
    cothetacmp = new int32_t [(int)coordslt.n_rows];
    z0cmp = new int32_t [(int)coordslt.n_rows];
#else
    cothetacmp = new double [(int)coordslt.n_rows];
    z0cmp = new double [(int)coordslt.n_rows];
#endif
    ptrs[PCA_COTTHETAIDX] = cothetacmp;
    ptrs[PCA_Z0IDX] = z0cmp;
  }
  else if (rphiplane)
  {
#ifdef INTBITEWISEFIT
    qoverptcmp = new int32_t [(int)coordslt.n_rows];
    phicmp = new int32_t [(int)coordslt.n_rows];
#else
    qoverptcmp = new double [(int)coordslt.n_rows];
    phicmp = new double [(int)coordslt.n_rows];
#endif
    ptrs[PCA_ONEOVERPTIDX] = qoverptcmp;
    ptrs[PCA_PHIIDX] = phicmp;
  }

  arma::rowvec chi2values;
  chi2values.resize(coordslt.n_rows);

  if (!fitter.compute_parameters (cmtx, q, amtx, vmtx, k, cm,
        coordslt, ptrs, fitter.get_paramdim(), 
        chi2values))
  {
    std::cerr << fitter.get_errmsg() << std::endl;
    return false;
  }

  //std::ofstream myfilechi2("chi2results.txt");
  //myfilechi2 << "chi2_value" << std::endl;
  //for (int i=0; i<(int) chi2values.n_cols; ++i)
  //  myfilechi2 << chi2values(i) << std::endl;
  //myfilechi2.close();

  delete [] ptrs; 

  std::ostringstream fname;
  fname << "results.txt";

  arma::running_stat<double> pcrelative[fitter.get_paramdim()];
  arma::running_stat<double> pcabsolute[fitter.get_paramdim()];
  arma::running_stat<double> chi2stat;
  std::ofstream myfile(fname.str().c_str());

  assert(chi2values.n_cols == coordslt.n_rows);

  if (rzplane)
  {
#ifdef INTBITEWISEFIT
    myfile << "pt cot0_orig cot0_fitt diff z0_orig z0_fitt diff chi2" << std::endl;
#else
    myfile << "pt eta_orig eta_fitt diff z0_orig z0_fitt diff chi2" << std::endl;
#endif
    
    arma::rowvec etadiffvct(coordslt.n_rows), 
      z0diffvct(coordslt.n_rows);
   
    for (int i=0; i<(int)coordslt.n_rows; ++i)
    {
      double thetacmp = atan(1.0e0 / cothetacmp[i]) ; 
      double etacmp = 0.0e0, tantheta2;
      tantheta2 = tan (thetacmp/2.0e0); 
      if (tantheta2 < 0.0)
        etacmp = 1.0e0 * log (-1.0e0 * tantheta2);
      else
        etacmp = -1.0e0 * log (tantheta2);

      double theta = atan(1.0e0 / paramslt(i, PCA_COTTHETAIDX));
      double etaorig = 0.0e0;
      tantheta2 = tan (theta/2.0e0);
      if (tantheta2 < 0.0)
        etaorig = 1.0e0 * log (-1.0e0 * tantheta2);
      else
        etaorig = -1.0e0 * log (tantheta2);

      double etadiff =  (etacmp - etaorig);
      etadiffvct(i) = etadiff;

      double z0cmps = z0cmp[i];
      double z0orig = paramslt(i, PCA_Z0IDX);
      double z0diff = z0cmp[i] - paramslt(i, PCA_Z0IDX);
      z0diffvct(i) = z0diff;

      pcrelative[PCA_COTTHETAIDX](etadiff/etaorig);
      pcrelative[PCA_Z0IDX](z0diff/z0orig);
      
      pcabsolute[PCA_COTTHETAIDX](etadiff);
      pcabsolute[PCA_Z0IDX](z0diff);

      chi2stat(chi2values(i));
        
#ifdef INTBITEWISEFIT
      myfile << ptvals(i) << "   " <<
	paramslt(i, PCA_COTTHETAIDX) << "   " << cothetacmp[i] << "   " <<
	(cothetacmp[i] - paramslt(i, PCA_COTTHETAIDX)) << " " <<
	z0orig << " " << z0cmps << " " <<
	(z0cmps - z0orig) << chi2values(i) << std::endl;
#else
      myfile << ptvals(i) << "   " <<
	paramslt(i, PCA_COTTHETAIDX) << "   " << cothetacmp[i] << "   " <<
	(cothetacmp[i] - paramslt(i, PCA_COTTHETAIDX)) << " " <<
	z0orig << " " << z0cmps << " " <<
	(z0cmps - z0orig) << " " << chi2values(i) << std::endl;
#endif
      
      if (verbose)
      {
        std::cout << "For track : " << i+1 << std::endl;
        std::cout << " cotheta      fitt " << cothetacmp[i] << std::endl;
        std::cout << " cotheta      orig " << paramslt(i, PCA_COTTHETAIDX) << std::endl;
        std::cout << " theta rad    fitt " << thetacmp << std::endl;
        std::cout << " theta rad    orig " << theta << std::endl;
        std::cout << " theta deg    fitt " << thetacmp*(180.0e0/M_PI) << std::endl;
        std::cout << " theta deg    orig " << theta*(180.0e0/M_PI) << std::endl;
        std::cout << " eta          fitt " << etacmp << std::endl;
        std::cout << " eta          orig " << etaorig << std::endl;
        std::cout << " z0           fitt " << z0cmps << std::endl;
        std::cout << " z0           orig " << z0orig << std::endl;
        std::cout << " chi2              " << chi2values(i) << std::endl;
      }
    }

    TH1D *hist_z0 = new TH1D("hist_diff_z0","z0 diff histogram",nbins, 
        z0diffvct.min(), z0diffvct.max());
    TH1D *hist_eta = new TH1D("hist_diff_eta","eta diff histogram",nbins, 
        etadiffvct.min(), etadiffvct.max());

    for (int i=0; i<(int)coordslt.n_rows; ++i)
    {
      hist_z0->Fill((Double_t) z0diffvct(i));
      hist_eta->Fill((double_t) etadiffvct(i));
    }

    hist_z0->Fit("gaus","","",z0diffvct.min(),z0diffvct.max());
    hist_eta->Fit("gaus","","",etadiffvct.min(),etadiffvct.max());

    TF1 *func_eta = (TF1*)hist_eta->GetFunction("gaus");
    TF1 *func_z0 = (TF1*)hist_z0->GetFunction("gaus");

    std::cout << 
      "Eta fitted mean: " << func_eta->GetParameter("Mean") << " +/- " << 
      func_eta->GetParError(1) << std::endl << 
      "Eta fitted sigma: " << func_eta->GetParameter("Sigma") << " +/- " <<
      func_eta->GetParError(2) << std::endl;

    std::cout << 
      "z0 fitted mean: " << func_z0->GetParameter("Mean") << " +/- " << 
      func_z0->GetParError(1) << std::endl << 
      "z0 fitted sigma: " << func_z0->GetParameter("Sigma") << " +/- " <<
      func_z0->GetParError(2) << std::endl;

  }
  else if (rphiplane)
  {
    std::ofstream myfile(fname.str().c_str());
      myfile << "pt q/pt_orig q/pt_fitt diff phi_orig phi_fitt diff chi2" << std::endl; 

    arma::rowvec qoverptdiffvct(coordslt.n_rows), 
      phidiffvct(coordslt.n_rows);
        
    for (int i=0; i<(int)coordslt.n_rows; ++i)
    {
      double qoverptorig = paramslt(i, PCA_ONEOVERPTIDX);
#ifdef INTBITEWISEFIT
      int32_t qoverptcmps = qoverptcmp[i];
#else
      double qoverptcmps = qoverptcmp[i];
#endif
      
      double diffqoverpt = qoverptcmps - qoverptorig;

      double phiorig = paramslt(i, PCA_PHIIDX);
#ifdef INTBITEWISEFIT
      int32_t phicmps = phicmp[i];
#else
      double phicmps = phicmp[i];      
#endif      
      double diffphi = phicmps - phiorig;

      pcrelative[PCA_PHIIDX](diffphi/phiorig);
      pcrelative[PCA_ONEOVERPTIDX](diffqoverpt/qoverptorig);
    
      pcabsolute[PCA_PHIIDX](diffphi);
      pcabsolute[PCA_ONEOVERPTIDX](diffqoverpt);

      chi2stat(chi2values(i));

      qoverptdiffvct(i) = diffqoverpt/qoverptorig;
      phidiffvct(i) = diffphi;
    
      myfile << ptvals(i) << " " <<
        qoverptorig << " " << qoverptcmps << " " << diffqoverpt << " " <<
        phiorig     << " " << phicmps     << " " << diffphi     << " " << 
        chi2values(i) << std::endl;
    
      if (verbose)
      {
        std::cout << "For track : " << i+1 << std::endl;
        std::cout << " q/pt         fitt " << qoverptcmps << std::endl;
        std::cout << " q/pt         orig " << qoverptorig << std::endl;
        std::cout << " phi          fitt " << phicmps << std::endl;
        std::cout << " phi          orig " << phiorig << std::endl;
        std::cout << " chi2              " << chi2values(i) << std::endl;
      }
    }

    TH1D *hist_qoverpt = new TH1D("hist_diff_qoverpt","q/pt diff histogram",nbins, 
        qoverptdiffvct.min(), qoverptdiffvct.max());
    TH1D *hist_phi = new TH1D("hist_diff_phi","phi diff histogram",nbins, 
        phidiffvct.min(), phidiffvct.max());

    for (int i=0; i<(int)coordslt.n_rows; ++i)
    {
      hist_qoverpt->Fill((Double_t) qoverptdiffvct(i));
      hist_phi->Fill((double_t) phidiffvct(i));
    }

    hist_qoverpt->Fit("gaus","","",qoverptdiffvct.min(),
        qoverptdiffvct.max());
    hist_phi->Fit("gaus","","",phidiffvct.min(),
        phidiffvct.max());

    TF1 *func_qoverpt = (TF1*)hist_qoverpt->GetFunction("gaus");
    TF1 *func_phi = (TF1*)hist_phi->GetFunction("gaus");

    std::cout << 
      "q/pt fitted mean: " << func_qoverpt->GetParameter("Mean")*100.0 << " +/- " << 
      func_qoverpt->GetParError(1)*100.0 << std::endl << 
      "p/pt fitted sigma: " << func_phi->GetParameter("Sigma")*100.0 << " +/- " <<
      func_phi->GetParError(2)*100.0 << std::endl;

    std::cout << 
      "Phi fitted mean: " << func_phi->GetParameter("Mean") << " +/- " << 
      func_phi->GetParError(1) << std::endl << 
      "Phi fitted sigma: " << func_phi->GetParameter("Sigma") << " +/- " <<
      func_phi->GetParError(2) << std::endl;
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

  std::cout << " " << std::endl;
  std::cout << "Chivalue mean " << chi2stat.mean() << " stdev " << 
    chi2stat.stddev() << std::endl;

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
  std::cerr << " -h, --help                       : display this help and exit" << std::endl;
  std::cerr << " -V, --verbose                    : verbose option on" << std::endl;
  std::cerr << " -v, --version                    : print version and exit" << std::endl;
  std::cerr << " -p, --dump-allcoords             : dump all stub coordinates to a file" << std::endl;
  std::cerr << " -c, --cmtx=[fillename]           : CMTX filename [default is c.[rz/rphi].bin]" << std::endl;
  std::cerr << " -q, --qvct=[fillename]           : QVCT filename [default is q.[rz/rphi].bin]" << std::endl;
  std::cerr << " -c, --amtx=[fillename]           : AMTX filename [default is a.[rz/rphi].bin]" << std::endl;
  std::cerr << " -y, --kvct=[fillename]           : KVCT filename [default is k.[rz/rphi].bin]" << std::endl;
  std::cerr << " -d, --cvct=[fillename]           : CVCT filename [default is cm.[rz/rphi].bin]" << std::endl;
  std::cerr << " -q, --vmtx=[fillename]           : VMTX filename [default is v.[rz/rphi].bin]" << std::endl;
  std::cerr << std::endl;                         
  std::cerr << " -z, --rz-plane                   : use rz plane view (fit eta and z0)" << std::endl;
  std::cerr << " -r, --rphi-plane                 : use r-phi plane view (fit ot and phi)" << std::endl;
  std::cerr << " -a, --relative                   : use relative coordinates (compute min values)" << std::endl;
  std::cerr << " -b, --relative-values=[v1;v2]    : use relative coordinates (using v1 (phi or z) and v2 (r) as min)" 
    << std::endl;
  std::cerr << " -f, --five-hits=[\"sequence\"]     : fit a specific 5 / 6 sequence, " << std::endl;
  std::cerr << " -l, --five-hits-lin=[\"sequence\"] : fit a specific the sequence using standard constat  " << std::endl;
  std::cerr << "                                    use linear interpolation to approximate the missed hit " << std::endl;
  std::cerr << std::endl;
  std::cerr << " -k, --check-layersids            : check exact layers sequence (is_a_valid_layers_seq for seq list)" 
    << std::endl;
  std::cerr << " -g, --charge-sign=[+/-]          : use only + particle or - paricle (again both planes)" << std::endl;
  std::cerr << " -t, --eta-range=\"etamin;etamax\"  : specify the eta range to use " << std::endl;
  std::cerr << " -n, --pt-range=\"ptmin;ptmax\"     : specify the pt range to use " << std::endl;
  std::cerr << " -m, --phi-range=\"phimin;phimax\"  : specify the phi range to use " << std::endl;
  std::cerr << " -o, --z0-range=\"z0min;z0max\"     : specify the z0 range to use " << std::endl;
  std::cerr << " -u, --d0-range=\"d0min;d0max\"     : specify the d0 range to use " << std::endl;
  std::cerr << std::endl;
  std::cerr << " -x, --exclude-s-module           : exclude S-module (last three layer) so 6 " << 
    "coordinates inseatd of 12 (rz)" << std::endl;

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
  std::string cmfname = "";
  bool verbose = false;
  bool rzplane = false;
  bool rphiplane = false;
  bool checklayersids = false;
  bool savecheckfiles = false;
  bool userelativecoord = false;
  bool lininterpolation = false;
  bool printallcoords = false;

  double etamin = -1.0e0 * INFINITY, etamax = +1.0e0 * INFINITY;
  double ptmin = -1.0e0 * INFINITY, ptmax = +1.0e0 * INFINITY;
  double phimin = -1.0e0 * INFINITY, phimax = +1.0e0 * INFINITY;
  double z0min = -1.0e0 * INFINITY, z0max = +1.0e0 * INFINITY;
  double d0min = -1.0e0 * INFINITY, d0max = +1.0e0 * INFINITY;
  double coord1min = std::numeric_limits<double>::infinity();
  double coord2min = std::numeric_limits<double>::infinity();

  std::string sequence;
  int numoflayers = 6;

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
      {"cvct", 1, NULL, 'd'},
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
      {"five-hits", 1, NULL, 'f'},
      {"five-hits-lin", 1, NULL, 'l'},
      {"dump-allcoords", 0, NULL, 'p'},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "pkxzrhaVl:f:d:y:b:A:B:t:g:c:q:n:s:m:o:u", 
        long_options, &option_index);

    if (c == -1)
      break;

    switch (c)
    {
      case 'p':
        printallcoords = true;
        break;
      case 'l':
        sequence = optarg;
        lininterpolation = true;
        break;
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
      case'd':
        cmfname = optarg;
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

  if (numoflayers == 5)
  {
    if (!pca::validate_barrel_sequence_5 (sequence))
    {
      std::cerr << "Wrong sequence" << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (lininterpolation)
  {
    if (!pca::validate_barrel_sequence_5 (sequence))
    {
      std::cerr << "Wrong sequence" << std::endl;
      return EXIT_FAILURE;
    }
  }

  if ((rzplane && rphiplane) ||
      (!rzplane && !rphiplane))
  {
    std::cerr << "r-phi or r-z plane ?" << std::endl;
    usage (argv[0]);
  }

  if (numoflayers == 5)
  {
    if (excludesmodule)
      fitter.set_coordim (2*2);
    else
      fitter.set_coordim (2*5);
  }
  else if (numoflayers == 6)
  {
    if (excludesmodule)
      fitter.set_coordim (2*3);
    else
      fitter.set_coordim (2*6);
  }
  else 
  {
    std::cerr << "Can use 5 or 6 layers" << std::endl;
    return EXIT_FAILURE;
  }

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
  arma::rowvec q, k, cm;

  // leggere file coordinate tracce simulate plus parametri
  if (!pca::file_exists(filename))
  {
    std::cerr << "Inout file does not exist" << std::endl;
    return EXIT_FAILURE;
  }

  if (rzplane)
  {
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

    if (cmfname == "")
      cmfname = "cm.rz.bin";
  }
  else if (rphiplane)
  {
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

    if (cmfname == "")
      cmfname = "cm.rphi.bin";
  }

  if (pca::file_exists(cmfname.c_str()))
  {
    std::cout << "Reading " << cmfname << std::endl;
    pca::read_armvct(cmfname.c_str(), cm);
  }
  else
  {
    std::cerr << cmfname << " does not exist" << std::endl;
    return 1;
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

  std::cout << "Reading data from " << filename << " file " << std::endl;

  arma::mat coord, param;
  arma::vec ptvals;

  pca::rootfilereader rootrdr;

  rootrdr.set_filename(filename);

  rootrdr.set_specificseq (sequence.c_str());
  rootrdr.set_maxnumoflayers(numoflayers);
 
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

  if (lininterpolation)
  {
    rootrdr.set_performlinearinterpolation (true);
    rootrdr.set_specificseq (sequence.c_str());
    rootrdr.set_maxnumoflayers (5);
  }

  if (!rootrdr.reading_from_root_file (fitter, param, coord, 
        ptvals))
  {
    std::cerr << rootrdr.get_errmsg() << std::endl;
    return EXIT_FAILURE;
  }

  if (userelativecoord)
    pca::global_to_relative(coord, coord1min, coord2min);

  if (printallcoords)
  {
    std::cout << "Printout coordinates " << std::endl;
    std::ofstream myfilect("allcoords_fit.txt");
    for (int i=0; i<(int)coord.n_rows; ++i)
      for (int j=0; j<fitter.get_coordim(); j=j+2)
        myfilect << coord(i, j) << " " << 
                    coord(i, j+1) << std::endl;
    myfilect.close();
  }

  std::cout << "Using " << param.n_rows << " tracks" << std::endl;

  if (rzplane)
  {
    if (savecheckfiles)
    {
      pca::write_to_file("tofit_cottheta.txt", param, PCA_COTTHETAIDX);
      pca::write_to_file("tofit_z0.txt", param, PCA_Z0IDX);
    }
  }
  else if (rphiplane)
  {
    if (savecheckfiles)
    {
      pca::write_to_file("tofit_phi.txt", param, PCA_PHIIDX);
      pca::write_to_file("tofit_oneoverpt.txt", param, PCA_ONEOVERPTIDX);
    }
  }

  if (!build_and_compare (param, coord, cmtx, q, amtx, vmtx, k, cm,
        verbose, fitter, rzplane, rphiplane, ptvals))
    return EXIT_FAILURE;

  std::cout << "Constants Used: C matrix: " << std::endl;
  std::cout << cmtx;
  std::cout << "Constants Used: q matrix: " << std::endl;
  std::cout << q;

#ifdef INTBITEWISEFIT
  std::cout << "Constants Used with precision:" << std::endl;
  std::cout << "Constants Used: C matrix: " << std::endl;
  for (int i=0; i<2; ++i)
    for (int j=0; j<12; ++j)
      std::cout << std::setprecision(9) << (double) cmtx(i, j) << std::endl;

  std::cout << "Constants Used: q matrix: " << std::endl;
  for (int i=0; i<2; ++i)
    std::cout << std::setprecision(9) << (double) q(i) << std::endl;
#endif
  
  return EXIT_SUCCESS;
}
#endif
