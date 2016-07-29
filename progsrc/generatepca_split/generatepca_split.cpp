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
  std::cerr << " -T, --int-bitewise              : Integer bitewise mode on" << std::endl;
  std::cerr << " -l, --correlation               : compute and print correlation" << std::endl;
  std::cerr << " -p, --dump-allcoords            : dump all stub coordinates to a file" << std::endl;
  std::cerr << " -d, --dump-bankfiles            : dump all coordinates files and more extracted from the rootfile" 
    << std::endl;
  std::cerr << " -X, --max-num-oftracks=[n]      : stop reading root file after n tracks" << std::endl;
  std::cerr << std::endl;
  std::cerr << " -z, --rz-plane                  : use rz plane view (fit eta and z0)" << std::endl;
  std::cerr << " -r, --rphi-plane                : use r-phi plane view (fit pt and phi)" << std::endl;
  std::cerr << " -a, --relative                  : use relative coordinates (compute min values)" << std::endl;
  std::cerr << " -b, --relative-values=[v1;v2]   : use relative coordinates (using v1 (phi or z) and v2 (r) as min)" 
    << std::endl;
  std::cerr << std::endl;
  std::cerr << " -k, --check-layersids           : check exact layers sequence (is_a_valid_layers_seq for seq list)" 
    << std::endl;
  std::cerr << std::endl;
  std::cerr << " -f, --five-hits=[\"sequence\"]    : build constants for 5 / 6, specify the sequence " << std::endl;
  std::cerr << "                                     it will use \"real 5 out of 6\" tracks " << std::endl;
  std::cerr << " -y, --fk-five-hits=[layerid]    : build constants for 5 / 6, specify the layr to be removed " 
    << std::endl;
  std::cerr << "                                   it will use 6 layers tracks, removing a layer " << std::endl;
  std::cerr << std::endl;
  std::cerr << " -g, --charge-sign=[+/-]         : use only + particle or - paricle (again both planes) " << std::endl;
  std::cerr << " -t, --eta-range=\"etamin;etamax\" : specify the eta range to use " << std::endl;
  std::cerr << " -n, --pt-range=\"ptmin;ptmax\"    : specify the pt range to use " << std::endl;
  std::cerr << " -m, --phi-range=\"phimin;phimax\" : specify the phi range to use " << std::endl;
  std::cerr << " -o, --z0-range=\"z0min;z0max\"    : specify the z0 range to use " << std::endl;
  std::cerr << " -u, --d0-range=\"d0min;d0max\"    : specify the d0 range to use " << std::endl;
  std::cerr << std::endl;
  std::cerr << " -x, --exclude-s-module          : exclude 2S-module (last three layer) so 6 coordinates " << 
    "instead of 12 (rz)" << std::endl;                                  
  std::cerr << " -c, --use-only-3-layers         : use three leyers ..." << std::endl;
  std::cerr << " -B, --write-binfiles            : will wite the PCA contants also as bin files " << std::endl;
  std::cerr << " -D, --towerid=[num]             : specify towid to be wriiten in the file " << std::endl;
  std::cerr << " -R, --region-type=[num]         : specify region-type 0=BARREL, 1=HYBRID, 2=ENDCAP " << std::endl;
  std::cerr << "                                   BARREL is the default " << std::endl;

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
    pca::pcafitter & fitter, 
    pca::rootfilereader & rootrdr,
    bool verbose, bool writebinfiles,
    bool intbitewise, int towerid,
    int regiontype)
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

  if (writebinfiles)
  {
    pca::write_armmat(cfname.c_str(), cmtx);
    pca::write_armvct(qfname.c_str(), q);
    pca::write_armmat(afname.c_str(), amtx);
    pca::write_armvct(kfname.c_str(), kivec);
    pca::write_armmat(vfname.c_str(), vmtx);
    pca::write_armvct(cmfname.c_str(), coordmvec);
  }

  if (intbitewise)
  {
    pca::matrixpcaconst<int32_t> 
      pcmtx(cmtx.n_rows, cmtx.n_cols), 
      pqvct(q.n_rows, q.n_cols), 
      pamtx(amtx.n_rows, amtx.n_cols), 
      pkvct(kivec.n_rows, kivec.n_cols);
    
    double ptmin, ptmax, etamin, etamax;
    rootrdr.get_ptlimits(ptmin, ptmax);
    rootrdr.get_etalimits(etamin, etamax);
    assert(rootrdr.get_rphiplane() != rootrdr.get_rzplane());
    
    pca::armamat_to_pcamat (cmtx, pcmtx);
    pcmtx.set_const_type (pca::CMTX);
    pcmtx.set_layersids (rootrdr.get_actualseq().c_str());
    
    if (regiontype == ISBARREL)
      pcmtx.set_sector_type (pca::BARREL);
    else if (regiontype == ISHYBRID)
      pcmtx.set_sector_type (pca::HYBRID);
    else if (regiontype == ISENDCAP)
      pcmtx.set_sector_type (pca::ENDCAP);

    pcmtx.set_towerid (towerid);
    pcmtx.set_ttype (pca::INTEGPT);
    pcmtx.set_chargesign(rootrdr.get_chargesign());
    if (rootrdr.get_rphiplane())
      pcmtx.set_plane_type (pca::RPHI);
    else if (rootrdr.get_rzplane())
      pcmtx.set_plane_type (pca::RZ);
    pcmtx.set_ptrange (ptmin, ptmax);
    pcmtx.set_etarange (etamin, etamax); 
    
    write_pcaconst_to_file (pcmtx, "pca_const.txt");
    
    pca::armamat_to_pcamat (q, pqvct);
    pqvct.set_const_type (pca::QVEC);
    pqvct.set_layersids (rootrdr.get_actualseq().c_str());
    if (regiontype == ISBARREL)
      pqvct.set_sector_type (pca::BARREL);
    else if (regiontype == ISHYBRID)
      pqvct.set_sector_type (pca::HYBRID);
    else if (regiontype == ISENDCAP)
      pqvct.set_sector_type (pca::ENDCAP);
    pqvct.set_towerid (towerid);
    pqvct.set_ttype (pca::INTEGPT);
    pqvct.set_chargesign(rootrdr.get_chargesign());
    if (rootrdr.get_rphiplane())
      pqvct.set_plane_type (pca::RPHI);
    else if (rootrdr.get_rzplane())
      pqvct.set_plane_type (pca::RZ);
    pqvct.set_ptrange (ptmin, ptmax);
    pqvct.set_etarange (etamin, etamax); 
    
    write_pcaconst_to_file (pqvct, "pca_const.txt");
    
    pca::armamat_to_pcamat (amtx, pamtx);
    pamtx.set_const_type (pca::AMTX);
    pamtx.set_layersids (rootrdr.get_actualseq().c_str());
    if (regiontype == ISBARREL)
      pamtx.set_sector_type (pca::BARREL);
    else if (regiontype == ISHYBRID)
      pamtx.set_sector_type (pca::HYBRID);
    else if (regiontype == ISENDCAP)
      pamtx.set_sector_type (pca::ENDCAP);
    pamtx.set_towerid (towerid);
    pamtx.set_ttype (pca::INTEGPT);
    pamtx.set_chargesign(rootrdr.get_chargesign());
    if (rootrdr.get_rphiplane())
      pamtx.set_plane_type (pca::RPHI);
    else if (rootrdr.get_rzplane())
      pamtx.set_plane_type (pca::RZ);
    pamtx.set_ptrange (ptmin, ptmax);
    pamtx.set_etarange (etamin, etamax); 
    
    write_pcaconst_to_file (pamtx, "pca_const.txt");
    
    pca::armamat_to_pcamat (kivec, pkvct);
    pkvct.set_const_type (pca::KVEC);
    pkvct.set_layersids (rootrdr.get_actualseq().c_str());
    if (regiontype == ISBARREL)
      pkvct.set_sector_type (pca::BARREL);
    else if (regiontype == ISHYBRID)
      pkvct.set_sector_type (pca::HYBRID);
    else if (regiontype == ISENDCAP)
      pkvct.set_sector_type (pca::ENDCAP);
    pkvct.set_towerid (towerid);
    pkvct.set_ttype (pca::INTEGPT);
    pkvct.set_chargesign(rootrdr.get_chargesign());
    if (rootrdr.get_rphiplane())
      pkvct.set_plane_type (pca::RPHI);
    else if (rootrdr.get_rzplane())
      pkvct.set_plane_type (pca::RZ);
    pkvct.set_ptrange (ptmin, ptmax);
    pkvct.set_etarange (etamin, etamax); 
    
    write_pcaconst_to_file (pkvct, "pca_const.txt");
  }
  else
  {
    pca::matrixpcaconst<double> 
      pcmtx(cmtx.n_rows, cmtx.n_cols), 
      pqvct(q.n_rows, q.n_cols), 
      pamtx(amtx.n_rows, amtx.n_cols), 
      pkvct(kivec.n_rows, kivec.n_cols);
    
    double ptmin, ptmax, etamin, etamax;
    rootrdr.get_ptlimits(ptmin, ptmax);
    rootrdr.get_etalimits(etamin, etamax);
    assert(rootrdr.get_rphiplane() != rootrdr.get_rzplane());
    
    pca::armamat_to_pcamat (cmtx, pcmtx);
    pcmtx.set_const_type (pca::CMTX);
    pcmtx.set_layersids (rootrdr.get_actualseq().c_str());
    if (regiontype == ISBARREL)
      pcmtx.set_sector_type (pca::BARREL);
    else if (regiontype == ISHYBRID)
      pcmtx.set_sector_type (pca::HYBRID);
    else if (regiontype == ISENDCAP)
      pcmtx.set_sector_type (pca::ENDCAP);
    pcmtx.set_towerid (towerid);
    pcmtx.set_ttype (pca::FLOATPT);
    pcmtx.set_chargesign(rootrdr.get_chargesign());
    if (rootrdr.get_rphiplane())
      pcmtx.set_plane_type (pca::RPHI);
    else if (rootrdr.get_rzplane())
      pcmtx.set_plane_type (pca::RZ);
    pcmtx.set_ptrange (ptmin, ptmax);
    pcmtx.set_etarange (etamin, etamax); 
    
    write_pcaconst_to_file (pcmtx, "pca_const.txt");
    
    pca::armamat_to_pcamat (q, pqvct);
    pqvct.set_const_type (pca::QVEC);
    pqvct.set_layersids (rootrdr.get_actualseq().c_str());
    if (regiontype == ISBARREL)
      pqvct.set_sector_type (pca::BARREL);
    else if (regiontype == ISHYBRID)
      pqvct.set_sector_type (pca::HYBRID);
    else if (regiontype == ISENDCAP)
      pqvct.set_sector_type (pca::ENDCAP);
    pqvct.set_towerid (towerid);
    pqvct.set_ttype (pca::FLOATPT);
    pqvct.set_chargesign(rootrdr.get_chargesign());
    if (rootrdr.get_rphiplane())
      pqvct.set_plane_type (pca::RPHI);
    else if (rootrdr.get_rzplane())
      pqvct.set_plane_type (pca::RZ);
    pqvct.set_ptrange (ptmin, ptmax);
    pqvct.set_etarange (etamin, etamax); 
    
    write_pcaconst_to_file (pqvct, "pca_const.txt");
    
    pca::armamat_to_pcamat (amtx, pamtx);
    pamtx.set_const_type (pca::AMTX);
    pamtx.set_layersids (rootrdr.get_actualseq().c_str());
    if (regiontype == ISBARREL)
      pamtx.set_sector_type (pca::BARREL);
    else if (regiontype == ISHYBRID)
      pamtx.set_sector_type (pca::HYBRID);
    else if (regiontype == ISENDCAP)
      pamtx.set_sector_type (pca::ENDCAP);
    pamtx.set_towerid (towerid);
    pamtx.set_ttype (pca::FLOATPT);
    pamtx.set_chargesign(rootrdr.get_chargesign());
    if (rootrdr.get_rphiplane())
      pamtx.set_plane_type (pca::RPHI);
    else if (rootrdr.get_rzplane())
      pamtx.set_plane_type (pca::RZ);
    pamtx.set_ptrange (ptmin, ptmax);
    pamtx.set_etarange (etamin, etamax); 
    
    write_pcaconst_to_file (pamtx, "pca_const.txt");
    
    pca::armamat_to_pcamat (kivec, pkvct);
    pkvct.set_const_type (pca::KVEC);
    pkvct.set_layersids (rootrdr.get_actualseq().c_str());
    if (regiontype == ISBARREL)
      pkvct.set_sector_type (pca::BARREL);
    else if (regiontype == ISHYBRID)
      pkvct.set_sector_type (pca::HYBRID);
    else if (regiontype == ISENDCAP)
      pkvct.set_sector_type (pca::ENDCAP);
    pkvct.set_towerid (towerid);
    pkvct.set_ttype (pca::FLOATPT);
    pkvct.set_chargesign(rootrdr.get_chargesign());
    if (rootrdr.get_rphiplane())
      pkvct.set_plane_type (pca::RPHI);
    else if (rootrdr.get_rzplane())
      pkvct.set_plane_type (pca::RZ);
    pkvct.set_ptrange (ptmin, ptmax);
    pkvct.set_etarange (etamin, etamax); 
    
    write_pcaconst_to_file (pkvct, "pca_const.txt");
  }
}

# ifndef __CINT__
int main (int argc, char ** argv)
{
  gROOT->ProcessLine("#include <vector>");

  pca::pcafitter fitter; 

  bool rzplane = false;
  bool rphiplane = false;
  bool correlation = false;
  bool writebinfiles = false;
  bool savecheckfiles = false;
  bool checklayersids = false;
  bool printallcoords = false;
  bool userelativecoord = false;
  bool usefakefiveoutofsix = false;

  int regiontype = 0;
  int chargesign = 0;
  int numoflayers = 6;
  int layeridtorm = -1;

  unsigned int maxnumoftracks = (unsigned int) INFINITY;

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

  bool intbitewise = false;
  bool use3layers = false;

  int towerid = -99;

  while (1)
  {
    int c, option_index;
    static struct option long_options[] = {
      {"help", 0, NULL, 'h'},
      {"version", 0, NULL, 'v'},
      {"verbose", 0, NULL, 'V'},
      {"int-bitewise", 0, NULL, 'T'}, 
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
      {"check-layersids", 0, NULL, 'k'},
      {"relative", 0, NULL, 'a'},
      {"five-hits", 1, NULL, 'f'},
      {"relative-values", 1, NULL, 'b'},
      {"dump-bankfiles", 0, NULL, 'd'},
      {"fk-five-hits", 1, NULL, 'y'},
      {"max-num-oftracks", 1, NULL, 'X'},
      {"write-binfiles", 1, NULL, 'B'},
      {"towerid", 1, NULL, 'D'},
      {"use-only-3-layers", 0, NULL, 'c'},
      {"region-type", 1, NULL, 'R'},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "hvVTlpg:zrxn:t:m:o:u:kaf:b:dy:X:B:D:cR:", 
        long_options, &option_index);

    if (c == -1)
      break;

    switch (c)
    {
      case 'h':
        usage (argv[0]);
        break;
      case 'v':
        std::cout << "Version: " << pca::pcafitter::get_version_string() << std::endl;
        exit(1);
        break;
      case 'V':
        verbose = true;
        break;
      case 'T':
        intbitewise = true;
        break;
      case 'l':
        correlation = true;
        break;
      case 'p':
        printallcoords = true;
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
      case 'k':
        checklayersids = true;
        break;
      case 'a':
        userelativecoord = true;
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
      case 'd':
        savecheckfiles = true;
        break;
      case 'y':
        usefakefiveoutofsix = true;
        layeridtorm = atoi(optarg);
        break;
      case 'X':
        maxnumoftracks = atoi(optarg);
        break;
      case 'B':
        writebinfiles = true;
        break;
      case 'D':
        towerid = atoi(optarg);
        break;
      case 'c':
        use3layers = true;
        break;
      case 'R':
        regiontype = atoi(optarg);
        if ((regiontype != 0) &&
            (regiontype != 1) &&
            (regiontype != 2))
        {
          std::cerr << "Regiontype is wrong " << std::endl;
          return EXIT_FAILURE;
        }
        break;
      default:
        usage (argv[0]);
        break;
    } 
  }

  fitter.set_useintbitewise(intbitewise);

  if (towerid == -99)
  {
    std::cerr << "Towid is mandatory for XY rotation" << std::endl;
    return EXIT_FAILURE;
  }

  if (numoflayers == 5)
  {
    if (usefakefiveoutofsix)
    {
      std::cerr << "Wrong options, cannot use both options together" << std::endl;
      return EXIT_FAILURE;
    }

    if (regiontype == ISBARREL)
    {
      if (!pca::validate_barrel_sequence_5 (sequence))
      {
        std::cerr << "Wrong sequence" << std::endl;
        return EXIT_FAILURE;
      }
    }
    else
    {
      std::cerr << "Not yet implemented" << std::endl;
      return EXIT_FAILURE;
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

  if (usefakefiveoutofsix)
  {
    if (use3layers)
    {
      std::cerr << "Not yet implemented" << std::endl;
      return EXIT_FAILURE; 
    }

    if (excludesmodule)
      fitter.set_coordim (2*2);
    else
      fitter.set_coordim (2*5);
  }
  else
  {
    if (numoflayers == 5)
    {
      if (use3layers)
      {
        std::cerr << "Not yet implemented" << std::endl;
        return EXIT_FAILURE; 
      }

      if (excludesmodule)
        fitter.set_coordim (2*2);
      else
        fitter.set_coordim (2*5);
    }
    else if (numoflayers == 6)
    {
      if (excludesmodule)
        fitter.set_coordim (2*3);
      else if (use3layers)
        fitter.set_coordim (2*3);
      else
        fitter.set_coordim (2*6);
    }
    else 
    {
      std::cerr << "Can use 5 or 6 layers" << std::endl;
      return EXIT_FAILURE;
    }
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

  // leggere file coordinate tracce simulate plus parametri
  if (!file_exists(filename))
  {
    std::cerr << "Inout file does not exist" << std::endl;
    return EXIT_FAILURE;
  }
                  
  arma::mat coordin, paramin;
  arma::vec ptvals, etavals;

  std::cout << "Reading data from " << filename << " file " << std::endl;

  pca::rootfilereader rootrdr;

  rootrdr.set_useintbitewise(intbitewise);

  rootrdr.set_specificseq (sequence.c_str());
  rootrdr.set_maxnumoflayers(numoflayers);

  rootrdr.set_filename(filename);
 
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
  rootrdr.set_region_type(regiontype);

  // should be removed 
  if (regiontype == ISHYBRID) // test single layers seq
  {
    rootrdr.set_specificseq("567181920");
  }
  //maxnumoftracks = 100000;
  rootrdr.set_maxnumoftracks(maxnumoftracks);
  if (use3layers)
  {
    if (regiontype != ISBARREL)
    {
      std::cerr << "Cannot be used in non BARREL regions " << std::endl;
      return EXIT_FAILURE;
    }

    std::set<int> layers;

    layers.insert(5);
    layers.insert(8);
    layers.insert(10);

    rootrdr.set_use3layers(layers);
  }


  // needed expecially for tow 19 20 27 28
  // only for the genration
  if ((towerid == 19) || (towerid == 20) || 
      (towerid == 27) || (towerid == 28))
  {
    rootrdr.apply_rotation_to_phi(true);
    rootrdr.apply_rotation_to_xy(true);
  }

  rootrdr.set_towid(towerid);

  rootrdr.set_fkfiveoutofsix(usefakefiveoutofsix, 
      layeridtorm);

  rootrdr.set_savecheckfiles(savecheckfiles);

  if (!rootrdr.reading_from_root_file (fitter, paramin, coordin, 
        ptvals, etavals))
  {
    std::cerr << rootrdr.get_errmsg() << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Extracted layers seq: " << rootrdr.get_actualseq() << std::endl;

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
      std::cout << "Correlation param " << i << " coord ";
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
      fitter, rootrdr, verbose, writebinfiles, 
      intbitewise, towerid, regiontype);

  return EXIT_SUCCESS;
}
#endif
