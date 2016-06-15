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

// lstorchi: basic code to fit tracks, using the PCA constants generated 
//           by the related generatepca

void usage (char * name)
{
  std::cerr << "usage: " << name << " [options] coordinatesfile " << std::endl;
  std::cerr << std::endl;
  std::cerr << " -h, --help                       : display this help and exit" << std::endl;
  std::cerr << " -V, --verbose                    : verbose option on" << std::endl;
  std::cerr << " -v, --version                    : print version and exit" << std::endl;
  std::cerr << " -p, --dump-allcoords             : dump all stub coordinates to a file" << std::endl;
  std::cerr << " -c, --pca-const-files=[file1;...;filen] " << std::endl;
  std::cerr << "                                  : PCA const txt filename [default is pca_const.txt]" << std::endl;
  std::cerr << std::endl;                         
  std::cerr << " -z, --rz-plane                   : use rz plane view (fit eta and z0)" << std::endl;
  std::cerr << " -r, --rphi-plane                 : use r-phi plane view (fit ot and phi)" << std::endl;
  std::cerr << " -a, --relative                   : use relative coordinates (compute min values)" << std::endl;
  std::cerr << " -b, --relative-values=[v1;v2]    : use relative coordinates (using v1 (phi or z) and v2 (r) as min)" 
    << std::endl;
  std::cerr << std::endl; 
  std::cerr << " -f, --five-hits=[\"sequence\"]     : fit a specific 5 / 6 sequence, it will use " << std::endl;
  std::cerr << "                                    \"real 5 out of 6\" tracks " << std::endl;
  std::cerr << " -l, --five-hits-lin=[\"sequence\"] : fit a specific the sequence using standard constat  " << std::endl;
  std::cerr << "                                      use linear interpolation to approximate the missed hit " 
    << std::endl;
  std::cerr << " -w, --fk-five-hits=[layerid]     : build constants for 5 / 6, specify the layr to be removed " 
    << std::endl;
  std::cerr << "                                   it will use 6 layers tracks, removing a layer " << std::endl;
  std::cerr << std::endl;
  std::cerr << " -k, --check-layersids            : check exact layers sequence (is_a_valid_layers_seq for seq list)" 
    << std::endl;
  std::cerr << std::endl;
  std::cerr << " -x, --exclude-s-module           : exclude S-module (last three layer) so 6 " << 
    "coordinates inseatd of 12 (rz)" << std::endl;
  std::cerr << " -D, --towerid=[num]              : MANDATORY: specify towid to be used for the XY rotation " 
    << std::endl;
  std::cerr << "                                    written in the file " << std::endl;
  std::cerr << " -N, --no-results                 : results file is not written, only mean and stdev " << std::endl;
  std::cerr << "                                    are computed and reported " << std::endl; 
  std::cerr << " -X, --max-num-oftracks=[n]       : stop reading root file after n tracks" << std::endl;
  std::cerr << std::endl;

  exit(1);
}


# ifndef __CINT__
int main (int argc, char ** argv)
{
  gROOT->ProcessLine("#include <vector>");

  pca::pcafitter fitter;

  bool rzplane = false, rphiplane = true, excludesmodule = false, 
       checklayersids = false, lininterpolation = false, 
       usefakefiveoutofsix = false, printallcoords = false, 
       writeresults = true, verbose = false, userelativecoord = false;
  double coord1min = std::numeric_limits<double>::infinity();
  double coord2min = std::numeric_limits<double>::infinity();

  unsigned int maxnumoftracks = (unsigned int) INFINITY;

  int layeridtorm = -1, towerid = -99, numoflayers = 6;

  std::string sequence = "";

  std::vector<std::string> cfnames, tokens;

  while (1)
  {
    int c, option_index;
    static struct option long_options[] = {
      {"help", 0, NULL, 'h'},
      {"pca-const-files", 1, NULL, 'c'},
      {"verbose", 0, NULL, 'V'},
      {"version", 0, NULL, 'v'},
      {"rz-plane", 0, NULL, 'z'},
      {"rphi-plane", 0, NULL, 'r'},
      {"exclude-s-module", 0, NULL, 'x'},
      {"check-layersids", 0, NULL, 'k'},
      {"relative", 0, NULL, 'a'},
      {"relative-values", 1, NULL, 'b'},
      {"five-hits", 1, NULL, 'f'},
      {"five-hits-lin", 1, NULL, 'l'},
      {"fk-five-hits", 1, NULL, 'w'},
      {"dump-allcoords", 0, NULL, 'p'},
      {"towerid", 1, NULL, 'D'},
      {"max-num-oftracks", 1, NULL, 'X'},
      {"no-results", 0, NULL, 'N'},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "hc:Vvzrxkab:f:l:w:pD:X:N", 
        long_options, &option_index);

    if (c == -1)
      break;

    switch (c)
    {
      case 'h':
        usage (argv[0]);
        break;
      case'c':
        cfnames.clear();
        pca::tokenize (optarg, cfnames, ";");
        break;
      case 'V':
        verbose = true;
        break;
      case 'v':
        std::cout << "Version: " << pca::pcafitter::get_version_string() << std::endl;
        exit(1);
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
      case 'k':
        checklayersids = true;
        break;
      case 'a':
        userelativecoord = true;
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
      case 'f':
        numoflayers = 5;
        sequence = optarg;
        break;
      case 'l':
        sequence = optarg;
        lininterpolation = true;
        break;
      case 'w':
        usefakefiveoutofsix = true;
        layeridtorm = atoi(optarg);
        break;
      case 'p':
        printallcoords = true;
        break;
      case 'D':
        towerid = atoi(optarg);
        break;
      case 'X':
        maxnumoftracks = atoi(optarg);
        break;
      case 'N':
        writeresults = false;
        break;
      default:
        usage (argv[0]);
        break;
    } 
  }

  if (optind >= argc) 
    usage (argv[0]);

  if (towerid == -99)
  {
    std::cerr << "Towid is mandatory for XY rotation" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
#endif
