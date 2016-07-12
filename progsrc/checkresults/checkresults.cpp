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
  std::cerr << "usage: " << name << " [options] resultsfile.txt " << std::endl;
  std::cerr << std::endl;
  std::cerr << " -h, --help                       : display this help and exit" << std::endl;
  std::cerr << " -v, --version                    : print version and exit" << std::endl;

  exit(1);
}


# ifndef __CINT__
int main (int argc, char ** argv)
{
  gROOT->ProcessLine("#include <vector>");

  while (1)
  {
    int c, option_index;
    static struct option long_options[] = {
      {"help", 0, NULL, 'h'},
      {"version", 0, NULL, 'v'},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "hv", 
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
      default:
        usage (argv[0]);
        break;
    } 
  }

  if (optind >= argc) 
    usage (argv[0]);

  char * filename = (char *) alloca (strlen(argv[optind]) + 1);
  strcpy (filename, argv[optind]);

  // leggere file coordinate tracce simulate plus parametri
  if (!pca::file_exists(filename))
  {
    std::cerr << "Inout file does not exist" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Reading data from " << filename << " file " << std::endl;

  std::ifstream fin (filename);

  fin.ignore (1024, '\n');

  double mmstdev, mpstdev, ptmin, ptmax;

  ptmin = 3.0;
  ptmax = 7.0;

  int nbins = 2000;
  mmstdev = -1.0;
  mpstdev = 1.0;
  TH1D *hist_qoverpt = new TH1D("hist_diff_qoverpt","q/pt diff histogram",nbins, 
        mmstdev, mpstdev);

  nbins = 100000;
  mmstdev = -5.0;
  mpstdev = 5.0;
  TH1D *hist_phi = new TH1D("hist_diff_phi","phi diff histogram",nbins, 
      mmstdev, mpstdev);

  while (!fin.eof())
  { 
    double pt, qpt_orig, qpt_fitt, pt_diff, phi_orig, 
           phi_fitt, phi_diff, fake;
    fin >> pt >> qpt_orig >> qpt_fitt >> pt_diff >> phi_orig >> 
      phi_fitt >> phi_diff >> fake;
 
    if ((pt >= ptmin) && (pt <= ptmax)) 
    {
      hist_qoverpt->Fill((Double_t) ((qpt_orig-qpt_fitt) / qpt_orig));
      hist_phi->Fill((double_t) phi_diff);
    }
  }

  mmstdev = -5.0;
  mpstdev = 5.0;
  hist_phi->Fit("gaus","","", mmstdev, mpstdev);
 
  mmstdev = -1.0;
  mpstdev = 1.0;
  hist_qoverpt->Fit("gaus","","",mmstdev, mpstdev);
  
  TF1 *func_qoverpt = (TF1*)hist_qoverpt->GetFunction("gaus");
  TF1 *func_phi = (TF1*)hist_phi->GetFunction("gaus");
  
  std::cout << 
    "q/pt fitted mean " 
    << ptmin << " " << ptmax << " " <<
    func_qoverpt->GetParameter("Mean")*100.0 << " +/- " << 
    func_qoverpt->GetParError(1)*100.0 << std::endl << 
    "q/pt fitted sigma " << 
    func_qoverpt->GetParameter("Sigma")*100.0 << " +/- " <<
    func_phi->GetParError(2)*100.0 << std::endl;


  fin.close();

  return EXIT_SUCCESS;
}
#endif
