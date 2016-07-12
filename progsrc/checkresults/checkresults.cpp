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

bool build_and_compare (arma::mat & paramslt, arma::mat & coordslt, 
     std::vector<pca::matrixpcaconst<double> > & allconst,
     bool verbose, pca::pcafitter & fitter, 
     bool rzplane, bool rphiplane, arma::vec & ptvals, 
     int towerid, double sec_phi, 
     bool writeresults, int layeridtorm, 
     double etamin, double etamax, 
     double ptmin, double ptmax, int chargesign, 
     const std::string & layersid, 
     const std::string & pslayersid, 
     bool coarsegrainpca, 
     std::vector<pca::matrixpcaconst<double> > & cgconst)
{
  int nbins = 100;

  double ** ptrs = NULL;

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

  arma::rowvec chi2values, chi2values_fake;
  chi2values.resize(coordslt.n_rows);
  chi2values_fake.resize(0);

  if (coarsegrainpca)
  {
    if (!fitter.compute_parameters_cgpca (cgconst, allconst, 
          coordslt, paramslt, layersid, pslayersid, 
          towerid, ptrs, fitter.get_paramdim(), 
          rphiplane, chi2values))
    {
      std::cerr << fitter.get_errmsg() << std::endl;
      return false;
    }
  }
  else
  {
    if (!fitter.compute_parameters (allconst, 
          coordslt, paramslt, layersid, pslayersid, 
          towerid, ptrs, fitter.get_paramdim(), 
          rphiplane, chi2values))
    {
      std::cerr << fitter.get_errmsg() << std::endl;
      return false;
    }
  }

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
    if (writeresults)
      myfile << "pt eta_orig eta_fitt diff z0_orig z0_fitt diff chi2 chi2" << std::endl;
  
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
      
      if (writeresults)
        myfile << ptvals(i) << "   " <<
          etaorig << "   " << etacmp << "   " <<
          (etacmp - etaorig) << " " <<
          z0orig << " " << z0cmps << " " <<
          (z0cmps - z0orig) << " " << chi2values(i) << std::endl;
    
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

    double mmstdev , mpstdev;

    nbins = 2000;
    //mmstdev = pcabsolute[PCA_Z0IDX].mean() - 1.5 * pcabsolute[PCA_Z0IDX].stddev();
    //mpstdev = pcabsolute[PCA_Z0IDX].mean() + 1.5 * pcabsolute[PCA_Z0IDX].stddev();
    mmstdev = -1.0;
    mpstdev = 1.0;
    //TH1D *hist_z0 = new TH1D("hist_diff_z0","z0 diff histogram",nbins, 
    //    z0diffvct.min(), z0diffvct.max());
    TH1D *hist_z0 = new TH1D("hist_diff_z0","z0 diff histogram",nbins, 
        mmstdev, mpstdev);
 
    nbins = 20000;
    //mmstdev = pcabsolute[PCA_COTTHETAIDX].mean() - 1.5 * pcabsolute[PCA_COTTHETAIDX].stddev();
    //mpstdev = pcabsolute[PCA_COTTHETAIDX].mean() + 1.5 * pcabsolute[PCA_COTTHETAIDX].stddev();
    mmstdev = -1.0;
    mpstdev = 1.0;
    //TH1D *hist_eta = new TH1D("hist_diff_eta","eta diff histogram",nbins, 
    //    etadiffvct.min(), etadiffvct.max());
    TH1D *hist_eta = new TH1D("hist_diff_eta","eta diff histogram",nbins, 
        mmstdev, mpstdev);
    
    for (int i=0; i<(int)coordslt.n_rows; ++i)
    {
      hist_z0->Fill((Double_t) z0diffvct(i));
      hist_eta->Fill((double_t) etadiffvct(i));
    }

    //mmstdev = pcabsolute[PCA_Z0IDX].mean() - 1.5 * pcabsolute[PCA_Z0IDX].stddev();
    //mpstdev = pcabsolute[PCA_Z0IDX].mean() + 1.5 * pcabsolute[PCA_Z0IDX].stddev();
    mmstdev = -1.0;
    mpstdev = 1.0;
    //hist_z0->Fit("gaus","","",z0diffvct.min(),z0diffvct.max());
    hist_z0->Fit("gaus","","",mmstdev,mpstdev);
    mmstdev = -1.0;
    mpstdev = 1.0;
    //mmstdev = pcabsolute[PCA_COTTHETAIDX].mean() - 1.5 * pcabsolute[PCA_COTTHETAIDX].stddev();
    //mpstdev = pcabsolute[PCA_COTTHETAIDX].mean() + 1.5 * pcabsolute[PCA_COTTHETAIDX].stddev();
    //hist_eta->Fit("gaus","","",etadiffvct.min(),etadiffvct.max());
    hist_eta->Fit("gaus","","",mmstdev,mpstdev);
    
    TF1 *func_eta = (TF1*)hist_eta->GetFunction("gaus");
    TF1 *func_z0 = (TF1*)hist_z0->GetFunction("gaus");
    
    std::cout << 
      "Eta fitted mean " 
      << etamin << " " << etamax << " "
      << layeridtorm << " " 
      << func_eta->GetParameter("Mean") << " +/- " << 
      func_eta->GetParError(1) << std::endl << 
      "Eta fitted sigma " << layeridtorm << " " 
      << func_eta->GetParameter("Sigma") << " +/- " <<
      func_eta->GetParError(2) << std::endl;
    
    std::cout << 
      "z0 fitted mean " 
      << etamin << " " << etamax << " "
      << layeridtorm << " "
      << func_z0->GetParameter("Mean") << " +/- " << 
      func_z0->GetParError(1) << std::endl << 
      "z0 fitted sigma " << layeridtorm << " " 
      << func_z0->GetParameter("Sigma") << " +/- " <<
      func_z0->GetParError(2) << std::endl;
  }
  else if (rphiplane)
  {
    std::ofstream myfile(fname.str().c_str());

    if (writeresults)
      myfile << "pt q/pt_orig q/pt_fitt diff phi_orig phi_fitt diff chi2 chi2" << std::endl; 
  
    arma::rowvec qoverptdiffvct(coordslt.n_rows), 
      phidiffvct(coordslt.n_rows);
        
    for (int i=0; i<(int)coordslt.n_rows; ++i)
    {
      double qoverptorig = paramslt(i, PCA_ONEOVERPTIDX);
      double qoverptcmps = qoverptcmp[i];
      
      double diffqoverpt = qoverptcmps - qoverptorig;
  
      double phiorig = paramslt(i, PCA_PHIIDX);
      double phicmps = phicmp[i];

      if ((towerid == 19) || (towerid == 20) || 
          (towerid == 27) || (towerid == 28))
      {
        phicmps = phicmps - sec_phi;
        phicmps = fmod(phicmps + M_PI, 2 * M_PI) - M_PI;
      }

      double diffphi = pca::delta_phi(phicmps, phiorig);
  
      pcrelative[PCA_PHIIDX](diffphi/phiorig);
      pcrelative[PCA_ONEOVERPTIDX](diffqoverpt/qoverptorig);
    
      pcabsolute[PCA_PHIIDX](diffphi);
      pcabsolute[PCA_ONEOVERPTIDX](diffqoverpt);
  
      chi2stat(chi2values(i));
  
      qoverptdiffvct(i) = diffqoverpt/qoverptorig;
      phidiffvct(i) = diffphi;
    
      if (writeresults) 
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

    double mmstdev, mpstdev;

    nbins = 2000;
    //double mmstdev = pcrelative[PCA_ONEOVERPTIDX].mean() - 1.5 * pcrelative[PCA_ONEOVERPTIDX].stddev();
    //double mpstdev = pcrelative[PCA_ONEOVERPTIDX].mean() + 1.5 * pcrelative[PCA_ONEOVERPTIDX].stddev();
    mmstdev = -1.0;
    mpstdev = 1.0;
    //TH1D *hist_qoverpt = new TH1D("hist_diff_qoverpt","q/pt diff histogram",nbins, 
    //    qoverptdiffvct.min(), qoverptdiffvct.max());
    TH1D *hist_qoverpt = new TH1D("hist_diff_qoverpt","q/pt diff histogram",nbins, 
          mmstdev, mpstdev);

    nbins = 100000;
    //mmstdev = pcabsolute[PCA_PHIIDX].mean() - 1.5 * pcabsolute[PCA_PHIIDX].stddev();
    //mpstdev = pcabsolute[PCA_PHIIDX].mean() + 1.5 * pcabsolute[PCA_PHIIDX].stddev();
    //TH1D *hist_phi = new TH1D("hist_diff_phi","phi diff histogram",nbins, 
    //    phidiffvct.min(), phidiffvct.max());
    mmstdev = -5.0;
    mpstdev = 5.0;
    TH1D *hist_phi = new TH1D("hist_diff_phi","phi diff histogram",nbins, 
        mmstdev, mpstdev);
  
    for (int i=0; i<(int)coordslt.n_rows; ++i)
    {
      hist_qoverpt->Fill((Double_t) qoverptdiffvct(i));
      hist_phi->Fill((double_t) phidiffvct(i));
    }

    mmstdev = -5.0;
    mpstdev = 5.0;
    //mmstdev = pcabsolute[PCA_PHIIDX].mean() - 1.5 * pcabsolute[PCA_PHIIDX].stddev();
    //mpstdev = pcabsolute[PCA_PHIIDX].mean() + 1.5 * pcabsolute[PCA_PHIIDX].stddev();
    //hist_phi->Fit("gaus","","",phidiffvct.min(),
    //    phidiffvct.max());
    hist_phi->Fit("gaus","","", mmstdev, mpstdev);
 
    mmstdev = -1.0;
    mpstdev = 1.0;
    //mmstdev = pcrelative[PCA_ONEOVERPTIDX].mean() - 1.5 * pcrelative[PCA_ONEOVERPTIDX].stddev();
    //mpstdev = pcrelative[PCA_ONEOVERPTIDX].mean() + 1.5 * pcrelative[PCA_ONEOVERPTIDX].stddev();
    //hist_qoverpt->Fit("gaus","","",qoverptdiffvct.min(),
    //    qoverptdiffvct.max());
    hist_qoverpt->Fit("gaus","","",mmstdev, mpstdev);
  
    TF1 *func_qoverpt = (TF1*)hist_qoverpt->GetFunction("gaus");
    TF1 *func_phi = (TF1*)hist_phi->GetFunction("gaus");
  
    std::cout << 
      "q/pt fitted mean " 
      << ptmin << " " << ptmax << " " << chargesign << " "
      << layeridtorm << " "  
      << func_qoverpt->GetParameter("Mean")*100.0 << " +/- " << 
      func_qoverpt->GetParError(1)*100.0 << std::endl << 
      "q/pt fitted sigma " << layeridtorm << " " 
      << func_qoverpt->GetParameter("Sigma")*100.0 << " +/- " <<
      func_phi->GetParError(2)*100.0 << std::endl;
  
    std::cout << 
      "Phi fitted mean " 
      << ptmin << " " << ptmax << " " << chargesign << " "
      << layeridtorm << " " 
      << func_phi->GetParameter("Mean") << " +/- " << 
      func_phi->GetParError(1) << std::endl << 
      "Phi fitted sigma: " << layeridtorm << " " 
      << func_phi->GetParameter("Sigma") << " +/- " <<
      func_phi->GetParError(2) << std::endl;
  }

  for (int i=0; i<fitter.get_paramdim(); ++i)
  {
    std::cout << "For " << fitter.paramidx_to_string(i) << " error " 
      << etamin << " " << etamax << " " << ptmin << " " << ptmax << " " 
      << chargesign << " " << layeridtorm << " "
      << pcabsolute[i].mean() << " " << pcabsolute[i].stddev() << std::endl;

    if (fitter.paramidx_to_string(i) == "q/pt")
      std::cout << "For " << fitter.paramidx_to_string(i) << " error " 
        << etamin << " " << etamax << " " << ptmin << " " << ptmax << " " 
        << chargesign << " " << layeridtorm << " "
        << 100.0 * pcrelative[i].mean() << " % " << 100.0 * pcrelative[i].stddev() << 
        " % " << std::endl;
  }

  std::cout << " " << std::endl;
  std::cout << "Chivalue mean " << layeridtorm << " " 
    << chi2stat.mean() << " stdev " << 
    chi2stat.stddev() << std::endl;
  
  myfile.close();

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
  while (!fin.eof())
  { 
    double pt, qpt_orig, qpt_fitt, diff, fake;
    fin >> pt >> qpt_orig >> qpt_fitt >> diff >> fake >> fake >> fake >> fake;

    double mmstdev, mpstdev;

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
  
    for (int i=0; i<(int)coordslt.n_rows; ++i)
    {
      hist_qoverpt->Fill((Double_t) diff);
      hist_phi->Fill((double_t) phidiffvct(i));
    }

    mmstdev = -5.0;
    mpstdev = 5.0;
    //mmstdev = pcabsolute[PCA_PHIIDX].mean() - 1.5 * pcabsolute[PCA_PHIIDX].stddev();
    //mpstdev = pcabsolute[PCA_PHIIDX].mean() + 1.5 * pcabsolute[PCA_PHIIDX].stddev();
    //hist_phi->Fit("gaus","","",phidiffvct.min(),
    //    phidiffvct.max());
    hist_phi->Fit("gaus","","", mmstdev, mpstdev);
 
    mmstdev = -1.0;
    mpstdev = 1.0;
    //mmstdev = pcrelative[PCA_ONEOVERPTIDX].mean() - 1.5 * pcrelative[PCA_ONEOVERPTIDX].stddev();
    //mpstdev = pcrelative[PCA_ONEOVERPTIDX].mean() + 1.5 * pcrelative[PCA_ONEOVERPTIDX].stddev();
    //hist_qoverpt->Fit("gaus","","",qoverptdiffvct.min(),
    //    qoverptdiffvct.max());
    hist_qoverpt->Fit("gaus","","",mmstdev, mpstdev);
  
    TF1 *func_qoverpt = (TF1*)hist_qoverpt->GetFunction("gaus");
    TF1 *func_phi = (TF1*)hist_phi->GetFunction("gaus");
  
    std::cout << 
      "q/pt fitted mean " 
      << ptmin << " " << ptmax << " " << chargesign << " "
      << layeridtorm << " "  
      << func_qoverpt->GetParameter("Mean")*100.0 << " +/- " << 
      func_qoverpt->GetParError(1)*100.0 << std::endl << 
      "q/pt fitted sigma " << layeridtorm << " " 
      << func_qoverpt->GetParameter("Sigma")*100.0 << " +/- " <<
      func_phi->GetParError(2)*100.0 << std::endl;
  
 
  }

  fin.close();

  return EXIT_SUCCESS;
}
#endif
