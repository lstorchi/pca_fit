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

bool import_pca_const_int (const std::string & cfname, 
    arma::mat & cmtx, arma::rowvec & qvec, 
    arma::mat & amtx, arma::rowvec & kvec, 
    bool rzplane, bool rphiplane, double etaminin, 
    double etamaxin, double ptminin, double ptmaxin, 
    int chargesignin,
    const std::string & layersid,
    const std::string & pslayersid)
{
  assert(rzplane != rphiplane); 

  std::vector<pca::matrixpcaconst<int32_t> > vct;
  if (read_pcaconst_from_file (vct, cfname.c_str()))
  {
    int hwmanygot = 0;
    std::vector<pca::matrixpcaconst<int32_t> >::const_iterator it = 
      vct.begin();
    for (; it != vct.end(); ++it)
    {
      double ptmin, ptmax, etamin, etamax;
      std::string actuallayids;
      int chargesign;

      it->get_ptrange(ptmin, ptmax);
      it->get_etarange(etamin, etamax);
      chargesign = it->get_chargesign();
      actuallayids = it->get_layersids();

      if (it->get_ttype() == pca::INTEGPT)
      {
        if (rzplane)
        {
          if (it->get_plane_type() == pca::RZ)
          {
            if (actuallayids == pslayersid)
            {
              if ((etaminin >= etamin) && (etaminin <= etamax) &&
                  (etamaxin >= etamin) && (etamaxin <= etamax))
              {
                switch(it->get_const_type())
                {
                  case pca::QVEC :
                    pcamat_to_armarowvec ((*it), qvec);
                    hwmanygot++;
                    break;
                  case pca::KVEC :
                    pcamat_to_armarowvec ((*it), kvec);
                    hwmanygot++;
                    break;
                  case pca::CMTX :
                    pcamat_to_armamat ((*it), cmtx);
                    hwmanygot++;
                    break;
                  case pca::AMTX :
                    pcamat_to_armamat ((*it), amtx);
                    hwmanygot++;
                    break;
                  default:
                    break;
                }
              } 
            }
          }
        }
        else if (rphiplane)
        {
          if (it->get_plane_type() == pca::RPHI)
          {
            if (actuallayids == layersid)
            {
              if (chargesignin == chargesign)
              {
                if ((ptminin >= ptmin) && (ptminin <= ptmax) &&
                    (ptmaxin >= ptmin) && (ptmaxin <= ptmax))
                {
                  switch(it->get_const_type())
                  {
                    case pca::QVEC : 
                      pcamat_to_armarowvec ((*it), qvec);
                      hwmanygot++;
                      break;
                    case pca::KVEC :
                      pcamat_to_armarowvec ((*it), kvec);
                      hwmanygot++;
                      break;
                    case pca::CMTX :
                      pcamat_to_armamat ((*it), cmtx);
                      hwmanygot++;
                      break;
                    case pca::AMTX :
                      pcamat_to_armamat ((*it), amtx);
                      hwmanygot++;
                      break;
                    default:
                      break;
                  }
                }
              } 
            }
          }
        }
      }
    }

    if (hwmanygot == 4)
      return true;
    else
    {
      std::cerr << "Found " << hwmanygot << " const instead of 4" << std::endl;
      return false;
    }
  }

  // TODO add consistency check for dims

  return false;
}


bool import_pca_const (const std::string & cfname, 
    arma::mat & cmtx, arma::rowvec & qvec, 
    arma::mat & amtx, arma::rowvec & kvec, 
    bool rzplane, bool rphiplane, double etaminin, 
    double etamaxin, double ptminin, double ptmaxin, 
    int chargesignin,
    const std::string & layersid,
    const std::string & pslayersid)
{
  assert(rzplane != rphiplane); 

  std::vector<pca::matrixpcaconst<double> > vct;
  if (read_pcaconst_from_file (vct, cfname.c_str()))
  {
    int hwmanygot = 0;
    std::vector<pca::matrixpcaconst<double> >::const_iterator it = 
      vct.begin();
    for (; it != vct.end(); ++it)
    {
      double ptmin, ptmax, etamin, etamax;
      std::string actuallayids;
      int chargesign;

      it->get_ptrange(ptmin, ptmax);
      it->get_etarange(etamin, etamax);
      chargesign = it->get_chargesign();
      actuallayids = it->get_layersids();

      if (it->get_ttype() == pca::FLOATPT)
      {
        if (rzplane)
        {
          if (it->get_plane_type() == pca::RZ)
          {
            if (actuallayids == pslayersid)
            {
              if ((etaminin >= etamin) && (etaminin <= etamax) &&
                  (etamaxin >= etamin) && (etamaxin <= etamax))
              {
                switch(it->get_const_type())
                {
                  case pca::QVEC :
                    pcamat_to_armarowvec ((*it), qvec);
                    hwmanygot++;
                    break;
                  case pca::KVEC :
                    pcamat_to_armarowvec ((*it), kvec);
                    hwmanygot++;
                    break;
                  case pca::CMTX :
                    pcamat_to_armamat ((*it), cmtx);
                    hwmanygot++;
                    break;
                  case pca::AMTX :
                    pcamat_to_armamat ((*it), amtx);
                    hwmanygot++;
                    break;
                  default:
                    break;
                }
              } 
            }
          }
        }
        else if (rphiplane)
        {
          if (it->get_plane_type() == pca::RPHI)
          {
            if (actuallayids == layersid)
            {
              if (chargesignin == chargesign)
              {
                if ((ptminin >= ptmin) && (ptminin <= ptmax) &&
                    (ptmaxin >= ptmin) && (ptmaxin <= ptmax))
                {
                  switch(it->get_const_type())
                  {
                    case pca::QVEC : 
                      pcamat_to_armarowvec ((*it), qvec);
                      hwmanygot++;
                      break;
                    case pca::KVEC :
                      pcamat_to_armarowvec ((*it), kvec);
                      hwmanygot++;
                      break;
                    case pca::CMTX :
                      pcamat_to_armamat ((*it), cmtx);
                      hwmanygot++;
                      break;
                    case pca::AMTX :
                      pcamat_to_armamat ((*it), amtx);
                      hwmanygot++;
                      break;
                    default:
                      break;
                  }
                }
              } 
            }
          }
        }
      }
    }

    if (hwmanygot == 4)
      return true;
    else
    {
      std::cerr << "Found " << hwmanygot << " const instead of 4" << std::endl;
      return false;
    }
  }

  // TODO add consistency check for dims

  return false;
}
 

bool build_and_compare (arma::mat & paramslt, arma::mat & coordslt, 
     arma::mat & cmtx, arma::rowvec & q, arma::mat & amtx, 
     arma::rowvec & k, bool verbose, pca::pcafitter & fitter, 
     bool rzplane, bool rphiplane, arma::vec & ptvals, 
     bool intbitewise, int towerid, double sec_phi, 
     bool writeresults, int layeridtorm, 
     double etamin, double etamax, 
     double ptmin, double ptmax, int chargesign)
{
  int nbins = 100;

  int32_t ** i_ptrs = NULL;
  double ** ptrs = NULL;

  if (intbitewise)
    i_ptrs = new int32_t* [fitter.get_paramdim()];
  else
    ptrs = new double* [fitter.get_paramdim()];

  int32_t * i_cothetacmp = NULL, * i_z0cmp = NULL, * i_qoverptcmp = NULL,
    * i_phicmp = NULL;
  double * cothetacmp = NULL, * z0cmp = NULL, * qoverptcmp = NULL,
    * phicmp = NULL;
  
  if (rzplane)
  {
    if (intbitewise)
    {
      i_cothetacmp = new int32_t [(int)coordslt.n_rows];
      i_z0cmp = new int32_t [(int)coordslt.n_rows];
      i_ptrs[PCA_COTTHETAIDX] = i_cothetacmp;
      i_ptrs[PCA_Z0IDX] = i_z0cmp;
    }
    else
    {
      cothetacmp = new double [(int)coordslt.n_rows];
      z0cmp = new double [(int)coordslt.n_rows];
      ptrs[PCA_COTTHETAIDX] = cothetacmp;
      ptrs[PCA_Z0IDX] = z0cmp;
    }  
  }
  else if (rphiplane)
  {
    if (intbitewise)
    {
      i_qoverptcmp = new int32_t [(int)coordslt.n_rows];
      i_phicmp = new int32_t [(int)coordslt.n_rows];
      i_ptrs[PCA_ONEOVERPTIDX] = i_qoverptcmp;
      i_ptrs[PCA_PHIIDX] = i_phicmp;
    }
    else
    {
      qoverptcmp = new double [(int)coordslt.n_rows];
      phicmp = new double [(int)coordslt.n_rows];
      ptrs[PCA_ONEOVERPTIDX] = qoverptcmp;
      ptrs[PCA_PHIIDX] = phicmp;
    }
  }

  arma::rowvec chi2values, chi2values_fake;
  chi2values.resize(coordslt.n_rows);
  chi2values_fake.resize(0);

  if (intbitewise)
  {
    if (!fitter.compute_parameters (cmtx, q, amtx, k, 
          coordslt, i_ptrs, fitter.get_paramdim(), 
          chi2values))
    {
      std::cerr << fitter.get_errmsg() << std::endl;
      return false;
    }
    delete [] i_ptrs;
  }
  else
  {
    if (!fitter.compute_parameters (cmtx, q, amtx, k, 
          coordslt, ptrs, fitter.get_paramdim(), 
          chi2values, chi2values_fake))
    {
      std::cerr << fitter.get_errmsg() << std::endl;
      return false;
    }
    delete [] ptrs; 
  }

  std::ostringstream fname;
  fname << "results.txt";

  arma::running_stat<double> pcrelative[fitter.get_paramdim()];
  arma::running_stat<double> pcabsolute[fitter.get_paramdim()];
  arma::running_stat<double> chi2stat;
  std::ofstream myfile(fname.str().c_str());

  assert(chi2values.n_cols == coordslt.n_rows);

  if (intbitewise)
  {
    if (rzplane)
    {
      if (writeresults)
        myfile << "pt eta_orig eta_fitt diff z0_orig z0_fitt diff chi2 chi2" << std::endl;
    
      arma::rowvec etadiffvct(coordslt.n_rows), 
        z0diffvct(coordslt.n_rows);
     
      for (int i=0; i<(int)coordslt.n_rows; ++i)
      {
        double thetacmp = atan(1.0e0 / i_cothetacmp[i]) ; 
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
     
        double etadiff = (etacmp - etaorig);
        etadiffvct(i) = etadiff;
     
        int32_t z0cmps = i_z0cmp[i];
        int32_t z0orig = paramslt(i, PCA_Z0IDX);
        int32_t z0diff = z0cmps - paramslt(i, PCA_Z0IDX);
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
      }
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
        int32_t qoverptorig = paramslt(i, PCA_ONEOVERPTIDX);
        int32_t qoverptcmps = i_qoverptcmp[i];
        
        int32_t diffqoverpt = qoverptcmps - qoverptorig;
    
        double phiorig = paramslt(i, PCA_PHIIDX);
        int32_t phicmps = i_phicmp[i];
        int32_t diffphi = phicmps - phiorig;
    
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
      }
    }
  }
  else
  {
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

  if (intbitewise)
  {
    if (rzplane)
    {
      delete [] i_cothetacmp;
      delete [] i_z0cmp;
    }
    else if (rphiplane)
    {
      delete [] i_qoverptcmp;
      delete [] i_phicmp;
    }
  }
  else
  {
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
  std::cerr << " -T, --int-bitewise               : Integer bitewise mode on" << std::endl;
  std::cerr << " -p, --dump-allcoords             : dump all stub coordinates to a file" << std::endl;
  std::cerr << " -c, --pca-const-file=[fillename] : PCA const txt filename [default is pca_const.txt]" << std::endl;
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
  std::cerr << "                                    use linear interpolation to approximate the missed hit " << std::endl;

  std::cerr << " -w, --fk-five-hits=[layerid]     : build constants for 5 / 6, specify the layr to be removed " 
    << std::endl;
  std::cerr << "                                   it will use 6 layers tracks, removing a layer " << std::endl;
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
  std::cerr << " -D, --towerid=[num]              : MANDATORY: specify towid to be used for the XY rotation " 
    << std::endl;
  std::cerr << "                                    written in the file " << std::endl;
  std::cerr << " -N, --no-results                 : results file is not written, only mean and stdev " << std::endl;
  std::cerr << "                                    are computed and reported " << std::endl; 
  std::cerr << " -X, --max-num-oftracks=[n]      : stop reading root file after n tracks" << std::endl;
  std::cerr << " -C, --use-only-3-layers         : use three leyers ..." << std::endl;
  std::cerr << std::endl;

  exit(1);
}


# ifndef __CINT__
int main (int argc, char ** argv)
{
  gROOT->ProcessLine("#include <vector>");

  pca::pcafitter fitter;

  std::string cfname = "pca_const.txt";
  bool verbose = false;
  bool rzplane = false;
  bool rphiplane = false;
  bool checklayersids = false;
  bool savecheckfiles = false;
  bool userelativecoord = false;
  bool lininterpolation = false;
  bool printallcoords = false;
  bool usefakefiveoutofsix = false;

  bool intbitewise = false;

  int layeridtorm = -1;

  unsigned int maxnumoftracks = (unsigned int) INFINITY;

  double etamin = -1.0e0 * INFINITY, etamax = +1.0e0 * INFINITY;
  double ptmin = -1.0e0 * INFINITY, ptmax = +1.0e0 * INFINITY;
  double phimin = -1.0e0 * INFINITY, phimax = +1.0e0 * INFINITY;
  double z0min = -1.0e0 * INFINITY, z0max = +1.0e0 * INFINITY;
  double d0min = -1.0e0 * INFINITY, d0max = +1.0e0 * INFINITY;
  double coord1min = std::numeric_limits<double>::infinity();
  double coord2min = std::numeric_limits<double>::infinity();

  std::string sequence = "";
  int numoflayers = 6;

  std::string layersid;
  std::string pslayersid;

  int chargesign = 0;

  std::vector<std::string> tokens;

  bool excludesmodule = false;

  int towerid = -99;

  bool writeresults = true;
  bool use3layers = false;

  while (1)
  {
    int c, option_index;
    static struct option long_options[] = {
      {"help", 0, NULL, 'h'},
      {"pca-const-file", 1, NULL, 'c'},
      {"verbose", 0, NULL, 'V'},
      {"version", 0, NULL, 'v'},
      {"int-bitewise", 0, NULL, 'T'}, 
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
      {"use-only-3-layers", 0, NULL, 'C'},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "CNTpkxzrhaVD:w:X:l:f:b:t:g:c:n:s:m:o:u:", 
        long_options, &option_index);

    if (c == -1)
      break;

    switch (c)
    {
      case 'N':
        writeresults = false;
        break;
      case 'D':
        towerid = atoi(optarg);
        break;
      case 'T':
        intbitewise = true;
        break;
      case 'w':
        usefakefiveoutofsix = true;
        layeridtorm = atoi(optarg);
        break;
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
      case 'X':
        maxnumoftracks = atoi(optarg);
        break;
      case 'C':
        use3layers = true;
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

  std::set<int> layers;
  if (use3layers)
  {
    layers.insert(5);
    layers.insert(8);
    layers.insert(10);
  }

  fitter.set_useintbitewise(intbitewise);

  // very quick 
  if (sequence == "")
  {
    std::ostringstream psosss, osss;
    std::cout << "Only for BARREL" << std::endl;
    for (int i =5; i<=10; ++i)
    {
      if (usefakefiveoutofsix)
        if (i == layeridtorm)
          continue;

      if (use3layers)
        if (layers.find(i) == layers.end())
          continue;
 
      osss << i << ":";
      if (i <= 7)
        psosss << i << ":";
    }

    layersid = osss.str();
    layersid.erase(layersid.end()-1);

    pslayersid = psosss.str();
    pslayersid.erase(pslayersid.end()-1);
  }
  else
  {
    std::cerr << "Import of constants to be implemented" << std::endl;
    return EXIT_FAILURE;
  }

  if (numoflayers == 5)
  {
    if (usefakefiveoutofsix)
    {
      std::cerr << "Wrong options, cannot use both options together" << std::endl;
      return EXIT_FAILURE;
    }

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
    // I am using cot(theta) internally
    if (!fitter.set_paramidx(PCA_COTTHETAIDX, "eta"))
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
  
  arma::mat cmtx, amtx;
  arma::rowvec qvec, kvec;

  // leggere file coordinate tracce simulate plus parametri
  if (!pca::file_exists(filename))
  {
    std::cerr << "Inout file does not exist" << std::endl;
    return EXIT_FAILURE;
  }

  if (pca::file_exists(cfname.c_str()))
  {
    std::cout << "Reading " << cfname << std::endl;

    if (intbitewise)
    {
      if (!import_pca_const_int (cfname, cmtx, qvec, amtx, kvec, 
            rzplane, rphiplane, etamin, etamax, ptmin, 
            ptmax, chargesign, layersid, pslayersid))
      {
        std::cerr << "Error in reading constants from file" << std::endl;
        return EXIT_FAILURE;
      }
    }
    else
    {
      if (!import_pca_const (cfname, cmtx, qvec, amtx, kvec, 
            rzplane, rphiplane, etamin, etamax, ptmin, 
            ptmax, chargesign, layersid, pslayersid))
      {
        std::cerr << "Error in reading constants from file" << std::endl;
        return EXIT_FAILURE;
      }
    }

    unsigned int coorddim = (unsigned int) fitter.get_coordim();
    unsigned int paramdim = (unsigned int) fitter.get_paramdim();

    std::cout << "using " << coorddim << " coordinates and " << paramdim <<
      " parameters " << std::endl;

    if (coorddim != cmtx.n_cols)
    {
      std::cerr << "Incompatible dimensions CMTX (2S-modules ?)" << std::endl;
      return EXIT_FAILURE;
    }
    if (paramdim != cmtx.n_rows)
    {
      std::cerr << "Incompatible dimensions CMTX (2S-modules ?)" << std::endl;
      return EXIT_FAILURE;
    }
    if (paramdim != qvec.n_elem)
    {
      std::cerr << "Incompatible dimensions QVEC (2S-modules ?)" << std::endl;
      return EXIT_FAILURE;
    }
  }
  else
  {
    std::cerr << cfname << " does not exist" << std::endl;
    return 1;
  }

  std::cout << "Reading data from " << filename << " file " << std::endl;

  arma::mat coord, param;
  arma::vec ptvals;

  pca::rootfilereader rootrdr;

  rootrdr.set_useintbitewise(intbitewise);
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
  //maxnumoftracks = 100000;
  rootrdr.set_maxnumoftracks(maxnumoftracks);
  if (use3layers)
  {
    rootrdr.set_use3layers(layers);
  }

  rootrdr.set_towid(towerid);

  if ((towerid == 19) || (towerid == 20) || 
      (towerid == 27) || (towerid == 28))
  {
    rootrdr.apply_rotation_to_xy(true);
  }

  rootrdr.set_fkfiveoutofsix(usefakefiveoutofsix, 
      layeridtorm);

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

  std::cout << "Rotation angle used: " << rootrdr.get_rotation_angle() << std::endl;

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

  if (!build_and_compare (param, coord, cmtx, qvec, amtx, kvec, 
        verbose, fitter, rzplane, rphiplane, ptvals, intbitewise,
        towerid, rootrdr.get_rotation_angle(), writeresults, 
        layeridtorm, etamin, etamax, ptmin, ptmax, chargesign))
    return EXIT_FAILURE;

  std::cout << "Constants Used: C matrix: " << std::endl;
  std::cout << cmtx;
  std::cout << "Constants Used: q matrix: " << std::endl;
  std::cout << qvec;

  if (intbitewise)
  {
    std::cout << "Constants Used with precision:" << std::endl;
    std::cout << "Constants Used: C matrix: " << std::endl;
    for (int i=0; i<2; ++i)
      for (int j=0; j<12; ++j)
        std::cout << std::setprecision(9) << (double) cmtx(i, j) << std::endl;
    
    std::cout << "Constants Used: q matrix: " << std::endl;
    for (int i=0; i<2; ++i)
      std::cout << std::setprecision(9) << (double) qvec(i) << std::endl;
  }
  
  return EXIT_SUCCESS;
}
#endif
