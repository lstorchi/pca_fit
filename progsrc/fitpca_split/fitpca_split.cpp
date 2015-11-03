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

// lstorchi : can be included in any case 
#ifdef INTBITEWISE
// lstorchi: should we use cstdint and -std=c++11 ? 
#include "stdint.h"
#endif

// lstorchi: basi code to fit tracks, using the PCA constants generated 
//           by the related generatepca


bool build_and_compare (arma::mat & paramslt, arma::mat & coordslt, 
     arma::mat & cmtx, arma::rowvec & q, arma::mat & amtx, arma::mat & vmtx, 
     bool verbose, pca::pcafitter & fitter, bool rzplane, bool rphiplane,
     bool usecharge, bool usealsod0, bool usex0y0, int singleparam,
     bool usealsox0, arma::vec & ptvals)
{

#ifdef INTBITEWISE
  int16_t ** ptrs;
  ptrs = new int16_t* [fitter.get_paramdim()];

  int16_t * cothetacmp = NULL, * z0cmp = NULL, * oneoverptcmp = NULL, 
    * phicmp = NULL, * d0cmp = NULL, * x0cmp = NULL, * y0cmp = NULL, 
    * singlep = NULL;
#else
  double ** ptrs;
  ptrs = new double* [fitter.get_paramdim()];

  double * cothetacmp = NULL, * z0cmp = NULL, * oneoverptcmp = NULL, 
    * phicmp = NULL, * d0cmp = NULL, * x0cmp = NULL, * y0cmp = NULL, 
    * singlep = NULL;
#endif
 
  if (usealsox0 && usex0y0)
    return false;

  if ((singleparam >= 1) && (singleparam <= 7))
  {
#ifdef INTBITEWISE    
    singlep = new int16_t [(int)coordslt.n_rows];
#else
    singlep = new double [(int)coordslt.n_rows];
#endif    
    ptrs[0] = singlep;
  }
  else
  {
    if (rzplane)
    {
      if (usex0y0)
      {
#ifdef INTBITEWISE        
        x0cmp = new int16_t [(int)coordslt.n_rows];
        y0cmp = new int16_t [(int)coordslt.n_rows];
#else
        x0cmp = new double [(int)coordslt.n_rows];
        y0cmp = new double [(int)coordslt.n_rows];
#endif
        ptrs[SPLIT_X0IDX] = x0cmp;
        ptrs[SPLIT_Y0IDX] = y0cmp;
      }
      else
      {
#ifdef INTBITEWISE        
        cothetacmp = new int16_t [(int)coordslt.n_rows];
        z0cmp = new int16_t [(int)coordslt.n_rows];
#else
        cothetacmp = new double [(int)coordslt.n_rows];
        z0cmp = new double [(int)coordslt.n_rows];
#endif        
        ptrs[SPLIT_COTTHETAIDX] = cothetacmp;
        ptrs[SPLIT_Z0IDX] = z0cmp;
      }
      
    }
    else if (rphiplane)
    {
      if (usex0y0)
      {
#ifdef INTBITEWISE        
        x0cmp = new int16_t [(int)coordslt.n_rows];
        y0cmp = new int16_t [(int)coordslt.n_rows];
#else
        cothetacmp = new double [(int)coordslt.n_rows];
        z0cmp = new double [(int)coordslt.n_rows];
#endif
        ptrs[SPLIT_X0IDX] = x0cmp;
        ptrs[SPLIT_Y0IDX] = y0cmp;
      }
      else 
      {
#ifdef INTBITEWISE        
        oneoverptcmp = new int16_t [(int)coordslt.n_rows];
        phicmp = new int16_t [(int)coordslt.n_rows];
#else
        oneoverptcmp = new double [(int)coordslt.n_rows];
        phicmp = new double [(int)coordslt.n_rows];
#endif
        ptrs[SPLIT_ONEOVERPTIDX] = oneoverptcmp;
        ptrs[SPLIT_PHIIDX] = phicmp;
      }
    
    }

#ifdef INTBITEWISE   
    if (usealsod0)
    {
      d0cmp = new int16_t [(int)coordslt.n_rows];
      ptrs[SPLIT_D0IDX] = d0cmp;
    }
    else if (usealsox0)
    {
      x0cmp = new int16_t [(int)coordslt.n_rows];
      ptrs[SPLIT_X0IDX_NS] = x0cmp;
    }
#else
    if (usealsod0)
    {
      d0cmp = new double [(int)coordslt.n_rows];
      ptrs[SPLIT_D0IDX] = d0cmp;
    }
    else if (usealsox0)
    {
      x0cmp = new double [(int)coordslt.n_rows];
      ptrs[SPLIT_X0IDX_NS] = x0cmp;
    }
#endif    

  }

  arma::rowvec chi2values;
  chi2values.resize(coordslt.n_rows);

  if (!fitter.compute_parameters (cmtx, q, amtx, vmtx, 
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

  if ((singleparam >= 1) && (singleparam <= 7))
  {
    myfile << "pt orig fitt diff " << std::endl;
    for (int i=0; i<(int)coordslt.n_rows; ++i)
    {
      pcrelative[0]((singlep[i] - paramslt(i, 0))/
          paramslt(i, 0));

      pcabsolute[0](singlep[i] - paramslt(i, 0));
    
      myfile << ptvals(i) << " " <<
        paramslt(i, 0) << " " << singlep[i] << " " <<
        singlep[i] - paramslt(i, 0) << std::endl;
    
      if (verbose)
      {
        std::cout << "For track : " << i+1 << std::endl;
        std::cout << " param        fitt " << singlep[i] << std::endl;
        std::cout << " param        orig " << paramslt(i, 0) << std::endl;
      }
    }
  }
  else
  {
    if (rzplane)
    {
      if (usex0y0)
      {
        if (usealsod0)
          myfile << "pt x0_orig x0_fitt diff y0_orig y0_fitt diff " <<
           " d0_orig d0_fitt diff" << std::endl;
        else  
          myfile << "pt x0_orig x0_fitt diff y0_orig y0_fitt diff" << std::endl; 
    
        for (int i=0; i<(int)coordslt.n_rows; ++i)
        {
          pcrelative[SPLIT_X0IDX]((x0cmp[i] - paramslt(i, SPLIT_X0IDX))/
              paramslt(i, SPLIT_X0IDX));
          pcrelative[SPLIT_Y0IDX]((y0cmp[i] - paramslt(i, SPLIT_Y0IDX))/
              paramslt(i, SPLIT_Y0IDX));
        
          pcabsolute[SPLIT_X0IDX](x0cmp[i] - paramslt(i, SPLIT_X0IDX));
          pcabsolute[SPLIT_Y0IDX](y0cmp[i] - paramslt(i, SPLIT_Y0IDX));
        
          if (usealsod0) 
          {
            pcrelative[SPLIT_D0IDX]((d0cmp[i] - paramslt(i, SPLIT_D0IDX))/
                paramslt(i, SPLIT_D0IDX));
            pcabsolute[SPLIT_D0IDX](d0cmp[i] - paramslt(i, SPLIT_D0IDX));
          }
          
          if (usealsod0)
            myfile << ptvals(i) << " " <<
              paramslt(i, SPLIT_X0IDX) << " " << x0cmp[i] << " " <<
              (x0cmp[i] - paramslt(i, SPLIT_X0IDX)) << " " <<
              paramslt(i, SPLIT_Y0IDX) << " " << y0cmp[i] << " " <<
              (y0cmp[i] - paramslt(i, SPLIT_Y0IDX)) << " " << 
              paramslt(i, SPLIT_D0IDX) << " " << d0cmp[i] << " " <<
              d0cmp[i] - paramslt(i, SPLIT_D0IDX) << std::endl;
          else
            myfile << ptvals(i) << " " <<
              paramslt(i, SPLIT_X0IDX) << " " << x0cmp[i] << " " <<
              (x0cmp[i] - paramslt(i, SPLIT_X0IDX)) << " " <<
              paramslt(i, SPLIT_Y0IDX) << " " << y0cmp[i] << " " <<
              (y0cmp[i] - paramslt(i, SPLIT_Y0IDX)) << std::endl;
        
          if (verbose)
          {
            std::cout << "For track : " << i+1 << std::endl;
            std::cout << " x0           fitt " << x0cmp[i] << std::endl;
            std::cout << " x0           orig " << paramslt(i, SPLIT_X0IDX) << std::endl;
            std::cout << " y0           fitt " << y0cmp[i] << std::endl;
            std::cout << " y0           orig " << paramslt(i, SPLIT_Y0IDX) << std::endl;
    
            if (usealsod0)
            { 
              std::cout << " d0           fitt " << d0cmp[i] << std::endl;
              std::cout << " d0           orig " << paramslt(i, SPLIT_D0IDX) << std::endl;
            }
          }
        }
    
      }
      else
      {
        if (usealsod0)
          myfile << "pt eta_orig eta_fitt diff z0_orig z0_fitt diff " <<
           " d0_orig d0_fitt diff" << std::endl;
        else if (usealsox0)
          myfile << "pt eta_orig eta_fitt diff z0_orig z0_fitt diff " <<
           " x0_orig x0_fitt diff" << std::endl;
        else  
#ifdef INTBITEWISE          
          myfile << "pt cot0_orig cot0_fitt diff z0_orig z0_fitt diff" << std::endl; 
#else
          myfile << "pt eta_orig eta_fitt diff z0_orig z0_fitt diff" << std::endl; 
#endif
    
        for (int i=0; i<(int)coordslt.n_rows; ++i)
        {
#ifdef INTBITEWISE          
          //It will not converge. Lets not try this now. Instead compare Cot(theta) itself (easier).
          //int16_t thetacmp = atan(1.0e0 / cothetacmp[i]) ; 
          //int16_t etacmps = 0.0e0, tantheta2;
#else
          double thetacmp = atan(1.0e0 / cothetacmp[i]) ; 
          double etacmps = 0.0e0, tantheta2;
#endif
          tantheta2 = tan (thetacmp/2.0e0); 
          if (tantheta2 < 0.0)
            etacmps = 1.0e0 * log (-1.0e0 * tantheta2);
          else
            etacmps = -1.0e0 * log (tantheta2);
#ifdef INTBITEWISE
          //int16_t theta = atan(1.0e0 / paramslt(i, SPLIT_COTTHETAIDX));
#else
          double theta = atan(1.0e0 / paramslt(i, SPLIT_COTTHETAIDX));
#endif
          double etaorig = 0.0e0;
          tantheta2 = tan (theta/2.0e0);
          if (tantheta2 < 0.0)
            etaorig = 1.0e0 * log (-1.0e0 * tantheta2);
          else
            etaorig = -1.0e0 * log (tantheta2);
        
          pcrelative[SPLIT_COTTHETAIDX]((etacmps - etaorig)/
              etaorig);
          pcrelative[SPLIT_Z0IDX]((z0cmp[i] - paramslt(i, SPLIT_Z0IDX))/
              paramslt(i, SPLIT_Z0IDX));
        
          pcabsolute[SPLIT_COTTHETAIDX](etacmps - etaorig);
          pcabsolute[SPLIT_Z0IDX](z0cmp[i] - paramslt(i, SPLIT_Z0IDX));
        
          if (usealsod0) 
          {
            pcrelative[SPLIT_D0IDX]((d0cmp[i] - paramslt(i, SPLIT_D0IDX))/
                paramslt(i, SPLIT_D0IDX));
            pcabsolute[SPLIT_D0IDX](d0cmp[i] - paramslt(i, SPLIT_D0IDX));
          }
          else if (usealsox0)
          {
            pcrelative[SPLIT_X0IDX_NS]((x0cmp[i] - paramslt(i, SPLIT_X0IDX_NS))/
                paramslt(i, SPLIT_X0IDX_NS));
            pcabsolute[SPLIT_X0IDX_NS](x0cmp[i] - paramslt(i, SPLIT_X0IDX_NS));
          }
          
          if (usealsod0)
            myfile << ptvals(i) << " " <<
              etaorig << " " << etacmps << " " <<
              (etacmps - etaorig) << " " <<
              paramslt(i, SPLIT_Z0IDX) << " " << z0cmp[i] << " " <<
              (z0cmp[i] - paramslt(i, SPLIT_Z0IDX)) << " " << 
              paramslt(i, SPLIT_D0IDX) << " " << d0cmp[i] << " " <<
              d0cmp[i] - paramslt(i, SPLIT_D0IDX) << std::endl;
          else if (usealsox0)
            myfile << ptvals(i) << " " <<
              etaorig << " " << etacmps << " " <<
              (etacmps - etaorig) << " " <<
              paramslt(i, SPLIT_Z0IDX) << " " << z0cmp[i] << " " <<
              (z0cmp[i] - paramslt(i, SPLIT_Z0IDX)) << " " << 
              paramslt(i, SPLIT_X0IDX_NS) << " " << x0cmp[i] << " " <<
              x0cmp[i] - paramslt(i, SPLIT_X0IDX_NS) << std::endl;
          else
#ifdef INTBITEWISE
            myfile << ptvals(i) << "   " <<
              paramslt(i, SPLIT_COTTHETAIDX) << "   " << cothetacmp[i] << "   " <<
              (cothetacmp[i] - paramslt(i, SPLIT_COTTHETAIDX)) << " " <<
              paramslt(i, SPLIT_Z0IDX) << " " << z0cmp[i] << " " <<
              (z0cmp[i] - paramslt(i, SPLIT_Z0IDX)) << std::endl;
#else
            myfile << ptvals(i) << " " <<
              etaorig << " " << etacmps << " " <<
              (etacmps - etaorig) << " " <<
              paramslt(i, SPLIT_Z0IDX) << " " << z0cmp[i] << " " <<
              (z0cmp[i] - paramslt(i, SPLIT_Z0IDX)) << std::endl;
#endif
        
          if (verbose)
          {
            std::cout << "For track : " << i+1 << std::endl;
            std::cout << " cotheta      fitt " << cothetacmp[i] << std::endl;
            std::cout << " cotheta      orig " << paramslt(i, SPLIT_COTTHETAIDX) << std::endl;
            std::cout << " theta rad    fitt " << thetacmp << std::endl;
            std::cout << " theta rad    orig " << theta << std::endl;
            std::cout << " theta deg    fitt " << thetacmp*(180.0e0/M_PI) << std::endl;
            std::cout << " theta deg    orig " << theta*(180.0e0/M_PI) << std::endl;
            std::cout << " eta          fitt " << etacmps << std::endl;
            std::cout << " eta          orig " << etaorig << std::endl;
            std::cout << " z0           fitt " << z0cmp[i] << std::endl;
            std::cout << " z0           orig " << paramslt(i, SPLIT_Z0IDX) << std::endl;
            if (usealsod0)
            { 
              std::cout << " d0           fitt " << d0cmp[i] << std::endl;
              std::cout << " d0           orig " << paramslt(i, SPLIT_D0IDX) << std::endl;
            }
            else if (usealsox0)
            {
              std::cout << " x0           fitt " << d0cmp[i] << std::endl;
              std::cout << " x0           orig " << paramslt(i, SPLIT_D0IDX) << std::endl;
            }
          }
        }
      }
    }
    else if (rphiplane)
    {
      std::ofstream myfile(fname.str().c_str());
      if (usex0y0)
      {
        if (usealsod0)
          myfile << "pt x0_orig x0_fitt diff y0_orig y0_fitt diff " <<
           " d0_orig d0_fitt diff" << std::endl;
        else  
          myfile << "pt x0_orig x0_fitt diff y0_orig y0_fitt diff" << std::endl; 
    
        for (int i=0; i<(int)coordslt.n_rows; ++i)
        {
          pcrelative[SPLIT_X0IDX]((x0cmp[i] - paramslt(i, SPLIT_X0IDX))/
              paramslt(i, SPLIT_X0IDX));
          pcrelative[SPLIT_Y0IDX]((y0cmp[i] - paramslt(i, SPLIT_Y0IDX))/
              paramslt(i, SPLIT_Y0IDX));
        
          pcabsolute[SPLIT_X0IDX](x0cmp[i] - paramslt(i, SPLIT_X0IDX));
          pcabsolute[SPLIT_Y0IDX](y0cmp[i] - paramslt(i, SPLIT_Y0IDX));
        
          if (usealsod0) 
          {
            pcrelative[SPLIT_D0IDX]((d0cmp[i] - paramslt(i, SPLIT_D0IDX))/
                paramslt(i, SPLIT_D0IDX));
            pcabsolute[SPLIT_D0IDX](d0cmp[i] - paramslt(i, SPLIT_D0IDX));
          }
          
          if (usealsod0)
            myfile << ptvals(i) << " " <<
              paramslt(i, SPLIT_X0IDX) << " " << x0cmp[i] << " " <<
              (x0cmp[i] - paramslt(i, SPLIT_X0IDX)) << " " <<
              paramslt(i, SPLIT_Y0IDX) << " " << y0cmp[i] << " " <<
              (y0cmp[i] - paramslt(i, SPLIT_Y0IDX)) << " " << 
              paramslt(i, SPLIT_D0IDX) << " " << d0cmp[i] << " " <<
              d0cmp[i] - paramslt(i, SPLIT_D0IDX) << std::endl;
          else
            myfile << ptvals(i) << " " <<
              paramslt(i, SPLIT_X0IDX) << " " << x0cmp[i] << " " <<
              (x0cmp[i] - paramslt(i, SPLIT_X0IDX)) << " " <<
              paramslt(i, SPLIT_Y0IDX) << " " << y0cmp[i] << " " <<
              (y0cmp[i] - paramslt(i, SPLIT_Y0IDX)) << std::endl;
        
          if (verbose)
          {
            std::cout << "For track : " << i+1 << std::endl;
            std::cout << " x0           fitt " << x0cmp[i] << std::endl;
            std::cout << " x0           orig " << paramslt(i, SPLIT_X0IDX) << std::endl;
            std::cout << " y0           fitt " << y0cmp[i] << std::endl;
            std::cout << " y0           orig " << paramslt(i, SPLIT_Y0IDX) << std::endl;
    
            if (usealsod0)
            { 
              std::cout << " d0           fitt " << d0cmp[i] << std::endl;
              std::cout << " d0           orig " << paramslt(i, SPLIT_D0IDX) << std::endl;
            }
          }
        } 
      }
      else
      {
        if (usecharge)
        {
          if (usealsod0)
            myfile << "pt q/pt_orig q/pt_fitt diff phi_orig phi_fitt diff " <<
             " d0_orig d0_fitt diff" << std::endl;
          else if (usealsox0)
            myfile << "pt q/pt_orig q/pt_fitt diff phi_orig phi_fitt diff " <<
             " x0_orig x0_fitt diff" << std::endl;
          else
            myfile << "pt q/pt_orig q/pt_fitt diff phi_orig phi_fitt diff" << std::endl; 
        
          for (int i=0; i<(int)coordslt.n_rows; ++i)
          {
            double qoverptorig = paramslt(i, SPLIT_ONEOVERPTIDX);
#ifdef INTBITEWISE            
            int16_t qoverptcmp = oneoverptcmp[i];
#else
            double qoverptcmp = oneoverptcmp[i];
#endif          
            pcrelative[SPLIT_PHIIDX]((phicmp[i] - paramslt(i, SPLIT_PHIIDX))/
                paramslt(i, SPLIT_PHIIDX));
            pcrelative[SPLIT_ONEOVERPTIDX]((qoverptcmp - qoverptorig)/
                qoverptorig);
          
            pcabsolute[SPLIT_PHIIDX](phicmp[i] - paramslt(i, SPLIT_PHIIDX));
            pcabsolute[SPLIT_ONEOVERPTIDX](qoverptcmp - qoverptorig);
        
            if (usealsod0) 
            {
              pcrelative[SPLIT_D0IDX]((d0cmp[i] - paramslt(i, SPLIT_D0IDX))/
                  paramslt(i, SPLIT_D0IDX));
              pcabsolute[SPLIT_D0IDX](d0cmp[i] - paramslt(i, SPLIT_D0IDX));
            }
            else if (usealsox0)
            {
              pcrelative[SPLIT_X0IDX_NS]((x0cmp[i] - paramslt(i, SPLIT_X0IDX_NS))/
                  paramslt(i, SPLIT_X0IDX_NS));
              pcabsolute[SPLIT_X0IDX_NS](x0cmp[i] - paramslt(i, SPLIT_X0IDX_NS));
            }
          
            if (usealsod0)
              myfile << ptvals(i) << " " <<
                qoverptorig << " " << qoverptcmp << " " <<
                (qoverptcmp - qoverptorig) <<  " " <<
                paramslt(i, SPLIT_PHIIDX) << " " << phicmp[i] << " " <<
                (phicmp[i] - paramslt(i, SPLIT_PHIIDX)) << " " << 
                paramslt(i, SPLIT_D0IDX) << " " << d0cmp[i] << " " <<
                d0cmp[i] - paramslt(i, SPLIT_D0IDX) << std::endl;
            else if (usealsox0)
              myfile << ptvals(i) << " " <<
                qoverptorig << " " << qoverptcmp << " " <<
                (qoverptcmp - qoverptorig) <<  " " <<
                paramslt(i, SPLIT_PHIIDX) << " " << phicmp[i] << " " <<
                (phicmp[i] - paramslt(i, SPLIT_PHIIDX)) << " " << 
                paramslt(i, SPLIT_X0IDX_NS) << " " << x0cmp[i] << " " <<
                x0cmp[i] - paramslt(i, SPLIT_X0IDX_NS) << std::endl;
            else
              myfile << ptvals(i) << " " <<
                qoverptorig << " " << qoverptcmp << " " <<
                (qoverptcmp - qoverptorig) <<  " " <<
                paramslt(i, SPLIT_PHIIDX) << " " << phicmp[i] << " " <<
                (phicmp[i] - paramslt(i, SPLIT_PHIIDX)) << std::endl;
          
            if (verbose)
            {
              std::cout << "For track : " << i+1 << std::endl;
              std::cout << " q/pt         fitt " << qoverptcmp << std::endl;
              std::cout << " q/pt         orig " << qoverptorig << std::endl;
              std::cout << " phi          fitt " << phicmp[i] << std::endl;
              std::cout << " phi          orig " << paramslt(i, SPLIT_PHIIDX) << std::endl;
              if (usealsod0)
              { 
                std::cout << " d0           fitt " << d0cmp[i] << std::endl;
                std::cout << " d0           orig " << paramslt(i, SPLIT_D0IDX) << std::endl;
              }
              else if (usealsox0)
              {
                std::cout << " x0           fitt " << x0cmp[i] << std::endl;
                std::cout << " x0           orig " << paramslt(i, SPLIT_X0IDX_NS) << std::endl;
              }
            }
          }
        }
        else 
        {
          if (usealsod0)
            myfile << "pt pt_orig pt_fitt diff phi_orig phi_fitt diff " <<
             " d0_orig d0_fitt diff" << std::endl;
          else if (usealsox0)
            myfile << "pt pt_orig pt_fitt diff phi_orig phi_fitt diff " <<
             " x0_orig x0_fitt diff" << std::endl;
          else
            myfile << "pt pt_orig pt_fitt diff phi_orig phi_fitt diff" << std::endl; 
        
          for (int i=0; i<(int)coordslt.n_rows; ++i)
          {
            // not sure can be double this one ...
            double ptorig = 1.0e0 / paramslt(i, SPLIT_ONEOVERPTIDX);
#ifdef INTBITEWISE            
            int16_t ptcmp = 1.0e0 / oneoverptcmp[i];
#else 
            double ptcmp = 1.0e0 / oneoverptcmp[i];
#endif
            pcrelative[SPLIT_PHIIDX]((phicmp[i] - paramslt(i, SPLIT_PHIIDX))/
                paramslt(i, SPLIT_PHIIDX));
            pcrelative[SPLIT_ONEOVERPTIDX]((ptcmp - ptorig)/
                ptorig);
          
            pcabsolute[SPLIT_PHIIDX](phicmp[i] - paramslt(i, SPLIT_PHIIDX));
            pcabsolute[SPLIT_ONEOVERPTIDX](ptcmp - ptorig);
        
            if (usealsod0) 
            {
              pcrelative[SPLIT_D0IDX]((d0cmp[i] - paramslt(i, SPLIT_D0IDX))/
                  paramslt(i, SPLIT_D0IDX));
              pcabsolute[SPLIT_D0IDX](d0cmp[i] - paramslt(i, SPLIT_D0IDX));
            }
            else if (usealsox0)
            {
              pcrelative[SPLIT_X0IDX_NS]((x0cmp[i] - paramslt(i, SPLIT_X0IDX_NS))/
                  paramslt(i, SPLIT_X0IDX_NS));
              pcabsolute[SPLIT_X0IDX_NS](x0cmp[i] - paramslt(i, SPLIT_X0IDX_NS));
            }
 
        
            if (usealsod0)
              myfile << ptvals(i) << " " <<
                ptorig << " " << ptcmp << " " <<
                (ptcmp - ptorig) <<  " " <<
                paramslt(i, SPLIT_PHIIDX) << " " << phicmp[i] << " " <<
                (phicmp[i] - paramslt(i, SPLIT_PHIIDX)) << " " << 
                paramslt(i, SPLIT_D0IDX) << " " << d0cmp[i] << " " <<  
                d0cmp[i] - paramslt(i, SPLIT_D0IDX) << std::endl;
            else if (usealsox0)
              myfile << ptvals(i) << " " <<
                ptorig << " " << ptcmp << " " <<
                (ptcmp - ptorig) <<  " " <<
                paramslt(i, SPLIT_PHIIDX) << " " << phicmp[i] << " " <<
                (phicmp[i] - paramslt(i, SPLIT_PHIIDX)) << " " << 
                paramslt(i, SPLIT_X0IDX_NS) << " " << x0cmp[i] << " " <<  
                x0cmp[i] - paramslt(i, SPLIT_X0IDX_NS) << std::endl;
            else
              myfile << ptvals(i) << " " <<
                ptorig << " " << ptcmp << " " <<
                (ptcmp - ptorig) <<  " " <<
                paramslt(i, SPLIT_PHIIDX) << " " << phicmp[i] << " " <<
                (phicmp[i] - paramslt(i, SPLIT_PHIIDX)) << std::endl;
          
            if (verbose)
            {
              std::cout << "For track : " << i+1 << std::endl;
              std::cout << " pt           fitt " << ptcmp << std::endl;
              std::cout << " pt           orig " << ptorig << std::endl;
              std::cout << " phi          fitt " << phicmp[i] << std::endl;
              std::cout << " phi          orig " << paramslt(i, SPLIT_PHIIDX) << std::endl;
              if (usealsod0)
              { 
                std::cout << " d0           fitt " << d0cmp[i] << std::endl;
                std::cout << " d0           orig " << paramslt(i, SPLIT_D0IDX) << std::endl;
              }
              else if (usealsox0)
              {
                std::cout << " x0           fitt " << x0cmp[i] << std::endl;
                std::cout << " x0           orig " << paramslt(i, SPLIT_X0IDX_NS) << std::endl;
              }
            }
          }
        }
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

  if ((singleparam >= 1) && (singleparam <= 7))
  {
    delete [] singlep;
  }
  else
  {
    if (usealsod0)
      delete [] d0cmp;

    if (usealsox0)
      delete [] x0cmp;

    if (rzplane)
    {
      if (usex0y0)
      {
        delete [] x0cmp;
        delete [] y0cmp;
      }
      else
      {
        delete [] cothetacmp;
        delete [] z0cmp;
      }
    }
    else if (rphiplane)
    {
      if (usex0y0)
      {
        delete [] x0cmp;
        delete [] y0cmp;
      }
      else
      {
        delete [] oneoverptcmp;
        delete [] phicmp;
      }
    }
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
  std::cerr << " -q, --vmtx=[fillename]          : VMTX filename [default is v.[rz/rphi].bin]" << std::endl;
  std::cerr << " -j, --jump-tracks               : perform the fittin only for odd tracks" << std::endl;
  std::cerr << " -e, --not-use-charge            : do not read charge from coordinatesfile " << std::endl;
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
  std::cerr << " -d, --use-d0                    : use also d0 param in both planes " << std::endl;
  std::cerr << " -X, --use-x0                    : use also x0 param in both planes " << std::endl;
  std::cerr << " -f, --fit-x0y0                  : use and fit x0 and y0 param in both planes instead of " << std::endl;
  std::cerr << "                                   eta, pt, z0, phi " << std::endl;
  std::cerr << " -s, --fit-single-param=[num]    : use and fit X param in both planes  " << std::endl;
  std::cerr << "                                   1=eta, 2=pt, 3=z0, 4=phi, 5=x0, 6=y0, 7=d0 " << std::endl;

  exit(1);
}

int main (int argc, char ** argv)
{
  pca::pcafitter fitter;

  std::string qfname = "";
  std::string cfname = "";
  std::string afname = "";
  std::string vfname = "";
  std::string subsec = "";
  std::string sublad = "";
  bool verbose = false;
  bool useonlyodd = false;
  bool rzplane = false;
  bool rphiplane = false;
  bool usecharge = true;
  bool usealsod0 = false;
  bool usex0y0 = false;
  bool usesingleparam = false;
  bool usealsox0 = false;
  bool checklayersids = false;
  bool userelativecoord = false;

  int singleparam=-1;

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
      {"verbose", 0, NULL, 'V'},
      {"version", 0, NULL, 'v'},
      {"jump-tracks", 0, NULL, 'j'},
      {"rz-plane", 0, NULL, 'z'},
      {"rphi-plane", 0, NULL, 'r'},
      {"not-use-charge", 0, NULL, 'e'},
      {"charge-sign", 1, NULL, 'g'},
      {"use-d0", 0, NULL, 'd'},
      {"use-x0", 0, NULL, 'X'},
      {"exclude-s-module", 0, NULL, 'x'},
      {"fit-x0y0", 0, NULL, 'f'},
      {"fit-single-param", 1, NULL, 's'},
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

    c = getopt_long (argc, argv, "XfkdxezrhaVjb:A:B:t:g:c:q:s:n:s:m:o:u", 
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
      case 'j':
        useonlyodd = true;
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
    if (usealsod0 || usealsox0)
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
        if (!fitter.set_paramidx(SPLIT_COTTHETAIDX, "cot(theta)"))
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

  // leggere file coordinate tracce e file costanti PCA
  // N righe di 9 double sono le coordinate
  // matrice C e vettore q sono le costanti
  
  arma::mat cmtx, amtx, vmtx;
  arma::rowvec q;

  // leggere file coordinate tracce simulate plus parametri
  if (!pca::file_exists(filename))
  {
    std::cerr << "Inout file does not exist" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Reading data from " << filename << " file " << std::endl;

  arma::mat coord, param;
  arma::vec ptvals;

  if (!pca::reading_from_file_split (fitter, filename, 
       param, coord, false, useonlyodd,
       rzplane, rphiplane, etamin, etamax, ptmin, ptmax, 
       usecharge, chargesign, excludesmodule, usealsod0,
       usex0y0, singleparam, phimin, phimax, z0min, z0max,
       d0min, d0max, usealsox0, verbose, ptvals, 
       checklayersids, 6))
    return EXIT_FAILURE;

  if (userelativecoord)
    pca::global_to_relative(coord, coord1min, coord2min);

  std::cout << "Using " << param.n_rows << " tracks" << std::endl;

  if (rzplane)
  {
    if (usesingleparam)
    {
      std::string name = pca::get_paramname_from_id(singleparam)+".txt"; 
      pca::write_to_file(("tofit_"+name).c_str(), param, 0);
    }
    else
    {
      if (usex0y0)
      {
        pca::write_to_file("tofit_x0.txt", param, SPLIT_X0IDX);
        pca::write_to_file("tofit_y0.txt", param, SPLIT_Y0IDX);
      }
      else
      {
        pca::write_to_file("tofit_cottheta.txt", param, SPLIT_COTTHETAIDX);
        pca::write_to_file("tofit_z0.txt", param, SPLIT_Z0IDX);
      }

      if (usealsod0)
        pca::write_to_file("tofit_d0.txt", param, SPLIT_D0IDX);
      else if (usealsox0)
        pca::write_to_file("tofit_x0.txt", param, SPLIT_X0IDX_NS);
    }

    if (cfname == "")
      cfname = "c.rz.bin";

    if (qfname == "")
      qfname = "q.rz.bin";

    if (afname == "")
      afname = "a.rz.bin";

    if (vfname == "")
      vfname = "v.rz.bin";
  }
  else if (rphiplane)
  {
    if (usesingleparam)
    {
      std::string name = pca::get_paramname_from_id(singleparam)+".txt"; 
      pca::write_to_file(("tofit_"+name).c_str(), param, 0);
    }
    else
    {
      if (usex0y0)
      {
        pca::write_to_file("tofit_x0.txt", param, SPLIT_X0IDX);
        pca::write_to_file("tofit_y0.txt", param, SPLIT_Y0IDX);
      }
      else
      {
        pca::write_to_file("tofit_phi.txt", param, SPLIT_PHIIDX);
        pca::write_to_file("tofit_oneoverpt.txt", param, SPLIT_ONEOVERPTIDX);
      }

      if (usealsod0)
        pca::write_to_file("tofit_d0.txt", param, SPLIT_D0IDX);
      else if (usealsox0)
        pca::write_to_file("tofit_x0.txt", param, SPLIT_X0IDX_NS);
    }

    if (cfname == "")
      cfname = "c.rphi.bin";

    if (qfname == "")
      qfname = "q.rphi.bin";

    if (afname == "")
      afname = "a.rphi.bin";

    if (vfname == "")
      vfname = "v.rphi.bin";
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

  if (!build_and_compare (param, coord, cmtx, q, amtx, vmtx, 
        verbose, fitter, rzplane, rphiplane, usecharge, 
        usealsod0, usex0y0, singleparam, usealsox0, ptvals))
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
