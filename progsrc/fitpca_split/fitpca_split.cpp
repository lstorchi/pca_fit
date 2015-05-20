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

// lstorchi: basi code to fit tracks, using the PCA constants generated 
//           by the related generatepca


bool build_and_compare (arma::mat & paramslt, arma::mat & coordslt, 
     arma::mat & cmtx, arma::rowvec & q, bool verbose, 
     pca::pcafitter & fitter, bool rzplane, bool rphiplane,
     bool usecharge, bool usealsod0, bool usex0y0, int singleparam)
{
  double ** ptrs;
  ptrs = new double* [fitter.get_paramdim()];

  double * cotethacmp, * z0cmp, * oneoverptcmp, * phicmp, 
         * d0cmp, * x0cmp, * y0cmp, *singlep;

  if ((singleparam >= 1) && (singleparam <= 7))
  {
    singlep = new double [(int)coordslt.n_rows];
    ptrs[0] = singlep;
  }
  else
  {
    if (rzplane)
    {
      if (usex0y0)
      {
        x0cmp = new double [(int)coordslt.n_rows];
        y0cmp = new double [(int)coordslt.n_rows];
        ptrs[SPLIT_X0IDX] = x0cmp;
        ptrs[SPLIT_Y0IDX] = y0cmp;
      }
      else
      {
        cotethacmp = new double [(int)coordslt.n_rows];
        z0cmp = new double [(int)coordslt.n_rows];
        ptrs[SPLIT_COTTETHAIDX] = cotethacmp;
        ptrs[SPLIT_Z0IDX] = z0cmp;
      }
      
      if (usealsod0)
      {
        d0cmp = new double [(int)coordslt.n_rows];
        ptrs[SPLIT_D0IDX] = d0cmp;
      }
    }
    else if (rphiplane)
    {
      if (usex0y0)
      {
        x0cmp = new double [(int)coordslt.n_rows];
        y0cmp = new double [(int)coordslt.n_rows];
        ptrs[SPLIT_X0IDX] = x0cmp;
        ptrs[SPLIT_Y0IDX] = y0cmp;
      }
      else 
      {
        oneoverptcmp = new double [(int)coordslt.n_rows];
        phicmp = new double [(int)coordslt.n_rows];
        ptrs[SPLIT_ONEOVERPTIDX] = oneoverptcmp;
        ptrs[SPLIT_PHIIDX] = phicmp;
      }
    
      if (usealsod0)
      {
        d0cmp = new double [(int)coordslt.n_rows];
        ptrs[SPLIT_D0IDX] = d0cmp;
      }
    }
  }

  if (!fitter.compute_parameters (cmtx, q, coordslt, ptrs, 
        fitter.get_paramdim()))
  {
    std::cerr << fitter.get_errmsg() << std::endl;
    return false;
  }

  delete [] ptrs; 

  std::ostringstream fname;
  fname << "results.txt";

  arma::running_stat<double> pc[fitter.get_paramdim()];
  arma::running_stat<double> pcabsolute[fitter.get_paramdim()];
  std::ofstream myfile(fname.str().c_str());

  if ((singleparam >= 1) && (singleparam <= 7))
  {
    myfile << "orig fitt diff " << std::endl;
    for (int i=0; i<(int)coordslt.n_rows; ++i)
    {
      pc[0](fabs(singlep[i] + paramslt(i, 0))/
          (fabs(singlep[i] + paramslt(i, 0))/2.0));
    
      pcabsolute[0](singlep[i] + paramslt(i, 0));
    
      myfile << 
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
          myfile << "x0_orig x0_fitt diff y0_orig y0_fitt diff " <<
           " d0_orig d0_fitt diff" << std::endl;
        else  
          myfile << "x0_orig x0_fitt diff y0_orig y0_fitt diff" << std::endl; 
    
        for (int i=0; i<(int)coordslt.n_rows; ++i)
        {
          pc[SPLIT_X0IDX](fabs(x0cmp[i] + paramslt(i, SPLIT_X0IDX))/
              (fabs(x0cmp[i] + paramslt(i, SPLIT_X0IDX))/2.0));
          pc[SPLIT_Y0IDX](fabs(y0cmp[i] - paramslt(i, SPLIT_Y0IDX))/
              (fabs(y0cmp[i] + paramslt(i, SPLIT_Y0IDX))/2.0));
        
          pcabsolute[SPLIT_X0IDX](x0cmp[i] + paramslt(i, SPLIT_X0IDX));
          pcabsolute[SPLIT_Y0IDX](y0cmp[i] - paramslt(i, SPLIT_Y0IDX));
        
          if (usealsod0) 
          {
            pc[SPLIT_D0IDX](fabs(d0cmp[i] - paramslt(i, SPLIT_D0IDX))/
                (fabs(d0cmp[i] + paramslt(i, SPLIT_D0IDX))/2.0));
            pcabsolute[SPLIT_D0IDX](d0cmp[i] - paramslt(i, SPLIT_D0IDX));
          }
          
          if (usealsod0)
            myfile << 
              paramslt(i, SPLIT_X0IDX) << " " << x0cmp[i] << " " <<
              (x0cmp[i] - paramslt(i, SPLIT_X0IDX)) << " " <<
              paramslt(i, SPLIT_Y0IDX) << " " << y0cmp[i] << " " <<
              (y0cmp[i] - paramslt(i, SPLIT_Y0IDX)) << " " << 
              paramslt(i, SPLIT_D0IDX) << " " << d0cmp[i] << " " <<
              d0cmp[i] - paramslt(i, SPLIT_D0IDX) << std::endl;
          else
            myfile << 
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
          myfile << "eta_orig eta_fitt diff z0_orig z0_fitt diff " <<
           " d0_orig d0_fitt diff" << std::endl;
        else  
          myfile << "eta_orig eta_fitt diff z0_orig z0_fitt diff" << std::endl; 
    
        for (int i=0; i<(int)coordslt.n_rows; ++i)
        {
          double tethacmp = atan(1.0e0 / cotethacmp[i]) ; 
          double etacmps = 0.0e0, tantetha2;
          tantetha2 = tan (tethacmp/2.0e0); 
          if (tantetha2 < 0.0)
            etacmps = 1.0e0 * log (-1.0e0 * tantetha2);
          else
            etacmps = -1.0e0 * log (tantetha2);
          double tetha = atan(1.0e0 / paramslt(i, SPLIT_COTTETHAIDX));
          double etaorig = 0.0e0;
          tantetha2 = tan (tetha/2.0e0);
          if (tantetha2 < 0.0)
            etaorig = 1.0e0 * log (-1.0e0 * tantetha2);
          else
            etaorig = -1.0e0 * log (tantetha2);
        
          pc[SPLIT_COTTETHAIDX](fabs(etacmps - etaorig)/
              (fabs(etacmps + etaorig)/2.0));
          pc[SPLIT_Z0IDX](fabs(z0cmp[i] - paramslt(i, SPLIT_Z0IDX))/
              (fabs(z0cmp[i] + paramslt(i, SPLIT_Z0IDX))/2.0));
        
          pcabsolute[SPLIT_COTTETHAIDX](etacmps - etaorig);
          pcabsolute[SPLIT_Z0IDX](z0cmp[i] - paramslt(i, SPLIT_Z0IDX));
        
          if (usealsod0) 
          {
            pc[SPLIT_D0IDX](fabs(d0cmp[i] - paramslt(i, SPLIT_D0IDX))/
                (fabs(d0cmp[i] + paramslt(i, SPLIT_D0IDX))/2.0));
            pcabsolute[SPLIT_D0IDX](d0cmp[i] - paramslt(i, SPLIT_D0IDX));
          }
          
          if (usealsod0)
            myfile << 
              etaorig << " " << etacmps << " " <<
              (etacmps - etaorig) << " " <<
              paramslt(i, SPLIT_Z0IDX) << " " << z0cmp[i] << " " <<
              (z0cmp[i] - paramslt(i, SPLIT_Z0IDX)) << " " << 
              paramslt(i, SPLIT_D0IDX) << " " << d0cmp[i] << " " <<
              d0cmp[i] - paramslt(i, SPLIT_D0IDX) << std::endl;
          else
            myfile << 
              etaorig << " " << etacmps << " " <<
              (etacmps - etaorig) << " " <<
              paramslt(i, SPLIT_Z0IDX) << " " << z0cmp[i] << " " <<
              (z0cmp[i] - paramslt(i, SPLIT_Z0IDX)) << std::endl;
        
          if (verbose)
          {
            std::cout << "For track : " << i+1 << std::endl;
            std::cout << " cotetha      fitt " << cotethacmp[i] << std::endl;
            std::cout << " cotetha      orig " << paramslt(i, SPLIT_COTTETHAIDX) << std::endl;
            std::cout << " tetha rad    fitt " << tethacmp << std::endl;
            std::cout << " tetha rad    orig " << tetha << std::endl;
            std::cout << " tetha deg    fitt " << tethacmp*(180.0e0/M_PI) << std::endl;
            std::cout << " tetha deg    orig " << tetha*(180.0e0/M_PI) << std::endl;
            std::cout << " eta          fitt " << etacmps << std::endl;
            std::cout << " eta          orig " << etaorig << std::endl;
            std::cout << " z0           fitt " << z0cmp[i] << std::endl;
            std::cout << " z0           orig " << paramslt(i, SPLIT_Z0IDX) << std::endl;
            if (usealsod0)
            { 
              std::cout << " d0           fitt " << d0cmp[i] << std::endl;
              std::cout << " d0           orig " << paramslt(i, SPLIT_D0IDX) << std::endl;
            }
          }
        }
      }
    }
    else if (rphiplane)
    {
      std::ofstream myfile(fname.str().c_str());
      arma::running_stat<double> pc[fitter.get_paramdim()];
    
      if (usex0y0)
      {
        if (usealsod0)
          myfile << "x0_orig x0_fitt diff y0_orig y0_fitt diff " <<
           " d0_orig d0_fitt diff" << std::endl;
        else  
          myfile << "x0_orig x0_fitt diff y0_orig y0_fitt diff" << std::endl; 
    
        for (int i=0; i<(int)coordslt.n_rows; ++i)
        {
          pc[SPLIT_X0IDX](fabs(x0cmp[i] + paramslt(i, SPLIT_X0IDX))/
              (fabs(x0cmp[i] + paramslt(i, SPLIT_X0IDX))/2.0));
          pc[SPLIT_Y0IDX](fabs(y0cmp[i] - paramslt(i, SPLIT_Y0IDX))/
              (fabs(y0cmp[i] + paramslt(i, SPLIT_Y0IDX))/2.0));
        
          pcabsolute[SPLIT_X0IDX](x0cmp[i] + paramslt(i, SPLIT_X0IDX));
          pcabsolute[SPLIT_Y0IDX](y0cmp[i] - paramslt(i, SPLIT_Y0IDX));
        
          if (usealsod0) 
          {
            pc[SPLIT_D0IDX](fabs(d0cmp[i] - paramslt(i, SPLIT_D0IDX))/
                (fabs(d0cmp[i] + paramslt(i, SPLIT_D0IDX))/2.0));
            pcabsolute[SPLIT_D0IDX](d0cmp[i] - paramslt(i, SPLIT_D0IDX));
          }
          
          if (usealsod0)
            myfile << 
              paramslt(i, SPLIT_X0IDX) << " " << x0cmp[i] << " " <<
              (x0cmp[i] - paramslt(i, SPLIT_X0IDX)) << " " <<
              paramslt(i, SPLIT_Y0IDX) << " " << y0cmp[i] << " " <<
              (y0cmp[i] - paramslt(i, SPLIT_Y0IDX)) << " " << 
              paramslt(i, SPLIT_D0IDX) << " " << d0cmp[i] << " " <<
              d0cmp[i] - paramslt(i, SPLIT_D0IDX) << std::endl;
          else
            myfile << 
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
            myfile << "q/pt_orig q/pt_fitt diff phi_orig phi_fitt diff " <<
             " d0_orig d0_fitt diff" << std::endl;
          else
            myfile << "q/pt_orig q/pt_fitt diff phi_orig phi_fitt diff" << std::endl; 
        
          for (int i=0; i<(int)coordslt.n_rows; ++i)
          {
            double qoverptorig = paramslt(i, SPLIT_ONEOVERPTIDX);
            double qoverptcmp = oneoverptcmp[i];
          
            pc[SPLIT_PHIIDX](fabs(phicmp[i] - paramslt(i, SPLIT_PHIIDX))/
                (fabs(phicmp[i] + paramslt(i, SPLIT_PHIIDX))/2.0));
            pc[SPLIT_ONEOVERPTIDX](fabs(qoverptcmp - qoverptorig)/
                (fabs(qoverptcmp + qoverptorig)/2.0));
          
            pcabsolute[SPLIT_PHIIDX](phicmp[i] - paramslt(i, SPLIT_PHIIDX));
            pcabsolute[SPLIT_ONEOVERPTIDX](qoverptcmp - qoverptorig);
        
            if (usealsod0) 
            {
              pc[SPLIT_D0IDX](fabs(d0cmp[i] - paramslt(i, SPLIT_D0IDX))/
                  (fabs(d0cmp[i] + paramslt(i, SPLIT_D0IDX))/2.0));
              pcabsolute[SPLIT_D0IDX](d0cmp[i] - paramslt(i, SPLIT_D0IDX));
            }
          
            if (usealsod0)
              myfile << 
                qoverptorig << " " << qoverptcmp << " " <<
                (qoverptcmp - qoverptorig) <<  " " <<
                paramslt(i, SPLIT_PHIIDX) << " " << phicmp[i] << " " <<
                (phicmp[i] - paramslt(i, SPLIT_PHIIDX)) << " " << 
                paramslt(i, SPLIT_D0IDX) << " " << d0cmp[i] << " " <<
                d0cmp[i] - paramslt(i, SPLIT_D0IDX) << std::endl;
            else
              myfile << 
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
            }
          }
        }
        else 
        {
          if (usealsod0)
            myfile << "pt_orig pt_fitt diff phi_orig phi_fitt diff " <<
             " d0_orig d0_fitt diff" << std::endl;
          else
            myfile << "pt_orig pt_fitt diff phi_orig phi_fitt diff" << std::endl; 
        
          for (int i=0; i<(int)coordslt.n_rows; ++i)
          {
            double ptorig = 1.0e0 / paramslt(i, SPLIT_ONEOVERPTIDX);
            double ptcmp = 1.0e0 / oneoverptcmp[i];
          
            pc[SPLIT_PHIIDX](fabs(phicmp[i] - paramslt(i, SPLIT_PHIIDX))/
                (fabs(phicmp[i] + paramslt(i, SPLIT_PHIIDX))/2.0));
            pc[SPLIT_ONEOVERPTIDX](fabs(ptcmp - ptorig)/
                (fabs(ptcmp + ptorig)/2.0));
          
            pcabsolute[SPLIT_PHIIDX](phicmp[i] - paramslt(i, SPLIT_PHIIDX));
            pcabsolute[SPLIT_ONEOVERPTIDX](ptcmp - ptorig);
        
            if (usealsod0) 
            {
              pc[SPLIT_D0IDX](fabs(d0cmp[i] - paramslt(i, SPLIT_D0IDX))/
                  (fabs(d0cmp[i] + paramslt(i, SPLIT_D0IDX))/2.0));
              pcabsolute[SPLIT_D0IDX](d0cmp[i] - paramslt(i, SPLIT_D0IDX));
            }
        
            if (usealsod0)
              myfile << 
                ptorig << " " << ptcmp << " " <<
                (ptcmp - ptorig) <<  " " <<
                paramslt(i, SPLIT_PHIIDX) << " " << phicmp[i] << " " <<
                (phicmp[i] - paramslt(i, SPLIT_PHIIDX)) << " " << 
                paramslt(i, SPLIT_D0IDX) << " " << d0cmp[i] << " " <<  
                d0cmp[i] - paramslt(i, SPLIT_D0IDX) << std::endl;
            else
              myfile << 
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
            }
          }
        }
      }
    }
  }

  myfile.close();

  for (int i=0; i<fitter.get_paramdim(); ++i)
     std::cout << "For " << fitter.paramidx_to_string(i) << " error " << 
       pcabsolute[i].mean() << " " << pcabsolute[i].stddev() << std::endl;

  if ((singleparam >= 1) && (singleparam <= 7))
  {
    delete [] singlep;
  }
  else
  {
    if (rzplane)
    {
      if (usex0y0)
      {
        delete [] x0cmp;
        delete [] y0cmp;
      }
      else
      {
        delete [] cotethacmp;
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
  std::cerr << " -j, --jump-tracks               : perform the fittin only for odd tracks" << std::endl;
  std::cerr << " -z, --rz-plane                  : use rz plane view (fit eta and z0)" << std::endl;
  std::cerr << " -r, --rphi-plane                : use r-phi plane view (fit ot and phi)" << std::endl;
  std::cerr << " -e, --not-use-charge            : do not read charge from coordinatesfile " << std::endl;
  std::cerr << " -g, --charge-sign=[+/-]         : use only + particle or - paricle (again both planes)" << std::endl;
  std::cerr << " -t, --eta-range=\"etamin;etamax\" : specify the eta range to use " << std::endl;
  std::cerr << " -n, --pt-range=\"ptmin;ptmax\"    : specify the pt range to use " << std::endl;
  std::cerr << " -x, --exclude-s-module          : exclude S-module (last three layer) so 6 " << 
    "coordinates inseatd of 12 " << std::endl;
  std::cerr << " -d, --use-d0                    : use also d0 param in both planes " << std::endl;
  std::cerr << " -f, --fit-x0y0                  : use and fit x0 and y0 param in both planes instead of " << std::endl;
  std::cerr << "                                   eta, pt, z0, phi " << std::endl;
  std::cerr << " -s, --fit-single-param=[num]    : use and fit X param in both planes  " << std::endl;
  std::cerr << "                                   1=eta, 2=pt, 3=z0, 4=phi, 5=x0, 6=y0, 7=d0 " << std::endl;

  exit(1);
}

int main (int argc, char ** argv)
{
  pca::pcafitter fitter;

  std::string qfname;
  std::string cfname;
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

  int singleparam=-1;

  double etamin = -1.0e0 * INFINITY, etamax = +1.0e0 * INFINITY;
  double ptmin = -1.0e0 * INFINITY, ptmax = +1.0e0 * INFINITY;

  int chargesign = 0;

  std::vector<std::string> tokens;

  bool excludesmodule = false;

  while (1)
  {
    int c, option_index;
    static struct option long_options[] = {
      {"help", 0, NULL, 'h'},
      {"cmtx", 1, NULL, 'c'},
      {"qvct", 1, NULL, 'q'},
      {"verbose", 0, NULL, 'V'},
      {"version", 0, NULL, 'v'},
      {"jump-tracks", 0, NULL, 'j'},
      {"rz-plane", 0, NULL, 'z'},
      {"rphi-plane", 0, NULL, 'r'},
      {"not-use-charge", 0, NULL, 'e'},
      {"charge-sign", 1, NULL, 'g'},
      {"eta-range", 1, NULL, 't'},
      {"pt-range", 1, NULL, 'n'},
      {"use-d0", 0, NULL, 'd'},
      {"exclude-s-module", 0, NULL, 'x'},
      {"fit-x0y0", 0, NULL, 'f'},
      {"fit-single-param", 1, NULL, 's'},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "fdxezrhVjt:g:c:q:s:n:s:", long_options, &option_index);

    if (c == -1)
      break;

    switch (c)
    {
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

  if (excludesmodule)
    fitter.set_coordim (2*3);
  else
    fitter.set_coordim (2*6);

  if (usesingleparam)
    fitter.set_paramdim(1);
  else
  {
    if (usealsod0)
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
  
  arma::mat cmtx;
  arma::rowvec q;

  // leggere file coordinate tracce simulate plus parametri
  if (!pca::file_exists(filename))
  {
    std::cerr << "Inout file does not exist" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Reading data from " << filename << " file " << std::endl;
  int num_of_line = pca::numofline(filename);
  std::cout << "file has " << num_of_line << " line " << std::endl;
  int num_of_ent_read = (num_of_line-1)/ENTDIM;

  int num_of_ent = num_of_ent_read;

  if (useonlyodd)
  {
    if (num_of_ent_read % 2)
      num_of_ent = (num_of_ent_read+1)/2;
    else
      num_of_ent = num_of_ent_read/2;
  }

  std::cout << "file has " << num_of_ent << " entries " << std::endl;

  arma::mat coord, param;
  coord.set_size(num_of_ent,fitter.get_coordim());
  param.set_size(num_of_ent,fitter.get_paramdim());

  std::cout << "Reading from file" << std::endl;
  if (!pca::reading_from_file_split (fitter, filename, 
       param, coord, num_of_ent, false, useonlyodd,
       rzplane, rphiplane, etamin, etamax, ptmin, ptmax, 
       usecharge, chargesign, excludesmodule, usealsod0,
       usex0y0, singleparam))
    return EXIT_FAILURE;
  std::cout << "Using " << param.n_rows << " tracks" << std::endl;

  if (rzplane)
  {
    cfname = "c.rz.bin";
    qfname = "q.rz.bin";
  }
  else if (rphiplane)
  {
    cfname = "c.rphi.bin";
    qfname = "q.rphi.bin";
  }

  pca::read_armmat(cfname.c_str(), cmtx);
  pca::read_armvct(qfname.c_str(), q);

  if (!build_and_compare (param, coord, cmtx, q, verbose, 
          fitter, rzplane, rphiplane, usecharge, usealsod0,
          usex0y0, singleparam))
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
