#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>

// Loriano: let's try Armadillo quick code 
#include <armadillo>
#include <cassert>

#include <getopt.h>
#include <unistd.h>
#include <alloca.h>

#include <pcafitter_private.hpp>

void usage (char * name)
{
  std::cerr << "usage: " << name << " [options] coordinatesfile " << std::endl;
  std::cerr << std::endl;
  std::cerr << " -h, --help               : display this help and exit" << std::endl;
  std::cerr << " -V, --verbose            : verbose option on" << std::endl;
  std::cerr << " -v, --version            : print version and exit" << std::endl;
  std::cerr << " -c, --cmtx=[fillename]   : CMTX filename [default is c.bin]" << std::endl;
  std::cerr << " -q, --qvct=[fillename]   : QVCT filename [default is q.bin]" << std::endl;
  std::cerr << " -s, --subsector=[subsec] : by default use values of the bigger subsector" << std::endl;
  std::cerr << "                            with this option you can speficy to perform " << std::endl;
  std::cerr << "                            prediction for subsector subsec " << std::endl;
  std::cerr << " -l, --subladder=[subld]  : by default use values of the bigger subladder " << std::endl;
  std::cerr << "                            with this option you can speficy to perform " << std::endl;
  std::cerr << "                            prediction for subladder subld " << std::endl;

  exit(1);
}

int main (int argc, char ** argv)
{
  std::string qfname = "q.bin";
  std::string cfname = "c.bin";
  std::string subsec = "";
  std::string sublad = "";
  bool verbose = false;

  while (1)
  {
    int c, option_index;
    static struct option long_options[] = {
      {"help", 0, NULL, 'h'},
      {"cmtx", 1, NULL, 'c'},
      {"qvct", 1, NULL, 'q'},
      {"subsector", 1, NULL, 's'},
      {"subladder", 1, NULL, 'l'},
      {"verbose", 0, NULL, 'V'},
      {"version", 0, NULL, 'v'},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "vhVc:q:s:l:", long_options, &option_index);

    if (c == -1)
      break;

    switch (c)
    {
      case 'V':
        verbose = true;
        break;
      case 'v':
        std::cout << "Version: " << pcafitter::get_version_string() << std::endl;
        exit(1);
        break;
      case 'h':
        usage (argv[0]);
        break;
      case 's':
        subsec = optarg;
        break;
      case 'l':
        sublad = optarg;
        break;
      case'c':
        cfname = optarg;
        break;
      case 'q':
        qfname = optarg;
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

  // leggere file coordinate tracce e file costanti PCA
  // N righe di 9 double sono le coordinate
  // matrice C e vettore q sono le costanti
  
  arma::mat cmtx;
  arma::rowvec q;

  std::cout << "Read constant from files (" << cfname << 
    " and " << qfname << ")" << std::endl;
  pcafitter::readarmmat(cfname.c_str(), cmtx);
  pcafitter::readarmvct(qfname.c_str(), q);

  std::cout << "Reading data from " << filename << " file " << std::endl;
  int num_of_line = pcafitter::numofline(filename);
  std::cout << "file has " << num_of_line << " line " << std::endl;
  int num_of_ent = (num_of_line-1)/ENTDIM;
  std::cout << "file has " << num_of_ent << " entries " << std::endl;

  arma::rowvec pt = arma::zeros<arma::rowvec>(num_of_ent);
  arma::rowvec phi = arma::zeros<arma::rowvec>(num_of_ent);
  arma::rowvec d0 = arma::zeros<arma::rowvec>(num_of_ent);
  arma::rowvec eta = arma::zeros<arma::rowvec>(num_of_ent);
  arma::rowvec z0 = arma::zeros<arma::rowvec>(num_of_ent);

  arma::mat layer, ladder, module, coord;
  layer.set_size(num_of_ent,COORDIM);
  ladder.set_size(num_of_ent,COORDIM);
  module.set_size(num_of_ent,COORDIM);
  coord.set_size(num_of_ent,3*COORDIM);

  std::map<std::string, int> subsectors, subladders;
  std::vector<std::string> subladderslist, subsectorslist;

  // leggere file coordinate tracce simulate plus parametri
  pcafitter::readingfromfile (filename, pt, phi, d0, eta, z0, coord, 
      layer, ladder, module, subsectors, subladders, 
      subsectorslist, subladderslist, num_of_ent);

  if ((subsec == "") && (sublad == ""))
  {
    int maxnumber;

    std::cout << "Looking for bigger subsector" << std::endl; 
    std::cout << "We  found " << subsectors.size() << 
      " subsectors " << std::endl;
    
    pcafitter::select_bigger_sub (subsectors, verbose, 
        maxnumber, subsec);
    
    std::cout << "Selected subsector " << subsec << 
      " numevt: " << maxnumber << std::endl;

    std::cout << "Looking for bigger subladder" << std::endl; 
    std::cout << "We  found " << subladders.size() << 
      " subladders " << std::endl;
    
    pcafitter::select_bigger_sub (subladders, verbose, 
        maxnumber, sublad);
    
    std::cout << "Selected subladder " << sublad << " numevt: " 
      << maxnumber << std::endl;
  }

  if (subsec != "")
  {
    std::cout << "Using subsector " << subsec << std::endl;

    arma::rowvec ptslt, phislt, d0slt, etaslt, z0slt;
    arma::mat coordslt;

    pcafitter::extract_sub (subsectorslist, 
        subsec, pt, phi, d0, eta, z0, coord, 
        ptslt, phislt, d0slt, etaslt, z0slt, 
        coordslt);

    double * ptcmp, * phicmp, * etacmp, * z0cmp, * d0cmp;
    ptcmp = new double [(int)coord.n_rows];
    phicmp = new double [(int)coord.n_rows];
    etacmp = new double [(int)coord.n_rows];
    z0cmp = new double [(int)coord.n_rows];
    d0cmp = new double [(int)coord.n_rows];
    pcafitter::computeparameters (cmtx, q, coordslt, ptcmp,
        phicmp, etacmp, z0cmp, d0cmp);

    arma::running_stat<double> pc[PARAMDIM];
    for (int i=0; i<(int)coordslt.n_rows; ++i)
    {
      pc[PTIDX](fabs(ptcmp[i] - ptslt(i))/
          (fabs(ptcmp[i] + ptslt(i))/2.0));
      pc[PHIIDX](fabs(phicmp[i] - phislt(i))/
          (fabs(phicmp[i] + phislt(i))/2.0));
      pc[ETAIDX](fabs(etacmp[i] - etaslt(i))/
          (fabs(etacmp[i] + etaslt(i))/2.0));
      pc[D0IDX](fabs(d0cmp[i] - d0slt(i))/
          (fabs(d0cmp[i] + d0slt(i))/2.0));
      pc[Z0IDX](fabs(z0cmp[i] - z0slt(i))/
          (fabs(z0cmp[i] + z0slt(i))/2.0));

      if (verbose)
      {
        std::cout << "For track : " << i+1 << std::endl;
        std::cout << " pt  cmpt " << ptcmp[i] << std::endl;
        std::cout << " pt  calc " << ptslt(i) << std::endl;
        std::cout << " phi cmpt " << phicmp[i] << std::endl;
        std::cout << " phi calc " << phislt(i) << std::endl;
        std::cout << " eta cmpt " << etacmp[i] << std::endl;
        std::cout << " eta calc " << etaslt(i) << std::endl;
        std::cout << " d0  cmpt " << d0cmp[i] << std::endl;
        std::cout << " d0  calc " << d0slt(i) << std::endl;
        std::cout << " z0  cmpt " << z0cmp[i] << std::endl;
        std::cout << " z0  calc " << z0slt(i) << std::endl;
      }
    }

    for (int i=0; i<PARAMDIM; ++i)
       std::cout << "For " << pcafitter::paramidxtostring(i) << " error " << 
         100.0*pc[i].mean() << " " << 100.0*pc[i].stddev() << std::endl;

    delete [] ptcmp;
    delete [] phicmp;
    delete [] etacmp;
    delete [] z0cmp;
    delete [] d0cmp;
  }
  
  if (sublad != "")
  {
    std::cout << "Using subladder " << sublad << std::endl;

    arma::rowvec ptslt, phislt, d0slt, etaslt, z0slt;
    arma::mat coordslt;

    pcafitter::extract_sub (subladderslist, 
        sublad, pt, phi, d0, eta, z0, coord, 
        ptslt, phislt, d0slt, etaslt, z0slt, 
        coordslt);

    double * ptcmp, * phicmp, * etacmp, * z0cmp, * d0cmp;
    ptcmp = new double [(int)coord.n_rows];
    phicmp = new double [(int)coord.n_rows];
    etacmp = new double [(int)coord.n_rows];
    z0cmp = new double [(int)coord.n_rows];
    d0cmp = new double [(int)coord.n_rows];
    pcafitter::computeparameters (cmtx, q, coordslt, ptcmp,
        phicmp, etacmp, z0cmp, d0cmp);

    arma::running_stat<double> pc[PARAMDIM];
    for (int i=0; i<(int)coordslt.n_rows; ++i)
    {
      pc[PTIDX](fabs(ptcmp[i] - ptslt(i))/
          (fabs(ptcmp[i] + ptslt(i))/2.0));
      pc[PHIIDX](fabs(phicmp[i] - phislt(i))/
          (fabs(phicmp[i] + phislt(i))/2.0));
      pc[ETAIDX](fabs(etacmp[i] - etaslt(i))/
          (fabs(etacmp[i] + etaslt(i))/2.0));
      pc[D0IDX](fabs(d0cmp[i] - d0slt(i))/
          (fabs(d0cmp[i] + d0slt(i))/2.0));
      pc[Z0IDX](fabs(z0cmp[i] - z0slt(i))/
          (fabs(z0cmp[i] + z0slt(i))/2.0));

      if (verbose)
      {
        std::cout << "For track : " << i+1 << std::endl;
        std::cout << " pt  cmpt " << ptcmp[i] << std::endl;
        std::cout << " pt  calc " << ptslt(i) << std::endl;
        std::cout << " phi cmpt " << phicmp[i] << std::endl;
        std::cout << " phi calc " << phislt(i) << std::endl;
        std::cout << " eta cmpt " << etacmp[i] << std::endl;
        std::cout << " eta calc " << etaslt(i) << std::endl;
        std::cout << " d0  cmpt " << d0cmp[i] << std::endl;
        std::cout << " d0  calc " << d0slt(i) << std::endl;
        std::cout << " z0  cmpt " << z0cmp[i] << std::endl;
        std::cout << " z0  calc " << z0slt(i) << std::endl;
      }
    }

    for (int i=0; i<PARAMDIM; ++i)
       std::cout << "For " << pcafitter::paramidxtostring(i) << " error " << 
         100.0*pc[i].mean() << " " << 100.0*pc[i].stddev() << std::endl;

    delete [] ptcmp;
    delete [] phicmp;
    delete [] etacmp;
    delete [] z0cmp;
    delete [] d0cmp;
  }

  return 0;
}
