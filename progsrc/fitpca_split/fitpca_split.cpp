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
     pca::pcafitter & fitter)
{
  double * cotethacmp, * z0cmp;
  cotethacmp = new double [(int)coordslt.n_rows];
  z0cmp = new double [(int)coordslt.n_rows];

  double ** ptrs;
  ptrs = new double* [fitter.get_paramdim()];
  ptrs[SPLIT_COTTETHAIDX] = cotethacmp;
  ptrs[SPLIT_Z0IDX] = z0cmp;

  if (!fitter.compute_parameters (cmtx, q, coordslt, ptrs, 
        fitter.get_paramdim()))
  {
    std::cerr << fitter.get_errmsg() << std::endl;
    return false;
  }

  delete [] ptrs; 

  std::ostringstream fname;
  fname << "results.txt";

  std::ofstream myfile(fname.str().c_str());
  myfile << "eta_orig eta_cmpt diff z0_orig z0_cmpt diff" << std::endl; 

  arma::running_stat<double> pc[fitter.get_paramdim()];
  for (int i=0; i<(int)coordslt.n_rows; ++i)
  {
    double tantetha = (1.0e0 / cotethacmp[i]) ; 
    double etacmps = -1.0e0 * log (tantetha);
    tantetha = (1.0e0 / paramslt(i, SPLIT_COTTETHAIDX));
    double etaorig = -1.0e0 * log (tantetha); 
  
    pc[SPLIT_COTTETHAIDX](fabs(etacmps - etaorig)/
        (fabs(etacmps + etaorig)/2.0));
    pc[SPLIT_Z0IDX](fabs(z0cmp[i] - paramslt(i, SPLIT_Z0IDX))/
        (fabs(z0cmp[i] + paramslt(i, SPLIT_Z0IDX))/2.0));

    myfile << 
      etaorig << " " << etacmps << " " <<
      (etacmps - etaorig) << " " <<
      paramslt(i, SPLIT_Z0IDX) << " " << z0cmp[i] << " " <<
      (z0cmp[i] - paramslt(i, SPLIT_Z0IDX)) << std::endl;

    if (verbose)
    {
      std::cout << "For track : " << i+1 << std::endl;
      std::cout << " eta          cmpt " << etacmps << std::endl;
      std::cout << " eta          calc " << etaorig << std::endl;
      std::cout << " z0           cmpt " << z0cmp[i] << std::endl;
      std::cout << " z0           calc " << paramslt(i, SPLIT_Z0IDX) << std::endl;
    }
  }

  myfile.close();

  for (int i=0; i<fitter.get_paramdim(); ++i)
     std::cout << "For " << fitter.paramidx_to_string(i) << " error " << 
       100.0*pc[i].mean() << " " << 100.0*pc[i].stddev() << std::endl;

  delete [] cotethacmp;
  delete [] z0cmp;

  return true;
} 

void usage (char * name)
{
  std::cerr << "usage: " << name << " [options] coordinatesfile " << std::endl;
  std::cerr << std::endl;
  std::cerr << " -h, --help               : display this help and exit" << std::endl;
  std::cerr << " -V, --verbose            : verbose option on" << std::endl;
  std::cerr << " -v, --version            : print version and exit" << std::endl;
  std::cerr << " -c, --cmtx=[fillename]   : CMTX filename [default is c.<selectedsubsecid>.bin]" << std::endl;
  std::cerr << " -q, --qvct=[fillename]   : QVCT filename [default is q.<selectedsubsecid>.bin]" << std::endl;
  std::cerr << " -j, --jump-tracks        : perform the fittin only for odd tracks" << std::endl;

  exit(1);
}

int main (int argc, char ** argv)
{
  pca::pcafitter fitter;

  std::string qfname = "q.bin";
  std::string cfname = "c.bin";
  std::string subsec = "";
  std::string sublad = "";
  bool verbose = false;
  bool useonlyodd = false;

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
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "hVc:q:s:j", long_options, &option_index);

    if (c == -1)
      break;

    switch (c)
    {
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
      default:
        usage (argv[0]);
        break;
    } 
  }

  // R-z
  fitter.set_coordim (2*6);

  fitter.set_paramdim(2);
  if (!fitter.set_paramidx(SPLIT_COTTETHAIDX, "cot(tetha/2)"))
  {
    std::cerr << fitter.get_errmsg() << std::endl;
    return EXIT_FAILURE;
  }
  if (!fitter.set_paramidx(SPLIT_Z0IDX, "z0"))
  {
    std::cerr << fitter.get_errmsg() << std::endl;
    return EXIT_FAILURE;
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

  pca::reading_from_file_split_rz (filename, 
     param, coord, num_of_ent, false, useonlyodd);

  pca::read_armmat(cfname.c_str(), cmtx);
  pca::read_armvct(qfname.c_str(), q);

  if (!build_and_compare (param, coord, cmtx, q, verbose, 
          fitter))
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
