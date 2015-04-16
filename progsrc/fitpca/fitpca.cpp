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
     const std::string & postfname, pca::pcafitter & fitter)
{
  double * oneoverptcmp, * phicmp, * cottethacmp, * z0cmp, * d0cmp;
  oneoverptcmp = new double [(int)coordslt.n_rows];
  phicmp = new double [(int)coordslt.n_rows];
  cottethacmp = new double [(int)coordslt.n_rows];
  z0cmp = new double [(int)coordslt.n_rows];
  d0cmp = new double [(int)coordslt.n_rows];

  double ** ptrs;
  ptrs = new double* [fitter.get_paramdim()];
  ptrs[ONEOVERPTIDX] = oneoverptcmp;
  ptrs[PHIIDX] = phicmp;
  ptrs[COTTETHAIDX] = cottethacmp;
  ptrs[Z0IDX] = z0cmp;
  ptrs[D0IDX] = d0cmp;

  if (!fitter.compute_parameters (cmtx, q, coordslt, ptrs, 
        fitter.get_paramdim()))
  {
    std::cerr << fitter.get_errmsg() << std::endl;
    return false;
  }

  delete [] ptrs; 

  std::ostringstream fname;
  fname << "results." << postfname << ".txt";

  std::ofstream myfile(fname.str().c_str());
  myfile << "pt_orig pt_cmpt diff phi_orig " <<
    "phi_cmpt diff eta_orig eta_cmpt diff " << 
    "d0_orig d0_cmpt diff z0_orig z0_cmpt diff" << std::endl; 

  arma::running_stat<double> pc[fitter.get_paramdim()];
  for (int i=0; i<(int)coordslt.n_rows; ++i)
  {
    double tetha = atan((1.0e0 / cottethacmp[i])); 
    double tantetha2 = tan (tetha/2.0e0);
    double etacmp = 0.0;
    if (tantetha2 < 0.0)
      etacmp = 1.0e0 * log (-1.0e0 * tantetha2);
    else
      etacmp = -1.0e0 * log (tantetha2);
    tetha = atan((1.0e0 / paramslt(i, COTTETHAIDX)));
    tantetha2 = tan (tetha/2.0e0);
    double etaorig = 0.0e0;
    if (tantetha2 < 0.0)
      etaorig = 1.0e0 * log (-1.0e0 * tantetha2);
    else
      etaorig = -1.0e0 * log (tantetha2);
    double deltaphi = pca::delta_phi(phicmp[i], paramslt(i, PHIIDX)); 
  
    pc[ONEOVERPTIDX](fabs(1.0e0/oneoverptcmp[i] - 1.0e0/paramslt(i, ONEOVERPTIDX))/
        (fabs(1.0e0/oneoverptcmp[i] + 1.0e0/paramslt(i, ONEOVERPTIDX))/2.0));
    pc[PHIIDX](fabs(deltaphi)/
        (fabs(phicmp[i] + paramslt(i, PHIIDX))/2.0));
    pc[COTTETHAIDX](fabs(etacmp - etaorig)/
        (fabs(etacmp + etaorig)/2.0));
    pc[Z0IDX](fabs(z0cmp[i] - paramslt(i, Z0IDX))/
        (fabs(z0cmp[i] + paramslt(i, Z0IDX))/2.0));
    pc[D0IDX](fabs(d0cmp[i] - paramslt(i, D0IDX))/
        (fabs(d0cmp[i] + paramslt(i, D0IDX))/2.0));

    myfile << 
      1.0e0/paramslt(i, ONEOVERPTIDX) << " " << 1.0e0/oneoverptcmp[i] << " " << 
      (1.0e0/oneoverptcmp[i] - 1.0e0/paramslt(i, ONEOVERPTIDX)) << " " << 
      paramslt(i, PHIIDX) << " " << phicmp[i] << " " <<
      deltaphi << " " << 
      etaorig << " " << etacmp << " " <<
      (etacmp - etaorig) << " " <<
      paramslt(i, D0IDX) << " " << d0cmp[i] << " " <<
      (d0cmp[i] - paramslt(i, D0IDX)) << " " <<
      paramslt(i, Z0IDX) << " " << z0cmp[i] << " " <<
      (z0cmp[i] - paramslt(i, Z0IDX)) << std::endl;

    if (verbose)
    {
      std::cout << "For track : " << i+1 << std::endl;
      std::cout << " pt           cmpt " << 1.0e0/oneoverptcmp[i] << std::endl;
      std::cout << " pt           calc " << 1.0e0/paramslt(i, ONEOVERPTIDX) << std::endl;
      std::cout << " phi          cmpt " << phicmp[i] << std::endl;
      std::cout << " phi          calc " << paramslt(i, PHIIDX) << std::endl;
      std::cout << " eta          cmpt " << etacmp << std::endl;
      std::cout << " eta          calc " << etaorig << std::endl;
      std::cout << " d0           cmpt " << d0cmp[i] << std::endl;
      std::cout << " d0           calc " << paramslt(i, D0IDX) << std::endl;
      std::cout << " z0           cmpt " << z0cmp[i] << std::endl;
      std::cout << " z0           calc " << paramslt(i, Z0IDX) << std::endl;
    }
  }

  myfile.close();

  for (int i=0; i<fitter.get_paramdim(); ++i)
     std::cout << "For " << fitter.paramidx_to_string(i) << " error " << 
       100.0*pc[i].mean() << " " << 100.0*pc[i].stddev() << std::endl;

  delete [] oneoverptcmp;
  delete [] cottethacmp;
  delete [] phicmp;
  delete [] z0cmp;
  delete [] d0cmp;

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
  std::cerr << " -i, --seg-id             : use SegId instead of coordinates" << std::endl;
  std::cerr << " -s, --subsector=[subsec] : by default use values of the bigger subsector" << std::endl;
  std::cerr << "                            with this option you can speficy to perform " << std::endl;
  std::cerr << "                            prediction for subsector subsec " << std::endl;
  std::cerr << " -l, --subladder=[subld]  : by default use values of the bigger subladder " << std::endl;
  std::cerr << "                            with this option you can speficy to perform " << std::endl;
  std::cerr << "                            prediction for subladder subld " << std::endl;
  std::cerr << " -a, --all-subsectors     : perform the fitting for all subsectors" << std::endl;
  std::cerr << " -r, --all-subladders     : perform the fitting for all subladders" << std::endl;
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
  bool useallsubsectors = false;
  bool useallsubladders = false;
  bool usesegid = false;
  bool useonlyodd = false;
  bool cfnameset = false;
  bool qfnameset = false;

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
      {"all-subsectors", 0, NULL, 'a'},
      {"all-subladders", 0, NULL, 'r'},
      {"seg-id", 0, NULL, 'i'},
      {"jump-tracks", 0, NULL, 'j'},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "iarjvhVc:q:s:l:", long_options, &option_index);

    if (c == -1)
      break;

    switch (c)
    {
      case 'j':
        useonlyodd = true;
        break;
      case 'i':
        usesegid = true;
        break;
      case 'a':
        useallsubsectors = true;
        break;
      case 'r':
        useallsubladders = true;
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
      case 's':
        subsec = optarg;
        break;
      case 'l':
        sublad = optarg;
        break;
      case'c':
        cfname = optarg;
        cfnameset = true;
        break;
      case 'q':
        qfname = optarg;
        qfnameset = true;
        break;
      default:
        usage (argv[0]);
        break;
    } 
  }

  if (usesegid)
    fitter.set_coordim(6);

  fitter.set_paramdim(5);
  if (!fitter.set_paramidx(ONEOVERPTIDX, "oneoverpt"))
  {
    std::cerr << fitter.get_errmsg() << std::endl;
    return EXIT_FAILURE;
  }
  if (!fitter.set_paramidx(PHIIDX, "phi"))
  {
    std::cerr << fitter.get_errmsg() << std::endl;
    return EXIT_FAILURE;
  }
  if (!fitter.set_paramidx(COTTETHAIDX, "cot(tetha/2)"))
  {
    std::cerr << fitter.get_errmsg() << std::endl;
    return EXIT_FAILURE;
  }
  if (!fitter.set_paramidx(Z0IDX, "z0"))
  {
    std::cerr << fitter.get_errmsg() << std::endl;
    return EXIT_FAILURE;
  }
  if (!fitter.set_paramidx(D0IDX, "d0"))
  {
    std::cerr << fitter.get_errmsg() << std::endl;
    return EXIT_FAILURE;
  }

  if (optind >= argc) 
    usage (argv[0]);

  if (useallsubladders && useallsubsectors)
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

  arma::mat layer, ladder, module, coord, param;
  layer.set_size(num_of_ent,fitter.get_coordim());
  ladder.set_size(num_of_ent,fitter.get_coordim());
  module.set_size(num_of_ent,fitter.get_coordim());
  coord.set_size(num_of_ent,fitter.get_coordim());
  param.set_size(num_of_ent,fitter.get_paramdim());

  std::map<std::string, int> subsectors, subladders;
  std::vector<std::string> subladderslist, subsectorslist;


  pca::reading_from_file (filename, param, coord, 
      layer, ladder, module, subsectors, subladders, 
      subsectorslist, subladderslist, num_of_ent_read, 
      usesegid, false, useonlyodd);

  if (!useallsubsectors && !useallsubladders)
  {
    if ((subsec == "") && (sublad == ""))
    {
      int maxnumber;
    
      std::cout << "Looking for bigger subsector" << std::endl; 
      std::cout << "We  found " << subsectors.size() << 
        " subsectors " << std::endl;
      
      fitter.select_bigger_sub (subsectors, verbose, 
          maxnumber, subsec);
      
      std::cout << "Selected subsector " << subsec << 
        " numevt: " << maxnumber << std::endl;
    
      std::cout << "Looking for bigger subladder" << std::endl; 
      std::cout << "We  found " << subladders.size() << 
        " subladders " << std::endl;
      
      fitter.select_bigger_sub (subladders, verbose, 
          maxnumber, sublad);
      
      std::cout << "Selected subladder " << sublad << " numevt: " 
        << maxnumber << std::endl;
    }
    
    if (subsec != "")
    {
      std::ostringstream cfnames, qfnames; 
      cfnames << "c." << subsec << ".bin";
      qfnames << "q." << subsec << ".bin";

      if (!cfnameset)
        cfname = cfnames.str();
      
      if (!qfnameset)
        qfname = qfnames.str();

      std::cout << "Read constant from files (" << cfname << 
        " and " << qfname << ")" << std::endl;
      if (!pca::file_exists(cfname) || !pca::file_exists(qfname))
      {
        std::cerr << "Constants file does not exist" << std::endl;
        return EXIT_FAILURE;
      }
      pca::read_armmat(cfname.c_str(), cmtx);
      pca::read_armvct(qfname.c_str(), q);

      arma::mat paramslt, coordslt;
    
      std::cout << "Using subsector " << subsec << std::endl;
    
      fitter.extract_sub (subsectorslist, 
          subsec, param, coord, paramslt,
          coordslt);

      std::cout << " numevt: " << coordslt.n_rows << std::endl;
    
      if (!build_and_compare (paramslt, coordslt, cmtx, q, verbose, 
            subsec, fitter))
        return EXIT_FAILURE;
    }
    
    if (sublad != "")
    {
      std::ostringstream cfnames, qfnames; 
      cfnames << "c." << sublad << ".bin";
      qfnames << "q." << sublad << ".bin";

      if (!cfnameset)
        cfname = cfnames.str();
      
      if (!qfnameset)
        qfname = qfnames.str();

      std::cout << "Read constant from files (" << cfname << 
        " and " << qfname << ")" << std::endl;
      if (!pca::file_exists(cfname) || !pca::file_exists(qfname))
      {
        std::cerr << "Constants file does not exist" << std::endl;
        return EXIT_FAILURE;
      }
      pca::read_armmat(cfname.c_str(), cmtx);
      pca::read_armvct(qfname.c_str(), q);

      arma::mat paramslt, coordslt;

      std::cout << "Using subladder " << sublad << std::endl;
    
      fitter.extract_sub (subladderslist, 
          sublad, param, coord, paramslt, 
          coordslt);

      std::cout << " numevt: " << coordslt.n_rows << std::endl;
     
      if (!build_and_compare (paramslt, coordslt, cmtx, q, verbose, 
            sublad, fitter))
        return EXIT_FAILURE;
    }
  }
  else
  {
    assert(useallsubladders != useallsubsectors);
    std::set<std::string> sectrorset;

    std::vector<std::string> * listtouse = NULL;

    if (useallsubsectors)
    {
      std::vector<std::string>::const_iterator selecteds = 
        subsectorslist.begin(); 
      for (; selecteds != subsectorslist.end(); ++selecteds)
        sectrorset.insert(*selecteds);

      listtouse = &subsectorslist;
    }
    else if (useallsubladders)
    {
      std::vector<std::string>::const_iterator selecteds = 
        subladderslist.begin(); 
      for (; selecteds != subladderslist.end(); ++selecteds)
        sectrorset.insert(*selecteds);

      listtouse = &subladderslist;
    }
    else 
      usage(argv[0]);

    std::set<std::string>::const_iterator selected = 
      sectrorset.begin(); 
    for (; selected != sectrorset.end(); ++selected)
    {
      std::ostringstream cfname, qfname; 
      cfname << "c." << *selected << ".bin";
      qfname << "q." << *selected << ".bin";

      if (pca::file_exists(cfname.str()) && pca::file_exists(qfname.str()))
      {
        std::cout << "Perfom fitting for " << *selected << std::endl;

        std::cout << "Read constants " << std::endl;
        pca::read_armmat(cfname.str().c_str(), cmtx);
        pca::read_armvct(qfname.str().c_str(), q);

        arma::mat paramslt, coordslt;
    
        fitter.extract_sub (*listtouse, 
            *selected, param, coord, paramslt,
            coordslt);
    
        if (!build_and_compare (paramslt, coordslt, cmtx, q, verbose, 
              *selected, fitter))
          return EXIT_FAILURE;
      }
    }
  }
 
  return EXIT_SUCCESS;
}
