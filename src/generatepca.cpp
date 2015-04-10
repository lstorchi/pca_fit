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

#define MINDIMLINIT 25

// lstorchi: basic quicj code to generate PCA constants

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
  std::cerr << "usage: " << name << " [options] coordinatesfile " << std::endl;
  std::cerr << std::endl;
  std::cerr << " -h, --help                 : display this help and exit" << std::endl;
  std::cerr << " -V, --verbose              : verbose option on" << std::endl;
  std::cerr << " -v, --version              : print version and exit" << std::endl;
  std::cerr << " -f, --fast                 : do not perfomr pca only diag matrix" << std::endl;
  std::cerr << " -i, --seg-id               : use SegId instead of coordinates" << std::endl;
  std::cerr << " -s, --bigger-sub-tower     : use values of the bigger sub-tower" << std::endl;
  std::cerr << "                              (connot be used with bigger-subladder)" << std::endl;
  std::cerr << " -l, --bigger-sub-sub-tower : use values of the bigger sub-sub-tower " << std::endl;
  std::cerr << "                              (connot be used with bigger-subsector)" << std::endl;
  std::cerr << " -a, --all-sub-tower        : generate the constants for all subsectors" << std::endl;
  std::cerr << " -r, --all-sub-sub-tower    : generate the constants for all subladders" << std::endl;
  std::cerr << " -j, --jump-tracks          : generate the constants using only even tracks" << std::endl;
  std::cerr << " -p, --dump-allcoords       : dump all stub coordinates to a file" << std::endl;

  exit(1);
}

void perform_main_computation (const bool fast, const bool verbose, 
    const arma::mat & coord, const arma::mat & param, 
    const std::string & cfname, const std::string & qfname,
    pcafitter & fitter)
{
  // ordered 
  arma::vec eigval;
  // by row or by column ?
  arma::mat eigvec;

  if (!fast)
  {
    std::cout << "Perform PCA " << std::endl;
    // projection 
    arma::mat score;
  
    arma::princomp(eigvec, score, eigval, coord);
    std::cout << score.n_rows << " " << score.n_cols << std::endl;
    
    std::cout << "Writing Scoreplot" << std::endl;
    std::ofstream myfilesc("scoreplot.txt");
    
    for (int i=0; i<(int)score.n_rows; ++i)
    { 
      myfilesc << score(i,0)  << " "
               << score(i,1) << " " << score(i,2) << std::endl;
    
      if (verbose)
      {
        double mainr = 0.0e0;
        for (int j=1; j<5; ++j)
          mainr += score(i,j) * score(i,j);
        
        double residue = 0.0;
        for (int j=5; j<fitter.get_coordim(); ++j)
          residue += score(i,j) * score(i,j);
        
        std::cout << "Track " << i+1 << " residue " << residue <<
                     " mainr " << mainr << std::endl;
      }
    }
    
    myfilesc.close();
  }
  else
  {
    std::cout << "Compute correlation mtx" << std::endl;
    arma::mat coordm = arma::zeros<arma::mat>(fitter.get_coordim());
    arma::mat hca = arma::zeros<arma::mat>(fitter.get_coordim(),
        fitter.get_coordim());
    arma::vec eigvaltmp = arma::zeros<arma::vec>(fitter.get_coordim());

    eigvec = arma::zeros<arma::mat>(fitter.get_coordim(),
        fitter.get_coordim());
    eigval = arma::zeros<arma::vec>(fitter.get_coordim());

    hca = arma::cov(coord);

    /* double check 
    double sum = 1.0e0;
    for (int l=0; l<(int)coord.n_rows; ++l) 
    {
      sum += 1.0e0;
      for (int i=0; i<fitter.get_coordim(); ++i)
        coordm(i) += (coord(l,i)-coordm(i))/sum;
    
      for (int i=0; i<fitter.get_coordim(); ++i)
      {
        for (int j=0; j<fitter.get_coordim(); ++j)
        {
          hca(i,j) += ((coord(l,i) - coordm(i))*
                       (coord(l,j) - coordm(j))-
                       (sum-1.0e0)*hca(i,j)/sum)/(sum-1.0e0);
        }
      }
    }
    */

    std::cout << "Eigensystem" << std::endl;
    arma::eig_sym(eigvaltmp, eigvec, hca);
 
    for (int i=0; i<fitter.get_coordim(); ++i)
      eigval(i) = eigvaltmp(fitter.get_coordim()-i-1);
  }


  double totval = 0.0e0;
  for (int i=0; i<fitter.get_coordim(); ++i)
    totval += eigval(i);

  std::cout << "Eigenvalues: " << std::endl;
  double totvar = 0.0e0; 
  for (int i=0; i<fitter.get_coordim(); ++i)
  {
    if (i < fitter.get_paramdim())
      totvar += 100.0e0*(eigval(i)/totval);

    std::cout << i+1 << " ==> " << 100.0e0*(eigval(i)/totval) 
              << "% value: " << eigval(i) <<  std::endl;
  }
  std::cout << "PARAMDIM eigenvalues: " << totvar << std::endl;

  std::cout << fitter.get_paramdim() << " X " << fitter.get_coordim() << std::endl;

  arma::mat cmtx = arma::zeros<arma::mat>(fitter.get_paramdim(),
      fitter.get_coordim());
  arma::rowvec q = arma::zeros<arma::rowvec>(fitter.get_paramdim());

  std::cout << "Compute PCA constants " << std::endl;
  fitter.compute_pca_constants (param,
      coord, cmtx, q);

  std::cout << "Write constant to file" << std::endl;
  fitter.write_armmat(cfname.c_str(), cmtx);
  fitter.write_armvct(qfname.c_str(), q);
}

int main (int argc, char ** argv)
{
  pcafitter fitter; 

  bool fast = false;
  bool verbose = false;
  bool selectsubsector = true;
  bool selectsubladder = false;
  bool useallsubsectors = false;
  bool useallsubladders = false;
  bool useonlyeven = false;
  bool usesegid = false;
  bool printallcoords = false;

  while (1)
  {
    int c, option_index;
    static struct option long_options[] = {
      {"help", 0, NULL, 'h'},
      {"verbose", 0, NULL, 'V'},
      {"bigger-sub-tower", 0, NULL, 's'},
      {"bigger-sub-sub-tower", 0, NULL, 'l'},
      {"fast", 0, NULL, 'f'},
      {"version", 0, NULL, 'v'},
      {"all-sub-tower", 0, NULL, 'a'},
      {"all-sub-sub-tower", 0, NULL, 'r'},
      {"seg-id", 0, NULL, 'i'},
      {"jump-tracks", 0, NULL, 'j'},
      {"dump-allcoords", 0, NULL, 'p'},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "iravhVslfjp", long_options, &option_index);

    if (c == -1)
      break;

    switch (c)
    {
      case 'p':
        printallcoords = true;
        break;
      case 'j':
        useonlyeven = true;
        break;
      case 'i':
        usesegid = true;
        break;
      case 'h':
        usage (argv[0]);
        break;
      case 'v':
        std::cout << "Version: " << pcafitter::get_version_string() << std::endl;
        exit(1);
        break;
      case 'a':
        useallsubsectors = true;
        break;
      case 'r':
        useallsubladders = true;
        break;
      case 'f':
        fast = true;
        break;
      case 'V':
        verbose = true;
        break;
      case's':
        selectsubsector = true;
        selectsubladder = false;
        break;
      case 'l':
        selectsubladder = true;
        selectsubsector = false;
        break;
      default:
        usage (argv[0]);
        break;
    } 
  }

  if (usesegid)
    fitter.set_coordim (6);

  if (optind >= argc) 
    usage (argv[0]);

  if (selectsubladder && selectsubsector)
    usage (argv[0]);

  if (useallsubsectors || useallsubladders)
  {
    selectsubladder = false;
    selectsubsector = false;
  }

  if (useallsubladders && useallsubsectors)
    usage (argv[0]);

  char * filename = (char *) alloca (strlen(argv[optind]) + 1);
  strcpy (filename, argv[optind]);

  // leggere file coordinate tracce simulate plus parametri
  if (!file_exists(filename))
  {
    std::cerr << "Inout file does not exist" << std::endl;
    return 1;
  }
                  
  int num_of_line = pcafitter::numofline(filename);
  std::cout << "file has " << num_of_line << " line " << std::endl;
  int num_of_ent_read = (num_of_line-1)/ENTDIM;

  int num_of_ent = num_of_ent_read;

  if (useonlyeven)
  {
    if (num_of_ent_read % 2)
      num_of_ent = (num_of_ent_read-1)/2;
    else
      num_of_ent = num_of_ent_read/2;
  }

  std::cout << "file has " << num_of_ent << " entries " << std::endl;

  // non perfomante ma easy to go
  arma::mat layer, ladder, module, coordin, paramin;
  layer.set_size(num_of_ent,fitter.get_coordim());
  ladder.set_size(num_of_ent,fitter.get_coordim());
  module.set_size(num_of_ent,fitter.get_coordim());
  coordin.set_size(num_of_ent,fitter.get_coordim());
  paramin.set_size(num_of_ent,fitter.get_paramdim());

  std::map<std::string, int> subsectors, subladders;
  std::vector<std::string> subladderslist, subsectorslist;

  // leggere file coordinate tracce simulate plus parametri
  std::cout << "Reading data from " << filename << " file " << std::endl;
  fitter.reading_from_file (filename, paramin, coordin, 
      layer, ladder, module, subsectors, subladders, 
      subsectorslist, subladderslist, num_of_ent_read, usesegid,
      useonlyeven, false);
  // write date to file 
 
  if (printallcoords)
  {
    std::cout << "Printout coordinates " << std::endl;
    std::ofstream myfilect("allcoords.txt");
    for (int i=0; i<(int)coordin.n_rows; ++i)
      for (int j=0; j<fitter.get_coordim(); j=j+3)
        myfilect << coordin(i, j) << " " << 
                    coordin(i, j+1) << " " <<
                    coordin(i, j+2) << std::endl;
    myfilect.close();
  }

  std::cout << "Writing parameters to files" << std::endl;
  fitter.write_to_file("oneoverpt.txt", paramin, PTIDX);
  fitter.write_to_file("phi.txt", paramin, PHIIDX);
  fitter.write_to_file("d0.txt", paramin, D0IDX);
  fitter.write_to_file("cotetha2.txt", paramin, TETHAIDX);
  fitter.write_to_file("z0.txt", paramin, Z0IDX);

  if (!useallsubsectors && !useallsubladders)
  {
    // values selection
    std::string slctsubsec = "";
    int maxnumber = 0;
  
    arma::mat coord, param;

    assert(selectsubsector != selectsubladder);
    
    if (selectsubladder || selectsubsector)
    {
      if (selectsubsector)
      {
        std::cout << "Looking for bigger subsector" << std::endl; 
        std::cout << "We  found " << subsectors.size() << 
          " subsectors " << std::endl;
      
        fitter.select_bigger_sub (subsectors, verbose, 
            maxnumber, slctsubsec);
        
        std::cout << "Selected subsector " << slctsubsec << 
          " numevt: " << maxnumber << std::endl;
    
        fitter.extract_sub (subsectorslist, 
            slctsubsec, paramin, coordin, param, 
            coord);
      }
      
      if (selectsubladder)
      {
        std::cout << "Looking for bigger subladder" << std::endl; 
        std::cout << "We  found " << subladders.size() << 
          " subladders " << std::endl;
      
        fitter.select_bigger_sub (subladders, verbose, 
            maxnumber, slctsubsec);
        
        std::cout << "Selected subladder " << slctsubsec << " numevt: " << maxnumber << std::endl;
    
        fitter.extract_sub (subladderslist, 
            slctsubsec, paramin, coordin, param, 
            coord);
      }
   
    }
    else
    {
      param = paramin;
      coord = coordin;
    }

    std::cout << "Printout selected coordinates " << std::endl;
    std::ofstream myfileslct("selectedcoords.txt");
    for (int i=0; i<(int)coord.n_rows; ++i)
      for (int j=0; j<fitter.get_coordim(); j=j+3)
        myfileslct << coord(i, j) << " " << 
                      coord(i, j+1) << " " <<
                      coord(i, j+2) << std::endl;
    myfileslct.close();
    
    // write date to file 
    std::cout << "Writing extracted parameters to files" << std::endl;
    fitter.write_to_file("oneoverpt_selected.txt", param, PTIDX);
    fitter.write_to_file("phi_selected.txt", param, PHIIDX);
    fitter.write_to_file("d0_selected.txt", param, D0IDX);
    fitter.write_to_file("cotetha2_selected.txt", param, TETHAIDX);
    fitter.write_to_file("z0_selected.txt", param, Z0IDX);

    std::ostringstream cfname, qfname; 
    cfname << "c." << slctsubsec << ".bin";
    qfname << "q." << slctsubsec << ".bin";
 
    perform_main_computation (fast, verbose, coord, param,
        cfname.str(), qfname.str(), fitter);
  }
  else 
  {
    assert(useallsubladders != useallsubsectors);

    arma::mat coord, param;
    unsigned int totalamnt = 0, usedamnt = 0;
    std::vector<std::string> * listtouse = NULL;

    std::set<std::string> sectrorset;

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
      fitter.extract_sub (*listtouse, 
          *selected, paramin, coordin, param, 
          coord);

      std::cout << "Selected " << *selected << " size " << 
        coord.n_rows << std::endl;

      totalamnt += coord.n_rows;

      if (coord.n_rows > MINDIMLINIT)
      {
        std::ostringstream cfname, qfname; 
        cfname << "c." << *selected << ".bin";
        qfname << "q." << *selected << ".bin";
        perform_main_computation (fast, verbose, coord, param,
            cfname.str(), qfname.str(), fitter);

        usedamnt += coord.n_rows;
      }
    }

    std::cout << "Total tracks: " << totalamnt << std::endl;
    std::cout << "Used tracks: " << usedamnt << std::endl;
  }

  return 0;
}
