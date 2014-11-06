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
  std::cerr << " -f, --fast               : do not perfomr pca only diag matrix" << std::endl;
  std::cerr << " -s, --bigger-subsector   : use values of the bigger subsector" << std::endl;
  std::cerr << "                            (connot be used with bigger-subladder)" << std::endl;
  std::cerr << " -l, --bigger-subladder   : use values of the bigger subladder " << std::endl;
  std::cerr << "                            (connot be used with bigger-subsector)" << std::endl;

  exit(1);
}

int main (int argc, char ** argv)
{
  bool fast = false;
  bool verbose = false;
  bool selectsubsector = true;
  bool selectsubladder = false;

  while (1)
  {
    int c, option_index;
    static struct option long_options[] = {
      {"help", 0, NULL, 'h'},
      {"verbose", 0, NULL, 'V'},
      {"bigger-subsector", 0, NULL, 's'},
      {"bigger-subladder", 0, NULL, 'l'},
      {"fast", 0, NULL, 'f'},
      {"version", 0, NULL, 'v'},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "vhVslf", long_options, &option_index);

    if (c == -1)
      break;

    switch (c)
    {
      case 'h':
        usage (argv[0]);
        break;
      case 'v':
        std::cout << "Version: " << pcafitter::get_version_string() << std::endl;
        exit(1);
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

  if (optind >= argc) 
    usage (argv[0]);

  if (selectsubladder && selectsubsector)
    usage (argv[0]);

  char * filename = (char *) alloca (strlen(argv[optind]) + 1);
  strcpy (filename, argv[optind]);
                  
  int num_of_line = pcafitter::numofline(filename);
  std::cout << "file has " << num_of_line << " line " << std::endl;
  int num_of_ent = (num_of_line-1)/ENTDIM;
  std::cout << "file has " << num_of_ent << " entries " << std::endl;

  // non perfomante ma easy to go
  arma::mat layer, ladder, module, coordin, paramin;
  layer.set_size(num_of_ent,COORDIM);
  ladder.set_size(num_of_ent,COORDIM);
  module.set_size(num_of_ent,COORDIM);
  coordin.set_size(num_of_ent,DIMPERCOORD*COORDIM);
  paramin.set_size(num_of_ent,PARAMDIM);

  std::map<std::string, int> subsectors, subladders;
  std::vector<std::string> subladderslist, subsectorslist;

  // leggere file coordinate tracce simulate plus parametri
  std::cout << "Reading data from " << filename << " file " << std::endl;
  pcafitter::readingfromfile (filename, paramin, coordin, 
      layer, ladder, module, subsectors, subladders, 
      subsectorslist, subladderslist, num_of_ent);
  // write date to file 
  std::cout << "Writing parameters to files" << std::endl;
  pcafitter::writetofile("oneoverpt.txt", paramin, PTIDX);
  pcafitter::writetofile("phi.txt", paramin, PHIIDX);
  pcafitter::writetofile("d0.txt", paramin, D0IDX);
  pcafitter::writetofile("cotetha.txt", paramin, TETHAIDX);
  pcafitter::writetofile("z0.txt", paramin, Z0IDX);

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
    
      pcafitter::select_bigger_sub (subsectors, verbose, 
          maxnumber, slctsubsec);
      
      std::cout << "Selected subsector " << slctsubsec << 
        " numevt: " << maxnumber << std::endl;

      pcafitter::extract_sub (subsectorslist, 
          slctsubsec, paramin, coordin, param, 
          coord);
    }
    
    if (selectsubladder)
    {
      std::cout << "Looking for bigger subladder" << std::endl; 
      std::cout << "We  found " << subladders.size() << 
        " subladders " << std::endl;
    
      pcafitter::select_bigger_sub (subladders, verbose, 
          maxnumber, slctsubsec);
      
      std::cout << "Selected subladder " << slctsubsec << " numevt: " << maxnumber << std::endl;

      pcafitter::extract_sub (subladderslist, 
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
    for (int j=0; j<(DIMPERCOORD*COORDIM); j=j+3)
      myfileslct << coord(i, j) << " " << 
                    coord(i, j+1) << " " <<
                    coord(i, j+2) << std::endl;
  myfileslct.close();

  // write date to file 
  std::cout << "Writing extracted parameters to files" << std::endl;
  pcafitter::writetofile("oneoverpt_selected.txt", param, PTIDX);
  pcafitter::writetofile("phi_selected.txt", param, PHIIDX);
  pcafitter::writetofile("d0_selected.txt", param, D0IDX);
  pcafitter::writetofile("cotetha_selected.txt", param, TETHAIDX);
  pcafitter::writetofile("z0_selected.txt", param, Z0IDX);

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
        for (int j=5; j<DIMPERCOORD*COORDIM; ++j)
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
    arma::mat coordm = arma::zeros<arma::mat>(DIMPERCOORD*COORDIM);
    arma::mat hca = arma::zeros<arma::mat>(DIMPERCOORD*COORDIM,DIMPERCOORD*COORDIM);
    arma::vec eigvaltmp = arma::zeros<arma::vec>(DIMPERCOORD*COORDIM);

    eigvec = arma::zeros<arma::mat>(DIMPERCOORD*COORDIM,DIMPERCOORD*COORDIM);
    eigval = arma::zeros<arma::vec>(DIMPERCOORD*COORDIM);

    hca = arma::cov(coord);

    /* double check 
    double sum = 1.0e0;
    for (int l=0; l<(int)coord.n_rows; ++l) 
    {
      sum += 1.0e0;
      for (int i=0; i<(DIMPERCOORD*COORDIM); ++i)
        coordm(i) += (coord(l,i)-coordm(i))/sum;
    
      for (int i=0; i<(DIMPERCOORD*COORDIM); ++i)
      {
        for (int j=0; j<(DIMPERCOORD*COORDIM); ++j)
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
 
    for (int i=0; i<(DIMPERCOORD*COORDIM); ++i)
      eigval(i) = eigvaltmp((DIMPERCOORD*COORDIM)-i-1);
  }


  double totval = 0.0e0;
  for (int i=0; i<(DIMPERCOORD*COORDIM); ++i)
    totval += eigval(i);

  std::cout << "Eigenvalues: " << std::endl;
  double totvar = 0.0e0; 
  for (int i=0; i<(DIMPERCOORD*COORDIM); ++i)
  {
    if (i < PARAMDIM)
      totvar += 100.0e0*(eigval(i)/totval);

    std::cout << i+1 << " ==> " << 100.0e0*(eigval(i)/totval) 
              << "% value: " << eigval(i) <<  std::endl;
  }
  std::cout << "PARAMDIM eigenvalues: " << totvar << std::endl;

  arma::mat cmtx = arma::zeros<arma::mat>(PARAMDIM,DIMPERCOORD*COORDIM);
  arma::rowvec q = arma::zeros<arma::rowvec>(PARAMDIM);

  std::cout << "Compute PCA constants " << std::endl;
  pcafitter::compute_pca_constants (param,
      coord, cmtx, q);

  std::cout << "Write constant to file" << std::endl;
  pcafitter::writearmmat("c.bin", cmtx);
  pcafitter::writearmvct("q.bin", q);

  return 0;
}
