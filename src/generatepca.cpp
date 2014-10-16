#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <string>
#include <map>

// Loriano: let's try Armadillo quick code 
#include <armadillo>

#include <getopt.h>
#include <unistd.h>
#include <alloca.h>

#define ENTDIM 8
#define COORDIM (ENTDIM-2)
#define PARAMDIM 5
#define THRSVAL 1000

namespace 
{
  int numofline (const char * fname) 
  { 
    int number_of_lines = 0;
    std::string line;
    std::ifstream myfile(fname);
    
    while (std::getline(myfile, line))
      ++number_of_lines;

    myfile.close();
                            
    return number_of_lines;
  }

  void writetofile (const char * fname, 
      const arma::rowvec & vec) 
  {
    std::ofstream myfile(fname);
      
    for (int i=0; i<(int)vec.n_cols; i++)
      myfile << i << " " << vec(i) << std::endl;
      
    myfile.close();
  }

  void readingfromfile (const char * filename, 
      arma::rowvec & ptin, arma::rowvec & phiin, 
      arma::rowvec & d0in, arma::rowvec & etain, 
      arma::rowvec & z0in, arma::mat & coordinp, 
      arma::mat & layer, arma::mat & ladder, 
      arma::mat & module, 
      std::map<std::string, int> & subsectors, 
      std::map<std::string, int> & subladders,
      int num_of_ent)
  {

    std::string line;
    std::ifstream mytfp;
    mytfp.open (filename, std::ios::in);
  
    std::getline (mytfp, line);
    for (int i = 0; i < num_of_ent; ++i)
    {
      int fake1, fake2;
      // valori aggiunti solo di controllo 
      mytfp >> fake1 >> fake2 ;
#ifdef DEBUG    
      std::cout << fake1 << " " << fake2 << std::endl;
#endif
      std::ostringstream osss, ossl;
      osss << std::setfill('0');
      ossl << std::setfill('0');
      
      for (int j = 0; j < COORDIM; ++j)
      {
        int a, b, c;
        mytfp >> coordinp(i, j*3) >> 
                 coordinp(i, j*3+1) >> 
                 coordinp(i, j*3+2) >> 
                 a >> b >> c; 
      
        layer(i, j) = a;
        ladder(i, j) = b;
        module(i, j) = c;
        
        osss << std::setw(2) << layer(i, j);
        osss << std::setw(2) << ladder(i, j);
        if (j != COORDIM-1)
          osss<<"-";

        ossl << std::setw(2) << layer(i, j);
        ossl << std::setw(2) << ladder(i, j);
        ossl << std::setw(2) << module(i, j);
        if (j != COORDIM-1)
          ossl<<"-";
      }
      
      std::map<std::string, int>::iterator its = subsectors.find(osss.str());
      if (its == subsectors.end())
        subsectors[osss.str()] = 1;
      else 
        subsectors[osss.str()] += 1;

      std::map<std::string, int>::iterator itl = subladders.find(ossl.str());
      if (itl == subladders.end())
        subladders[ossl.str()] = 1;
      else 
        subladders[ossl.str()] += 1;
      
      mytfp >> ptin(i) >> 
               phiin(i) >> 
               d0in(i) >> 
               etain(i) >> 
               z0in(i);
    }

    mytfp.close();
  }
}

void usage (char * name)
{
  std::cerr << "usage: " << name << " [options] coordinatesfile " << std::endl;
  std::cerr << std::endl;
  std::cerr << " -h, --help               : display this help and exit" << std::endl;
  std::cerr << " -v, --verbose            : verbose option on" << std::endl;
  std::cerr << " -s, --bigger-subsector   : use values of the bigger subsector" << std::endl;
  std::cerr << "                            (connot be used with bigger-subladder)" << std::endl;
  std::cerr << " -l, --bigger-subladder   : use values of the bigger subladder " << std::endl;
  std::cerr << "                            (connot be used with bigger-subsector)" << std::endl;

  exit(1);
}

int main (int argc, char ** argv)
{
  bool verbose = false;
  bool selectsubsecor = false;
  bool selectsubladder = false;

  while (1)
  {
    int c, option_index;
    static struct option long_options[] = {
      {"help", 0, NULL, 'h'},
      {"verbose", 0, NULL, 'V'},
      {"bigger-subsector", 0, NULL, 's'},
      {"bigger-subladder", 0, NULL, 'l'},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "hVsl", long_options, &option_index);

    if (c == -1)
      break;

    switch (c)
    {
      case 'h':
        usage (argv[0]);
        break;
      case 'v':
        verbose = true;
        break;
      case's':
        selectsubsecor = true;
        if (selectsubladder)
          usage (argv[0]);
        break;
      case 'l':
        selectsubladder = true;
        if (selectsubsecor)
          usage (argv[0]);
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
                  
  int num_of_line = numofline(filename);
  std::cout << "file has " << num_of_line << " line " << std::endl;
  int num_of_ent = (num_of_line-1)/ENTDIM;
  std::cout << "file has " << num_of_ent << " entries " << std::endl;

  // non perfomante ma easy to go
  arma::rowvec ptin = arma::zeros<arma::rowvec>(num_of_ent);
  arma::rowvec phiin = arma::zeros<arma::rowvec>(num_of_ent);
  arma::rowvec d0in = arma::zeros<arma::rowvec>(num_of_ent);
  arma::rowvec etain = arma::zeros<arma::rowvec>(num_of_ent);
  arma::rowvec z0in = arma::zeros<arma::rowvec>(num_of_ent);

  arma::mat layer, ladder, module, coordinp;
  layer.set_size(num_of_ent,COORDIM);
  ladder.set_size(num_of_ent,COORDIM);
  module.set_size(num_of_ent,COORDIM);
  coordinp.set_size(num_of_ent,3*COORDIM);

  std::map<std::string, int> subsectors;
  std::map<std::string, int> subladders;

  // leggere file coordinate tracce simulate plus parametri
  std::cout << "Reading data from " << filename << " file " << std::endl;
  readingfromfile (filename, ptin, phiin, d0in, etain, z0in, coordinp, 
      layer, ladder, module, subsectors, subladders, num_of_ent);
 
  // write date to file 
  std::cout << "Writing parameters to files" << std::endl;
  writetofile("pt.txt", ptin);
  writetofile("phi.txt", phiin);
  writetofile("d0.txt", d0in);
  writetofile("eta.txt", etain);
  writetofile("z0t.txt", z0in);

  // selection subsector 
  std::string slctsubsec = "";
  int maxnumber = 0;
  
  if (selectsubsecor)
  {
    std::cout << "Looking for bigger subsector" << std::endl; 
    std::cout << "We  found " << subsectors.size() << 
      " subsectors " << std::endl;
    
    std::map<std::string, int>::iterator it = subsectors.begin();
    for (; it != subsectors.end(); ++it)
    {
      if (it->second > maxnumber)
      {
        maxnumber = it->second;
        slctsubsec = it->first;
      }
    
      if (verbose)
        std::cout << "Subsector " << it->first << " has " << 
            it->second << " tracks " << std::endl;
    } 
    
    std::cout << "Selected subsector " << slctsubsec << " numevt: " << maxnumber << std::endl;
  }

  return 0;
}
