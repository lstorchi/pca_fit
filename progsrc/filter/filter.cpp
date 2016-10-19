#include <iostream>
#include <fstream>
#include <iomanip>

// Internal includes

#include "filter.h"
#include "jobparams.h"
#include "TROOT.h"

using namespace std;

///////////////////////////////////
//
//
// Base code for the AM filtering tool
//
//
///////////////////////////////////

int main(int argc, char** argv) {

  // Necessary lines to make branches containing vectors
  gROOT->ProcessLine(".L Loader.C+");

  // read jobParams
  jobparams params(argc,argv);

  // Depending on the option chosen, process the information

  // Option 1: basic filtering stage for banks

  //
  // How to call it :
  //
  // ./AM_filter -c filter -i INPUT -n SECNUM -l HITLIM -s SECFILE -o OUTPUT
  //
  // where:
  //
  // INPUT  : directory and name of the input root file to filter. If you have more than one file, 
  //          the code can also process a txt file containing the adresses of all the files
  //
  // SECNUM : the sector number (between 0 and 47 for 8x6)
  //
  // HITLIM : the minimum number of layer/disks with stubs in the sector (never use less than 4 here)
  //
  // SECFILE: directory and name of the csv file defining the sectors
  //
  // OUTPUT : name of the ouput root file you want to create
  //


  if (params.option()=="filter")
  {
    filter* my_test = new filter(params.inputfile(),params.sector(),
				 params.outfile(),params.nsec(),params.lim());
    delete my_test;
  }
  
  return 0;
}
