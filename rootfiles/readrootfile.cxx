#include <iostream>
#include <string>

#include "TFile.h"

void readandtest (const std::string & fname)
{
  TFile* inputFile = new TFile(fname.c_str(),"READ");
  
  inputFile->Print();

  inputFile->Close();
}

# ifndef __CINT__
int main(int argc, char ** argv) 
{
  if (argc != 2) 
  {
    std::cerr << "usage: " << argv[0] << " rootfilename " << std::endl;
    return 1;
  }
   
  readandtest(argv[1]);

  return 0;
}
# endif
