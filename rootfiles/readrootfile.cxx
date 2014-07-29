#include <iostream>
#include <string>

#include "TFile.h"
#include "TNtuple.h"

void readandtest (const std::string & fname)
{
  TFile* inputFile = new TFile(fname.c_str(),"READ");
  
  std::cout << "Print file info: " << std::endl;
  inputFile->Print();
  std::cout << std::endl;

  std::cout << "List file TTree: " << std::endl;
  inputFile->ls();
  std::cout << std::endl;

  std::cout << "Print TkStubs info: " << std::endl;
  TTree* t1 = (TTree*) inputFile->Get("TkStubs");
  t1->Print();
  std::cout << std::endl;

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
