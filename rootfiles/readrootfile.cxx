#include <iostream>
#include <string>
#include <cassert>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h" 
#include "TBasket.h"

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
  TChain* TT = (TChain*) inputFile->Get("TkStubs");
  t1->Print();
  std::cout << std::endl;

  /*
  TT->SetBranchAddress("STUB_n", &m_stub);
  TT->SetBranchAddress("STUB_modid", &p_m_stub_modid);
  TT->SetBranchAddress("STUB_ptGEN", &p_m_stub_ptGEN);
  TT->SetBranchAddress("STUB_phiGEN", &p_m_stub_phiGEN);
  TT->SetBranchAddress("STUB_etaGEN", &p_m_stub_etaGEN);
  TT->SetBranchAddress("STUB_x", &p_m_stub_x);
  TT->SetBranchAddress("STUB_y", &p_m_stub_y);
  TT->SetBranchAddress("STUB_z", &p_m_stub_z);
  TT->SetBranchAddress("STUB_pdgGEN", &p_m_stub_pdg);
  TT->SetBranchAddress("STUB_xGEN", &p_m_stub_x0);
  TT->SetBranchAddress("STUB_yGEN", &p_m_stub_y0);
  TT->SetBranchAddress("STUB_zGEN", &p_m_stub_z0);
  */

  std::vector<int> layerid, * p_layerid, moduleid, * p_moduleid, 
    ladderid, * p_ladderid;
  p_layerid = &layerid;
  p_ladderid = &ladderid;
  p_moduleid = &moduleid;
  TT->SetBranchAddress("L1TkSTUB_layer", &p_layerid);
  TT->SetBranchAddress("L1TkSTUB_ladder", &p_ladderid);
  TT->SetBranchAddress("L1TkSTUB_module", &p_moduleid);

  unsigned int countevt = 0;
  Int_t nevent = t1->GetEntries(); 
  std::cout << "We got " << nevent << " events " << std::endl;
  for (Int_t i=0; i<nevent; ++i) 
  { 
     t1->GetEvent(i);
     assert (layerid.size() == ladderid.size());
     assert (layerid.size() == moduleid.size());

     if (layerid.size() == 6)
     {
       std::cout << "Event: " << i+1 << std::endl;
       
       for (int j=0; j<(int)layerid.size(); ++j)
       {
         std::cout << layerid[j] << " " << ladderid[j] << " " << 
           moduleid[j] << std::endl;
       }

       countevt++;
     }

     //t1->Show(i);
  }

  std::cout << "Event with 6 layer " << countevt << std::endl;

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
