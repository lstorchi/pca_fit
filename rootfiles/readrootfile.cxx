#include <iostream>
#include <string>
#include <cassert>
#include <vector>

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

  std::vector<int> layerid, * p_layerid, moduleid, * p_moduleid, 
    ladderid, * p_ladderid, tp, * p_tp;
  p_layerid = &layerid;
  p_ladderid = &ladderid;
  p_moduleid = &moduleid;
  p_tp = &tp;

  TT->SetBranchAddress("L1TkSTUB_layer", &p_layerid);
  TT->SetBranchAddress("L1TkSTUB_ladder", &p_ladderid);
  TT->SetBranchAddress("L1TkSTUB_module", &p_moduleid);
  TT->SetBranchAddress("L1TkSTUB_tp", &p_tp);

  std::vector<float> stubx, * p_stubx, stuby, * p_stuby, stubz, * p_stubz,
    px, * p_px, py, * p_py, x0, * p_x0, y0, * p_y0, z0, * p_z0, eta, * p_eta,
    phi, * p_phi;
  p_stubx = &stubx;
  p_stuby = &stuby;
  p_stubz = &stubz;
  p_px = &px;
  p_py = &py;
  p_x0 = &x0;
  p_y0 = &y0;
  p_z0 = &z0;
  p_eta = &eta;
  p_phi = &phi;

  TT->SetBranchAddress("L1TkSTUB_x", &p_stubx);
  TT->SetBranchAddress("L1TkSTUB_y", &p_stuby);
  TT->SetBranchAddress("L1TkSTUB_z", &p_stubz);

  TT->SetBranchAddress("L1TkSTUB_pxGEN", &p_px);
  TT->SetBranchAddress("L1TkSTUB_pyGEN", &p_py);
  TT->SetBranchAddress("L1TkSTUB_X0", &p_x0);
  TT->SetBranchAddress("L1TkSTUB_Y0", &p_y0);
  TT->SetBranchAddress("L1TkSTUB_Z0", &p_z0);
  TT->SetBranchAddress("L1TkSTUB_etaGEN", &p_eta);
  TT->SetBranchAddress("L1TkSTUB_PHI0", &p_phi);


  unsigned int countevt = 0;
  Int_t nevent = t1->GetEntries(); 
  std::cout << "We got " << nevent << " events " << std::endl;
  for (Int_t i=0; i<nevent; ++i) 
  { 
     //L1TkSTUB_tp tutti gli stub appartenti alla stessa traccia hanno stesso tp 
     //  posso avere eventi con molte tracce devo solezionarle usando tp
     //
     //StubExtractor e' utile punto di partenza visto che contiene il codice 
     //   che scrive il Tree L1Tk

     t1->GetEvent(i);
     assert (layerid.size() == ladderid.size());
     assert (layerid.size() == moduleid.size());
     assert (layerid.size() == tp.size());
     
     assert (layerid.size() == stubx.size());
     assert (layerid.size() == stuby.size());
     assert (layerid.size() == stubz.size());
     assert (layerid.size() == px.size());
     assert (layerid.size() == py.size());
     assert (layerid.size() == x0.size());
     assert (layerid.size() == y0.size());
     assert (layerid.size() == z0.size());
     assert (layerid.size() == eta.size());
     assert (layerid.size() == phi.size());

     /*
      *
      * pT = sqrt(pxGEN^2 + pyGEN^2)
      * Phi = PHI0 
      * d0 = sqrt(X0^2 + Y0^2)
      * Eta = etaGEN
      * z0 = Z0
      *
      */

     if (layerid.size() == 6)
     {
       std::cout << "Event: " << i+1 << " size " << tp.size() << std::endl;
       
       for (int j=0; j<(int)layerid.size(); ++j)
       {
         std::cout << stubx[j] << " " << stuby[j] << " " <<
           stubz[j];
         std::cout << layerid[j] << " " << ladderid[j] << " " << 
           moduleid[j] << " "  << tp[j] << std::endl;
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
