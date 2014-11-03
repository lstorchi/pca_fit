#include <algorithm>
#include <iostream>
#include <cassert>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h" 
#include "TBasket.h"

void print_bankstub (TFile * inputFile)
{
  TChain* TT = (TChain*) inputFile->Get("BankStubs");

  std::vector<int> moduleid, * p_moduleid; 
  p_moduleid = &moduleid;

  TT->SetBranchAddress("STUB_modid", &p_moduleid); // QA come determino layerid e altro ? 
                                                   //    devo caricare la geometria ?

  std::vector<float> stubx, * p_stubx, stuby, * p_stuby, stubz, * p_stubz,
    pt, * p_pt, x0, * p_x0, y0, * p_y0, z0, * p_z0, eta, * p_eta,
    phi, * p_phi;
  p_stubx = &stubx;
  p_stuby = &stuby;
  p_stubz = &stubz;
  p_pt = &pt;
  p_x0 = &x0;
  p_y0 = &y0;
  p_z0 = &z0;
  p_eta = &eta;
  p_phi = &phi;

  TT->SetBranchAddress("STUB_x", &p_stubx);
  TT->SetBranchAddress("STUB_y", &p_stuby);
  TT->SetBranchAddress("STUB_z", &p_stubz);

  TT->SetBranchAddress("STUB_ptGEN", &p_pt);
  TT->SetBranchAddress("STUB_xGEN", &p_x0);
  TT->SetBranchAddress("STUB_yGEN", &p_y0);
  TT->SetBranchAddress("STUB_zGEN", &p_z0);
  TT->SetBranchAddress("STUB_etaGEN", &p_eta);
  TT->SetBranchAddress("STUB_phiGEN", &p_phi);

  unsigned int countevt = 0;
  Int_t nevent = TT->GetEntries(); 
  std::cout << "We got " << nevent << " events in BankStubs" << std::endl; 
  // QA perche' il numero di eventi qui e' molto maggiore ?

  std::ofstream ptfile("pt_BankStubs.txt");
  std::ofstream phifile("phi_BankStubs.txt");
  std::ofstream d0file("d0_BankStubs.txt");
  std::ofstream etafile("eta_BankStubs.txt");
  std::ofstream z0file("z0_BankStubs.txt");
  for (Int_t i=0; i<nevent; ++i) 
  { 
     TT->GetEntry(i);
     
     assert (moduleid.size() == stubx.size());
     assert (moduleid.size() == stuby.size());
     assert (moduleid.size() == stubz.size());
     assert (moduleid.size() == pt.size());
     assert (moduleid.size() == x0.size());
     assert (moduleid.size() == y0.size());
     assert (moduleid.size() == z0.size());
     assert (moduleid.size() == eta.size());
     assert (moduleid.size() == phi.size());

     bool allAreEqual = ((std::find_if(z0.begin() + 1, z0.end(), 
        std::bind1st(std::not_equal_to<int>(), z0.front())) == z0.end()) &&
                        (std::find_if(x0.begin() + 1, x0.end(), 
        std::bind1st(std::not_equal_to<int>(), x0.front())) == x0.end()) &&
                        (std::find_if(y0.begin() + 1, y0.end(), 
        std::bind1st(std::not_equal_to<int>(), y0.front())) == y0.end()) &&
                        (std::find_if(pt.begin() + 1, pt.end(), 
        std::bind1st(std::not_equal_to<int>(), pt.front())) == pt.end()) &&
                        (std::find_if(eta.begin() + 1, eta.end(), 
        std::bind1st(std::not_equal_to<int>(), eta.front())) == eta.end()) &&
                        (std::find_if(phi.begin() + 1, phi.end(), 
        std::bind1st(std::not_equal_to<int>(), phi.front())) == phi.end()));


     if ((moduleid.size() == 6)  && allAreEqual) // QA nel caso dei BankStubs questo check e' utile ?
     {
       ptfile << pt[0] << std::endl;
       phifile << phi[0] << std::endl;
       d0file << sqrt(pow(x0[0],2.0) + pow(y0[0],2.0)) 
         << std::endl;
       etafile << eta[0] << std::endl;
       z0file << z0[0] << std::endl;

       std::cerr << i+1 << " " << moduleid.size() << std::endl;

       int j = 0;
       for (; j<(int)moduleid.size(); ++j)
       {
        std::cerr << stubx[j] << " " << stuby[j] << " " <<
           stubz[j] << " ";
        std::cerr << moduleid[j] << " " << moduleid[j] << " " << 
          moduleid[j] << " " << std::endl;
       }
       --j;

       std::cerr << pt[j]<< " "  <<
         phi[j] << " " << sqrt(pow(x0[j],2.0) + pow(y0[j],2.0)) << " " 
         << eta[j] << " " << z0[j] << std::endl;
     }
  }
  ptfile.close();
  phifile.close();
  d0file.close();
  etafile.close();
  z0file.close();
}

void readandtest (const std::string & fname)
{
  TFile* inputFile = new TFile(fname.c_str(),"READ");

#if 0
  std::cout << "Print file info: " << std::endl;
  inputFile->Print();
  std::cout << std::endl;

  std::cout << "List file TTree: " << std::endl;
  inputFile->ls();
  std::cout << std::endl;

  std::cout << "Print BankStubs info: " << std::endl;
  TTree* tbs = (TTree*) inputFile->Get("BankStubs");
  tbs->Print(); 
  std::cout << std::endl;

  std::cout << "Print Pixels info: " << std::endl;
  TTree* t0 = (TTree*) inputFile->Get("Pixels");
  t0->Print(); 
  std::cout << std::endl;

  std::cout << "Print MC info: " << std::endl;
  TTree* t01 = (TTree*) inputFile->Get("MC");
  t01->Print(); 
  std::cout << std::endl;

  std::cout << "Print TkStubs info: " << std::endl;
  TTree* t1 = (TTree*) inputFile->Get("TkStubs");
  t1->Print();
  std::cout << std::endl;
#endif

  TChain* TT = (TChain*) inputFile->Get("TkStubs");

  TTree* tmc = (TTree*) inputFile->Get("MC");
  //tmc->Print();
  TChain* TTmc = (TChain*) inputFile->Get("MC");

  std::vector<float> genpx, * p_genpx;
  p_genpx = &genpx;
  TTmc->SetBranchAddress("gen_px", &p_genpx);
  
  Int_t neventmc = tmc->GetEntries();
  for (Int_t i=0; i<neventmc; ++i) 
  { 
     tmc->GetEvent(i);

     int j = 0;
     for (; j<(int)genpx.size(); ++j)
     {
       //std::cout << genpx[j] << std::endl;
     }
  }

  print_bankstub (inputFile);
  // QA le distribuzioni dei vari parametri sono simili ma 
  //    ancora in bankstub ci sono tantssime tracce in piu'.
  //    ed il module id ha un valore numeri che non "capisco" immagino
  //    sia una specie di identificativo univoco del modulo

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
  Int_t nevent = TT->GetEntries(); 
  std::cout << "We got " << nevent << " events " << std::endl;
  std::ofstream ptfile("pt_L1TkSTUB.txt");
  std::ofstream phifile("phi_L1TkSTUB.txt");
  std::ofstream d0file("d0_L1TkSTUB.txt");
  std::ofstream etafile("eta_L1TkSTUB.txt");
  std::ofstream z0file("z0_L1TkSTUB.txt");
  for (Int_t i=0; i<nevent; ++i) 
  { 
     //L1TkSTUB_tp tutti gli stub appartenti alla stessa traccia hanno stesso tp 
     //  posso avere eventi con molte tracce devo solezionarle usando tp
     //
     //StubExtractor e' utile punto di partenza visto che contiene il codice 
     //   che scrive il Tree L1Tk

     TT->GetEntry(i);
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
       bool allAreEqual = ((std::find_if(z0.begin() + 1, z0.end(), 
          std::bind1st(std::not_equal_to<int>(), z0.front())) == z0.end()) &&
                          (std::find_if(x0.begin() + 1, x0.end(), 
          std::bind1st(std::not_equal_to<int>(), x0.front())) == x0.end()) &&
                          (std::find_if(y0.begin() + 1, y0.end(), 
          std::bind1st(std::not_equal_to<int>(), y0.front())) == y0.end()) &&
                          (std::find_if(px.begin() + 1, px.end(), 
          std::bind1st(std::not_equal_to<int>(), px.front())) == px.end()) &&
                          (std::find_if(py.begin() + 1, py.end(), 
          std::bind1st(std::not_equal_to<int>(), py.front())) == py.end()) &&
                          (std::find_if(eta.begin() + 1, eta.end(), 
          std::bind1st(std::not_equal_to<int>(), eta.front())) == eta.end()) &&
                          (std::find_if(phi.begin() + 1, phi.end(), 
          std::bind1st(std::not_equal_to<int>(), phi.front())) == phi.end()));

       //if ((z0[0] <= 15.0) && (z0[0] >= -15) && allAreEqual && (tp[0] == 0))
       if ((tp[0] == 0) && allAreEqual) // QA perche' devo fare un controllo anche su allAreEqual 
                                        // non dovrebbe bazter il controllo su tp ? 
                                        // tp non indica le traccie primarie ? 
       //if (tp[0] == 0)
       {
         std::cout << i+1 << " " << tp.size() << std::endl;

         int j = 0;
         for (; j<(int)layerid.size(); ++j)
         {
          std::cout << stubx[j] << " " << stuby[j] << " " <<
             stubz[j] << " ";
           
#if 0     
           std::cout << sqrt(pow(px[j],2.0) + pow(py[j],2.0)) << " "  <<
             phi[j] << " " << sqrt(pow(x0[j],2.0) + pow(y0[j],2.0)) << " " 
             << eta[j] << " " << z0[j] << " ";
#endif   
          
           std::cout << layerid[j] << " " << ladderid[j] << " " << 
             moduleid[j] << " "
#if 0             
           << tp[j] << std::endl;
#else     
           << std::endl;
#endif    
          
         }
         --j;

         ptfile << sqrt(pow(px[j],2.0) + pow(py[j],2.0))  
           << std::endl;
         phifile << phi[j] << std::endl;
         d0file << sqrt(pow(x0[j],2.0) + pow(y0[j],2.0)) 
           << std::endl;
         etafile << eta[j] << std::endl;
         z0file << z0[j] << std::endl;

         std::cout << sqrt(pow(px[j],2.0) + pow(py[j],2.0)) << " "  <<
           phi[j] << " " << sqrt(pow(x0[j],2.0) + pow(y0[j],2.0)) << " " 
           << eta[j] << " " << z0[j] << std::endl;

         countevt++;
       }
     }

     //t1->Show(i);
  }
  ptfile.close();
  phifile.close();
  d0file.close();
  etafile.close();
  z0file.close();

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
