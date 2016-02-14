#include <algorithm>
#include <iostream>
#include <cassert>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#include <getopt.h>
#include <unistd.h>
#include <alloca.h>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h" 
#include "TBasket.h"

// can be included in nay case if -std=c++11
#ifdef INTBITEWISE
#include "stdint.h"
#endif

#define STOPAFTERMAXEVT INFINITY

// lstorchi: as all the code here is a very quick implementation
//    of a simil stubs extractors

void usage (char * name)
{
  std::cerr << "usage: " << name << " [options] rootfile " << std::endl;
  std::cerr << std::endl;
  std::cerr << " -h, --help               : display this help and exit" << std::endl;
  std::cerr << " -l, --l1tk-stubs         : extract l1tkstubs" << std::endl;
  std::cerr << " -m, --max-tracks[=value] : max value of tracks to be extracted" << std::endl;

  exit(1);
}

void print_bankstub_new (TFile * inputFile, std::ostream& ss, 
    std::ostream& ssext, std::ostream& ssext1, unsigned int maxtracks)
{
  TChain* TT = (TChain*) inputFile->Get("BankStubs");

  std::vector<int> moduleid, * p_moduleid; 
  p_moduleid = &moduleid;

  TT->SetBranchAddress("STUB_modid", &p_moduleid); // QA come determino layerid e altro ? 
                                                   //    devo caricare la geometria ?
  std::vector<float> stubx, * p_stubx, stuby, * p_stuby, stubz, * p_stubz,
    pt, * p_pt, x0, * p_x0, y0, * p_y0, z0, * p_z0, eta, * p_eta,
    phi, * p_phi;
  std::vector<float> pdg, * p_pdg;

  p_stubx = &stubx;
  p_stuby = &stuby;
  p_stubz = &stubz;
  p_pt = &pt;
  p_x0 = &x0;
  p_y0 = &y0;
  p_z0 = &z0;
  p_eta = &eta;
  p_phi = &phi;
  p_pdg = &pdg;

  TT->SetBranchAddress("STUB_x", &p_stubx);
  TT->SetBranchAddress("STUB_y", &p_stuby);
  TT->SetBranchAddress("STUB_z", &p_stubz);

  TT->SetBranchAddress("STUB_ptGEN", &p_pt);
  TT->SetBranchAddress("STUB_X0", &p_x0);
  TT->SetBranchAddress("STUB_Y0", &p_y0);
  TT->SetBranchAddress("STUB_Z0", &p_z0);
  TT->SetBranchAddress("STUB_etaGEN", &p_eta);
  TT->SetBranchAddress("STUB_PHI0", &p_phi);
  TT->SetBranchAddress("STUB_pdg", &p_pdg);

  unsigned int countevt = 0;
  Int_t nevent = TT->GetEntries(); 
  ss << "We got " << nevent << " events in BankStubs" << std::endl; 
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
     assert (moduleid.size() == pdg.size());

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
       double d0val;
       //d0val = (y0[0]-(tan(phi[0])*x0[0]))*cos(phi[0]);
       d0val = y0[0]*cos(phi[0])-x0[0]*sin(phi[0]);
       //double d0val = x0[0];

       ptfile << pt[0] << std::endl;
       phifile << phi[0] << std::endl;
       d0file << d0val << std::endl;
       etafile << eta[0] << std::endl;
       z0file << z0[0] << std::endl;

       ss << i+1 << " " << moduleid.size() << std::endl;

       int j = 0;
       for (; j<(int)moduleid.size(); ++j)
       {
#ifdef INTBITEWISE         
        //Can we provide these scale factors from outside
        int16_t stubX = stubx[j]*10;
        int16_t stubY = stuby[j]*10;
        int16_t stubZ = stubz[j]*10;

        ss << stubX << " " << stubY << " " <<
           stubZ << " ";
#else
        ss << stubx[j] << " " << stuby[j] << " " <<
          stubz[j] << " ";
#endif

        int value = moduleid[j];
        int layer = value/1000000;
        value = value-layer*1000000;
        int ladder = value/10000;
        value = value-ladder*10000;
        int module = value/100;
        value = value-module*100;
        int segid = value; // QA is just this ? from the source code seems so, I need to / by 10 ?

        ss << layer << " " << ladder << " " << 
          module << " " << segid << " " << pdg[j] << std::endl;
       }
       --j;

       ss << pt[j]<< " "  <<
         phi[j] << " " << d0val << " " 
         << eta[j] << " " << z0[j] << " " <<
         x0[j] << " " << y0[j] << std::endl;

       countevt++;
     }
     else if ((moduleid.size() > 6) && allAreEqual)
     {
       double d0val;
       //d0val = (y0[0]-(tan(phi[0])*x0[0]))*cos(phi[0]);
       d0val = y0[0]*cos(phi[0])-x0[0]*sin(phi[0]);
       //double d0val = x0[0];

       ptfile << pt[0] << std::endl;
       phifile << phi[0] << std::endl;
       d0file << d0val << std::endl;
       etafile << eta[0] << std::endl;
       z0file << z0[0] << std::endl;

       ss << i+1 << " " << moduleid.size() << std::endl;

       int j = 0;
       for (; j<(int)moduleid.size(); ++j)
       {
#ifdef INTBITEWISE
        int16_t stubX = stubx[j]*10;
        int16_t stubY = stuby[j]*10;
        int16_t stubZ = stubz[j]*10;

        ss << stubX << " " << stubY << " " <<
           stubZ << " ";
#else
        ss << stubx[j] << " " << stuby[j] << " " <<
           stubz[j] << " ";
#endif
        int value = moduleid[j];
        int layer = value/1000000;
        value = value-layer*1000000;
        int ladder = value/10000;
        value = value-ladder*10000;
        int module = value/100;
        value = value-module*100;
        int segid = value; // QA is just this ? from the source code seems so, I need to / by 10 ?

        ss << layer << " " << ladder << " " << 
          module << " " << segid << " " << pdg[j] << std::endl;
       }
       --j;

       ss << pt[j]<< " "  <<
         phi[j] << " " << d0val << " " 
         << eta[j] << " " << z0[j] << " " <<
         x0[j] << " " << y0[j] << std::endl;

       countevt++;
     }
     else if ((moduleid.size() < 6) && allAreEqual)
     {
       double d0val;
       //d0val = (y0[0]-(tan(phi[0])*x0[0]))*cos(phi[0]);
       d0val = y0[0]*cos(phi[0])-x0[0]*sin(phi[0]);
       //double d0val = x0[0];

       ptfile << pt[0] << std::endl;
       phifile << phi[0] << std::endl;
       d0file << d0val << std::endl;
       etafile << eta[0] << std::endl;
       z0file << z0[0] << std::endl;

       ssext1 << i+1 << " " << moduleid.size() << std::endl;

       int j = 0;
       for (; j<(int)moduleid.size(); ++j)
       {
#ifdef INTBITEWISE
        int16_t stubX = stubx[j]*10;
        int16_t stubY = stuby[j]*10;
        int16_t stubZ = stubz[j]*10;

        ssext1 << stubX << " " << stubY << " " <<
           stubZ << " ";
#else
        ssext1 << stubx[j] << " " << stuby[j] << " " <<
           stubz[j] << " ";
#endif
        int value = moduleid[j];
        int layer = value/1000000;
        value = value-layer*1000000;
        int ladder = value/10000;
        value = value-ladder*10000;
        int module = value/100;
        value = value-module*100;
        int segid = value; // QA is just this ? from the source code seems so, I need to / by 10 ?

        ssext1 << layer << " " << ladder << " " << 
          module << " " << segid << " " << pdg[j] << std::endl;
       }
       --j;

       ssext1 << pt[j]<< " "  <<
         phi[j] << " " << d0val << " " 
         << eta[j] << " " << z0[j] << " " <<
         x0[j] << " " << y0[j] << std::endl;

       countevt++;
     }
     else
     {
       double d0val;
       //d0val = (y0[0]-(tan(phi[0])*x0[0]))*cos(phi[0]);
       d0val = y0[0]*cos(phi[0])-x0[0]*sin(phi[0]);
       //double d0val = x0[0];

       ptfile << pt[0] << std::endl;
       phifile << phi[0] << std::endl;
       d0file << d0val << std::endl;
       etafile << eta[0] << std::endl;
       z0file << z0[0] << std::endl;

       ssext << i+1 << " " << moduleid.size() << std::endl;

       int j = 0;
       for (; j<(int)moduleid.size(); ++j)
       {
#ifdef INTBITEWISE
        int16_t stubX = stubx[j]*10;
        int16_t stubY = stuby[j]*10;
        int16_t stubZ = stubz[j]*10;

        ssext << stubX << " " << stubY << " " <<
           stubZ << " ";
#else
        ssext << stubx[j] << " " << stuby[j] << " " <<
           stubz[j] << " ";
#endif
        int value = moduleid[j];
        int layer = value/1000000;
        value = value-layer*1000000;
        int ladder = value/10000;
        value = value-ladder*10000;
        int module = value/100;
        value = value-module*100;
        int segid = value; // QA is just this ? from the source code seems so, I need to / by 10 ?

        ssext << layer << " " << ladder << " " << 
          module << " " << segid << " " << pdg[j] << std::endl;
       }
       --j;

       ssext << pt[j]<< " "  <<
         phi[j] << " " << d0val << " " 
         << eta[j] << " " << z0[j] << " " <<
         x0[j] << " " << y0[j] << std::endl;

       countevt++;
     }
 

     if (countevt >= maxtracks)
       break;
  }
  ptfile.close();
  phifile.close();
  d0file.close();
  etafile.close();
  z0file.close();
}


void print_l1tkstub (TFile * inputFile, std::ostream & ss, unsigned int maxtracks)
{
  TChain* TT = (TChain*) inputFile->Get("TkStubs");
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
  ss << "We got " << nevent << " events " << std::endl;
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
                                        // I should use bankstubs, so useless 
       //if (tp[0] == 0)
       {
         ss << i+1 << " " << tp.size() << std::endl;

         int j = 0;
         for (; j<(int)layerid.size(); ++j)
         {
          ss << stubx[j] << " " << stuby[j] << " " <<
             stubz[j] << " ";
           
#if 0     
           ss << sqrt(pow(px[j],2.0) + pow(py[j],2.0)) << " "  <<
             phi[j] << " " << sqrt(pow(x0[j],2.0) + pow(y0[j],2.0)) << " " 
             << eta[j] << " " << z0[j] << " ";
#endif   
          
           ss << layerid[j] << " " << ladderid[j] << " " << 
             moduleid[j] << " "
           << tp[j] << std::endl;
          
         }
         --j;

         ptfile << sqrt(pow(px[j],2.0) + pow(py[j],2.0))  
           << std::endl;
         phifile << phi[j] << std::endl;
         d0file << sqrt(pow(x0[j],2.0) + pow(y0[j],2.0)) 
           << std::endl;
         etafile << eta[j] << std::endl;
         z0file << z0[j] << std::endl;

         ss << sqrt(pow(px[j],2.0) + pow(py[j],2.0)) << " "  <<
           phi[j] << " " << sqrt(pow(x0[j],2.0) + pow(y0[j],2.0)) << " " 
           << eta[j] << " " << z0[j] << std::endl;

         countevt++;

       }
     }

     //t1->Show(i);
     
     if (countevt >= maxtracks)
       break;
  }
  ptfile.close();
  phifile.close();
  d0file.close();
  etafile.close();
  z0file.close();

  ss << "Event with 6 layer " << countevt << std::endl;

}

void readandtest (const std::string & fname, bool tkstubs, 
    bool bkstubs, int maxtracks)
{
  //TFile* inputFile = new TFile(fname.c_str(),"READ");
  // use xrootd as suggested 
  TFile* inputFile = TFile::Open(fname.c_str());

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

  /*
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
  */

  if (bkstubs)
  {
    std::ofstream bankstbfile("bankstub.txt");
    std::ofstream bankstbfileex("bankstub_notequal.txt");
    std::ofstream bankstbfileex1("bankstub_lesst6layers.txt");
    print_bankstub_new (inputFile, bankstbfile, bankstbfileex, 
        bankstbfileex1, (unsigned int)maxtracks);
    bankstbfile.close();
    bankstbfileex.close();
    bankstbfileex1.close();
  }

  if (tkstubs)
  {
    std::ofstream l1tkstubfile("l1tkstub.txt");
    //print_l1tkstub (inputFile, std::cout);
    print_l1tkstub (inputFile, l1tkstubfile, (unsigned int)maxtracks);
    l1tkstubfile.close();
  }

  inputFile->Close();
}

# ifndef __CINT__
int main(int argc, char ** argv) 
{
  bool tkstubs = false;
  bool bkstubs = true;
  int maxtracks = (int) STOPAFTERMAXEVT;

  while (1)
  {
    int c, option_index;
    static struct option long_options[] = {
      {"help", 0, NULL, 'h'},
      {"l1tk-stubs", 0, NULL, 'l'},
      {"max-track", 0, NULL, 'm'}, 
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "hlm:", long_options, &option_index);

    if (c == -1)
      break;

    switch (c)
    {
      case 'h':
        usage (argv[0]);
        break;
      case 'l':
        tkstubs = true;
        break;
      case 'm':
        maxtracks = atoi(optarg);
        break;
      default:
        usage (argv[0]);
        break;
    } 
  }

  if (optind >= argc) 
    usage (argv[0]);

  readandtest(argv[optind], tkstubs, bkstubs, maxtracks);

  return 0;
}
# endif
