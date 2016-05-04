#include <iostream>
#include <fstream>
#include <string>

#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TAxis.h"


TCanvas* canvas[20];
TLegend* leg[20];
TGraphErrors* g_phi_bias[16];
TGraphErrors* g_phi_sigma[16];
TGraphErrors* g_phi_bias_miss[6][16];
TGraphErrors* g_phi_sigma_miss[6][16];

int colors[16] = {1, 3, 2, 4, kOrange-3, 6, 7, 8,
		  kBlue-7, kGreen+4, kRed+3, kMagenta+3, kYellow+3, kViolet+9, kCyan-2, kBlue-6};

float pt_bin[7] = {4.,9.5,15.,21.5,37.5,75,150.};
float pt_bin_width[7] = {2.,2.5,3.,3.5,12.5,25.,50.};

TLatex* t = new TLatex(0.,0.,"#mu^{-}"); 


void do_pt_plots_1(unsigned int flag=0){

  TString filename = "pt_fitted_results.txt";

  ifstream infile;
  infile.open(filename.Data());

  string name, charge;
  float pt_bin_min[2][16][100], pt_bin_max[2][16][1000], phi_dummy[2][16][1000],
    phi_bias[2][16][1000], phi_bias_e[2][16][1000],
    phi_sigma[2][16][1000], phi_sigma_e[2][16][1000],
    phi_bias_m5[2][16][1000], phi_bias_m6[2][16][1000], phi_bias_m7[2][16][1000], 
    phi_bias_m8[2][16][1000], phi_bias_m9[2][16][1000], phi_bias_m10[2][16][1000], 
    phi_sigma_m5[2][16][1000], phi_sigma_m6[2][16][1000], phi_sigma_m7[2][16][1000], 
    phi_sigma_m8[2][16][1000], phi_sigma_m9[2][16][1000], phi_sigma_m10[2][16][1000]; 

  int nlines = 0;

  if (infile.is_open()){
    while ( !infile.eof() ) {

      int itow = nlines/14;
      int ibin = (nlines-14*itow)/2;

      //cout << nlines << " " << itow << " " << ibin << endl;
      
      phi_dummy[0][itow][ibin] = 0.;

      infile >> name >> pt_bin_min[0][itow][ibin] >> pt_bin_max[0][itow][ibin] >> charge
	     >> phi_bias[0][itow][ibin] >> phi_bias_e[0][itow][ibin]
	     >> phi_sigma[0][itow][ibin] >> phi_sigma_e[0][itow][ibin]
	     >> phi_bias_m5[0][itow][ibin] >> phi_bias_m6[0][itow][ibin] >> phi_bias_m7[0][itow][ibin] 
	     >> phi_bias_m8[0][itow][ibin] >> phi_bias_m9[0][itow][ibin] >> phi_bias_m10[0][itow][ibin] 
	     >> phi_sigma_m5[0][itow][ibin] >> phi_sigma_m6[0][itow][ibin] >> phi_sigma_m7[0][itow][ibin]
	     >> phi_sigma_m8[0][itow][ibin] >> phi_sigma_m9[0][itow][ibin] >> phi_sigma_m10[0][itow][ibin];
      
      nlines++;

      phi_dummy[1][itow][ibin] = 0.;

      infile >> name >> pt_bin_min[1][itow][ibin] >> pt_bin_max[1][itow][ibin] >> charge
	     >> phi_bias[1][itow][ibin] >> phi_bias_e[1][itow][ibin]
	     >> phi_sigma[1][itow][ibin] >> phi_sigma_e[1][itow][ibin]
	     >> phi_bias_m5[1][itow][ibin] >> phi_bias_m6[1][itow][ibin] >> phi_bias_m7[1][itow][ibin] 
	     >> phi_bias_m8[1][itow][ibin] >> phi_bias_m9[1][itow][ibin] >> phi_bias_m10[1][itow][ibin] 
	     >> phi_sigma_m5[1][itow][ibin] >> phi_sigma_m6[1][itow][ibin] >> phi_sigma_m7[1][itow][ibin]
	     >> phi_sigma_m8[1][itow][ibin] >> phi_sigma_m9[1][itow][ibin] >> phi_sigma_m10[1][itow][ibin];
      
      nlines++;
      
    }
    infile.close();
  }
  else std::cout << "Unable to open file"; 


  t->SetNDC(kTRUE);
  t->SetX(0.4);
  t->SetY(0.8);
  t->SetTextSize(0.06);

  //====================================================================================================================
  //====================================================================================================================
  if ( flag == 1 ) {

    canvas[0] = new TCanvas("c0","c0",800,600);
    canvas[0]->SetGridy(kTRUE);
    canvas[0]->SetLogx(kTRUE);
    canvas[0]->SetFillColor(kWhite);

    leg[0] = new TLegend(0.136935, 0.447552, 0.346734, 0.868881);
    
    for (int itow=0; itow<16; ++itow) {
      g_phi_bias[itow] = new TGraphErrors(7,pt_bin,phi_bias[1][itow],pt_bin_width,phi_bias_e[1][itow]);
      g_phi_bias[itow]->SetTitle("");
      g_phi_bias[itow]->SetMarkerStyle(20);
      g_phi_bias[itow]->SetMarkerColor(colors[itow]);
      g_phi_bias[itow]->SetLineColor(colors[itow]);
      
      if ( itow == 0 ){
	g_phi_bias[itow]->Draw("AP");
	g_phi_bias[itow]->SetMinimum(-0.055);
	g_phi_bias[itow]->SetMaximum(0.4);
	g_phi_bias[itow]->GetXaxis()->SetTitle("p_{T}^{bin} [GeV/c]");
	g_phi_bias[itow]->GetYaxis()->SetTitle("#LTp_{T}^{gen}-p_{T}^{fit}#GT/p_{T}^{gen}");
	g_phi_bias[itow]->GetYaxis()->SetTitleOffset(1.2);
      }
      else {
	g_phi_bias[itow]->Draw("P");
      }

      leg[0]->AddEntry(g_phi_bias[itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[0]->Draw();
    t->Draw();
    canvas[0]->Update();
    canvas[0]->SaveAs("pt_bias_mum.png");


    //---------------------------------------------------------------------------------------------------------------


    canvas[1] = new TCanvas("c1","c1",800,600);
    canvas[1]->SetGridy(kTRUE);
    canvas[1]->SetLogx(kTRUE);
    canvas[1]->SetFillColor(kWhite);

    leg[1] = new TLegend(0.136935, 0.447552, 0.346734, 0.868881);

    for (int itow=0; itow<16; ++itow) {
      g_phi_sigma[itow] = new TGraphErrors(7,pt_bin,phi_sigma[1][itow],pt_bin_width,phi_sigma_e[1][itow]);
      g_phi_sigma[itow]->SetTitle("");
      g_phi_sigma[itow]->SetMarkerStyle(20);
      g_phi_sigma[itow]->SetMarkerColor(colors[itow]);
      g_phi_sigma[itow]->SetLineColor(colors[itow]);
      
      if ( itow == 0 ){
	g_phi_sigma[itow]->Draw("AP");
	g_phi_sigma[itow]->SetMinimum(0.5);
	g_phi_sigma[itow]->SetMaximum(4.4);
	g_phi_sigma[itow]->GetXaxis()->SetTitle("p_{T}^{bin} [GeV/c]");
	g_phi_sigma[itow]->GetYaxis()->SetTitle("#sigma_{p_{T}^{gen}-p_{T}^{fit}}/p_{T}^{gen}");
	g_phi_sigma[itow]->GetYaxis()->SetTitleOffset(1.2);
      }
      else {
	g_phi_sigma[itow]->Draw("P");
      }

      leg[1]->AddEntry(g_phi_sigma[itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[1]->Draw();
    t->Draw();
    canvas[1]->Update();
    canvas[1]->SaveAs("pt_sigma_mum.png");

    
  }
  //====================================================================================================================
  //====================================================================================================================
  else if ( flag == 2 ) {
    
    //for (int itow=2;itow<3;++itow){
    //  cout << " ----> " << itow << endl;
    //  for ( int i=0; i<7; ++i){
    //	cout << "    bin " << i << "  -  " <<  phi_bias_m5[1][itow][i] << " " <<  phi_bias_m6[1][itow][i]
    //	     << " " <<  phi_bias_m7[1][itow][i]<< " " <<  phi_bias_m8[1][itow][i]<< " " <<  phi_bias_m9[1][itow][i]<< " " <<  phi_bias_m10[1][itow][i]
    //	     << endl;
    //  }
    //
    //}
  
  //return;

  

    canvas[3] = new TCanvas("c3","c3",800,600);
    canvas[3]->SetGridy(kTRUE);
    canvas[3]->SetLogx(kTRUE);
    canvas[3]->SetFillColor(kWhite);

    leg[3] = new TLegend(0.449749,0.438811,0.659548,0.858392);

    //for (int itow=0; itow<16; ++itow) {
      for (int itow=2; itow<3; ++itow) {
      g_phi_bias_miss[0][itow] = new TGraphErrors(7,pt_bin,phi_bias_m5[1][itow],pt_bin_width,phi_dummy[1][itow]);
      //g_phi_bias_miss[0][itow]->Print();
      //g_phi_bias_miss[0][itow]->SetTitle("Missing layer 5");
      g_phi_bias_miss[0][itow]->SetTitle("Trigger tower 18");
      g_phi_bias_miss[0][itow]->SetMarkerStyle(20);
      g_phi_bias_miss[0][itow]->SetMarkerColor(colors[itow]);
      g_phi_bias_miss[0][itow]->SetLineColor(colors[itow]);
      
      if ( itow == 2 ){
	g_phi_bias_miss[0][itow]->Draw("AP");
	//g_phi_bias_miss[0][itow]->SetMinimum(-5.e-6);
	//g_phi_bias_miss[0][itow]->SetMaximum(10.e-6);
	g_phi_bias_miss[0][itow]->GetXaxis()->SetTitle("p_{T}^{bin} [GeV/c]");
	g_phi_bias_miss[0][itow]->GetYaxis()->SetTitle("#LT#phi_{gen} - #phi_{fit}#GT [rad]");
	g_phi_bias_miss[0][itow]->GetYaxis()->SetTitleOffset(1.2);
      }
      else {
	g_phi_bias_miss[0][itow]->Draw("P");
      }

      leg[3]->AddEntry(g_phi_bias_miss[0][itow],Form("Tower %d", itow+16),"P");
    
    }

      //leg[3]->Draw();
    canvas[3]->Update();
    canvas[3]->SaveAs("phi_bias_miss5.png");


    //---------------------------------------------------------------------------------------------------------------


    canvas[4] = new TCanvas("c4","c4",800,600);
    canvas[4]->SetGridy(kTRUE);
    canvas[4]->SetLogx(kTRUE);
    canvas[4]->SetFillColor(kWhite);

    leg[4] = new TLegend(0.449749,0.438811,0.659548,0.858392);

    //for (int itow=0; itow<16; ++itow) {
      for (int itow=2; itow<3; ++itow) {
      g_phi_bias_miss[1][itow] = new TGraphErrors(7,pt_bin,phi_bias_m6[1][itow],pt_bin_width,phi_dummy[1][itow]);
      g_phi_bias_miss[1][itow]->SetTitle("Missing layer 6");
      g_phi_bias_miss[1][itow]->SetMarkerStyle(20);
      g_phi_bias_miss[1][itow]->SetMarkerColor(colors[itow]);
      g_phi_bias_miss[1][itow]->SetLineColor(colors[itow]);
      
      if ( itow == 2 ){
	g_phi_bias_miss[1][itow]->Draw("AP");
	//g_phi_bias_miss[1][itow]->SetMinimum(-5.e-6);
	//g_phi_bias_miss[1][itow]->SetMaximum(10.e-6);
	g_phi_bias_miss[1][itow]->GetXaxis()->SetTitle("p_{T}^{bin} [GeV/c]");
	g_phi_bias_miss[1][itow]->GetYaxis()->SetTitle("#LT#phi_{gen} - #phi_{fit}#GT [rad]");
	g_phi_bias_miss[1][itow]->GetYaxis()->SetTitleOffset(1.2);
      }
      else {
	g_phi_bias_miss[1][itow]->Draw("P");
      }

      leg[4]->AddEntry(g_phi_bias_miss[1][itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[4]->Draw();
    canvas[4]->Update();
    canvas[4]->SaveAs("phi_bias_miss6.png");
    

    //---------------------------------------------------------------------------------------------------------------


    canvas[5] = new TCanvas("c5","c5",800,600);
    canvas[5]->SetGridy(kTRUE);
    canvas[5]->SetLogx(kTRUE);
    canvas[5]->SetFillColor(kWhite);

    leg[5] = new TLegend(0.449749,0.438811,0.659548,0.858392);

    //for (int itow=0; itow<16; ++itow) {
      for (int itow=2; itow<3; ++itow) {
      g_phi_bias_miss[2][itow] = new TGraphErrors(7,pt_bin,phi_bias_m7[1][itow],pt_bin_width,phi_dummy[1][itow]);
      g_phi_bias_miss[2][itow]->SetTitle("Missing layer 7");
      g_phi_bias_miss[2][itow]->SetMarkerStyle(20);
      g_phi_bias_miss[2][itow]->SetMarkerColor(colors[itow]);
      g_phi_bias_miss[2][itow]->SetLineColor(colors[itow]);
      
      if ( itow == 2 ){
	g_phi_bias_miss[2][itow]->Draw("AP");
	//g_phi_bias_miss[2][itow]->SetMinimum(-5.e-6);
	//g_phi_bias_miss[2][itow]->SetMaximum(10.e-6);
	g_phi_bias_miss[2][itow]->GetXaxis()->SetTitle("p_{T}^{bin} [GeV/c]");
	g_phi_bias_miss[2][itow]->GetYaxis()->SetTitle("#LT#phi_{gen} - #phi_{fit}#GT [rad]");
	g_phi_bias_miss[2][itow]->GetYaxis()->SetTitleOffset(1.2);
      }
      else {
	g_phi_bias_miss[2][itow]->Draw("P");
      }

      leg[5]->AddEntry(g_phi_bias_miss[2][itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[5]->Draw();
    canvas[5]->Update();
    canvas[5]->SaveAs("phi_bias_miss7.png");


    //---------------------------------------------------------------------------------------------------------------


    canvas[6] = new TCanvas("c6","c6",800,600);
    canvas[6]->SetGridy(kTRUE);
    canvas[6]->SetLogx(kTRUE);
    canvas[6]->SetFillColor(kWhite);

    leg[6] = new TLegend(0.449749,0.438811,0.659548,0.858392);

    //for (int itow=0; itow<16; ++itow) {
      for (int itow=2; itow<3; ++itow) {
      g_phi_bias_miss[3][itow] = new TGraphErrors(7,pt_bin,phi_bias_m8[1][itow],pt_bin_width,phi_dummy[1][itow]);
      g_phi_bias_miss[3][itow]->SetTitle("Missing layer 8");
      g_phi_bias_miss[3][itow]->SetMarkerStyle(20);
      g_phi_bias_miss[3][itow]->SetMarkerColor(colors[itow]);
      g_phi_bias_miss[3][itow]->SetLineColor(colors[itow]);
      
      if ( itow == 2 ){
	g_phi_bias_miss[3][itow]->Draw("AP");
	//g_phi_bias_miss[3][itow]->SetMinimum(-5.e-6);
	//g_phi_bias_miss[3][itow]->SetMaximum(10.e-6);
	g_phi_bias_miss[3][itow]->GetXaxis()->SetTitle("p_{T}^{bin} [GeV/c]");
	g_phi_bias_miss[3][itow]->GetYaxis()->SetTitle("#LT#phi_{gen} - #phi_{fit}#GT [rad]");
	g_phi_bias_miss[3][itow]->GetYaxis()->SetTitleOffset(1.2);
      }
      else {
	g_phi_bias_miss[3][itow]->Draw("P");
      }

      leg[6]->AddEntry(g_phi_bias_miss[3][itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[6]->Draw();
    canvas[6]->Update();
    canvas[6]->SaveAs("phi_bias_miss8.png");


    //---------------------------------------------------------------------------------------------------------------


    canvas[7] = new TCanvas("c7","c7",800,600);
    canvas[7]->SetGridy(kTRUE);
    canvas[7]->SetLogx(kTRUE);
    canvas[7]->SetFillColor(kWhite);

    leg[7] = new TLegend(0.449749,0.438811,0.659548,0.858392);

    //for (int itow=0; itow<16; ++itow) {
      for (int itow=2; itow<3; ++itow) {
      g_phi_bias_miss[4][itow] = new TGraphErrors(7,pt_bin,phi_bias_m9[1][itow],pt_bin_width,phi_dummy[1][itow]);
      g_phi_bias_miss[4][itow]->SetTitle("Missing layer 9");
      g_phi_bias_miss[4][itow]->SetMarkerStyle(20);
      g_phi_bias_miss[4][itow]->SetMarkerColor(colors[itow]);
      g_phi_bias_miss[4][itow]->SetLineColor(colors[itow]);
      
      if ( itow == 2 ){
	g_phi_bias_miss[4][itow]->Draw("AP");
	//g_phi_bias_miss[4][itow]->SetMinimum(-5.e-6);
	//g_phi_bias_miss[4][itow]->SetMaximum(10.e-6);
	g_phi_bias_miss[4][itow]->GetXaxis()->SetTitle("p_{T}^{bin} [GeV/c]");
	g_phi_bias_miss[4][itow]->GetYaxis()->SetTitle("#LT#phi_{gen} - #phi_{fit}#GT [rad]");
	g_phi_bias_miss[4][itow]->GetYaxis()->SetTitleOffset(1.2);
      }
      else {
	g_phi_bias_miss[4][itow]->Draw("P");
      }

      leg[7]->AddEntry(g_phi_bias_miss[4][itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[7]->Draw();
    canvas[7]->Update();
    canvas[7]->SaveAs("phi_bias_miss9.png");

    
    //---------------------------------------------------------------------------------------------------------------


    canvas[8] = new TCanvas("c8","c8",800,600);
    canvas[8]->SetGridy(kTRUE);
    canvas[8]->SetLogx(kTRUE);
    canvas[8]->SetFillColor(kWhite);

    leg[8] = new TLegend(0.449749,0.438811,0.659548,0.858392);

    //for (int itow=0; itow<16; ++itow) {
      for (int itow=2; itow<3; ++itow) {
      g_phi_bias_miss[5][itow] = new TGraphErrors(7,pt_bin,phi_bias_m10[1][itow],pt_bin_width,phi_dummy[1][itow]);
      g_phi_bias_miss[5][itow]->SetTitle("Missing layer 10");
      g_phi_bias_miss[5][itow]->SetMarkerStyle(20);
      g_phi_bias_miss[5][itow]->SetMarkerColor(colors[itow]);
      g_phi_bias_miss[5][itow]->SetLineColor(colors[itow]);
      
      if ( itow == 2 ){
	g_phi_bias_miss[5][itow]->Draw("AP");
	//g_phi_bias_miss[5][itow]->SetMinimum(-5.e-6);
	//g_phi_bias_miss[5][itow]->SetMaximum(10.e-6);
	g_phi_bias_miss[5][itow]->GetXaxis()->SetTitle("p_{T}^{bin} [GeV/c]");
	g_phi_bias_miss[5][itow]->GetYaxis()->SetTitle("#LT#phi_{gen} - #phi_{fit}#GT [rad]");
	g_phi_bias_miss[5][itow]->GetYaxis()->SetTitleOffset(1.2);
      }
      else {
	g_phi_bias_miss[5][itow]->Draw("P");
      }

      leg[8]->AddEntry(g_phi_bias_miss[5][itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[8]->Draw();
    canvas[8]->Update();
    canvas[8]->SaveAs("phi_bias_miss10.png");
    


    canvas[18] = new TCanvas("c18","c18",800,600);
    canvas[18]->SetGridy(kTRUE);
    canvas[18]->SetLogx(kTRUE);
    canvas[18]->SetFillColor(kWhite);

    leg[18] = new TLegend(0.136935, 0.447552, 0.346734, 0.868881);

    for (int i=0; i<6; ++i){
      
      int itow = 2;
      
      g_phi_bias_miss[i][itow]->SetMarkerColor(colors[i]);
      g_phi_bias_miss[i][itow]->SetLineColor(colors[i]);

      if ( i == 0 ){
	g_phi_bias_miss[i][itow]->Draw("AP");
	g_phi_bias_miss[i][itow]->SetMinimum(-0.08);
	g_phi_bias_miss[i][itow]->SetMaximum(0.24);
	g_phi_bias_miss[i][itow]->GetXaxis()->SetTitle("p_{T}^{bin} [GeV/c]");
	g_phi_bias_miss[i][itow]->GetYaxis()->SetTitle("#LTp_{T}^{gen}-p_{T}^{fit}#GT/p_{T}^{gen}");
	g_phi_bias_miss[i][itow]->GetYaxis()->SetTitleOffset(1.2);
      }
      else
	g_phi_bias_miss[i][itow]->Draw("P");


      
      leg[18]->AddEntry(g_phi_bias_miss[i][itow],Form("Missing layer %d", i+5),"P");
    
    }

    leg[18]->Draw();
    t->Draw();
    canvas[18]->Update();
    canvas[18]->SaveAs("pt_bias_miss_mum.png");






  }
  //====================================================================================================================
  //====================================================================================================================
  else if ( flag == 3 ) {

    canvas[9] = new TCanvas("c9","c9",800,600);
    canvas[9]->SetGridy(kTRUE);
    canvas[9]->SetLogx(kTRUE);
    canvas[9]->SetFillColor(kWhite);

    leg[9] = new TLegend(0.449749,0.438811,0.659548,0.858392);

    //for (int itow=0; itow<16; ++itow) {
      for (int itow=2; itow<3; ++itow) {
      g_phi_sigma_miss[0][itow] = new TGraphErrors(7,pt_bin,phi_sigma_m5[1][itow],pt_bin_width,phi_dummy[1][itow]);
      //g_phi_sigma_miss[0][itow]->SetTitle("Missing layer 5");
      g_phi_sigma_miss[0][itow]->SetTitle("Trigger tower 18");
      g_phi_sigma_miss[0][itow]->SetMarkerStyle(20);
      g_phi_sigma_miss[0][itow]->SetMarkerColor(colors[itow]);
      g_phi_sigma_miss[0][itow]->SetLineColor(colors[itow]);
      
      if ( itow == 2 ){
	g_phi_sigma_miss[0][itow]->Draw("AP");
	//g_phi_sigma_miss[0][itow]->SetMinimum(-5.e-6);
	//g_phi_sigma_miss[0][itow]->SetMaximum(10.e-6);
	g_phi_sigma_miss[0][itow]->GetXaxis()->SetTitle("p_{T}^{bin} [GeV/c]");
	g_phi_sigma_miss[0][itow]->GetYaxis()->SetTitle("#LT#phi_{gen} - #phi_{fit}#GT [rad]");
	g_phi_sigma_miss[0][itow]->GetYaxis()->SetTitleOffset(1.2);
      }
      else {
	g_phi_sigma_miss[0][itow]->Draw("P");
      }

      leg[9]->AddEntry(g_phi_sigma_miss[0][itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[9]->Draw();
    canvas[9]->Update();
    canvas[9]->SaveAs("phi_sigma_miss5.png");


    //---------------------------------------------------------------------------------------------------------------


    canvas[10] = new TCanvas("c10","c10",800,600);
    canvas[10]->SetGridy(kTRUE);
    //canvas[10]->SetLogx(kTRUE);
    canvas[10]->SetFillColor(kWhite);

    leg[10] = new TLegend(0.449749,0.438811,0.659548,0.858392);

    //for (int itow=0; itow<16; ++itow) {
      for (int itow=2; itow<3; ++itow) {
      g_phi_sigma_miss[1][itow] = new TGraphErrors(7,pt_bin,phi_sigma_m6[1][itow],pt_bin_width,phi_dummy[1][itow]);
      g_phi_sigma_miss[1][itow]->SetTitle("Missing layer 6");
      g_phi_sigma_miss[1][itow]->SetMarkerStyle(20);
      g_phi_sigma_miss[1][itow]->SetMarkerColor(colors[itow]);
      g_phi_sigma_miss[1][itow]->SetLineColor(colors[itow]);
      
      if ( itow == 2 ){
	g_phi_sigma_miss[1][itow]->Draw("AP");
	g_phi_sigma_miss[1][itow]->SetMinimum(-5.e-6);
	g_phi_sigma_miss[1][itow]->SetMaximum(10.e-6);
	g_phi_sigma_miss[1][itow]->GetXaxis()->SetTitle("p_{T}^{bin} [GeV/c]");
	g_phi_sigma_miss[1][itow]->GetYaxis()->SetTitle("#LT#phi_{gen} - #phi_{fit}#GT [rad]");
	g_phi_sigma_miss[1][itow]->GetYaxis()->SetTitleOffset(1.2);
      }
      else {
	g_phi_sigma_miss[1][itow]->Draw("P");
      }

      leg[10]->AddEntry(g_phi_sigma_miss[1][itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[10]->Draw();
    canvas[10]->Update();
    canvas[10]->SaveAs("phi_sigma_miss6.png");
    

    //---------------------------------------------------------------------------------------------------------------


    canvas[11] = new TCanvas("c11","c11",800,600);
    canvas[11]->SetGridy(kTRUE);
    //canvas[11]->SetLogx(kTRUE);
    canvas[11]->SetFillColor(kWhite);
 
    leg[11] = new TLegend(0.449749,0.438811,0.659548,0.858392);
 
    //for (int itow=0; itow<16; ++itow) {
      for (int itow=2; itow<3; ++itow) {
      g_phi_sigma_miss[2][itow] = new TGraphErrors(7,pt_bin,phi_sigma_m7[1][itow],pt_bin_width,phi_dummy[1][itow]);
      g_phi_sigma_miss[2][itow]->SetTitle("Missing layer 7");
      g_phi_sigma_miss[2][itow]->SetMarkerStyle(20);
      g_phi_sigma_miss[2][itow]->SetMarkerColor(colors[itow]);
      g_phi_sigma_miss[2][itow]->SetLineColor(colors[itow]);
      
      if ( itow == 2 ){
 	g_phi_sigma_miss[2][itow]->Draw("AP");
 	//g_phi_sigma_miss[2][itow]->SetMinimum(-5.e-6);
 	//g_phi_sigma_miss[2][itow]->SetMaximum(10.e-6);
 	g_phi_sigma_miss[2][itow]->GetXaxis()->SetTitle("p_{T}^{bin} [GeV/c]");
 	g_phi_sigma_miss[2][itow]->GetYaxis()->SetTitle("#LT#phi_{gen} - #phi_{fit}#GT [rad]");
 	g_phi_sigma_miss[2][itow]->GetYaxis()->SetTitleOffset(1.2);
      }
      else {
 	g_phi_sigma_miss[2][itow]->Draw("P");
      }
 
      leg[11]->AddEntry(g_phi_sigma_miss[2][itow],Form("Tower %d", itow+16),"P");
    
    }
 
    leg[11]->Draw();
    canvas[11]->Update();
    canvas[11]->SaveAs("phi_sigma_miss7.png");


    //---------------------------------------------------------------------------------------------------------------


    canvas[12] = new TCanvas("c12","c12",800,600);
    canvas[12]->SetGridy(kTRUE);
    //canvas[12]->SetLogx(kTRUE);
    canvas[12]->SetFillColor(kWhite);

    leg[12] = new TLegend(0.449749,0.438811,0.659548,0.858392);

    //for (int itow=0; itow<16; ++itow) {
      for (int itow=2; itow<3; ++itow) {
      g_phi_sigma_miss[3][itow] = new TGraphErrors(7,pt_bin,phi_sigma_m8[1][itow],pt_bin_width,phi_dummy[1][itow]);
      g_phi_sigma_miss[3][itow]->SetTitle("Missing layer 8");
      g_phi_sigma_miss[3][itow]->SetMarkerStyle(20);
      g_phi_sigma_miss[3][itow]->SetMarkerColor(colors[itow]);
      g_phi_sigma_miss[3][itow]->SetLineColor(colors[itow]);
      
      if ( itow == 2 ){
	g_phi_sigma_miss[3][itow]->Draw("AP");
	//g_phi_sigma_miss[3][itow]->SetMinimum(-5.e-6);
	//g_phi_sigma_miss[3][itow]->SetMaximum(10.e-6);
	g_phi_sigma_miss[3][itow]->GetXaxis()->SetTitle("p_{T}^{bin} [GeV/c]");
	g_phi_sigma_miss[3][itow]->GetYaxis()->SetTitle("#LT#phi_{gen} - #phi_{fit}#GT [rad]");
	g_phi_sigma_miss[3][itow]->GetYaxis()->SetTitleOffset(1.2);
      }
      else {
	g_phi_sigma_miss[3][itow]->Draw("P");
      }

      leg[12]->AddEntry(g_phi_sigma_miss[3][itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[12]->Draw();
    canvas[12]->Update();
    canvas[12]->SaveAs("phi_sigma_miss8.png");


    //---------------------------------------------------------------------------------------------------------------


    canvas[13] = new TCanvas("c13","c13",800,600);
    canvas[13]->SetGridy(kTRUE);
    //canvas[13]->SetLogx(kTRUE);
    canvas[13]->SetFillColor(kWhite);

    leg[13] = new TLegend(0.449749,0.438811,0.659548,0.858392);

    //for (int itow=0; itow<16; ++itow) {
      for (int itow=2; itow<3; ++itow) {
      g_phi_sigma_miss[4][itow] = new TGraphErrors(7,pt_bin,phi_sigma_m9[1][itow],pt_bin_width,phi_dummy[1][itow]);
      g_phi_sigma_miss[4][itow]->SetTitle("Missing layer 9");
      g_phi_sigma_miss[4][itow]->SetMarkerStyle(20);
      g_phi_sigma_miss[4][itow]->SetMarkerColor(colors[itow]);
      g_phi_sigma_miss[4][itow]->SetLineColor(colors[itow]);
      
      if ( itow == 2 ){
	g_phi_sigma_miss[4][itow]->Draw("AP");
	//g_phi_sigma_miss[4][itow]->SetMinimum(-5.e-6);
	//g_phi_sigma_miss[4][itow]->SetMaximum(10.e-6);
	g_phi_sigma_miss[4][itow]->GetXaxis()->SetTitle("p_{T}^{bin} [GeV/c]");
	g_phi_sigma_miss[4][itow]->GetYaxis()->SetTitle("#LT#phi_{gen} - #phi_{fit}#GT [rad]");
	g_phi_sigma_miss[4][itow]->GetYaxis()->SetTitleOffset(1.2);
      }
      else {
	g_phi_sigma_miss[4][itow]->Draw("P");
      }

      leg[13]->AddEntry(g_phi_sigma_miss[4][itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[13]->Draw();
    canvas[13]->Update();
    canvas[13]->SaveAs("phi_sigma_miss9.png");

    
    //---------------------------------------------------------------------------------------------------------------


    canvas[14] = new TCanvas("c14","c14",800,600);
    canvas[14]->SetGridy(kTRUE);
    //canvas[14]->SetLogx(kTRUE);
    canvas[14]->SetFillColor(kWhite);

    leg[14] = new TLegend(0.449749,0.438811,0.659548,0.858392);

    //for (int itow=0; itow<16; ++itow) {
      for (int itow=2; itow<3; ++itow) {
      g_phi_sigma_miss[5][itow] = new TGraphErrors(7,pt_bin,phi_sigma_m10[1][itow],pt_bin_width,phi_dummy[1][itow]);
      g_phi_sigma_miss[5][itow]->SetTitle("Missing layer 10");
      g_phi_sigma_miss[5][itow]->SetMarkerStyle(20);
      g_phi_sigma_miss[5][itow]->SetMarkerColor(colors[itow]);
      g_phi_sigma_miss[5][itow]->SetLineColor(colors[itow]);
      
      if ( itow == 2 ){
	g_phi_sigma_miss[5][itow]->Draw("AP");
	//g_phi_sigma_miss[5][itow]->SetMinimum(-5.e-6);
	//g_phi_sigma_miss[5][itow]->SetMaximum(10.e-6);
	g_phi_sigma_miss[5][itow]->GetXaxis()->SetTitle("p_{T}^{bin} [GeV/c]");
	g_phi_sigma_miss[5][itow]->GetYaxis()->SetTitle("#LT#phi_{gen} - #phi_{fit}#GT [rad]");
	g_phi_sigma_miss[5][itow]->GetYaxis()->SetTitleOffset(1.2);
      }
      else {
	g_phi_sigma_miss[5][itow]->Draw("P");
      }

      leg[14]->AddEntry(g_phi_sigma_miss[5][itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[14]->Draw();
    canvas[14]->Update();
    canvas[14]->SaveAs("phi_sigma_miss10.png");



    canvas[19] = new TCanvas("c19","c19",800,600);
    canvas[19]->SetGridy(kTRUE);
    canvas[19]->SetLogx(kTRUE);
    canvas[19]->SetFillColor(kWhite);

    leg[19] = new TLegend(0.136935, 0.447552, 0.346734, 0.868881);

    for (int i=0; i<6; ++i){
      
      int itow = 2;
      
      g_phi_sigma_miss[i][itow]->SetMarkerColor(colors[i]);
      g_phi_sigma_miss[i][itow]->SetLineColor(colors[i]);

      if ( i == 0 ){
	g_phi_sigma_miss[i][itow]->Draw("AP");
	g_phi_sigma_miss[i][itow]->SetMinimum(0.5);
	g_phi_sigma_miss[i][itow]->SetMaximum(6.3);
	g_phi_sigma_miss[i][itow]->GetXaxis()->SetTitle("p_{T}^{bin} [GeV/c]");
	g_phi_sigma_miss[i][itow]->GetYaxis()->SetTitle("#sigma_{p_{T}^{gen}-p_{T}^{fit}}/p_{T}^{gen}");
	g_phi_sigma_miss[i][itow]->GetYaxis()->SetTitleOffset(1.2);
      }
      else
	g_phi_sigma_miss[i][itow]->Draw("P");


      
      leg[19]->AddEntry(g_phi_sigma_miss[i][itow],Form("Missing layer %d", i+5),"P");
    
    }

    leg[19]->Draw();
    t->Draw();
    canvas[19]->Update();
    canvas[19]->SaveAs("pt_sigma_miss_mum.png");



  }

}
