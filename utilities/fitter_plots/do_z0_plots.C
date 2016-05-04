#include <iostream>
#include <fstream>
#include <string>

#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TAxis.h"


TCanvas* canvas[17];
TLegend* leg[17];
TGraphErrors* g_eta_bias[16];
TGraphErrors* g_eta_sigma[16];
TGraphErrors* g_eta_bias_miss[3][16];
TGraphErrors* g_eta_sigma_miss[3][16];

int colors[16] = {1, 2, 3, 4, kOrange-3, 6, 7, 8,
		  kBlue-7, kGreen+4, kRed+3, kMagenta+3, kYellow+3, kViolet+9, kCyan-2, kBlue-6};


void do_z0_plots(unsigned int flag=0){

  TString filename = "z0_fitted_results.txt";

  ifstream infile;
  infile.open(filename.Data());

  string name;
  float eta_bin[16][1000], eta_bin_width[16][1000], eta_dummy[16][1000],
    eta_bias[16][1000], eta_bias_e[16][1000],
    eta_sigma[16][1000], eta_sigma_e[16][1000],
    eta_bias_m5[16][1000], eta_bias_m6[16][1000], eta_bias_m7[16][1000], 
    eta_sigma_m5[16][1000], eta_sigma_m6[16][1000], eta_sigma_m7[16][1000]; 

  int nlines = 0;

  if (infile.is_open()){
    while ( !infile.eof() ) {

      int itow = nlines/20;
      int ibin = nlines-20*itow;

      eta_bin_width[itow][ibin] = 0.025;
      eta_dummy[itow][ibin] = 0.;
      infile >> name >> eta_bin[itow][ibin] 
	     >> eta_bias[itow][ibin] >> eta_bias_e[itow][ibin]
	     >> eta_sigma[itow][ibin] >> eta_sigma_e[itow][ibin]
	     >> eta_bias_m5[itow][ibin] >> eta_bias_m6[itow][ibin] >> eta_bias_m7[itow][ibin] 
	     >> eta_sigma_m5[itow][ibin] >> eta_sigma_m6[itow][ibin] >> eta_sigma_m7[itow][ibin]; 
      
      //cout << nlines << " --> " << itow << " " << ibin << endl;
      //cout << nlines << " " << name << " " <<  eta_bin[itow][ibin] << " " <<  eta_bias[itow][ibin] << " " <<  eta_bias_e[itow][ibin] << endl; 

      nlines++;

    }
    infile.close();
  }
  else std::cout << "Unable to open file"; 

  
  //====================================================================================================================
  //====================================================================================================================
  if ( flag == 1 ) {

    canvas[0] = new TCanvas("c0","c0",800,600);
    canvas[0]->SetGridy(kTRUE);
    canvas[0]->SetFillColor(kWhite);

    leg[0] = new TLegend(0.768844,0.543706,0.977387,0.963287);

    for (int itow=0; itow<8; ++itow) {
      g_eta_bias[itow] = new TGraphErrors(20,eta_bin[itow],eta_bias[itow],eta_bin_width[itow],eta_bias_e[itow]);
      g_eta_bias[itow]->SetTitle("");
      g_eta_bias[itow]->SetMarkerStyle(20);
      g_eta_bias[itow]->SetMarkerColor(colors[itow]);
      g_eta_bias[itow]->SetLineColor(colors[itow]);
      
      if ( itow == 0 ){
	g_eta_bias[itow]->Draw("AP");
	g_eta_bias[itow]->SetMinimum(-1e-2);
	g_eta_bias[itow]->SetMaximum( 1e-2);
	g_eta_bias[itow]->GetXaxis()->SetTitle("#eta_{bin}");
	g_eta_bias[itow]->GetYaxis()->SetTitle("#LTz_{0}^{gen}-z_{0}^{fit}#GT [cm]");
	g_eta_bias[itow]->GetYaxis()->SetTitleOffset(1.2);
      }
      else {
	g_eta_bias[itow]->Draw("P");
      }

      leg[0]->AddEntry(g_eta_bias[itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[0]->Draw();
    canvas[0]->Update();
    canvas[0]->SaveAs("z0_bias_tow16-23.png");

    // -----------------------------------------------------------------------------------------------------------------

    canvas[1] = new TCanvas("c1","c1",800,600);
    canvas[1]->SetGridy(kTRUE);
    canvas[1]->SetFillColor(kWhite);

    leg[1] = new TLegend(0.768844,0.543706,0.977387,0.963287);

    for (int itow=8; itow<16; ++itow) {
      g_eta_bias[itow] = new TGraphErrors(20,eta_bin[itow],eta_bias[itow],eta_bin_width[itow],eta_bias_e[itow]);
      g_eta_bias[itow]->SetTitle("");
      g_eta_bias[itow]->SetMarkerStyle(20);
      g_eta_bias[itow]->SetMarkerColor(colors[itow-8]);
      g_eta_bias[itow]->SetLineColor(colors[itow-8]);

      if ( itow == 8 ){
	g_eta_bias[itow]->Draw("AP");
	g_eta_bias[itow]->SetMinimum(-1e-2);
	g_eta_bias[itow]->SetMaximum( 1e-2);
	g_eta_bias[itow]->GetXaxis()->SetTitle("#eta_{bin}");
	g_eta_bias[itow]->GetYaxis()->SetTitle("#LTz_{0}^{gen}-z_{0}^{fit}#GT [cm]");
	g_eta_bias[itow]->GetYaxis()->SetTitleOffset(1.2);
      }
      else {
	g_eta_bias[itow]->Draw("P");
      }

      leg[1]->AddEntry(g_eta_bias[itow],Form("Tower %d", itow+16),"P");

    }
  
    leg[1]->Draw();
    canvas[1]->Update();
    canvas[1]->SaveAs("z0_bias_tow24-31.png");
  
  }
  //====================================================================================================================
  //====================================================================================================================
  else if ( flag == 2 ) {

    canvas[3] = new TCanvas("c3","c3",800,600);
    canvas[3]->SetGridy(kTRUE);
    canvas[3]->SetFillColor(kWhite);

    leg[3] = new TLegend(0.449749,0.438811,0.659548,0.858392);

    for (int itow=0; itow<8; ++itow) {
      g_eta_sigma[itow] = new TGraphErrors(20,eta_bin[itow],eta_sigma[itow],eta_bin_width[itow],eta_sigma_e[itow]);
      g_eta_sigma[itow]->SetTitle("");
      g_eta_sigma[itow]->SetMarkerStyle(20);
      g_eta_sigma[itow]->SetMarkerColor(colors[itow]);
      g_eta_sigma[itow]->SetLineColor(colors[itow]);
      
      if ( itow == 0 ){
	g_eta_sigma[itow]->Draw("AP");
	g_eta_sigma[itow]->SetMinimum(0.0);
	g_eta_sigma[itow]->SetMaximum(0.5);
	g_eta_sigma[itow]->GetXaxis()->SetTitle("#eta_{bin}");
	g_eta_sigma[itow]->GetYaxis()->SetTitle("#sigma_{z_{0}^{gen}-z_{0}^{fit}} [cm]");
	g_eta_sigma[itow]->GetYaxis()->SetTitleOffset(1.2);
      }
      else {
	g_eta_sigma[itow]->Draw("P");
      }

      leg[3]->AddEntry(g_eta_sigma[itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[3]->Draw();
    canvas[3]->Update();
    canvas[3]->SaveAs("z0_sigma_tow16-23.png");

    // -----------------------------------------------------------------------------------------------------------------

    canvas[4] = new TCanvas("c4","c4",800,600);
    canvas[4]->SetGridy(kTRUE);
    canvas[4]->SetFillColor(kWhite);

    leg[4] = new TLegend(0.449749,0.438811,0.659548,0.858392);

    for (int itow=8; itow<16; ++itow) {
      g_eta_sigma[itow] = new TGraphErrors(20,eta_bin[itow],eta_sigma[itow],eta_bin_width[itow],eta_sigma_e[itow]);
      g_eta_sigma[itow]->SetTitle("");
      g_eta_sigma[itow]->SetMarkerStyle(20);
      g_eta_sigma[itow]->SetMarkerColor(colors[itow-8]);
      g_eta_sigma[itow]->SetLineColor(colors[itow-8]);
      
      if ( itow == 8 ){
	g_eta_sigma[itow]->Draw("AP");
	g_eta_sigma[itow]->SetMinimum(0.);
	g_eta_sigma[itow]->SetMaximum(0.5);
	g_eta_sigma[itow]->GetXaxis()->SetTitle("#eta_{bin}");
	g_eta_sigma[itow]->GetYaxis()->SetTitle("#sigma_{z_{0}^{gen}-z_{0}^{fit}} [cm]");
	g_eta_sigma[itow]->GetYaxis()->SetTitleOffset(1.2);
      }
      else {
	g_eta_sigma[itow]->Draw("P");
      }

      leg[4]->AddEntry(g_eta_sigma[itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[4]->Draw();
    canvas[4]->Update();
    canvas[4]->SaveAs("z0_sigma_tow24-31.png");

  }
  //====================================================================================================================
  //====================================================================================================================
  else if ( flag == 3 ) {

    canvas[5] = new TCanvas("c5","c5",800,600);
    canvas[5]->SetGridy(kTRUE);
    canvas[5]->SetFillColor(kWhite);

    leg[5] = new TLegend(0.768844,0.543706,0.977387,0.963287);

    for (int itow=0; itow<8; ++itow) {
      g_eta_bias_miss[0][itow] = new TGraphErrors(20,eta_bin[itow],eta_bias_m5[itow],eta_bin_width[itow],eta_dummy[itow]);
      g_eta_bias_miss[0][itow]->SetTitle("Missing layer 5");
      g_eta_bias_miss[0][itow]->SetMarkerStyle(20);
      g_eta_bias_miss[0][itow]->SetMarkerColor(colors[itow]);
      g_eta_bias_miss[0][itow]->SetLineColor(colors[itow]);

      if ( itow == 0 ){
	g_eta_bias_miss[0][itow]->Draw("AP");
	g_eta_bias_miss[0][itow]->SetMinimum(-1.e-2);
	g_eta_bias_miss[0][itow]->SetMaximum(1.e-2);
	g_eta_bias_miss[0][itow]->GetXaxis()->SetTitle("#eta_{bin}");
	g_eta_bias_miss[0][itow]->GetYaxis()->SetTitle("#LTz_{0}^{gen}-z_{0}^{fit}#GT [cm]");
	g_eta_bias_miss[0][itow]->GetYaxis()->SetTitleOffset(1.2);

      }
      else {
	g_eta_bias_miss[0][itow]->Draw("P");
      }

      leg[5]->AddEntry(g_eta_bias_miss[0][itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[5]->Draw();
    canvas[5]->Update();
    canvas[5]->SaveAs("z0_bias_missing5_tow16-23.png");


    // -----------------------------------------------------------------------------------------------------------------


    canvas[6] = new TCanvas("c6","c6",800,600);
    canvas[6]->SetGridy(kTRUE);
    canvas[6]->SetFillColor(kWhite);
    
    leg[6] = new TLegend(0.768844,0.543706,0.977387,0.963287);

    for (int itow=0; itow<8; ++itow) {
      g_eta_bias_miss[1][itow] = new TGraphErrors(20,eta_bin[itow],eta_bias_m6[itow],eta_bin_width[itow],eta_dummy[itow]);
      g_eta_bias_miss[1][itow]->SetTitle("Missing layer 6");
      g_eta_bias_miss[1][itow]->SetMarkerStyle(20);
      g_eta_bias_miss[1][itow]->SetMarkerColor(colors[itow]);
      g_eta_bias_miss[1][itow]->SetLineColor(colors[itow]);

      if ( itow == 0 ){
	g_eta_bias_miss[1][itow]->Draw("AP");
	g_eta_bias_miss[1][itow]->SetMinimum(-1.e-2);
	g_eta_bias_miss[1][itow]->SetMaximum( 1.e-2);
	g_eta_bias_miss[1][itow]->GetXaxis()->SetTitle("#eta_{bin}");
	g_eta_bias_miss[1][itow]->GetYaxis()->SetTitle("#LTz_{0}^{gen}-z_{0}^{fit}#GT [cm]");
	g_eta_bias_miss[1][itow]->GetYaxis()->SetTitleOffset(1.2);

      }
      else {
	g_eta_bias_miss[1][itow]->Draw("P");
      }

      leg[6]->AddEntry(g_eta_bias_miss[0][itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[6]->Draw();
    canvas[6]->Update();
    canvas[6]->SaveAs("z0_bias_missing6_tow16-23.png");


    // -----------------------------------------------------------------------------------------------------------------


    canvas[7] = new TCanvas("c7","c7",800,600);
    canvas[7]->SetGridy(kTRUE);
    canvas[7]->SetFillColor(kWhite);

    leg[7] = new TLegend(0.768844,0.543706,0.977387,0.963287);

    for (int itow=0; itow<8; ++itow) {
      g_eta_bias_miss[2][itow] = new TGraphErrors(20,eta_bin[itow],eta_bias_m7[itow],eta_bin_width[itow],eta_dummy[itow]);
      g_eta_bias_miss[2][itow]->SetTitle("Missing layer 7");
      g_eta_bias_miss[2][itow]->SetMarkerStyle(20);
      g_eta_bias_miss[2][itow]->SetMarkerColor(colors[itow]);
      g_eta_bias_miss[2][itow]->SetLineColor(colors[itow]);

      if ( itow == 0 ){
	g_eta_bias_miss[2][itow]->Draw("AP");
	g_eta_bias_miss[2][itow]->SetMinimum(-1.e-2);
	g_eta_bias_miss[2][itow]->SetMaximum(1.e-2);
	g_eta_bias_miss[2][itow]->GetXaxis()->SetTitle("#eta_{bin}");
	g_eta_bias_miss[2][itow]->GetYaxis()->SetTitle("#LTz_{0}^{gen}-z_{0}^{fit}#GT [cm]");
	g_eta_bias_miss[2][itow]->GetYaxis()->SetTitleOffset(1.2);

      }
      else {
	g_eta_bias_miss[2][itow]->Draw("P");
      }

      leg[7]->AddEntry(g_eta_bias_miss[0][itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[7]->Draw();
    canvas[7]->Update();
    canvas[7]->SaveAs("z0_bias_missing7_tow16-23.png");


    // -----------------------------------------------------------------------------------------------------------------


    canvas[8] = new TCanvas("c8","c8",800,600);
    canvas[8]->SetGridy(kTRUE);
    canvas[8]->SetFillColor(kWhite);

    leg[8] = new TLegend(0.768844,0.543706,0.977387,0.963287);

    for (int itow=8; itow<16; ++itow) {
      g_eta_bias_miss[0][itow] = new TGraphErrors(20,eta_bin[itow],eta_bias_m5[itow],eta_bin_width[itow],eta_dummy[itow]);
      g_eta_bias_miss[0][itow]->SetTitle("Missing layer 5");
      g_eta_bias_miss[0][itow]->SetMarkerStyle(20);
      g_eta_bias_miss[0][itow]->SetMarkerColor(colors[itow-8]);
      g_eta_bias_miss[0][itow]->SetLineColor(colors[itow-8]);

      if ( itow == 8 ){
	g_eta_bias_miss[0][itow]->Draw("AP");
	g_eta_bias_miss[0][itow]->SetMinimum(-1.e-2);
	g_eta_bias_miss[0][itow]->SetMaximum(1.e-2);
	g_eta_bias_miss[0][itow]->GetXaxis()->SetTitle("#eta_{bin}");
	g_eta_bias_miss[0][itow]->GetYaxis()->SetTitle("#LTz_{0}^{gen}-z_{0}^{fit}#GT [cm]");
	g_eta_bias_miss[0][itow]->GetYaxis()->SetTitleOffset(1.2);

      }
      else {
	g_eta_bias_miss[0][itow]->Draw("P");
      }

      leg[8]->AddEntry(g_eta_bias_miss[0][itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[8]->Draw();
    canvas[8]->Update();
    canvas[8]->SaveAs("z0_bias_missing5_tow24-31.png");


    // -----------------------------------------------------------------------------------------------------------------


    canvas[9] = new TCanvas("c9","c9",800,600);
    canvas[9]->SetGridy(kTRUE);
    canvas[9]->SetFillColor(kWhite);
    
    leg[9] = new TLegend(0.768844,0.543706,0.977387,0.963287);
 
    for (int itow=8; itow<16; ++itow) {
      g_eta_bias_miss[1][itow] = new TGraphErrors(20,eta_bin[itow],eta_bias_m6[itow],eta_bin_width[itow],eta_dummy[itow]);
      g_eta_bias_miss[1][itow]->SetTitle("Missing layer 6");
      g_eta_bias_miss[1][itow]->SetMarkerStyle(20);
      g_eta_bias_miss[1][itow]->SetMarkerColor(colors[itow-8]);
      g_eta_bias_miss[1][itow]->SetLineColor(colors[itow-8]);
 
      if ( itow == 8 ){
 	g_eta_bias_miss[1][itow]->Draw("AP");
	g_eta_bias_miss[1][itow]->SetMinimum(-1.e-2);
	g_eta_bias_miss[1][itow]->SetMaximum(1.e-2);
 	g_eta_bias_miss[1][itow]->GetXaxis()->SetTitle("#eta_{bin}");
 	g_eta_bias_miss[1][itow]->GetYaxis()->SetTitle("#LTz_{0}^{gen}-z_{0}^{fit}#GT [cm]");
 	g_eta_bias_miss[1][itow]->GetYaxis()->SetTitleOffset(1.2);
 
      }
      else {
 	g_eta_bias_miss[1][itow]->Draw("P");
      }
 
      leg[9]->AddEntry(g_eta_bias_miss[0][itow],Form("Tower %d", itow+16),"P");
    
    }
 
    leg[9]->Draw();
    canvas[9]->Update();
    canvas[9]->SaveAs("z0_bias_missing6_tow24-31.png");
 
 
    // -----------------------------------------------------------------------------------------------------------------


    canvas[10] = new TCanvas("c10","c10",800,600);
    canvas[10]->SetGridy(kTRUE);
    canvas[10]->SetFillColor(kWhite);
 
    leg[10] = new TLegend(0.768844,0.543706,0.977387,0.963287);
 
    for (int itow=8; itow<16; ++itow) {
      g_eta_bias_miss[2][itow] = new TGraphErrors(20,eta_bin[itow],eta_bias_m7[itow],eta_bin_width[itow],eta_dummy[itow]);
      g_eta_bias_miss[2][itow]->SetTitle("Missing layer 7");
      g_eta_bias_miss[2][itow]->SetMarkerStyle(20);
      g_eta_bias_miss[2][itow]->SetMarkerColor(colors[itow-8]);
      g_eta_bias_miss[2][itow]->SetLineColor(colors[itow-8]);
 
      if ( itow == 8 ){
 	g_eta_bias_miss[2][itow]->Draw("AP");
	g_eta_bias_miss[2][itow]->SetMinimum(-1.e-2);
	g_eta_bias_miss[2][itow]->SetMaximum(1.e-2);
 	g_eta_bias_miss[2][itow]->GetXaxis()->SetTitle("#eta_{bin}");
 	g_eta_bias_miss[2][itow]->GetYaxis()->SetTitle("#LTz_{0}^{gen}-z_{0}^{fit}#GT [cm]");
 	g_eta_bias_miss[2][itow]->GetYaxis()->SetTitleOffset(1.2);
 
      }
      else {
 	g_eta_bias_miss[2][itow]->Draw("P");
      }
 
      leg[10]->AddEntry(g_eta_bias_miss[0][itow],Form("Tower %d", itow+16),"P");
    
    }
 
    leg[10]->Draw();
    canvas[10]->Update();
    canvas[10]->SaveAs("z0_bias_missing7_tow24-31.png");

  }
  //====================================================================================================================
  //====================================================================================================================
  else if ( flag == 4 ) {

    canvas[11] = new TCanvas("c11","c11",800,600);
    canvas[11]->SetGridy(kTRUE);
    canvas[11]->SetFillColor(kWhite);

    leg[11] = new TLegend(0.66206, 0.449301,0.871859,0.872378);

    for (int itow=0; itow<8; ++itow) {
      g_eta_sigma_miss[0][itow] = new TGraphErrors(20,eta_bin[itow],eta_sigma_m5[itow],eta_bin_width[itow],eta_dummy[itow]);
      g_eta_sigma_miss[0][itow]->SetTitle("Missing layer 5");
      g_eta_sigma_miss[0][itow]->SetMarkerStyle(20);
      g_eta_sigma_miss[0][itow]->SetMarkerColor(colors[itow]);
      g_eta_sigma_miss[0][itow]->SetLineColor(colors[itow]);

      if ( itow == 0 ){
	g_eta_sigma_miss[0][itow]->Draw("AP");
	g_eta_sigma_miss[0][itow]->SetMinimum(0.);
	g_eta_sigma_miss[0][itow]->SetMaximum(0.5);
	g_eta_sigma_miss[0][itow]->GetXaxis()->SetTitle("#eta_{bin}");
	g_eta_sigma_miss[0][itow]->GetYaxis()->SetTitle("#sigma_{z_{0}^{gen}-z_{0}^{fit}} [cm]");
	g_eta_sigma_miss[0][itow]->GetYaxis()->SetTitleOffset(1.2);

      }
      else {
	g_eta_sigma_miss[0][itow]->Draw("P");
      }

      leg[11]->AddEntry(g_eta_sigma_miss[0][itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[11]->Draw();
    canvas[11]->Update();
    canvas[11]->SaveAs("z0_sigma_missing5_tow16-23.png");


    // -----------------------------------------------------------------------------------------------------------------


    canvas[12] = new TCanvas("c12","c12",800,600);
    canvas[12]->SetGridy(kTRUE);
    canvas[12]->SetFillColor(kWhite);

    leg[12] = new TLegend(0.66206,0.449301,0.871859,0.872378);

    for (int itow=0; itow<8; ++itow) {
      g_eta_sigma_miss[0][itow] = new TGraphErrors(20,eta_bin[itow],eta_sigma_m6[itow],eta_bin_width[itow],eta_dummy[itow]);
      g_eta_sigma_miss[0][itow]->SetTitle("Missing layer 6");
      g_eta_sigma_miss[0][itow]->SetMarkerStyle(20);
      g_eta_sigma_miss[0][itow]->SetMarkerColor(colors[itow]);
      g_eta_sigma_miss[0][itow]->SetLineColor(colors[itow]);

      if ( itow == 0 ){
	g_eta_sigma_miss[0][itow]->Draw("AP");
	g_eta_sigma_miss[0][itow]->SetMinimum(0.);
	g_eta_sigma_miss[0][itow]->SetMaximum(0.5);
	g_eta_sigma_miss[0][itow]->GetXaxis()->SetTitle("#eta_{bin}");
	g_eta_sigma_miss[0][itow]->GetYaxis()->SetTitle("#sigma_{z_{0}^{gen}-z_{0}^{fit}} [cm]");
	g_eta_sigma_miss[0][itow]->GetYaxis()->SetTitleOffset(1.2);

      }
      else {
	g_eta_sigma_miss[0][itow]->Draw("P");
      }

      leg[12]->AddEntry(g_eta_sigma_miss[0][itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[12]->Draw();
    canvas[12]->Update();
    canvas[12]->SaveAs("z0_sigma_missing6_tow16-23.png");


    // -----------------------------------------------------------------------------------------------------------------


    canvas[13] = new TCanvas("c13","c13",800,600);
    canvas[13]->SetGridy(kTRUE);
    canvas[13]->SetFillColor(kWhite);

    leg[13] = new TLegend(0.66206,0.449301,0.871859,0.872378);

    for (int itow=0; itow<8; ++itow) {
      g_eta_sigma_miss[0][itow] = new TGraphErrors(20,eta_bin[itow],eta_sigma_m7[itow],eta_bin_width[itow],eta_dummy[itow]);
      g_eta_sigma_miss[0][itow]->SetTitle("Missing layer 7");
      g_eta_sigma_miss[0][itow]->SetMarkerStyle(20);
      g_eta_sigma_miss[0][itow]->SetMarkerColor(colors[itow]);
      g_eta_sigma_miss[0][itow]->SetLineColor(colors[itow]);

      if ( itow == 0 ){
	g_eta_sigma_miss[0][itow]->Draw("AP");
	g_eta_sigma_miss[0][itow]->SetMinimum(0.);
	g_eta_sigma_miss[0][itow]->SetMaximum(0.5);
	g_eta_sigma_miss[0][itow]->GetXaxis()->SetTitle("#eta_{bin}");
	g_eta_sigma_miss[0][itow]->GetYaxis()->SetTitle("#sigma_{z_{0}^{gen}-z_{0}^{fit}} [cm]");
	g_eta_sigma_miss[0][itow]->GetYaxis()->SetTitleOffset(1.2);

      }
      else {
	g_eta_sigma_miss[0][itow]->Draw("P");
      }

      leg[13]->AddEntry(g_eta_sigma_miss[0][itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[13]->Draw();
    canvas[13]->Update();
    canvas[13]->SaveAs("z0_sigma_missing7_tow16-23.png");


    // -----------------------------------------------------------------------------------------------------------------


    canvas[14] = new TCanvas("c14","c14",800,600);
    canvas[14]->SetGridy(kTRUE);
    canvas[14]->SetFillColor(kWhite);

    leg[14] = new TLegend(0.66206,0.449301,0.871859,0.872378);

    for (int itow=8; itow<16; ++itow) {
      g_eta_sigma_miss[0][itow] = new TGraphErrors(20,eta_bin[itow],eta_sigma_m5[itow],eta_bin_width[itow],eta_dummy[itow]);
      g_eta_sigma_miss[0][itow]->SetTitle("Missing layer 5");
      g_eta_sigma_miss[0][itow]->SetMarkerStyle(20);
      g_eta_sigma_miss[0][itow]->SetMarkerColor(colors[itow-8]);
      g_eta_sigma_miss[0][itow]->SetLineColor(colors[itow-8]);

      if ( itow == 8 ){
	g_eta_sigma_miss[0][itow]->Draw("AP");
	g_eta_sigma_miss[0][itow]->SetMinimum(0.);
	g_eta_sigma_miss[0][itow]->SetMaximum(0.5);
	g_eta_sigma_miss[0][itow]->GetXaxis()->SetTitle("#eta_{bin}");
	g_eta_sigma_miss[0][itow]->GetYaxis()->SetTitle("#sigma_{z_{0}^{gen}-z_{0}^{fit}} [cm]");
	g_eta_sigma_miss[0][itow]->GetYaxis()->SetTitleOffset(1.2);

      }
      else {
	g_eta_sigma_miss[0][itow]->Draw("P");
      }

      leg[14]->AddEntry(g_eta_sigma_miss[0][itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[14]->Draw();
    canvas[14]->Update();
    canvas[14]->SaveAs("z0_sigma_missing5_tow24-31.png");


    // -----------------------------------------------------------------------------------------------------------------


    canvas[15] = new TCanvas("c15","c15",800,600);
    canvas[15]->SetGridy(kTRUE);
    canvas[15]->SetFillColor(kWhite);

    leg[15] = new TLegend(0.66206,0.449301,0.871859,0.872378);

    for (int itow=8; itow<16; ++itow) {
      g_eta_sigma_miss[0][itow] = new TGraphErrors(20,eta_bin[itow],eta_sigma_m6[itow],eta_bin_width[itow],eta_dummy[itow]);
      g_eta_sigma_miss[0][itow]->SetTitle("Missing layer 6");
      g_eta_sigma_miss[0][itow]->SetMarkerStyle(20);
      g_eta_sigma_miss[0][itow]->SetMarkerColor(colors[itow-8]);
      g_eta_sigma_miss[0][itow]->SetLineColor(colors[itow-8]);

      if ( itow == 8 ){
	g_eta_sigma_miss[0][itow]->Draw("AP");
	g_eta_sigma_miss[0][itow]->SetMinimum(0.);
	g_eta_sigma_miss[0][itow]->SetMaximum(0.5);
	g_eta_sigma_miss[0][itow]->GetXaxis()->SetTitle("#eta_{bin}");
	g_eta_sigma_miss[0][itow]->GetYaxis()->SetTitle("#sigma_{z_{0}^{gen}-z_{0}^{fit}} [cm]");
	g_eta_sigma_miss[0][itow]->GetYaxis()->SetTitleOffset(1.2);

      }
      else {
	g_eta_sigma_miss[0][itow]->Draw("P");
      }

      leg[15]->AddEntry(g_eta_sigma_miss[0][itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[15]->Draw();
    canvas[15]->Update();
    canvas[15]->SaveAs("z0_sigma_missing6_tow24-31.png");


    // -----------------------------------------------------------------------------------------------------------------


    canvas[16] = new TCanvas("c16","c16",800,600);
    canvas[16]->SetGridy(kTRUE);
    canvas[16]->SetFillColor(kWhite);

    leg[16] = new TLegend(0.66206,0.449301,0.871859,0.872378);

    for (int itow=8; itow<16; ++itow) {
      g_eta_sigma_miss[0][itow] = new TGraphErrors(20,eta_bin[itow],eta_sigma_m7[itow],eta_bin_width[itow],eta_dummy[itow]);
      g_eta_sigma_miss[0][itow]->SetTitle("Missing layer 7");
      g_eta_sigma_miss[0][itow]->SetMarkerStyle(20);
      g_eta_sigma_miss[0][itow]->SetMarkerColor(colors[itow-8]);
      g_eta_sigma_miss[0][itow]->SetLineColor(colors[itow-8]);

      if ( itow == 8 ){
	g_eta_sigma_miss[0][itow]->Draw("AP");
	g_eta_sigma_miss[0][itow]->SetMinimum(0.);
	g_eta_sigma_miss[0][itow]->SetMaximum(0.5);
	g_eta_sigma_miss[0][itow]->GetXaxis()->SetTitle("#eta_{bin}");
	g_eta_sigma_miss[0][itow]->GetYaxis()->SetTitle("#sigma_{z_{0}^{gen}-z_{0}^{fit}} [cm]");
	g_eta_sigma_miss[0][itow]->GetYaxis()->SetTitleOffset(1.2);

      }
      else {
	g_eta_sigma_miss[0][itow]->Draw("P");
      }

      leg[16]->AddEntry(g_eta_sigma_miss[0][itow],Form("Tower %d", itow+16),"P");
    
    }

    leg[16]->Draw();
    canvas[16]->Update();
    canvas[16]->SaveAs("z0_sigma_missing7_tow24-31.png");



  }

}
