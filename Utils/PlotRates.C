#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <iomanip>
#include <cmath>
#include <fstream>


/*
  ROOT macro for stub/cluster rates visualization, to be used with a file
  created with ./AM_ana -c rates 

  Use:

  root[1]-> .L MakeRatesGraphs.C

  Then you can choose between the following methods 
  knowing that filename denotes the extracted ROOT file 
  containing the sectors definition (created by the 
  sectorMaker program) 


  //
  // Rate plots and maps
  //

  Solution 1 do all the plots with complete rates
  root[2]-> do_rates_plots(std::string filename) 

  // Then you can make rate plots for only a type of stubs.
  // Stubs could indeed be put in two categories
  // 
  // Primaries
  // Secondaries
  // Fakes
  // 

  In order to do rate plots only for primaries do:
  root[2]-> do_the_plots(std::string filename,1,0,0) 

  In order to do rate plots only for secondaries do:
  root[2]-> do_the_plots(std::string filename,0,1,0) 

  In order to do rate plots only for fakes do:
  root[2]-> do_the_plots(std::string filename,0,1,0) 

  // You can also make plots only for a layer/disk only

  Solution 2 do the plot for barrel layer i
  root[2]-> do_layer_map(std::string filename,int layer,p,s,f)

  Solution 3 do the plot for endcap disk i
  root[2]-> do_disk_map(std::string filename,int disk,p,s,f)

  Finally you can do 1D average plots, with the rates averaged over rings/ladders:

  root[2]-> do_layers_1D_map(std::string filename, bool inhz)
  root[2]-> do_disks_1D_map(std::string filename, bool inhz)

  if inhz is set to 0, the rate is given in stubs/module/BX, if set to 1 the rate is 
  given in MHz/cm²

  The rates obtained on this plots can be explicitely using:

  root[2]-> print_ENDCAP_rate(std::string filename) // Ring in the disks
  root[2]-> print_TIB_rate(std::string filename)    // Ladder in TIB
  root[2]-> print_TOB_rate(std::string filename)    // Ladder in TOB


  //
  // Module/chip plots
  //

  // You can also plot the rate per chip for a given module

  root[2]-> do_module_map(std::string filename,int layer, int ladder, int module)

  with layer numbered as:

  5-10 : Barrel
  11-17: Endcap +Z
  18-24: Endcap -Z

  Information about this macro may be found in the following page:

  http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto (STEP III)

  Author: viret@in2p3_dot_fr
  Date:        23/05/2013
  Last update: 17/01/2014
*/


void fill_histo(TH2F *hist,float val[], double nlad, double dec, int disk, int ring,float max)
{
  float scale = 8*atan(1.)/nlad;

  float phi_sec = 0.;

  for (int i=0;i<int(nlad);++i)  
  {
    phi_sec = (i+.5-dec)*scale;

    if (phi_sec<0.) phi_sec+=8*atan(1.); 

    int idx = 10000*disk + 100*ring + i;
    hist->Fill(phi_sec,0.5+ring,val[idx]);
  }
 
  hist->Fill(0.5*scale,25.5,0);
  hist->Fill(1.5*scale,25.5,max);
}

// Rate plot for the barrel layer (1D official plot)

void do_layers_1D_map(std::string filename, bool inhz, int pu, int tilt)
{
  // First get the data
  // by merging all the available files

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // Module dims are extracted from TkLayout info:
  // http://cms-tklayout.web.cern.ch/cms-tklayout/layouts/reference/OT613_200_IT4025/layout.html

  // Active surface of modules (in cm²)

  double PS_surf = 45.07;
  double SS_surf = 91.90;

  double n_mod_surf_barrel[6];

  for (int i=0;i<3;++i) n_mod_surf_barrel[i] = PS_surf;
  for (int i=3;i<6;++i) n_mod_surf_barrel[i] = SS_surf;
 
  // Bit number per stub

  float bit_per_s = 1;

  float count_s[600];
  float count_c[600];

 
  double n_lad_flat[6]  = {16,24,34,48,62,76};
  double n_mod_flat[6]  = {63,55,54,24,24,24};
  
  double n_lad_tilt[6]  = {18,26,36,48,60,78};
  double n_mod_tilt[6]  = {31,35,39,24,24,24};
  
  double n_lad_barrel[6];
  double n_mod_barrel[6];
  
  if (tilt==1)
  {
    for (int i=0;i<6;++i) n_lad_barrel[i] = n_lad_tilt[i];
    for (int i=0;i<6;++i) n_mod_barrel[i] = n_mod_tilt[i];
  }
  else
  {
    for (int i=0;i<6;++i) n_lad_barrel[i] = n_lad_flat[i];
    for (int i=0;i<6;++i) n_mod_barrel[i] = n_mod_flat[i];
  }

  TChain *newtree = new TChain("L1Rates");
    
  newtree->Add(filename.c_str());
  newtree->SetBranchAddress("STUB_b_l_rates",   &count_s);
  newtree->SetBranchAddress("CLUS_b_l_rates",   &count_c);
  newtree->GetEntry(0);

  int idx;
  float maxsval = 0;
  float maxcval = 0;
  float maxrval = 0;
  float rates;
  float ratec;
  float rater;

  for (int i=0;i<6;++i)
  { 
    for (int j=i*100;j<(i+1)*100;++j)
    { 
      rates = count_s[j]/n_lad_barrel[i];
      ratec = count_c[j]/n_lad_barrel[i];
      rater = 0;
      if (rates!=0) rater = ratec/rates;

      if (inhz==true) 
      {	
	rates /= n_mod_surf_barrel[i];
	ratec /= n_mod_surf_barrel[i];
	rates *= bit_per_s*40.; // Give the value in MHz/cm2
	ratec *= bit_per_s*40.; // Give the value in MHz/cm2
      }

      if (rates>1000) continue;	      

      if (rates>maxsval) maxsval = rates;
      if (ratec>maxcval) maxcval = ratec;
      if (rater>maxrval) maxrval = rater;
    }
  }


  TH2F *cadres  = new TH2F("c_s","c_s",1000,-100.,100.,100,0.,1.1*maxsval);
  TH2F *cadrec  = new TH2F("c_c","c_c",1000,-100.,100.,100,0.,1.1*maxcval);
  TH2F *cadrer  = new TH2F("c_r","c_r",1000,-100.,100.,100,0.,1.1*maxrval);

  std::vector<TH2F*> stub_plots; 
  std::vector<TH2F*> clus_plots; 
  std::vector<TH2F*> rat_plots; 

  for (int i=0;i<6;++i)
  {
    TH2F *cadre_stub = new TH2F("","",20*n_mod_barrel[i],-100.,100.,100,0.,1.1*maxsval);
    TH2F *cadre_clus = new TH2F("","",20*n_mod_barrel[i],-100.,100.,100,0.,1.1*maxcval);
    TH2F *cadre_rati = new TH2F("","",20*n_mod_barrel[i],-100.,100.,100,0.,1.1*maxrval);

    stub_plots.push_back(cadre_stub);
    clus_plots.push_back(cadre_clus);
    rat_plots.push_back(cadre_rati);
  }

  float bin_w;

  for (int i=0;i<6;++i)
  {
    bin_w = 200./(20.*n_mod_barrel[i]);

    for (int j=0;j<n_mod_barrel[i];++j)
    {
      rates = count_s[100*i+j]/n_lad_barrel[i];
      ratec = count_c[100*i+j]/n_lad_barrel[i];
      if (rates!=0) rater = ratec/rates;
      if (rates==0) continue;
        
      if (inhz) 
      {	
	rates /= n_mod_surf_barrel[i];
	rates *= bit_per_s*40.; // Give the value in MHz/cm2
	
	ratec /= n_mod_surf_barrel[i];
	ratec *= bit_per_s*40.; // Give the value in MHz/cm2
      }

      stub_plots[i]->Fill(-100+20*bin_w*(j+0.5),rates);

      clus_plots[i]->Fill(-100+20*bin_w*(j+0.5),ratec);

      rat_plots[i]->Fill(-100+20*bin_w*(j+0.5),rater);
    }
  }


  TCanvas *c1 = new TCanvas("c1","Layer map",0,0,640,480);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c1->Range(-1.,-1.,1.,1.);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetLeftMargin(0.07692308);
  c1->SetTopMargin(0.07124352);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderMode(0);
  
  
  cadres->GetXaxis()->SetTitle("Module index along Z");
  cadres->GetYaxis()->SetTitle("Stub rate (number of stubs/module/BX)"); 

  if (inhz)  cadres->GetYaxis()->SetTitle("Stub rate (in MHz/cm^{2})"); 

  cadres->GetXaxis()->SetLabelFont(42);
  cadres->GetXaxis()->SetLabelSize(0.035);
  cadres->GetXaxis()->SetTitleSize(0.035);
  cadres->GetXaxis()->SetTitleFont(42);
  cadres->GetYaxis()->SetLabelFont(42);
  cadres->GetYaxis()->SetLabelOffset(0.004);
  cadres->GetYaxis()->SetLabelSize(0.035);
  cadres->GetYaxis()->SetTitleSize(0.035);
  cadres->GetYaxis()->SetTitleOffset(0.83);
  cadres->Draw();

  int color[6]  = {24,20,27,34,28,25};

  for (int i=0;i<6;++i)
  {
    stub_plots[i]->SetMarkerStyle(color[i]);
    stub_plots[i]->SetMarkerSize(1.4);
    stub_plots[i]->Draw("same");
  }

  TLegend *leg = new TLegend(0.63,0.75,0.81,0.98);
  char buffer[100];

  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(stub_plots[0],"TBPS layer 1","p");
  leg->AddEntry(stub_plots[1],"TBPS layer 2","p");
  leg->AddEntry(stub_plots[2],"TBPS layer 3","p");
  leg->AddEntry(stub_plots[3],"TB2S layer 1","p");
  leg->AddEntry(stub_plots[4],"TB2S layer 2","p");
  leg->AddEntry(stub_plots[5],"TB2S layer 3","p");
  leg->Draw();

  sprintf (buffer, "#sqrt{s}=14TeV, %d PU", pu);

  TLatex Tl;
  Tl.SetTextSize(0.04);
  Tl.DrawLatex(-101, 1.11*maxsval, "CMS Phase-2 Simulation");

  TLatex Tl2;
  Tl2.SetTextSize(0.03);
  Tl2.SetTextFont(52);
  Tl2.DrawLatex(-15, 1.11*maxsval, buffer);

  c1->Modified();
  c1->Update();
      
  sprintf (buffer, "Barrel_%d_%d.pdf", pu, inhz);
  c1->Print(buffer);
  sprintf (buffer, "Barrel_%d_%d.C", pu, inhz);
  c1->SaveSource(buffer);

  
  TCanvas *c2 = new TCanvas("c2","Layer map",0,0,640,480);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c2->Range(-1.,-1.,1.,1.);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetBorderSize(2);
  c2->SetGridx();
  c2->SetGridy();
  c2->SetLeftMargin(0.07692308);
  c2->SetTopMargin(0.07124352);
  c2->SetFrameBorderMode(0);
  c2->SetFrameBorderMode(0);
  
  cadrec->GetXaxis()->SetTitle("Module index along Z");
  cadrec->GetYaxis()->SetTitle("Cluster rate (number of clusters/module/BX)"); 

  if (inhz)  cadrec->GetYaxis()->SetTitle("Cluster rate (in MHz/cm^{2})"); 

  cadrec->GetXaxis()->SetLabelFont(42);
  cadrec->GetXaxis()->SetLabelSize(0.035);
  cadrec->GetXaxis()->SetTitleSize(0.035);
  cadrec->GetXaxis()->SetTitleFont(42);
  cadrec->GetYaxis()->SetLabelFont(42);
  cadrec->GetYaxis()->SetLabelOffset(0.004);
  cadrec->GetYaxis()->SetLabelSize(0.035);
  cadrec->GetYaxis()->SetTitleSize(0.035);
  cadrec->GetYaxis()->SetTitleOffset(0.83);
  cadrec->Draw();

  for (int i=0;i<6;++i)
  {
    clus_plots[i]->SetMarkerStyle(color[i]);
    clus_plots[i]->SetMarkerSize(1.4);
    clus_plots[i]->Draw("same");
  }


  leg = new TLegend(0.63,0.75,0.81,0.98);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(clus_plots[0],"TBPS layer 1","p");
  leg->AddEntry(clus_plots[1],"TBPS layer 2","p");
  leg->AddEntry(clus_plots[2],"TBPS layer 3","p");
  leg->AddEntry(clus_plots[3],"TB2S layer 1","p");
  leg->AddEntry(clus_plots[4],"TB2S layer 2","p");
  leg->AddEntry(clus_plots[5],"TB2S layer 3","p");
  leg->Draw();


  sprintf (buffer, "#sqrt{s}=14TeV, %d PU", pu);

  Tl.SetTextSize(0.04);
  Tl.DrawLatex(-101, 1.11*maxcval, "CMS Phase-2 Simulation");

  Tl2.SetTextSize(0.03);
  Tl2.SetTextFont(52);
  Tl2.DrawLatex(-15, 1.11*maxcval, buffer);

  c2->Modified();
  c2->Update();
      
  sprintf (buffer, "Barrel_clus_%d_%d.pdf", pu, inhz);
  c2->Print(buffer);
  sprintf (buffer, "Barrel_clus_%d_%d.C", pu, inhz);
  c2->SaveSource(buffer);

  TCanvas *c3 = new TCanvas("c3","Layer map",0,0,640,480);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c3->Range(-1.,-1.,1.,1.);
  c3->SetFillColor(0);
  c3->SetBorderMode(0);
  c3->SetBorderSize(2);
  c3->SetGridx();
  c3->SetGridy();
  c3->SetLeftMargin(0.07692308);
  c3->SetTopMargin(0.07124352);
  c3->SetFrameBorderMode(0);
  c3->SetFrameBorderMode(0);
  
  
  cadrer->GetXaxis()->SetTitle("Module index along Z");
  cadrer->GetYaxis()->SetTitle("Cluster over stub ratio"); 

  cadrer->GetXaxis()->SetLabelFont(42);
  cadrer->GetXaxis()->SetLabelSize(0.035);
  cadrer->GetXaxis()->SetTitleSize(0.035);
  cadrer->GetXaxis()->SetTitleFont(42);
  cadrer->GetYaxis()->SetLabelFont(42);
  cadrer->GetYaxis()->SetLabelOffset(0.004);
  cadrer->GetYaxis()->SetLabelSize(0.035);
  cadrer->GetYaxis()->SetTitleSize(0.035);
  cadrer->GetYaxis()->SetTitleOffset(0.83);
  cadrer->Draw();

  for (int i=0;i<6;++i)
  {
    rat_plots[i]->SetMarkerStyle(color[i]);
    rat_plots[i]->SetMarkerSize(1.4);
    rat_plots[i]->Draw("same");
  }


  leg = new TLegend(0.63,0.75,0.81,0.98);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(rat_plots[0],"TBPS layer 1","p");
  leg->AddEntry(rat_plots[1],"TBPS layer 2","p");
  leg->AddEntry(rat_plots[2],"TBPS layer 3","p");
  leg->AddEntry(rat_plots[3],"TB2S layer 1","p");
  leg->AddEntry(rat_plots[4],"TB2S layer 2","p");
  leg->AddEntry(rat_plots[5],"TB2S layer 3","p");
  leg->Draw();


  sprintf (buffer, "#sqrt{s}=14TeV, %d PU", pu);

  Tl.SetTextSize(0.04);
  Tl.DrawLatex(-101, 1.11*maxrval, "CMS Phase-2 Simulation");

  Tl2.SetTextSize(0.03);
  Tl2.SetTextFont(52);
  Tl2.DrawLatex(-15, 1.11*maxrval, buffer);

  c3->Modified();
  c3->Update();
      
  sprintf (buffer, "Barrel_ratio_%d_%d.pdf", pu, inhz);
  c3->Print(buffer);
  sprintf (buffer, "Barrel_ratio_%d_%d.C", pu, inhz);
  c3->SaveSource(buffer);
}

// Rate plot for the disks (1D commercial plot)

void do_disks_1D_map(std::string filename, bool inhz, int pu, int tilt)
{
  // First get the data
  // by merging all the available files

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // Module dims are extracted from TkLayout info:
  // http://cms-tklayout.web.cern.ch/cms-tklayout/layouts/reference/OT613_200_IT4025/layout.html

  // Active surface of modules (in cm²)

  double PS_surf = 45.07;
  double SS_surf = 91.90;

  double n_mod_surf_endcap[15];

  for (int i=0;i<9;++i) n_mod_surf_endcap[i] = PS_surf;
  for (int i=9;i<15;++i) n_mod_surf_endcap[i] = SS_surf;
 
  // Bit number per stub

  float bit_per_s = 1;

  float count_s[1500];
  float count_c[1500];

  double n_mod_endcap[15]   = {20,24,24,28,32,32,36,40,40,44,52,60,64,72,76};
  float radii[15]           = {25,30,33,38,41,46,48,53,56,63,71,81,88,98,104};

 TChain *newtree = new TChain("L1Rates");
    
  newtree->Add(filename.c_str());

  newtree->SetBranchAddress("STUB_e_l_rates",   &count_s);
  newtree->SetBranchAddress("CLUS_e_l_rates",   &count_c);
  newtree->GetEntry(0);

  int idx;
  float maxsval = 0;
  float maxcval = 0;
  float maxrval = 0;
  float rates;
  float ratec;
  float rater;  

  for (int i=7;i<12;++i)
  {   
    for (int j=i*100;j<i*100+15;++j)
    { 
      rates = count_s[j]/n_mod_endcap[j-100*i];
      ratec = count_c[j]/n_mod_endcap[j-100*i];
      rater = 0;
      if (rates !=0) rater = ratec/rates;

      if (inhz==true) 
      {	
	rates /= n_mod_surf_endcap[j-100*i];
	rates *= bit_per_s*40.; // Give the value in MHz/cm2

	ratec /= n_mod_surf_endcap[j-100*i];
	ratec *= bit_per_s*40.; // Give the value in MHz/cm2
      }
      
      if (rates>maxsval) maxsval = rates;
      if (ratec>maxcval) maxcval = ratec;
      if (rater>maxrval) maxrval = rater;
    }
  }

  TH2F *cadres  = new TH2F("c_s","c_s",200,20.,120.,100,0.,1.1*maxsval);
  TH2F *cadrer  = new TH2F("c_r","c_r",200,20.,120.,100,0.,1.1*maxrval); 
  TH2F *cadrec  = new TH2F("c_c","c_c",200,20.,120.,100,0.,1.1*maxcval);

  std::vector<TH2F*> stub_plots; 
  std::vector<TH2F*> clus_plots; 
  std::vector<TH2F*> rat_plots; 

  for (int i=0;i<5;++i)
  {
    TH2F *cadre_stub = new TH2F("","",200,20.,120.,100,0.,1.1*maxsval);
    TH2F *cadre_clus = new TH2F("","",200,20.,120.,100,0.,1.1*maxcval);
    TH2F *cadre_rati = new TH2F("","",200,20.,120.,100,0.,1.1*maxrval);

    stub_plots.push_back(cadre_stub);
    clus_plots.push_back(cadre_clus);
    rat_plots.push_back(cadre_rati);
  }

  float bin_w;

  int mv;

  for (int i=7;i<12;++i)
  {
    for (int j=i*100;j<i*100+15;++j)
    {
      mv = j%100;
      if (i>=9) mv+=3;
      if (mv>=15) continue;

      rates = count_s[j]/n_mod_endcap[mv];
      ratec = count_c[j]/n_mod_endcap[mv];
      if (rates!=0) rater = count_c[j]/count_s[j];

      if (inhz) 
      {	
	rates /= n_mod_surf_endcap[mv];
	rates *= bit_per_s*40.; // Give the value in MHz/cm2

	ratec /= n_mod_surf_endcap[mv];
	ratec *= bit_per_s*40.; // Give the value in MHz/cm2
      }

      if (rates !=0) stub_plots[i-7]->Fill(radii[mv],rates);
      if (ratec !=0) clus_plots[i-7]->Fill(radii[mv],ratec);
      if (rater !=0) rat_plots[i-7]->Fill(radii[mv],rater);
    }
  }
  
  TCanvas *c1 = new TCanvas("c1","Layer map",0,0,640,480);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c1->Range(-1.,-1.,1.,1.);

  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetLeftMargin(0.07692308);
  c1->SetTopMargin(0.07124352);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderMode(0);
  
  cadres->GetXaxis()->SetTitle("r (cm)");
  cadres->GetYaxis()->SetTitle("Stub rate (number of stubs/module/BX)"); 

  if (inhz)  cadres->GetYaxis()->SetTitle("Stub rate (in MHz/cm^{2})"); 

  cadres->GetXaxis()->SetLabelFont(42);
  cadres->GetXaxis()->SetLabelSize(0.035);
  cadres->GetXaxis()->SetTitleSize(0.035);
  cadres->GetXaxis()->SetTitleFont(42);
  cadres->GetYaxis()->SetLabelFont(42);
  cadres->GetYaxis()->SetLabelOffset(0.004);
  cadres->GetYaxis()->SetLabelSize(0.035);
  cadres->GetYaxis()->SetTitleSize(0.035);
  cadres->GetYaxis()->SetTitleOffset(0.83);
  cadres->Draw();
    
  int color[6]  = {24,20,27,34,28,25};
  
  for (int i=0;i<5;++i)
  {
    stub_plots[i]->SetMarkerStyle(color[i]);
    stub_plots[i]->SetMarkerSize(1.4);
    stub_plots[i]->Draw("same");
  }
  
  TLegend *leg = new TLegend(0.45,0.6,0.75,0.86);
  char buffer[100];

  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->AddEntry(stub_plots[0],"TEDD double-disc 1","p");
  leg->AddEntry(stub_plots[1],"TEDD double-disc 2","p");
  leg->AddEntry(stub_plots[2],"TEDD double-disc 3","p");
  leg->AddEntry(stub_plots[3],"TEDD double-disc 4","p");
  leg->AddEntry(stub_plots[4],"TEDD double-disc 5","p");
  leg->Draw();

  sprintf (buffer, "#sqrt{s}=14TeV, %d PU", pu);

  TLatex Tl;
  Tl.SetTextSize(0.04);
  Tl.DrawLatex(20, 1.11*maxsval, "CMS Phase-2 Simulation");

  TLatex Tl2;
  Tl2.SetTextSize(0.03);
  Tl2.SetTextFont(52);
  Tl2.DrawLatex(100, 1.11*maxsval, buffer);

  c1->Modified();
  c1->Update();
      
  sprintf (buffer, "Endcap_%d_%d.pdf", pu, inhz);
  c1->Print(buffer);
  sprintf (buffer, "Endcap_%d_%d.C", pu, inhz);
  c1->SaveSource(buffer);

 
  TCanvas *c2 = new TCanvas("c2","Layer map",0,0,640,480);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c2->Range(-1.,-1.,1.,1.);
    
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetBorderSize(2);
  c2->SetGridx();
  c2->SetGridy();
  c2->SetLeftMargin(0.07692308);
  c2->SetTopMargin(0.07124352);
  c2->SetFrameBorderMode(0);
  c2->SetFrameBorderMode(0);
  
  cadrec->GetXaxis()->SetTitle("r (cm)");
  cadrec->GetYaxis()->SetTitle("Cluster rate (number of clusters/module/BX)"); 

  if (inhz)  cadrec->GetYaxis()->SetTitle("Cluster rate (in MHz/cm^{2})"); 

  cadrec->GetXaxis()->SetLabelFont(42);
  cadrec->GetXaxis()->SetLabelSize(0.035);
  cadrec->GetXaxis()->SetTitleSize(0.035);
  cadrec->GetXaxis()->SetTitleFont(42);
  cadrec->GetYaxis()->SetLabelFont(42);
  cadrec->GetYaxis()->SetLabelOffset(0.004);
  cadrec->GetYaxis()->SetLabelSize(0.035);
  cadrec->GetYaxis()->SetTitleSize(0.035);
  cadrec->GetYaxis()->SetTitleOffset(0.83);
  cadrec->Draw();
  
  for (int i=0;i<5;++i)
  {
    clus_plots[i]->SetMarkerStyle(color[i]);
    clus_plots[i]->SetMarkerSize(1.4);
    clus_plots[i]->Draw("same");
  }
  
  leg = new TLegend(0.45,0.6,0.75,0.86);

  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->AddEntry(clus_plots[0],"TEDD double-disc 1","p");
  leg->AddEntry(clus_plots[1],"TEDD double-disc 2","p");
  leg->AddEntry(clus_plots[2],"TEDD double-disc 3","p");
  leg->AddEntry(clus_plots[3],"TEDD double-disc 4","p");
  leg->AddEntry(clus_plots[4],"TEDD double-disc 5","p");
  leg->Draw();

  sprintf (buffer, "#sqrt{s}=14TeV, %d PU", pu);

  Tl.SetTextSize(0.04);
  Tl.DrawLatex(20, 1.11*maxcval, "CMS Phase-2 Simulation");

  Tl2.SetTextSize(0.03);
  Tl2.SetTextFont(52);
  Tl2.DrawLatex(100, 1.11*maxcval, buffer);

  c2->Modified();
  c2->Update();
      
  sprintf (buffer, "Endcap_clus_%d_%d.pdf", pu, inhz);
  c2->Print(buffer);
  sprintf (buffer, "Endcap_clus_%d_%d.C", pu, inhz);
  c2->SaveSource(buffer);

  TCanvas *c3 = new TCanvas("c3","Layer map",0,0,640,480);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c3->Range(-1.,-1.,1.,1.);
  c3->SetFillColor(0);
  c3->SetBorderMode(0);
  c3->SetBorderSize(2);
  c3->SetGridx();
  c3->SetGridy();
  c3->SetLeftMargin(0.07692308);
  c3->SetTopMargin(0.07124352);
  c3->SetFrameBorderMode(0);
  c3->SetFrameBorderMode(0);
  
  cadrer->GetXaxis()->SetTitle("r (cm)");
  cadrer->GetYaxis()->SetTitle("Cluster over stub ratio"); 
  cadrer->GetXaxis()->SetLabelFont(42);
  cadrer->GetXaxis()->SetLabelSize(0.035);
  cadrer->GetXaxis()->SetTitleSize(0.035);
  cadrer->GetXaxis()->SetTitleFont(42);
  cadrer->GetYaxis()->SetLabelFont(42);
  cadrer->GetYaxis()->SetLabelOffset(0.004);
  cadrer->GetYaxis()->SetLabelSize(0.035);
  cadrer->GetYaxis()->SetTitleSize(0.035);
  cadrer->GetYaxis()->SetTitleOffset(0.83);
  cadrer->Draw();
  
  for (int i=0;i<5;++i)
  {
    rat_plots[i]->SetMarkerStyle(color[i]);
    rat_plots[i]->SetMarkerSize(1.4);
    rat_plots[i]->Draw("same");
  }
  
  leg = new TLegend(0.45,0.6,0.75,0.86);

  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->AddEntry(rat_plots[0],"TEDD double-disc 1","p");
  leg->AddEntry(rat_plots[1],"TEDD double-disc 2","p");
  leg->AddEntry(rat_plots[2],"TEDD double-disc 3","p");
  leg->AddEntry(rat_plots[3],"TEDD double-disc 4","p");
  leg->AddEntry(rat_plots[4],"TEDD double-disc 5","p");
  leg->Draw();

  sprintf (buffer, "#sqrt{s}=14TeV, %d PU", pu);

  Tl.SetTextSize(0.04);
  Tl.DrawLatex(20, 1.11*maxrval, "CMS Phase-2 Simulation");

  Tl2.SetTextSize(0.03);
  Tl2.SetTextFont(52);
  Tl2.DrawLatex(100, 1.11*maxrval, buffer);

  c3->Modified();
  c3->Update();
      
  sprintf (buffer, "Endcap_ratio_%d_%d.pdf", pu, inhz);
  c3->Print(buffer);
  sprintf (buffer, "Endcap_ratio_%d_%d.C", pu, inhz);
  c3->SaveSource(buffer);
}

void do_TDR_plots(std::string filename, int pu, int tilt)
{
  do_disks_1D_map(filename, 0, pu, tilt);
  do_layers_1D_map(filename, 0, pu, tilt);
}


// Plot everything
/*
void do_rates_plots(std::string filename)
{
  do_the_plots(filename,1,1,1);
}



// Same as before, but you can plot only one part of the info
//
// pri = primaries
// sec = secondaries
// fax = fake

void do_the_plots(std::string filename, float pri, float sec, float fak)
{
  for (int i=0;i<5;++i)
  { 
    do_disk_map(filename,i,pri,sec,fak,1);
    do_disk_map(filename,i,pri,sec,fak,2);
    do_disk_map(filename,i,pri,sec,fak,3);

    do_disk_map(filename,i+7,pri,sec,fak,1);
    do_disk_map(filename,i+7,pri,sec,fak,2);
    do_disk_map(filename,i+7,pri,sec,fak,3);
  }

  
  do_disks_1D_map(filename, 0);
  do_disks_1D_map(filename, 1);

  for (int i=5;i<11;++i)
  { 
    do_layer_map(filename,i,pri,sec,fak);
  }

  do_layers_1D_map(filename, 0);
  do_layers_1D_map(filename, 1);
}
*/

// Rate plot for a given barrel layer

void do_layer_map(std::string filename,int layer, float pri, float sec, float fak, int tilt)
{
  // First get the data
  // by merging all the available files

  if (layer<5) return;

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  float count_p[58000];
  float count_f[58000];
  float count_s[58000];
  float count_c[58000];
    

  double n_lad_flat[6]  = {16,24,34,48,62,76};
  double n_mod_flat[6]  = {63,55,54,24,24,24};

  double n_lad_tilt[6]  = {18,26,36,48,60,78};
  double n_mod_tilt[6]  = {31,35,39,24,24,24};
  
  double n_lad_barrel[6];
  double n_mod_barrel[6];
  
  if (tilt==1)
  {
    for (int i=0;i<6;++i) n_lad_barrel[i] = n_lad_tilt[i];
    for (int i=0;i<6;++i) n_mod_barrel[i] = n_mod_tilt[i];
  }
  else
  {
    for (int i=0;i<6;++i) n_lad_barrel[i] = n_lad_flat[i];
    for (int i=0;i<6;++i) n_mod_barrel[i] = n_mod_flat[i];
  }


  TChain *newtree = new TChain("L1Rates");
    
  newtree->Add(filename.c_str());

  newtree->SetBranchAddress("STUB_b_rates_prim",  &count_p);
  newtree->SetBranchAddress("STUB_b_rates_sec",   &count_f);
  newtree->SetBranchAddress("STUB_b_rates_f",     &count_s);
  newtree->SetBranchAddress("CLUS_b_rates",       &count_c);
  newtree->GetEntry(0);

  int idx;
  float maxval  = 0;
  float maxcval = 0;
  float rate = 0;
  int maxmod = 0;
  float ratio = 0;
  float maxratio = 0;

  for (int j=(layer-5)*10000;j<(layer-5)*10000+8000;++j)
  { 
    rate = pri*count_p[j]+fak*count_f[j]+sec*count_s[j];
   
    if (rate!=0.)
    {
      ratio = count_c[j]/rate;
      if (ratio>maxratio) maxratio=ratio;
    }

    if (rate>maxval) 
    {
      maxval = rate;
      maxmod = j;
    }

    if (count_c[j]>maxcval) 
    {
      maxcval = count_c[j];
    }
  }
  
  cout << "max module for layer " << layer << " is " << maxmod << " with "
       << maxval << " stubs/mod/BX" << endl;  
  
  TH2F *cadre  = new TH2F("zz","zz",n_mod_barrel[layer-5],0.,n_mod_barrel[layer-5],n_lad_barrel[layer-5],0.,n_lad_barrel[layer-5]);
  TH2F *cadre2 = new TH2F("zz2","zz2",n_mod_barrel[layer-5]+1,0.,n_mod_barrel[layer-5]+1,n_lad_barrel[layer-5],0.,n_lad_barrel[layer-5]);
  TH2F *cadre3 = new TH2F("zz3","zz3",n_mod_barrel[layer-5]+1,0.,n_mod_barrel[layer-5]+1,n_lad_barrel[layer-5],0.,n_lad_barrel[layer-5]);
  TH2F *cadre4 = new TH2F("zz4","zz4",n_mod_barrel[layer-5]+1,0.,n_mod_barrel[layer-5]+1,n_lad_barrel[layer-5],0.,n_lad_barrel[layer-5]);

  cadre2->Fill(n_mod_barrel[layer-5]+0.5,0.5,0.);
  cadre2->Fill(n_mod_barrel[layer-5]+0.5,0.5,maxval);

  cadre3->Fill(n_mod_barrel[layer-5]+0.5,0.5,0.);
  cadre3->Fill(n_mod_barrel[layer-5]+0.5,0.5,maxcval);

  cadre4->Fill(n_mod_barrel[layer-5]+0.5,0.5,0.);
  cadre4->Fill(n_mod_barrel[layer-5]+0.5,0.5,maxratio);
  
  for (int j=0;j<n_mod_barrel[layer-5];++j)
  {
    for (int i=0;i<n_lad_barrel[layer-5];++i)
    {
      idx = 10000*(layer-5) + 100*i + j;
      rate= pri*count_p[idx]+fak*count_f[idx]+sec*count_s[idx];

      cadre2->Fill(j+0.5,i+0.5,rate);
      cadre3->Fill(j+0.5,i+0.5,count_c[idx]);

      if (rate!=0) cadre4->Fill(j+0.5,i+0.5,count_c[idx]/rate);
    }
  }
 
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  
  TCanvas *c1 = new TCanvas("c1","Layer map",201,77,1470,858);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c1->Range(-5.887851,-1.930603,70.65421,17.37543);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetLeftMargin(0.07692308);
  c1->SetTopMargin(0.07124352);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderMode(0);
  

  cadre->GetXaxis()->SetTitle("Module Z index");
  cadre->GetXaxis()->SetNdivisions(n_mod_barrel[layer-5]);
  cadre->GetXaxis()->SetLabelFont(42);
  cadre->GetXaxis()->SetLabelSize(0.02);
  cadre->GetXaxis()->SetTitleSize(0.035);
  cadre->GetXaxis()->SetTickLength(1);
  cadre->GetXaxis()->SetTitleFont(42);
  cadre->GetYaxis()->SetTitle("Module #phi index");
  cadre->GetYaxis()->SetNdivisions(n_lad_barrel[layer-5]);
  cadre->GetYaxis()->SetLabelFont(42);
  cadre->GetYaxis()->SetLabelOffset(0.004);
  cadre->GetYaxis()->SetLabelSize(0.02);
  cadre->GetYaxis()->SetTitleSize(0.03);
  cadre->GetYaxis()->SetTickLength(0.99);
  cadre->GetYaxis()->SetTitleFont(42);
  cadre->GetZaxis()->SetLabelFont(42);
  cadre->GetZaxis()->SetLabelSize(0.035);
  cadre->GetZaxis()->SetTitleSize(0.035);
  cadre->GetZaxis()->SetTitleFont(42);  

  cadre->Draw();
  cadre2->Draw("colzsame");
   
  TPaveText *pt = new TPaveText(0.,n_lad_barrel[layer-5],18./55.*n_mod_barrel[layer-5],17./16.*n_lad_barrel[layer-5],"br");

  char buffer[50];
  sprintf (buffer, "Barrel Layer %d stub rate (in Stub/Module/BX)",layer);

  TLatex Tl;
  Tl.SetTextSize(0.02);

  Tl.DrawLatex(40./55.*n_mod_barrel[layer-5], 16.5/16.*n_lad_barrel[layer-5], "CMS Preliminary Simulation");  

  pt->SetFillColor(0);
  pt->SetTextSize(0.02);
  pt->AddText(buffer);
  pt->Draw();
  c1->Modified();

  c1->Update();

  sprintf (buffer, "Barrel_%d_rate_140.eps", layer);
  c1->Print(buffer);

  sprintf (buffer, "Barrel_%d_rate_140.png", layer);
  c1->Print(buffer);

  c1->Update();


  TCanvas *c2 = new TCanvas("c2","Cluster map",201,77,1470,858);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c2->Range(-5.887851,-1.930603,70.65421,17.37543);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetBorderSize(2);
  c2->SetGridx();
  c2->SetGridy();
  c2->SetLeftMargin(0.07692308);
  c2->SetTopMargin(0.07124352);
  c2->SetFrameBorderMode(0);
  c2->SetFrameBorderMode(0);
  

  cadre->GetXaxis()->SetTitle("Module Z index");
  cadre->GetXaxis()->SetNdivisions(n_mod_barrel[layer-5]);
  cadre->GetXaxis()->SetLabelFont(42);
  cadre->GetXaxis()->SetLabelSize(0.02);
  cadre->GetXaxis()->SetTitleSize(0.035);
  cadre->GetXaxis()->SetTickLength(1);
  cadre->GetXaxis()->SetTitleFont(42);
  cadre->GetYaxis()->SetTitle("Module #phi index");
  cadre->GetYaxis()->SetNdivisions(n_lad_barrel[layer-5]);
  cadre->GetYaxis()->SetLabelFont(42);
  cadre->GetYaxis()->SetLabelOffset(0.004);
  cadre->GetYaxis()->SetLabelSize(0.02);
  cadre->GetYaxis()->SetTitleSize(0.03);
  cadre->GetYaxis()->SetTickLength(0.99);
  cadre->GetYaxis()->SetTitleFont(42);
  cadre->GetZaxis()->SetLabelFont(42);
  cadre->GetZaxis()->SetLabelSize(0.035);
  cadre->GetZaxis()->SetTitleSize(0.035);
  cadre->GetZaxis()->SetTitleFont(42);  

  cadre->Draw();
  cadre3->Draw("colzsame");
   
  TPaveText *pt2 = new TPaveText(0.,n_lad_barrel[layer-5],18./55.*n_mod_barrel[layer-5],17./16.*n_lad_barrel[layer-5],"br");

  sprintf (buffer, "Barrel Layer %d cluster rate (in Cluster/Module/BX)",layer);


  Tl.DrawLatex(40./55.*n_mod_barrel[layer-5], 16.5/16.*n_lad_barrel[layer-5], "CMS Preliminary Simulation");  
  pt2->SetFillColor(0);
  pt2->SetTextSize(0.02);
  pt2->AddText(buffer);
  pt2->Draw();
  c2->Modified();

  c2->Update();

  sprintf (buffer, "Barrel_%d_clus_rate_140.eps", layer);
  c2->Print(buffer);

  sprintf (buffer, "Barrel_%d_clus_rate_140.png", layer);
  c2->Print(buffer);

  c2->Update();


  TCanvas *c3 = new TCanvas("c3","Cluster over Stub map",201,77,1470,858);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c3->Range(-5.887851,-1.930603,70.65421,17.37543);
  c3->SetFillColor(0);
  c3->SetBorderMode(0);
  c3->SetBorderSize(2);
  c3->SetGridx();
  c3->SetGridy();
  c3->SetLeftMargin(0.07692308);
  c3->SetTopMargin(0.07124352);
  c3->SetFrameBorderMode(0);
  c3->SetFrameBorderMode(0);
  
  cadre->GetXaxis()->SetTitle("Module Z index");
  cadre->GetXaxis()->SetNdivisions(n_mod_barrel[layer-5]);
  cadre->GetXaxis()->SetLabelFont(42);
  cadre->GetXaxis()->SetLabelSize(0.02);
  cadre->GetXaxis()->SetTitleSize(0.035);
  cadre->GetXaxis()->SetTickLength(1);
  cadre->GetXaxis()->SetTitleFont(42);
  cadre->GetYaxis()->SetTitle("Module #phi index");
  cadre->GetYaxis()->SetNdivisions(n_lad_barrel[layer-5]);
  cadre->GetYaxis()->SetLabelFont(42);
  cadre->GetYaxis()->SetLabelOffset(0.004);
  cadre->GetYaxis()->SetLabelSize(0.02);
  cadre->GetYaxis()->SetTitleSize(0.03);
  cadre->GetYaxis()->SetTickLength(0.99);
  cadre->GetYaxis()->SetTitleFont(42);
  cadre->GetZaxis()->SetLabelFont(42);
  cadre->GetZaxis()->SetLabelSize(0.035);
  cadre->GetZaxis()->SetTitleSize(0.035);
  cadre->GetZaxis()->SetTitleFont(42);  

  cadre->Draw();
  cadre4->Draw("colzsame");


  TPaveText *pt3 = new TPaveText(0.,n_lad_barrel[layer-5],18./55.*n_mod_barrel[layer-5],17./16.*n_lad_barrel[layer-5],"br");

  sprintf (buffer, "Barrel Layer %d cluster/stub ratio",layer);

  Tl.DrawLatex(40./55.*n_mod_barrel[layer-5], 16.5/16.*n_lad_barrel[layer-5], "CMS Preliminary Simulation");  

  pt3->SetFillColor(0);
  pt3->SetTextSize(0.02);
  pt3->AddText(buffer);
  pt3->Draw();
  c3->Modified();

  c3->Update();

  sprintf (buffer, "Barrel_%d_ratio_140.eps", layer);
  c3->Print(buffer);

  sprintf (buffer, "Barrel_%d_ratio_140.png", layer);
  c3->Print(buffer);

  c3->Update();
 

}


void do_layer_maps(std::string filename,int layer, int tilt)
{
  if (layer<5) return;

  // First get the data
  // by merging all the available files
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  float count_p[58000];
  float count_f[58000];
  float count_s[58000];
  float count_c[58000];
  
  double n_lad_flat[6]  = {16,24,34,48,62,76};
  double n_mod_flat[6]  = {63,55,54,24,24,24};
  
  double n_lad_tilt[6]  = {18,26,36,48,60,78};
  double n_mod_tilt[6]  = {31,35,39,24,24,24};
  
  double n_lad_barrel[6];
  double n_mod_barrel[6];
  
  if (tilt==1)
  {
    for (int i=0;i<6;++i) n_lad_barrel[i] = n_lad_tilt[i];
    for (int i=0;i<6;++i) n_mod_barrel[i] = n_mod_tilt[i];
  }
  else
  {
    for (int i=0;i<6;++i) n_lad_barrel[i] = n_lad_flat[i];
    for (int i=0;i<6;++i) n_mod_barrel[i] = n_mod_flat[i];
  }

  TChain *newtree = new TChain("L1Rates");
    
  newtree->Add(filename.c_str());
    
  newtree->SetBranchAddress("STUB_b_rates_prim",  &count_p);
  newtree->SetBranchAddress("STUB_b_rates_sec",   &count_f);
  newtree->SetBranchAddress("STUB_b_rates_f",     &count_s);
  newtree->SetBranchAddress("CLUS_b_rates",       &count_c);
  newtree->GetEntry(0);
    
  int idx;
  float maxval  = 0;
  float maxcval = 0;
  float rate = 0;
  int maxmod = 0;
  float ratio = 0;
  float maxratio = 0;

  for (int j=(layer-5)*10000;j<(layer-5)*10000+8000;++j)
  {
    rate = count_p[j]+count_f[j]+count_s[j];
    
    if (rate>maxval)
    {
      maxval = rate;
      maxmod = j;
    }
  }
    
  TH2F *cadre = new TH2F("zz","zz",n_mod_barrel[layer-5],0.,n_mod_barrel[layer-5],n_lad_barrel[layer-5],0.,n_lad_barrel[layer-5]);
  TH2F *cadre1 = new TH2F("zz2","zz2",n_mod_barrel[layer-5]+1,0.,n_mod_barrel[layer-5]+1,n_lad_barrel[layer-5],0.,n_lad_barrel[layer-5]);
  TH2F *cadre2 = new TH2F("zz3","zz3",n_mod_barrel[layer-5]+1,0.,n_mod_barrel[layer-5]+1,n_lad_barrel[layer-5],0.,n_lad_barrel[layer-5]);
  TH2F *cadre3 = new TH2F("zz4","zz4",n_mod_barrel[layer-5]+1,0.,n_mod_barrel[layer-5]+1,n_lad_barrel[layer-5],0.,n_lad_barrel[layer-5]);
  TH2F *cadre4 = new TH2F("zz5","zz5",n_mod_barrel[layer-5]+1,0.,n_mod_barrel[layer-5]+1,n_lad_barrel[layer-5],0.,n_lad_barrel[layer-5]);
  
  cadre1->Fill(n_mod_barrel[layer-5]+0.5,0.5,0.);
  cadre1->Fill(n_mod_barrel[layer-5]+0.5,0.5,1.5);
  cadre2->Fill(n_mod_barrel[layer-5]+0.5,0.5,0.);
  cadre2->Fill(n_mod_barrel[layer-5]+0.5,0.5,1.);
  cadre3->Fill(n_mod_barrel[layer-5]+0.5,0.5,0.);
  cadre3->Fill(n_mod_barrel[layer-5]+0.5,0.5,2.5);
  cadre4->Fill(n_mod_barrel[layer-5]+0.5,0.5,0.);
  cadre4->Fill(n_mod_barrel[layer-5]+0.5,0.5,4);
  
  for (int j=0;j<n_mod_barrel[layer-5];++j)
  {
    for (int i=0;i<n_lad_barrel[layer-5];++i)
    {
      idx = 10000*(layer-5) + 100*i + j;
      rate= count_p[idx]+count_f[idx]+count_s[idx];
      
      cadre1->Fill(j+0.5,i+0.5,count_p[idx]);
      cadre2->Fill(j+0.5,i+0.5,count_f[idx]);
      cadre3->Fill(j+0.5,i+0.5,count_s[idx]);
      cadre4->Fill(j+0.5,i+0.5,rate);
    }
  }
    
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  
  TCanvas *c1 = new TCanvas("c1","Layer map",201,77,1470,858);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c1->Range(-5.887851,-1.930603,70.65421,17.37543);
  c1->SetFillColor(0);
  c1->Divide(2,2);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetLeftMargin(0.07692308);
  c1->SetTopMargin(0.07124352);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderMode(0);
  
  c1->cd(1);
  cadre1->GetXaxis()->SetTitle("Module Z index");
  cadre1->GetXaxis()->SetNdivisions(n_mod_barrel[layer-5]);
  cadre1->GetXaxis()->SetLabelFont(42);
  cadre1->GetXaxis()->SetLabelSize(0.02);
  cadre1->GetXaxis()->SetTitleSize(0.035);
  cadre1->GetXaxis()->SetTickLength(1);
  cadre1->GetXaxis()->SetTitleFont(42);
  cadre1->GetYaxis()->SetTitle("Module #phi index");
  cadre1->GetYaxis()->SetNdivisions(n_lad_barrel[layer-5]);
  cadre1->GetYaxis()->SetLabelFont(42);
  cadre1->GetYaxis()->SetLabelOffset(0.004);
  cadre1->GetYaxis()->SetLabelSize(0.02);
  cadre1->GetYaxis()->SetTitleSize(0.03);
  cadre1->GetYaxis()->SetTickLength(0.99);
  cadre1->GetYaxis()->SetTitleFont(42);
  cadre1->GetZaxis()->SetLabelFont(42);
  cadre1->GetZaxis()->SetLabelSize(0.035);
  cadre1->GetZaxis()->SetTitleSize(0.035);
  cadre1->GetZaxis()->SetTitleFont(42);
  
  cadre->Draw("");
  cadre1->Draw("colzsame");
  
  TPaveText *pt = new TPaveText(0.,n_lad_barrel[layer-5],28./55.*n_mod_barrel[layer-5],17./16.*n_lad_barrel[layer-5],"br");
  
  char buffer[50];
  sprintf (buffer, "Barrel Layer %d primary stub rate (in Stub/Module/BX)",layer);
    
  TLatex Tl;
  Tl.SetTextSize(0.03);
    
  Tl.DrawLatex(40./55.*n_mod_barrel[layer-5], 16.5/16.*n_lad_barrel[layer-5], "CMS Preliminary Simulation");

  pt->SetFillColor(0);
  pt->SetTextSize(0.03);
  pt->AddText(buffer);
  pt->Draw();
  
  c1->cd(2);
  cadre2->GetXaxis()->SetTitle("Module Z index");
  cadre2->GetXaxis()->SetNdivisions(n_mod_barrel[layer-5]);
  cadre2->GetXaxis()->SetLabelFont(42);
  cadre2->GetXaxis()->SetLabelSize(0.02);
  cadre2->GetXaxis()->SetTitleSize(0.035);
  cadre2->GetXaxis()->SetTickLength(1);
  cadre2->GetXaxis()->SetTitleFont(42);
  cadre2->GetYaxis()->SetTitle("Module #phi index");
  cadre2->GetYaxis()->SetNdivisions(n_lad_barrel[layer-5]);
  cadre2->GetYaxis()->SetLabelFont(42);
  cadre2->GetYaxis()->SetLabelOffset(0.004);
  cadre2->GetYaxis()->SetLabelSize(0.02);
  cadre2->GetYaxis()->SetTitleSize(0.03);
  cadre2->GetYaxis()->SetTickLength(0.99);
  cadre2->GetYaxis()->SetTitleFont(42);
  cadre2->GetZaxis()->SetLabelFont(42);
  cadre2->GetZaxis()->SetLabelSize(0.035);
  cadre2->GetZaxis()->SetTitleSize(0.035);
  cadre2->GetZaxis()->SetTitleFont(42);
  
  cadre->Draw("");
  cadre2->Draw("colzsame");
  
  TPaveText *pt2 = new TPaveText(0.,n_lad_barrel[layer-5],28./55.*n_mod_barrel[layer-5],17./16.*n_lad_barrel[layer-5],"br");
    
  sprintf (buffer, "Barrel Layer %d secondary stub rate (in Stub/Module/BX)",layer);
    
  TLatex Tl2;
  Tl2.SetTextSize(0.03);
    
  Tl2.DrawLatex(40./55.*n_mod_barrel[layer-5], 16.5/16.*n_lad_barrel[layer-5], "CMS Preliminary Simulation");

  pt2->SetFillColor(0);
  pt2->SetTextSize(0.03);
  pt2->AddText(buffer);
  pt2->Draw();
        
  c1->cd(3);
  cadre3->GetXaxis()->SetTitle("Module Z index");
  cadre3->GetXaxis()->SetNdivisions(n_mod_barrel[layer-5]);
  cadre3->GetXaxis()->SetLabelFont(42);
  cadre3->GetXaxis()->SetLabelSize(0.02);
  cadre3->GetXaxis()->SetTitleSize(0.035);
  cadre3->GetXaxis()->SetTickLength(1);
  cadre3->GetXaxis()->SetTitleFont(42);
  cadre3->GetYaxis()->SetTitle("Module #phi index");
  cadre3->GetYaxis()->SetNdivisions(n_lad_barrel[layer-5]);
  cadre3->GetYaxis()->SetLabelFont(42);
  cadre3->GetYaxis()->SetLabelOffset(0.004);
  cadre3->GetYaxis()->SetLabelSize(0.02);
  cadre3->GetYaxis()->SetTitleSize(0.03);
  cadre3->GetYaxis()->SetTickLength(0.99);
  cadre3->GetYaxis()->SetTitleFont(42);
  cadre3->GetZaxis()->SetLabelFont(42);
  cadre3->GetZaxis()->SetLabelSize(0.035);
  cadre3->GetZaxis()->SetTitleSize(0.035);
  cadre3->GetZaxis()->SetTitleFont(42);
  
  cadre->Draw("");
  cadre3->Draw("colzsame");
  
  TPaveText *pt3 = new TPaveText(0.,n_lad_barrel[layer-5],28./55.*n_mod_barrel[layer-5],17./16.*n_lad_barrel[layer-5],"br");
    
  sprintf (buffer, "Barrel Layer %d fake stub rate (in Stub/Module/BX)",layer);
    
  TLatex Tl3;
  Tl3.SetTextSize(0.03);
  
  Tl3.DrawLatex(40./55.*n_mod_barrel[layer-5], 16.5/16.*n_lad_barrel[layer-5], "CMS Preliminary Simulation");

  pt3->SetFillColor(0);
  pt3->SetTextSize(0.03);
  pt3->AddText(buffer);
  pt3->Draw();
      
  c1->cd(4);
  cadre4->GetXaxis()->SetTitle("Module Z index");
  cadre4->GetXaxis()->SetNdivisions(n_mod_barrel[layer-5]);
  cadre4->GetXaxis()->SetLabelFont(42);
  cadre4->GetXaxis()->SetLabelSize(0.02);
  cadre4->GetXaxis()->SetTitleSize(0.035);
  cadre4->GetXaxis()->SetTickLength(1);
  cadre4->GetXaxis()->SetTitleFont(42);
  cadre4->GetYaxis()->SetTitle("Module #phi index");
  cadre4->GetYaxis()->SetNdivisions(n_lad_barrel[layer-5]);
  cadre4->GetYaxis()->SetLabelFont(42);
  cadre4->GetYaxis()->SetLabelOffset(0.004);
  cadre4->GetYaxis()->SetLabelSize(0.02);
  cadre4->GetYaxis()->SetTitleSize(0.03);
  cadre4->GetYaxis()->SetTickLength(0.99);
  cadre4->GetYaxis()->SetTitleFont(42);
  cadre4->GetZaxis()->SetLabelFont(42);
  cadre4->GetZaxis()->SetLabelSize(0.035);
  cadre4->GetZaxis()->SetTitleSize(0.035);
  cadre4->GetZaxis()->SetTitleFont(42);
  
  cadre->Draw("");
  cadre4->Draw("colzsame");
  
  
  TPaveText *pt4 = new TPaveText(0.,n_lad_barrel[layer-5],28./55.*n_mod_barrel[layer-5],17./16.*n_lad_barrel[layer-5],"br");
  
  sprintf (buffer, "Barrel Layer %d TOTAL stub rate (in Stub/Module/BX)",layer);
  
  TLatex Tl4;
  Tl4.SetTextSize(0.03);
    
  Tl4.DrawLatex(40./55.*n_mod_barrel[layer-5], 16.5/16.*n_lad_barrel[layer-5], "CMS Preliminary Simulation");
  

  pt4->SetFillColor(0);
  pt4->SetTextSize(0.03);
  pt4->AddText(buffer);
  pt4->Draw();
  
  c1->Modified();
  c1->Update();
  
  sprintf (buffer, "Barrel_%d_rate_140.eps", layer);
  c1->Print(buffer);
  
  sprintf (buffer, "Barrel_%d_rate_140.png", layer);
  c1->Print(buffer);
  
  c1->Update();
}



// Rate plot for a given disk

void do_disk_map(std::string filename,int disk, float pri, float sec, float fak,int type)
{

  // First get the data
  // by merging all the available files

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  float count_e[142000];

  float count_p[142000];
  float count_f[142000];
  float count_s[142000];
  float count_c[142000];

  double n_lad_endcap[15]   = {20,24,24,28,32,32,36,40,40,44,52,60,64,72,76};
  double n_start_endcap[15] = {0,0,0,0,0,0,1,0,0,0,1,0,0,0,1};

  if (disk>=7) // -Z case
  {
    for (int i=0;i<15;++i) n_start_endcap[i] = 0; // There is no problem in this case
  }

  TChain *newtree = new TChain("L1Rates");
    
  newtree->Add(filename.c_str()); 
  newtree->SetBranchAddress("STUB_e_rates_prim",  &count_p);
  newtree->SetBranchAddress("STUB_e_rates_sec",   &count_f);
  newtree->SetBranchAddress("STUB_e_rates_f",     &count_s);
  newtree->SetBranchAddress("CLUS_e_rates",       &count_c);

  // Initialize some params
  newtree->GetEntry(0);

  float maxval = 0;
  float rate   = 0;

 
  for (int i=0;i<142000;++i) count_e[i]=0.; 
   
  for (int i=disk;i<disk+1;++i) 
  { 
    for (int j=10000*i;j<10000*i+2000;++j)
    {
      rate = pri*count_p[j]+fak*count_f[j]+sec*count_s[j];

      if (type==1)            count_e[j]=rate; 
      if (type==2)            count_e[j]=count_c[j]; 
      if (type==3 && rate!=0) count_e[j]=count_c[j]/rate; 

      if (count_e[j]>maxval) maxval = count_e[j];
    }
  }

  TH2F *cadre  = new TH2F("zz","zz",30,-15.,15.,30,-15.,15.);
  TH2F *ring1  = new TH2F("r1","r1",n_lad_endcap[0],0.,8*atan(1.),26,0.,26.);
  TH2F *ring2  = new TH2F("r2","r2",n_lad_endcap[1],0.,8*atan(1.),26,0.,26.);
  TH2F *ring3  = new TH2F("r3","r3",n_lad_endcap[2],0.,8*atan(1.),26,0.,26.);
  TH2F *ring4  = new TH2F("r4","r4",n_lad_endcap[3],0.,8*atan(1.),26,0.,26.);
  TH2F *ring5  = new TH2F("r5","r5",n_lad_endcap[4],0.,8*atan(1.),26,0.,26.);
  TH2F *ring6  = new TH2F("r6","r6",n_lad_endcap[5],0.,8*atan(1.),26,0.,26.);
  TH2F *ring7  = new TH2F("r7","r7",n_lad_endcap[6],0.,8*atan(1.),26,0.,26.);
  TH2F *ring8  = new TH2F("r8","r8",n_lad_endcap[7],0.,8*atan(1.),26,0.,26.);
  TH2F *ring9  = new TH2F("r9","r9",n_lad_endcap[8],0.,8*atan(1.),26,0.,26.);
  TH2F *ring10 = new TH2F("r10","r10",n_lad_endcap[9],0.,8*atan(1.),26,0.,26.);
  TH2F *ring11 = new TH2F("r11","r11",n_lad_endcap[10],0.,8*atan(1.),26,0.,26.);
  TH2F *ring12 = new TH2F("r12","r12",n_lad_endcap[11],0.,8*atan(1.),26,0.,26.);
  TH2F *ring13 = new TH2F("r13","r13",n_lad_endcap[12],0.,8*atan(1.),26,0.,26.);
  TH2F *ring14 = new TH2F("r14","r14",n_lad_endcap[13],0.,8*atan(1.),26,0.,26.);
  TH2F *ring15 = new TH2F("r15","r15",n_lad_endcap[14],0.,8*atan(1.),26,0.,26.);

  fill_histo(ring1,count_e,n_lad_endcap[0],n_start_endcap[0],disk,0,maxval);
  fill_histo(ring2,count_e,n_lad_endcap[1],n_start_endcap[1],disk,1,maxval);
  fill_histo(ring3,count_e,n_lad_endcap[2],n_start_endcap[2],disk,2,maxval);
  fill_histo(ring4,count_e,n_lad_endcap[3],n_start_endcap[3],disk,3,maxval);
  fill_histo(ring5,count_e,n_lad_endcap[4],n_start_endcap[4],disk,4,maxval);
  fill_histo(ring6,count_e,n_lad_endcap[5],n_start_endcap[5],disk,5,maxval);
  fill_histo(ring7,count_e,n_lad_endcap[6],n_start_endcap[6],disk,6,maxval);
  fill_histo(ring8,count_e,n_lad_endcap[7],n_start_endcap[7],disk,7,maxval);
  fill_histo(ring9,count_e,n_lad_endcap[8],n_start_endcap[8],disk,8,maxval);
  fill_histo(ring10,count_e,n_lad_endcap[9],n_start_endcap[9],disk,9,maxval);
  fill_histo(ring11,count_e,n_lad_endcap[10],n_start_endcap[10],disk,10,maxval);
  fill_histo(ring12,count_e,n_lad_endcap[11],n_start_endcap[11],disk,11,maxval);
  fill_histo(ring13,count_e,n_lad_endcap[12],n_start_endcap[12],disk,12,maxval);
  fill_histo(ring14,count_e,n_lad_endcap[13],n_start_endcap[13],disk,13,maxval);
  fill_histo(ring15,count_e,n_lad_endcap[14],n_start_endcap[14],disk,14,maxval);
  
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  
  TCanvas *c1 = new TCanvas("c1","Layer map",166,77,947,880);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c1->Range(-15.23757,-14.89125,18.52486,15.55968);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetLeftMargin(0.03665521);
  c1->SetRightMargin(0.1340206);
  c1->SetTopMargin(0.05121951);
  c1->SetBottomMargin(0.02926829);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderMode(0);

  cadre->GetXaxis()->SetLabelSize(0.);
  cadre->GetXaxis()->SetLabelOffset(999);
  cadre->GetYaxis()->SetLabelSize(0.);
  cadre->GetYaxis()->SetLabelOffset(999);
  
  cadre->Draw("");
  ring1->Draw("polcolsame");
  ring2->Draw("polcolsame"); 
  ring3->Draw("polcolsame");
  ring4->Draw("polcolsame");
  ring5->Draw("polcolsame");
  ring6->Draw("polcolsame");
  ring7->Draw("polcolsame");
  ring8->Draw("polcolsame");
  ring9->Draw("polcolsame");
  ring10->Draw("polcolsame");
  ring11->Draw("polcolsame");
  ring12->Draw("polcolsame");
  ring13->Draw("polcolsame");
  ring14->Draw("polcolzsame");
  ring15->Draw("polcolzsame");
  
  TPaveText *pt = new TPaveText(-15.,14.9,3.,16.,"br");

  char buffer[100];

  if (disk>=7)
  {
    if (type==1) sprintf(buffer, "Endcap Disk %d (%sZ side) stub rate (in Stub/Module/BX)", (disk-7)%7+1,"-");
    if (type==2) sprintf(buffer, "Endcap Disk %d (%sZ side) cluster rate (in Cluster/Module/BX)", (disk-7)%7+1,"-");
    if (type==3) sprintf(buffer, "Endcap Disk %d (%sZ side) cluster/stub ratio", (disk-7)%7+1,"-");
  }
  else
  {
    if (type==1) sprintf(buffer, "Endcap Disk %d (%sZ side) stub rate (in Stub/Module/BX)", disk%7+1,"+");
    if (type==2) sprintf(buffer, "Endcap Disk %d (%sZ side) cluster rate (in Cluster/Module/BX)", disk%7+1,"+");
    if (type==3) sprintf(buffer, "Endcap Disk %d (%sZ side) cluster/stub ratio", disk%7+1,"+");
  }

  pt->SetFillColor(0);
  pt->SetTextSize(0.02);

  TLatex Tl;
  Tl.SetTextSize(0.02);
  Tl.DrawLatex(9., 15.3, "CMS Preliminary Simulation");     

  TText *text = pt->AddText(buffer);
  pt->Draw();
  c1->Modified();
 
  c1->Update();
  
  if (type==1) sprintf (buffer, "Disk_%d_rate_140.eps", disk+1);
  if (type==2) sprintf (buffer, "Disk_%d_crate_140.eps", disk+1);
  if (type==3) sprintf (buffer, "Disk_%d_ratio_140.eps", disk+1);
  c1->Print(buffer);

  if (type==1) sprintf (buffer, "Disk_%d_rate_140.png", disk+1);
  if (type==2) sprintf (buffer, "Disk_%d_crate_140.png", disk+1);
  if (type==3) sprintf (buffer, "Disk_%d_ratio_140.png", disk+1);
  c1->Print(buffer);
   
}



//
// Print average rates over a given ladder/ring (in stub/bx/module)
// 


void print_TIB_rate(std::string filename)
{
  double n_lad_barrel[6]  = {18,26,36};
  double n_mod_barrel[6]  = {31,35,39};
    
  float count_p[58000];
  float count_pp[58000];
  float count_f[58000];
  float count_s[58000];

  TChain *newtree = new TChain("L1Rates");
    
  newtree->Add(filename.c_str()); 
  newtree->SetBranchAddress("STUB_b_rates_prim",  &count_p);
  newtree->SetBranchAddress("STUB_b_rates_prim2", &count_pp);
  newtree->SetBranchAddress("STUB_b_rates_sec",   &count_f);
  newtree->SetBranchAddress("STUB_b_rates_f",     &count_s);
  newtree->GetEntry(0);
  
  double rate   = 0;
  double rate_c = 0;
  double rate_s = 0;
  double rate_p = 0;
  double rate_pp= 0;
  int idx;
  
  std::cout.precision(3);

  for (int i=0;i<3;++i)
  { 
    cout << endl;
    cout << "Layer " << i+1 << endl;
    cout << "        /  TOT  / PRIM>2/ PRIM<2/  SEC  /  COMB" << endl;
    
    
    for (int j=0;j<n_mod_barrel[i];++j)
    {
      cout << "Ring " << j+1 << "  / ";

      rate = 0;
      rate_c = 0;
      rate_p = 0;
      rate_pp= 0;
      rate_s = 0;
	
      for (int k=0;k<n_lad_barrel[i];++k)
      {
	idx = 10000*(i) + 100*k + j;

	rate_c += count_f[idx];
	rate_p += count_p[idx]-count_pp[idx];
	rate_pp+= count_pp[idx];
	rate_s += count_s[idx];
      }
 
      float tot_rate = (rate_pp+rate_p+rate_s+rate_c);
      float factor   = 100./tot_rate;
        
      cout << std::fixed << std::setprecision(1) << std::setw(5) << tot_rate/n_lad_barrel[i]
	   << " / "<< std::setw(5) << rate_pp*factor
	   << " / "<< std::setw(5) << rate_p*factor
	   << " / "<< std::setw(5) << rate_s*factor
	   << " / "<< std::setw(5) << rate_c*factor << endl;
    }
  }
}

void print_TOB_rate(std::string filename)
{
  double n_lad_barrel[6]  = {48,60,78};
  double n_mod_barrel[6]  = {24,24,24};
  double n_off_barrel[6]  = {12,12,12};

  float count_p[58000];
  float count_pp[58000];
  float count_f[58000];
  float count_s[58000];

  TChain *newtree = new TChain("L1Rates");
    
  newtree->Add(filename.c_str()); 
  newtree->SetBranchAddress("STUB_b_rates_prim",  &count_p);
  newtree->SetBranchAddress("STUB_b_rates_prim2", &count_pp);
  newtree->SetBranchAddress("STUB_b_rates_sec",   &count_f);
  newtree->SetBranchAddress("STUB_b_rates_f",     &count_s);
  newtree->GetEntry(0);
  
  double rate   = 0;
  double rate_c = 0;
  double rate_s = 0;
  double rate_p = 0;
  double rate_pp= 0;
  int idx;
  
  std::cout.precision(3);

  for (int i=0;i<3;++i)
  { 
    cout << endl;
    cout << "Layer " << i+1 << endl;
    cout << "        /  TOT  / PRIM>2/ PRIM<2/  SEC  /  COMB" << endl;
    
    
    for (int j=n_off_barrel[i];j<n_mod_barrel[i];++j)
    {
      if (int(j+1-n_off_barrel[i])<10)
	cout << "Ring " << int(j+1-n_off_barrel[i]) << "  / ";
      if (int(j+1-n_off_barrel[i])>=10)
	cout << "Ring " << int(j+1-n_off_barrel[i]) << " / ";

      rate = 0;
      rate_c = 0;
      rate_p = 0;
      rate_pp= 0;
      rate_s = 0;
	
      for (int k=0;k<n_lad_barrel[i];++k)
      {
	idx = 10000*(i+3) + 100*k + j;

	rate_c += count_f[idx];
	rate_p += count_p[idx]-count_pp[idx];
	rate_pp+= count_pp[idx];
	rate_s += count_s[idx];
      }

        float tot_rate = (rate_pp+rate_p+rate_s+rate_c);
        float factor   = 100./tot_rate;
        
        cout << std::fixed << std::setprecision(1) << std::setw(5) << tot_rate/n_lad_barrel[i]
        << " / "<< std::setw(5) << rate_pp*factor
        << " / "<< std::setw(5) << rate_p*factor
        << " / "<< std::setw(5) << rate_s*factor
        << " / "<< std::setw(5) << rate_c*factor << endl;

    }
  }
}

void print_ENDCAP_rate(std::string filename)
{

  double n_lad_endcap[15]   = {20,24,24,28,32,32,36,40,40,44,52,60,64,72,76};

  float count_p[142000];
  float count_pp[142000];
  float count_f[142000];
  float count_s[142000];

  TChain *newtree = new TChain("L1Rates");
    
  newtree->Add(filename.c_str()); 
  newtree->SetBranchAddress("STUB_e_rates_prim",  &count_p);
  newtree->SetBranchAddress("STUB_e_rates_prim2", &count_pp);
  newtree->SetBranchAddress("STUB_e_rates_sec",   &count_f);
  newtree->SetBranchAddress("STUB_e_rates_f",     &count_s);
  newtree->GetEntry(0);
  
  double rate   = 0;
  double rate_c = 0;
  double rate_s = 0;
  double rate_p = 0;
  double rate_pp= 0;
  int idx;
  
  std::cout.precision(3);

  for (int i=0;i<5;++i)
  { 
    cout << endl;
    cout << "DISK " << i+1 << endl;
    cout << "        /  TOT  / PRIM>2/ PRIM<2/  SEC  /  COMB" << endl;
    
    
    for (int j=0;j<15;++j)
    {
      if (j<10)
	cout << "Ring " << j << "  / ";
      if (j>=10)
	cout << "Ring " << j << " / ";

      rate = 0;
      rate_c = 0;
      rate_p = 0;
      rate_pp= 0;
      rate_s = 0;
	
      for (int k=0;k<n_lad_endcap[j];++k)
      {
	idx = 10000*i + 100*j + k;

	rate_c += count_f[idx];
	rate_p += count_p[idx]-count_pp[idx];
	rate_pp+= count_pp[idx];
	rate_s += count_s[idx];
      }
        
        float tot_rate = (rate_pp+rate_p+rate_s+rate_c);
        float factor   = 100./tot_rate;
        
        cout << std::fixed << std::setprecision(1) << std::setw(5) << tot_rate/n_lad_endcap[j]
        << " / "<< std::setw(5) << rate_pp*factor
        << " / "<< std::setw(5) << rate_p*factor
        << " / "<< std::setw(5) << rate_s*factor
        << " / "<< std::setw(5) << rate_c*factor << endl;
        
    }
  }
}

