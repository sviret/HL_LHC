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
  ROOT macro for stub efficiencies visualization

  Use:

  root[1]-> .L MakeEffGraphs.C

  Then you can choose between the following methods 
  knowing that filename denotes the extracted ROOT file 
  containing the stub efficiencies (created by the 
  sectorMaker program).

  For infos about the efficiency definition, have a look at the following presentation:

  https://indico.cern.ch/getFile.py/access?contribId=0&resId=0&materialId=slides&confId=263068

  Solution 1 do all the plots (lot of infos)
  root[2]-> do_efficiency_plots(std::string filename) 

  Solution 2 just plot the comparicon pT plots for the barrel
  root[2]-> do_pt_barrel_summary(std::string filename,int type) 

  where type=0 is for the module efficiency, and type=1 for the 
  layer/disk efficiency

  Solution 3 just plot the comparison pT plots for the endcap
  root[2]-> do_pt_endcap_summary(std::string filename,int type) 

  Solution 4 do the plot for barrel layer i
  root[2]-> do_layer_map(std::string filename,int layer,int type)

  Solution 5 do the plot for endcap disk i
  root[2]-> do_disk_map(std::string filename,int disk,int type)

  Information about this macro may be found in the following page:

  http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto (STEP III part 2.4)

  Author: viret@in2p3_dot_fr
  Date: 20/07/2013
*/

#include <iomanip>      // for std::setprecision
/*
void do_efficiency_plots(std::string filename)
{
  for (int i=0;i<20;++i)
  { 
    do_pt_plots(filename,i,0);
    do_pt_plots(filename,i,1);
  }

  do_pt_barrel_summary(filename,0);
  do_pt_barrel_summary(filename,1);
  do_pt_endcap_summary(filename,0);
  do_pt_endcap_summary(filename,1);


  for (int i=0;i<13;++i)
  { 
    do_eta_plots(filename,i,0);
    do_eta_plots(filename,i,1);
  }
}
*/
void do_pt_plots(std::string filename,int layer,int per_lay,float pmin)
{
  // First get the data
  // by merging all the available files

  if (layer<5) return;

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
 
  float pt_val[100];

  float digi_pt[30][100];

  float clus_pt[30][100];
  float stub_pt[30][100];

  int disk =0;

  if (layer>=11  && layer<18) disk=int(layer-11)%8;
  if (layer>=18 && layer<25)  disk=-int(layer-18)%8;

  TChain *newtree = new TChain("Efficiencies");
  newtree->Add(filename.c_str());
  newtree->SetBranchAddress("pt_val",&pt_val);
  newtree->SetBranchAddress("digi_pt",&digi_pt);
  newtree->SetBranchAddress("clus_pt",&clus_pt);

  (per_lay==1)
    ? newtree->SetBranchAddress("stub_pt_lay",&stub_pt)
    : newtree->SetBranchAddress("stub_pt",&stub_pt);

  newtree->GetEntry(0);

  TH2F *cadre    = new TH2F("zz","zz",300,0.,pmin,200,0.,1.02);
  TH2F *digi_eff = new TH2F("deff","deff",300,0.,20.,200,0.,1.02);
  TH2F *coff_eff = new TH2F("coeff","coeff",300,0.,20.,200,0.,1.02);
  TH2F *soff_eff = new TH2F("soeff","soeff",300,0.,20.,200,0.,1.02);

  for (int i=0;i<100;++i)
  { 
    digi_eff->Fill(pt_val[i],digi_pt[layer-5][i]);
    coff_eff->Fill(pt_val[i],clus_pt[layer-5][i]);
    soff_eff->Fill(pt_val[i],stub_pt[layer-5][i]);
  }

  TCanvas *c2 = new TCanvas("c2","Official effs",5,75,670,660);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c2->Range(-1.2,-0.1,10.5,1.1);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetBorderSize(2);
  c2->SetGridx();
  c2->SetGridy();
  c2->SetLeftMargin(0.04);
  c2->SetRightMargin(0.05);
  c2->SetFrameBorderMode(0);
  c2->SetFrameBorderMode(0);

  cadre->GetXaxis()->SetTitle("Particle p_{T} (GeV)");
  cadre->GetXaxis()->SetNdivisions(518);
  cadre->GetXaxis()->SetLabelFont(42);
  cadre->GetXaxis()->SetLabelSize(0.03);
  cadre->GetXaxis()->SetTitleFont(42);
  cadre->GetYaxis()->SetTitle("Efficiency");
  cadre->GetYaxis()->SetNdivisions(511);
  cadre->GetYaxis()->SetLabelFont(42);
  cadre->GetYaxis()->SetLabelSize(0.03);
  cadre->GetYaxis()->SetTitleOffset(1.25);
  cadre->GetYaxis()->SetTitleFont(42);
  cadre->Draw();
  
  digi_eff->SetMarkerStyle(3);
  coff_eff->SetMarkerStyle(4);
  soff_eff->SetMarkerStyle(20);
  digi_eff->Draw("same");
  coff_eff->Draw("same");
  soff_eff->Draw("same");


  char buffer[80];

  if (disk==0)
  {
    sprintf(buffer, "Layer %d",layer);
  }
  else
  {
    (disk>0)
      ? sprintf(buffer, "Disk %d",disk)
      : sprintf(buffer, "Disk %d",-disk);
  }

  TLegend *leg = new TLegend(0.7,0.28,0.9,0.45);

  leg->SetTextSize(0.03);
  leg->SetHeader(buffer); 
  leg->AddEntry(digi_eff,"Digis","p");
  leg->AddEntry(coff_eff,"Clusters","p");
  leg->AddEntry(soff_eff,"Stubs","p");
  leg->Draw();
    
  sprintf (buffer, "#sqrt{s}=14TeV, 0 PU");

  TLatex Tl;
  Tl.SetTextSize(0.04);
  Tl.DrawLatex(0., 1.03, "CMS Phase-2 Simulation");

  TLatex Tl2;
  Tl2.SetTextSize(0.03);
  Tl2.SetTextFont(52);
  Tl2.DrawLatex(0.6*pmin, 1.03, buffer);

  c2->Modified();
  c2->Update();

  if (per_lay==1)
  {
    if (disk==0) sprintf (buffer, "Barrel_%d_layer_eff.png", layer+5);
    if (disk>0)  sprintf (buffer, "Endcap_p%d_layer_eff.png", disk); 
    if (disk<0)  sprintf (buffer, "Endcap_m%d_layer_eff.png", abs(disk)); 
  }
  else
  {
    if (disk==0) sprintf (buffer, "Barrel_%d_module_eff.png", layer+5);
    if (disk>0)  sprintf (buffer, "Endcap_p%d_module_eff.png", disk); 
    if (disk<0)  sprintf (buffer, "Endcap_m%d_module_eff.png", abs(disk)); 
  }

  c2->Print(buffer);
}


void do_pt_barrel_summary(std::string filename, int per_lay)
{
  // First get the data
  // by merging all the available files

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
 
  float pt_val[100];
  float stub_pt[30][100];

  TChain *newtree = new TChain("Efficiencies");
  newtree->Add(filename.c_str());
  newtree->SetBranchAddress("pt_val",&pt_val);

  (per_lay==1)
    ? newtree->SetBranchAddress("stub_pt_lay",&stub_pt)
    : newtree->SetBranchAddress("stub_pt",&stub_pt);

  newtree->GetEntry(0);

  TH2F *cadre    = new TH2F("zz","zz",300,0.,10.,100,0.,1.15);
  TH2F *L1_off   = new TH2F("L1off","L1off",300,0.,20.,100,0.,1.02);
  TH2F *L2_off   = new TH2F("L2off","L2off",300,0.,20.,100,0.,1.02);
  TH2F *L3_off   = new TH2F("L3off","L3off",300,0.,20.,100,0.,1.02);
  TH2F *L4_off   = new TH2F("L4off","L4off",300,0.,20.,100,0.,1.02);
  TH2F *L5_off   = new TH2F("L5off","L5off",300,0.,20.,100,0.,1.02);
  TH2F *L6_off   = new TH2F("L6off","L6off",300,0.,20.,100,0.,1.02);

  for (int i=0;i<100;++i)
  { 
    if (stub_pt[0][i]!=0) L1_off->Fill(pt_val[i],stub_pt[0][i]);
    if (stub_pt[1][i]!=0) L2_off->Fill(pt_val[i],stub_pt[1][i]);
    if (stub_pt[2][i]!=0) L3_off->Fill(pt_val[i],stub_pt[2][i]);
    if (stub_pt[3][i]!=0) L4_off->Fill(pt_val[i],stub_pt[3][i]);
    if (stub_pt[4][i]!=0) L5_off->Fill(pt_val[i],stub_pt[4][i]);
    if (stub_pt[5][i]!=0) L6_off->Fill(pt_val[i],stub_pt[5][i]);
  }

  TCanvas *c2 = new TCanvas("c1","Official effs",0,0,640,640);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c2->Range(-1.2,-0.1,10.6,1.24);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetBorderSize(2);
  c2->SetGridx();
  c2->SetGridy();
  c2->SetRightMargin(0.07);
  c2->SetTopMargin(0.065);
  c2->SetFrameBorderMode(0);
  c2->SetFrameBorderMode(0);

    
  cadre->GetXaxis()->SetTitle("Particle p_{T} (GeV)");
  cadre->GetXaxis()->SetNdivisions(518);
  cadre->GetXaxis()->SetLabelFont(42);
  cadre->GetXaxis()->SetLabelSize(0.035);
  cadre->GetXaxis()->SetTitleFont(42);
  cadre->GetXaxis()->SetTitleSize(0.045);
  cadre->GetYaxis()->SetTitleSize(0.035);
  cadre->GetYaxis()->SetTitleOffset(0.88);
  cadre->GetXaxis()->SetTitleOffset(0.9);
  cadre->GetYaxis()->SetTitle("Stub efficiency");
  cadre->GetYaxis()->SetNdivisions(511);
  cadre->GetYaxis()->SetLabelFont(42);
  cadre->GetYaxis()->SetLabelSize(0.035);
  cadre->GetYaxis()->SetTitleFont(42);
  cadre->Draw();
  
  L1_off->SetMarkerStyle(24);
  L2_off->SetMarkerStyle(20);
  L3_off->SetMarkerStyle(27);
  L4_off->SetMarkerStyle(34);
  L5_off->SetMarkerStyle(28);
  L6_off->SetMarkerStyle(25);
    
  L1_off->SetMarkerSize(1.4);
  L2_off->SetMarkerSize(1.4);
  L3_off->SetMarkerSize(1.4);
  L4_off->SetMarkerSize(1.4);
  L5_off->SetMarkerSize(1.4);
  L6_off->SetMarkerSize(1.4);

  L1_off->Draw("same");
  L2_off->Draw("same");
  L3_off->Draw("same");
  L4_off->Draw("same");
  L5_off->Draw("same");
  L6_off->Draw("same");

  TLegend *leg = new TLegend(0.6,0.2,0.86,0.45);

  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->AddEntry(L1_off,"TBPS layer 1","p");
  leg->AddEntry(L2_off,"TBPS layer 2","p");
  leg->AddEntry(L3_off,"TBPS layer 3","p");
  leg->AddEntry(L4_off,"TB2S layer 1","p");
  leg->AddEntry(L5_off,"TB2S layer 2","p");
  leg->AddEntry(L6_off,"TB2S layer 3","p");
  leg->Draw();

  char buffer[80];
       
  TLatex Tl;
  Tl.SetTextSize(0.03);
  Tl.DrawLatex(-0.02, 1.18, "CMS Phase-2 simulation");
    
  sprintf (buffer, "#sqrt{s}=14TeV, 0 PU");

  TLatex Tl2;
  Tl2.SetTextSize(0.03);
  Tl2.SetTextFont(52);
  Tl2.DrawLatex(5., 1.18, buffer);
    
  c2->Modified();
  c2->Update();
    
  if (per_lay==1)
  {
    sprintf (buffer, "Barrel_layer_eff.pdf");
    
    c2->Print(buffer);
    sprintf (buffer, "Barrel_layer_eff.C");
    c2->SaveSource(buffer);
  }
  else
  {
    sprintf (buffer, "Barrel_module_eff.pdf");
    c2->Print(buffer);
  }

}


void do_pt_endcap_summary(std::string filename, int per_lay)
{
  // First get the data
  // by merging all the available files

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
 
  float pt_val[100];
  float stub_pt[30][100];

  TChain *newtree = new TChain("Efficiencies");
  newtree->Add(filename.c_str());
  newtree->SetBranchAddress("pt_val",&pt_val);

  (per_lay==1)
    ? newtree->SetBranchAddress("stub_pt_lay",&stub_pt)
    : newtree->SetBranchAddress("stub_pt",&stub_pt);

  newtree->GetEntry(0);

  TH2F *cadre    = new TH2F("zz","zz",300,0.,10.,100,0.,1.15);
  TH2F *D1_off   = new TH2F("L1off","L1off",300,0.,20.,100,0.,1.05);
  TH2F *D2_off   = new TH2F("L2off","L2off",300,0.,20.,100,0.,1.05);
  TH2F *D3_off   = new TH2F("L3off","L3off",300,0.,20.,100,0.,1.05);
  TH2F *D4_off   = new TH2F("L4off","L4off",300,0.,20.,100,0.,1.05);
  TH2F *D5_off   = new TH2F("L5off","L5off",300,0.,20.,100,0.,1.05);

  for (int i=0;i<100;++i)
  { 
    if (stub_pt[6][i]!=0) D1_off->Fill(pt_val[i],stub_pt[6][i]);
    if (stub_pt[7][i]!=0) D2_off->Fill(pt_val[i],stub_pt[7][i]);
    if (stub_pt[8][i]!=0) D3_off->Fill(pt_val[i],stub_pt[8][i]);
    if (stub_pt[9][i]!=0) D4_off->Fill(pt_val[i],stub_pt[9][i]);
    if (stub_pt[10][i]!=0) D5_off->Fill(pt_val[i],stub_pt[10][i]);
  }

  TCanvas *c2 = new TCanvas("c2","Official effs",0,0,640,640);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c2->Range(-1.2,-0.1,10.6,1.24);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetBorderSize(2);
  c2->SetGridx();
  c2->SetGridy();
  c2->SetRightMargin(0.07);
  c2->SetTopMargin(0.065);
  c2->SetFrameBorderMode(0);
  c2->SetFrameBorderMode(0);

    
  cadre->GetXaxis()->SetTitle("Particle p_{T} (GeV)");
  cadre->GetXaxis()->SetNdivisions(518);
  cadre->GetXaxis()->SetLabelFont(42);
  cadre->GetXaxis()->SetLabelSize(0.035);
  cadre->GetXaxis()->SetTitleFont(42);
  cadre->GetXaxis()->SetTitleSize(0.045);
  cadre->GetYaxis()->SetTitleSize(0.035);
  cadre->GetYaxis()->SetTitleOffset(0.88);
  cadre->GetXaxis()->SetTitleOffset(0.9);
  cadre->GetYaxis()->SetTitle("Stub efficiency");
  cadre->GetYaxis()->SetNdivisions(511);
  cadre->GetYaxis()->SetLabelFont(42);
  cadre->GetYaxis()->SetLabelSize(0.035);
  cadre->GetYaxis()->SetTitleFont(42);
  cadre->Draw();

  
  D1_off->SetMarkerStyle(24);
  D2_off->SetMarkerStyle(20);
  D3_off->SetMarkerStyle(27);
  D4_off->SetMarkerStyle(34);
  D5_off->SetMarkerStyle(28);
    
  D1_off->SetMarkerSize(1.4);
  D2_off->SetMarkerSize(1.4);
  D3_off->SetMarkerSize(1.4);
  D4_off->SetMarkerSize(1.4);
  D5_off->SetMarkerSize(1.4);

  D1_off->Draw("same");
  D2_off->Draw("same");
  D3_off->Draw("same");
  D4_off->Draw("same");
  D5_off->Draw("same");

  TLegend *leg = new TLegend(0.55,0.2,0.9,0.45);

  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(D1_off,"TEDD double-disc 1","p");
  leg->AddEntry(D2_off,"TEDD double-disc 2","p");
  leg->AddEntry(D3_off,"TEDD double-disc 3","p");
  leg->AddEntry(D4_off,"TEDD double-disc 4","p");
  leg->AddEntry(D5_off,"TEDD double-disc 5","p");
  leg->Draw();
    
  char buffer[80];
       
  TLatex Tl;
  Tl.SetTextSize(0.03);
  Tl.DrawLatex(-0.02, 1.18, "CMS Phase-2 simulation");
    
  sprintf (buffer, "#sqrt{s}=14TeV, 0 PU");

  TLatex Tl2;
  Tl2.SetTextSize(0.03);
  Tl2.SetTextFont(52);
  Tl2.DrawLatex(5., 1.18, buffer);


  c2->Modified();
  c2->Update();
    
  if (per_lay==1)
  {
    sprintf (buffer, "Endcap_layer_eff.pdf");
    c2->Print(buffer);
    sprintf (buffer, "Endcap_layer_eff.C");
    c2->SaveSource(buffer);
  }
  else
  {
    sprintf (buffer, "Endcap_module_eff.pdf");
    c2->Print(buffer);
  }

}
/*

void do_eta_comp(std::string flat, std::string tilt, int layer)
{
  // First get the data
  // by merging all the available files

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
 
  float eta_val[50];

  float stub_flat_eta[20][50];
  float stub_tilt_eta[20][50];

  TFile *flatf    = TFile::Open(flat.c_str());
  TFile *tiltf    = TFile::Open(tilt.c_str());

  TTree *newflat = (TTree*)flatf->Get("Efficiencies");
  TTree *newtilt = (TTree*)tiltf->Get("Efficiencies");
    
  newflat->SetBranchAddress("eta_val",&eta_val);

  newflat->SetBranchAddress("stub_eta_lay",&stub_flat_eta);
  newtilt->SetBranchAddress("stub_eta_lay",&stub_tilt_eta);

  newflat->GetEntry(0);
  newtilt->GetEntry(0);

  TH2F *cadre    = new TH2F("zz","zz",300,-2.5,2.5,100,0.4,1.1);
  TH2F *stil_eff = new TH2F("steff","steff",300,-2.5,2.5,100,0.,1.02);
  TH2F *sfla_eff = new TH2F("sfeff","sfeff",300,-2.5,2.5,100,0.,1.02);

  for (int i=0;i<50;++i)
  { 
    if (stub_flat_eta[layer][i]!=0.) sfla_eff->Fill(eta_val[i],stub_flat_eta[layer][i]);
    if (stub_tilt_eta[layer][i]!=0.) stil_eff->Fill(eta_val[i],stub_tilt_eta[layer][i]);
  }

  c3 = new TCanvas("c3","Compare effs",0,0,640,480);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c3->Range(-2.9,-0.1,2.7,1.1);
  c3->SetFillColor(0);
  c3->SetBorderMode(0);
  c3->SetBorderSize(2);
  c3->SetGridx();
  c3->SetGridy();
  c3->SetLeftMargin(0.08);
  c3->SetRightMargin(0.05);
  c3->SetTopMargin(0.04);
  c3->SetBottomMargin(0.1);
  c3->SetFrameBorderMode(0);
  c3->SetFrameBorderMode(0);

  cadre->GetXaxis()->SetTitle("Tracking particle #eta");
  cadre->GetXaxis()->SetNdivisions(522);
    cadre->GetXaxis()->SetNdivisions(518);
    cadre->GetXaxis()->SetLabelFont(42);
    cadre->GetXaxis()->SetLabelSize(0.035);
    cadre->GetXaxis()->SetTitleFont(42);
    cadre->GetXaxis()->SetTitleSize(0.045);
    cadre->GetYaxis()->SetTitleSize(0.045);
    cadre->GetYaxis()->SetTitleOffset(0.64);
    cadre->GetXaxis()->SetTitleOffset(0.9);
    cadre->GetYaxis()->SetTitle("Stub efficiency");
    cadre->GetYaxis()->SetNdivisions(511);
    cadre->GetYaxis()->SetLabelFont(42);
    cadre->GetYaxis()->SetLabelSize(0.035);
    cadre->GetYaxis()->SetTitleFont(42);
    cadre->Draw();
  
  sfla_eff->SetMarkerStyle(20);
  stil_eff->SetMarkerStyle(4);
  sfla_eff->SetMarkerSize(1.25);
  stil_eff->SetMarkerSize(1.25);
  sfla_eff->Draw("same");
  stil_eff->Draw("same");

  leg = new TLegend(0.4,0.33,0.65,0.5);

  leg->SetTextSize(0.04);
  leg->AddEntry(sfla_eff,"Flat geometry","p");
  leg->AddEntry(stil_eff,"Tilted geometry","p");
  leg->Draw();
    
  TLatex Tl;
  Tl.SetTextSize(0.04);
  Tl.DrawLatex(-2.2, 1.06, "#it{CMS phase II simulation}");

    
  c3->Modified();
  c3->Update();

    char buffer[80];
  sprintf (buffer, "Plots/Barrel_layer_%d_effcomp.pdf", layer+5);
 
  c3->Print(buffer);
    
    sprintf (buffer, "Plots/Barrel_layer_%d_effcomp.C", layer+5);
    c3->SaveSource(buffer);
}
*/
void do_eta_plots(std::string filename,int layer,int per_lay)
{
  // First get the data
  // by merging all the available files
  
  if (layer<5) return;
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  float eta_val[50];
  
  float digi_et[30][50];    
  float clus_eta[30][50];
  float stub_eta[30][50];
   
  for (int i=0;i<30;++i)
  { 
    for (int j=0;j<50;++j)
    { 
      digi_et[i][j]  = 0;
      clus_eta[i][j] = 0;
      stub_eta[i][j] = 0;
    }
  }
 
  int disk =0;

  if (layer>=11  && layer<18) disk=int(layer-11)%8;
  if (layer>=18 && layer<25)  disk=-int(layer-18)%8;

  TChain *newtree = new TChain("Efficiencies");
  newtree->Add(filename.c_str());
  newtree->SetBranchAddress("eta_val",&eta_val);
  newtree->SetBranchAddress("digi_eta",&digi_et);
  newtree->SetBranchAddress("clus_eta",&clus_eta);

  (per_lay==1)
    ? newtree->SetBranchAddress("stub_eta_lay",&stub_eta)
    : newtree->SetBranchAddress("stub_eta",&stub_eta);  
    
  newtree->GetEntry(0);
    
  TH2F *cadre    = new TH2F("zz","zz",300,-2.5,2.5,200,0.,1.02);
  TH2F *digi_eff = new TH2F("deff","deff",300,-2.5,2.5,200,0.,1.02);
  TH2F *ceff = new TH2F("coeff","coeff",300,-2.5,2.5,200,0.,1.02);
  TH2F *seff = new TH2F("soeff","soeff",300,-2.5,2.5,200,0.,1.02);
    
  for (int i=0;i<50;++i)
  { 
    digi_eff->Fill(eta_val[i],digi_et[layer-5][i]);
    ceff->Fill(eta_val[i],clus_eta[layer-5][i]);
    seff->Fill(eta_val[i],stub_eta[layer-5][i]);
  }

  TCanvas *c2 = new TCanvas("c2","Official effs",5,75,670,660);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c2->Range(-1.833856,-0.1286626,21.1442,1.157964);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetBorderSize(2);
  c2->SetGridx();
  c2->SetGridy();
  c2->SetLeftMargin(0.08);
  c2->SetRightMargin(0.05);
  c2->SetFrameBorderMode(0);
  c2->SetFrameBorderMode(0);
    
  cadre->GetXaxis()->SetTitle("Particle #eta");
  cadre->GetXaxis()->SetNdivisions(522);
  cadre->GetXaxis()->SetLabelFont(42);
  cadre->GetXaxis()->SetLabelSize(0.03);
  cadre->GetXaxis()->SetTitleFont(42);
  cadre->GetYaxis()->SetTitle("Efficiency");
  cadre->GetYaxis()->SetNdivisions(511);
  cadre->GetYaxis()->SetLabelFont(42);
  cadre->GetYaxis()->SetLabelSize(0.03);
  cadre->GetYaxis()->SetTitleOffset(0.8);
  cadre->GetYaxis()->SetTitleFont(42);
  cadre->Draw();
    
  digi_eff->SetMarkerStyle(3);
  ceff->SetMarkerStyle(4);
  seff->SetMarkerStyle(20);
  digi_eff->Draw("same");
  ceff->Draw("same");
  seff->Draw("same");
    

  char buffer[80];

  if (disk==0)
  {
    sprintf(buffer, "Layer %d",layer);
  }
  else
  {
    (disk>0)
      ? sprintf(buffer, "Disk %d",disk)
      : sprintf(buffer, "Disk %d",-disk);
  }

  TLegend *leg = new TLegend(0.7,0.28,0.9,0.45);

  leg->SetTextSize(0.03);
  leg->SetHeader(buffer); 
  leg->AddEntry(digi_eff,"Digis","p");
  leg->AddEntry(ceff,"Clusters","p");
  leg->AddEntry(seff,"Stubs","p");
  leg->Draw();
    
  sprintf (buffer, "#sqrt{s}=14TeV, 0 PU");

  TLatex Tl;
  Tl.SetTextSize(0.04);
  Tl.DrawLatex(-2.5, 1.03, "CMS Phase-2 Simulation");

  TLatex Tl2;
  Tl2.SetTextSize(0.03);
  Tl2.SetTextFont(52);
  Tl2.DrawLatex(0., 1.03, buffer);

  c2->Modified();
  c2->Update();

  if (per_lay==1)
  {
    if (disk==0) sprintf (buffer, "Barrel_%d_layer_eta.png", layer);
    if (disk>0)  sprintf (buffer, "Endcap_p%d_layer_eta.png", disk); 
    if (disk<0)  sprintf (buffer, "Endcap_m%d_layer_eta.png", abs(disk)); 
  }
  else
  {
    if (disk==0) sprintf (buffer, "Barrel_%d_module_eta.png", layer+5);
    if (disk>0)  sprintf (buffer, "Endcap_p%d_module_eta.png", disk); 
    if (disk<0)  sprintf (buffer, "Endcap_m%d_module_eta.png", abs(disk)); 
  }
    
  c2->Print(buffer);  
}

