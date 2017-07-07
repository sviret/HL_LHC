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

#include <vector>
#include <list>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

void detector_acceptance(std::string filename, int nhits, float ptmin, float ptmax, float rhomax, float domax, float etamax, int pdid, int nevt)
{
    // First get the data
    // and link the file
    
    TChain *FullI = new TChain("FullInfo");
    
    FullI->Add(filename.c_str());
    
    //
    //1/ FullInfo TREE content:
    //
    // https://github.com/sviret/HL_LHC/blob/master/Utils/AM_ana/track_eff.h
    //
    
    
    int n_part;                        // The total number of particles inducing at least one stub in the event
    
    std::vector<int>     *part_pdg;       // PDG id of the particles
    std::vector<int>     *part_nstubs;    // How many stubs are induced by the particle in the tracker?
    std::vector<int>     *part_nhits;     // How many different layers/disks are hit by the particle?
    std::vector<int>     *part_nhits_fe;  // How many different layers/disks are hit by the particle and pass the FE cuts?
    
    std::vector<float>   *part_pt;     // pt of the particles
    std::vector<float>   *part_rho;    // rho0 of the particles
    std::vector<float>   *part_d0;     // d0 of the particles
    std::vector<float>   *part_z0;     // z0 of the particles
    std::vector<float>   *part_eta;    // eta of the particles
    std::vector<float>   *part_phi;    // phi of the particles
    std::vector<float>   *part_dist;
   
    part_pdg=0;
    part_nstubs=0;
    part_nhits=0;
    part_nhits_fe=0;
    
    part_pt=0;
    part_rho=0;
    part_d0=0;
    part_z0=0;
    part_eta=0;
    part_phi=0;
    part_dist=0;
    
    
    TEfficiency *myEff = new TEfficiency();
    
    FullI->SetBranchAddress("n_part",       &n_part);
    FullI->SetBranchAddress("part_pdg",     &part_pdg);
    FullI->SetBranchAddress("part_nstubs",  &part_nstubs);
    FullI->SetBranchAddress("part_nhits",   &part_nhits);
    FullI->SetBranchAddress("part_nhits_fe",&part_nhits_fe);

    FullI->SetBranchAddress("part_pt",      &part_pt);
    FullI->SetBranchAddress("part_rho",     &part_rho);
    FullI->SetBranchAddress("part_d0",      &part_d0);
    FullI->SetBranchAddress("part_z0",      &part_z0);
    FullI->SetBranchAddress("part_eta",     &part_eta);
    FullI->SetBranchAddress("part_phi",     &part_phi);
    FullI->SetBranchAddress("part_dist",    &part_dist);

    int n_entries = FullI->GetEntries();

    
    const int nbins = 25;
    
    bool insec,inpatt,intc,intrack;
    
    float min_pt  = ptmin;
    float min_eta = -etamax;
    float min_z0  = -15.;
    float max_z0  = 15.;
    float min_d0  = -5.;
    float max_d0  = 5.;

    float max_pt  = ptmax;
    float max_eta = etamax;

    
    
    float bin_pt  = (max_pt-min_pt)/nbins;
    float bin_eta = (max_eta-min_eta)/nbins;
    float bin_z0  = (max_z0-min_z0)/nbins;
    float bin_d0  = (max_d0-min_d0)/nbins;
    
    
    TH2F *eff_pt  = new TH2F("eff_pt","eff_pt",nbins,min_pt,max_pt,200,0.1,1.05);
    TH2F *eff_eta = new TH2F("eff_eta","eff_eta",nbins,min_eta,max_eta,200,0.1,1.05);
    TH2F *eff_z   = new TH2F("eff_z","eff_z",nbins,min_z0,max_z0,200,0.1,1.05);
    TH2F *eff_d   = new TH2F("eff_d","eff_d",nbins,min_d0,max_d0,200,0.1,1.05);
    
    float eff_val[5][20][nbins];
    float sig_val[5][20][nbins];
    float count[5][20][nbins];
    
    float eff_vals[5][20];
    
    float abscissa[5][2][nbins];
    
    for (int k=0;k<5;++k)
    {
        for (int l=0;l<20;++l)
        {
            eff_vals[k][l]  = 0.;
        }
    }
    
    for (int j=0;j<nbins;++j)
    {
        abscissa[0][0][j] = min_pt+(j+0.5)*bin_pt;
        abscissa[1][0][j] = min_eta+(j+0.5)*bin_eta;
        abscissa[3][0][j] = min_z0+(j+0.5)*bin_z0;
        abscissa[4][0][j] = min_d0+(j+0.5)*bin_d0;
        
        for (int k=0;k<5;++k)
        {
            abscissa[k][1][j]  = 0.;
            
            for (int l=0;l<20;++l)
            {
                eff_val[k][l][j]  = 0.;
                sig_val[k][l][j]  = 0.;
                count[k][l][j]    = 0.;
            }
        }
    }
    
    // Definition of eff_val values
    
    // eff_val[][0] : all the particles in the acceptance with at least 1 stub
    // eff_val[][1] : all the particles in [0] with at least nhits stubs
    // eff_val[][2] : all the particles in [1] with stubs in at least nhits layer/disks
    // eff_val[][3] : all the particles in [2] with stubs passing FE in at least nhits layer/disks
    //
    //
    
    
    // Loop over events
    
    int bin[5];
    
    float npart_0=0;
    float npart_1=0;
    float npart_2=0;
    float npart_3=0;
    float npart_4=0;
    
    for (int j=0;j<std::min(n_entries,nevt);++j)
    {
        if (j%1000==0) cout << j << endl;
        
        int n_comb_tot = 0;
        
        FullI->GetEntry(j); // Get the tree info
        
        // First of all count the number of patterns in the sector
        
        int npatt=0;
        int npatt_full=0;
        int ptbin;
        int inacc;
        int sectype;
        
        for (int k=0;k<n_part;++k) // Loop over all the particles
        {
            if (abs(part_pdg->at(k))!=pdid && pdid!=-1) continue;
	    if (part_nstubs->at(k)<1) continue;
            if (part_pt->at(k)<ptmin) continue;
	    if (part_pt->at(k)>ptmax) continue;
            if (part_rho->at(k)>rhomax) continue;
            if (fabs(part_d0->at(k))>domax) continue;

            if (part_z0->at(k)>max_z0) continue;
	    if (part_z0->at(k)<min_z0) continue;
            if (part_d0->at(k)>max_d0) continue;
            if (part_d0->at(k)<min_d0) continue;
            if (part_eta->at(k)>max_eta) continue;
            if (part_eta->at(k)<min_eta) continue;

            bin[0] = std::min(int((part_pt->at(k)-min_pt)/bin_pt),nbins-1);
            bin[1] = std::min(int((part_eta->at(k)-min_eta)/bin_eta),nbins-1);
	    bin[2] = 0;
            bin[3] = std::min(int((part_z0->at(k)-min_z0)/bin_z0),nbins-1);
            bin[4] = std::min(int((part_d0->at(k)-min_d0)/bin_d0),nbins-1);
            
            for (int l=0;l<5;++l) ++eff_val[l][0][bin[l]];
            
            if (part_nstubs->at(k)<nhits) continue;
            
            for (int l=0;l<5;++l) ++eff_val[l][1][bin[l]];

            if (part_nhits->at(k)<nhits) continue;            

            for (int l=0;l<5;++l) ++eff_val[l][2][bin[l]];

            if (part_nhits_fe->at(k)<nhits) continue;            

            for (int l=0;l<5;++l) ++eff_val[l][3][bin[l]];

        } // End of loop over particles
    } // End of loop on events
    
    for (int j=0;j<nbins;++j)
    {
        for (int k=0;k<5;++k)
        {
            for (int l=0;l<20;++l)
            {
                eff_vals[k][l] += eff_val[k][l][j];
            }
            
            if (eff_val[k][0][j]>=1) eff_val[k][6][j] = eff_val[k][1][j]/eff_val[k][0][j]; // SR
            if (eff_val[k][1][j]>=1) eff_val[k][7][j] = eff_val[k][2][j]/eff_val[k][1][j]; // SRNHITS
            if (eff_val[k][2][j]>=1) eff_val[k][8][j] = eff_val[k][3][j]/eff_val[k][2][j]; // SRNHITS_FE
	    if (eff_val[k][0][j]>=1) eff_val[k][9][j] = eff_val[k][3][j]/eff_val[k][0][j]; // ALL
            
            if (eff_val[k][0][j]>=1) sig_val[k][1][j] = sqrt(eff_val[k][6][j]*(1-eff_val[k][6][j])/eff_val[k][0][j]);
            if (eff_val[k][1][j]>=1) sig_val[k][2][j] = sqrt(eff_val[k][7][j]*(1-eff_val[k][7][j])/eff_val[k][1][j]);
            if (eff_val[k][2][j]>=1) sig_val[k][3][j] = sqrt(eff_val[k][8][j]*(1-eff_val[k][8][j])/eff_val[k][2][j]);            
            if (eff_val[k][0][j]>=1) sig_val[k][4][j] = sqrt(eff_val[k][9][j]*(1-eff_val[k][9][j])/eff_val[k][0][j]);            
        }
    }
    
    cout << eff_vals[0][0] << " particles pass the kinematic cuts and have at least one stub" << endl;
    cout << eff_vals[0][1] << " of them have at least " << nhits << " stubs" << endl;
    cout << eff_vals[0][2] << " of them have stubs in at least " << nhits << " layer disks" << endl;
    cout << eff_vals[0][3] << " of them have stubs passing FE cuts in at least " << nhits << " layer disks" << endl;
    
    
    TGraphErrors  *pt_stubs_eff     = new TGraphErrors(nbins,abscissa[0][0],eff_val[0][6],abscissa[0][1],sig_val[0][1]);
    TGraphErrors  *pt_tower_eff     = new TGraphErrors(nbins,abscissa[0][0],eff_val[0][7],abscissa[0][1],sig_val[0][2]);
    TGraphErrors  *pt_AM_eff        = new TGraphErrors(nbins,abscissa[0][0],eff_val[0][8],abscissa[0][1],sig_val[0][3]);
    
    TGraphErrors  *eta_stubs_eff     = new TGraphErrors(nbins,abscissa[1][0],eff_val[1][6],abscissa[1][1],sig_val[1][1]);
    TGraphErrors  *eta_tower_eff     = new TGraphErrors(nbins,abscissa[1][0],eff_val[1][7],abscissa[1][1],sig_val[1][2]);
    TGraphErrors  *eta_AM_eff        = new TGraphErrors(nbins,abscissa[1][0],eff_val[1][8],abscissa[1][1],sig_val[1][3]);

    TGraphErrors  *z_stubs_eff     = new TGraphErrors(nbins,abscissa[3][0],eff_val[3][6],abscissa[3][1],sig_val[3][1]);
    TGraphErrors  *z_tower_eff     = new TGraphErrors(nbins,abscissa[3][0],eff_val[3][7],abscissa[3][1],sig_val[3][2]);
    TGraphErrors  *z_AM_eff        = new TGraphErrors(nbins,abscissa[3][0],eff_val[3][8],abscissa[3][1],sig_val[3][3]);

    TGraphErrors  *d_stubs_eff     = new TGraphErrors(nbins,abscissa[4][0],eff_val[4][6],abscissa[4][1],sig_val[4][1]);
    TGraphErrors  *d_tower_eff     = new TGraphErrors(nbins,abscissa[4][0],eff_val[4][7],abscissa[4][1],sig_val[4][2]);
    TGraphErrors  *d_AM_eff        = new TGraphErrors(nbins,abscissa[4][0],eff_val[4][8],abscissa[4][1],sig_val[4][3]);

    TGraphErrors  *pt_inclusive_eff  = new TGraphErrors(nbins,abscissa[0][0],eff_val[0][9],abscissa[0][1],sig_val[0][4]);    
    TGraphErrors  *eta_inclusive_eff = new TGraphErrors(nbins,abscissa[1][0],eff_val[1][9],abscissa[1][1],sig_val[1][4]);    
    TGraphErrors  *z_inclusive_eff   = new TGraphErrors(nbins,abscissa[3][0],eff_val[3][9],abscissa[3][1],sig_val[3][4]);    
    TGraphErrors  *d_inclusive_eff   = new TGraphErrors(nbins,abscissa[4][0],eff_val[4][9],abscissa[4][1],sig_val[4][4]);    


    char buffer[80];
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    float low,up;
    
    
    TCanvas *c1 = new TCanvas("c1","pT efficiencies",201,77,1170,608);
    c1->Divide(3,1);
    
    for (int i=0;i<3;++i)
    {
        c1->cd(i+1);
        c1->cd(i+1)->Range(-5.887851,-1.930603,70.65421,17.37543);
        c1->cd(i+1)->SetFillColor(0);
        c1->cd(i+1)->SetBorderMode(0);
        c1->cd(i+1)->SetBorderSize(2);
        c1->cd(i+1)->SetGridx();
        c1->cd(i+1)->SetGridy();
        c1->cd(i+1)->SetLeftMargin(0.07692308);
        c1->cd(i+1)->SetTopMargin(0.07124352);
        c1->cd(i+1)->SetFrameBorderMode(0);
        c1->cd(i+1)->SetFrameBorderMode(0);
        eff_pt->GetXaxis()->SetTitle("Tracking particle p_{T} (in GeV/c)");
        eff_pt->GetYaxis()->SetTitle("Efficiency");
        eff_pt->Draw();
        
        if (i==0)
        {
            pt_stubs_eff->SetMarkerStyle(20);
            pt_stubs_eff->Draw("Psame");
            
            low = 100*(myEff->ClopperPearson(eff_vals[0][0],eff_vals[0][1],0.68,0)-eff_vals[0][1]/eff_vals[0][0]);
            up  = 100*(myEff->ClopperPearson(eff_vals[0][0],eff_vals[0][1],0.68,1)-eff_vals[0][1]/eff_vals[0][0]);
            
            TPaveText *pt1 = new TPaveText(min_pt+0.1*(max_pt-min_pt),0.26,min_pt+0.7*(max_pt-min_pt),0.34,"br");
            sprintf(buffer,"#epsilon_{stubs}=%.2f_{%.2f}^{+%.2f}%%",100*eff_vals[0][1]/eff_vals[0][0],low,up);
            pt1->AddText(buffer);
            pt1->SetTextFont(102);
            pt1->SetTextSize(0.04);
            pt1->Draw();                                    
        }
        
        if (i==1)
        {
            pt_tower_eff->SetMarkerStyle(20);
            pt_tower_eff->Draw("Psame");
            
            low = 100*(myEff->ClopperPearson(eff_vals[0][1],eff_vals[0][2],0.68,0)-eff_vals[0][2]/eff_vals[0][1]);
            up  = 100*(myEff->ClopperPearson(eff_vals[0][1],eff_vals[0][2],0.68,1)-eff_vals[0][2]/eff_vals[0][1]);
            
            TPaveText *pt1 = new TPaveText(min_pt+0.1*(max_pt-min_pt),0.26,min_pt+0.7*(max_pt-min_pt),0.34,"br");
            sprintf(buffer,"#epsilon_{layers}=%.2f_{%.2f}^{+%.2f}%%",100*eff_vals[0][2]/eff_vals[0][1],low,up);
            pt1->AddText(buffer);
            pt1->SetTextFont(102);
            pt1->SetTextSize(0.04);
            pt1->Draw();
        }
        
        if (i==2)
        {
            pt_AM_eff->SetMarkerStyle(20);
            pt_AM_eff->Draw("Psame");
            pt_inclusive_eff->SetMarkerStyle(24);
            pt_inclusive_eff->SetMarkerColor(4);
            pt_inclusive_eff->Draw("Psame");
            
            low = 100*(myEff->ClopperPearson(eff_vals[0][2],eff_vals[0][3],0.68,0)-eff_vals[0][3]/eff_vals[0][2]);
            up  = 100*(myEff->ClopperPearson(eff_vals[0][2],eff_vals[0][3],0.68,1)-eff_vals[0][3]/eff_vals[0][2]);
            
            TPaveText *pt1 = new TPaveText(min_pt+0.1*(max_pt-min_pt),0.26,min_pt+0.7*(max_pt-min_pt),0.34,"br");
            sprintf(buffer,"#epsilon_{layers_fe}=%.2f_{%.2f}^{+%.2f}%%",100*eff_vals[0][3]/eff_vals[0][2],low,up);
            pt1->AddText(buffer);
            pt1->SetTextFont(102);
            pt1->SetTextSize(0.04);
            pt1->Draw();
       
            low = 100*(myEff->ClopperPearson(eff_vals[0][0],eff_vals[0][4],0.68,0)-eff_vals[0][4]/eff_vals[0][0]);
            up  = 100*(myEff->ClopperPearson(eff_vals[0][0],eff_vals[0][4],0.68,1)-eff_vals[0][4]/eff_vals[0][0]);
            
            TPaveText *pt = new TPaveText(min_pt+0.1*(max_pt-min_pt),0.15,min_pt+0.9*(max_pt-min_pt),0.23,"br");
            sprintf(buffer,"#epsilon_{incl}^{N#geq%d,p_{T}#geq%.1f}=%.2f_{%.2f}^{+%.2f}%%",nhits,ptmin,100*eff_vals[0][3]/eff_vals[0][0],low,up);
            pt->AddText(buffer);
            pt->SetFillColor(1);
            pt->SetTextColor(0);
            pt->SetTextFont(112);
            pt->SetTextSize(0.05);
            pt->Draw();
        }    
    }
    
    c1->Update();
    

    TCanvas *c2 = new TCanvas("c2","eta efficiencies",201,77,1170,608);
    c2->Divide(3,1);
    
    for (int i=0;i<3;++i)
    {
        c2->cd(i+1);
        c2->cd(i+1)->Range(-5.887851,-1.930603,70.65421,17.37543);
        c2->cd(i+1)->SetFillColor(0);
        c2->cd(i+1)->SetBorderMode(0);
        c2->cd(i+1)->SetBorderSize(2);
        c2->cd(i+1)->SetGridx();
        c2->cd(i+1)->SetGridy();
        c2->cd(i+1)->SetLeftMargin(0.07692308);
        c2->cd(i+1)->SetTopMargin(0.07124352);
        c2->cd(i+1)->SetFrameBorderMode(0);
        c2->cd(i+1)->SetFrameBorderMode(0);
        eff_eta->GetXaxis()->SetTitle("Tracking particle #eta");
        eff_eta->GetYaxis()->SetTitle("Efficiency");
        eff_eta->Draw();
        
        if (i==0)
        {
            eta_stubs_eff->SetMarkerStyle(20);
            eta_stubs_eff->Draw("Psame");
            
            low = 100*(myEff->ClopperPearson(eff_vals[0][0],eff_vals[0][1],0.68,0)-eff_vals[0][1]/eff_vals[0][0]);
            up  = 100*(myEff->ClopperPearson(eff_vals[0][0],eff_vals[0][1],0.68,1)-eff_vals[0][1]/eff_vals[0][0]);
            
            TPaveText *eta1 = new TPaveText(min_eta+0.1*(max_eta-min_eta),0.26,min_eta+0.7*(max_eta-min_eta),0.34,"br");
            sprintf(buffer,"#epsilon_{stubs}=%.2f_{%.2f}^{+%.2f}%%",100*eff_vals[0][1]/eff_vals[0][0],low,up);
            eta1->AddText(buffer);
            eta1->SetTextFont(102);
            eta1->SetTextSize(0.04);
            eta1->Draw(); 
        }
        
        if (i==1)
        {
            eta_tower_eff->SetMarkerStyle(20);
            eta_tower_eff->Draw("Psame");

            low = 100*(myEff->ClopperPearson(eff_vals[0][1],eff_vals[0][2],0.68,0)-eff_vals[0][2]/eff_vals[0][1]);
            up  = 100*(myEff->ClopperPearson(eff_vals[0][1],eff_vals[0][2],0.68,1)-eff_vals[0][2]/eff_vals[0][1]);
            
            TPaveText *eta1 = new TPaveText(min_eta+0.1*(max_eta-min_eta),0.26,min_eta+0.7*(max_eta-min_eta),0.34,"br");
            sprintf(buffer,"#epsilon_{layers}=%.2f_{%.2f}^{+%.2f}%%",100*eff_vals[0][2]/eff_vals[0][1],low,up);
            eta1->AddText(buffer);
            eta1->SetTextFont(102);
            eta1->SetTextSize(0.04);
            eta1->Draw();
        }
        
        if (i==2)
        {
            eta_AM_eff->SetMarkerStyle(20);
            eta_AM_eff->Draw("Psame");
            eta_inclusive_eff->SetMarkerStyle(24);
            eta_inclusive_eff->SetMarkerColor(4);
            eta_inclusive_eff->Draw("Psame");

            
            low = 100*(myEff->ClopperPearson(eff_vals[0][2],eff_vals[0][3],0.68,0)-eff_vals[0][3]/eff_vals[0][2]);
            up  = 100*(myEff->ClopperPearson(eff_vals[0][2],eff_vals[0][3],0.68,1)-eff_vals[0][3]/eff_vals[0][2]);
            
            TPaveText *eta1 = new TPaveText(min_eta+0.1*(max_eta-min_eta),0.26,min_eta+0.7*(max_eta-min_eta),0.34,"br");
            sprintf(buffer,"#epsilon_{layers_fe}=%.2f_{%.2f}^{+%.2f}%%",100*eff_vals[0][3]/eff_vals[0][2],low,up);
            eta1->AddText(buffer);
            eta1->SetTextFont(102);
            eta1->SetTextSize(0.04);
            eta1->Draw();
       
            low = 100*(myEff->ClopperPearson(eff_vals[0][0],eff_vals[0][4],0.68,0)-eff_vals[0][4]/eff_vals[0][0]);
            up  = 100*(myEff->ClopperPearson(eff_vals[0][0],eff_vals[0][4],0.68,1)-eff_vals[0][4]/eff_vals[0][0]);
            
            TPaveText *eta = new TPaveText(min_eta+0.1*(max_eta-min_eta),0.15,min_eta+0.9*(max_eta-min_eta),0.23,"br");
            sprintf(buffer,"#epsilon_{incl}^{N#geq%d,p_{T}#geq%.1f}=%.2f_{%.2f}^{+%.2f}%%",nhits,ptmin,100*eff_vals[0][3]/eff_vals[0][0],low,up);
            eta->AddText(buffer);
            eta->SetFillColor(1);
            eta->SetTextColor(0);
            eta->SetTextFont(112);
            eta->SetTextSize(0.05);
            eta->Draw();

        }
    }
    
    c2->Update();
    
    
    
    TCanvas *c3 = new TCanvas("c3","z efficiencies",201,77,1170,608);
    c3->Divide(3,1);
    
    for (int i=0;i<3;++i)
    {
        c3->cd(i+1);
        c3->cd(i+1)->Range(-5.887851,-1.930603,70.65421,17.37543);
        c3->cd(i+1)->SetFillColor(0);
        c3->cd(i+1)->SetBorderMode(0);
        c3->cd(i+1)->SetBorderSize(2);
        c3->cd(i+1)->SetGridx();
        c3->cd(i+1)->SetGridy();
        c3->cd(i+1)->SetLeftMargin(0.07692308);
        c3->cd(i+1)->SetTopMargin(0.07124352);
        c3->cd(i+1)->SetFrameBorderMode(0);
        c3->cd(i+1)->SetFrameBorderMode(0);
        eff_z->GetXaxis()->SetTitle("Tracking particle #z");
        eff_z->GetYaxis()->SetTitle("Efficiency");
        eff_z->Draw();
        
        if (i==0)
        {
            z_stubs_eff->SetMarkerStyle(20);
            z_stubs_eff->Draw("Psame");
        }
        
        if (i==1)
        {
            z_tower_eff->SetMarkerStyle(20);
            z_tower_eff->Draw("Psame");
        }
        
        if (i==2)
        {
            z_AM_eff->SetMarkerStyle(20);
            z_AM_eff->Draw("Psame");
	    z_inclusive_eff->SetMarkerStyle(24);
	    z_inclusive_eff->SetMarkerColor(4);
	    z_inclusive_eff->Draw("Psame");
        }
    }
    
    c3->Update();
    
    
    
    
    TCanvas *c4 = new TCanvas("c4","d0 efficiencies",201,77,1170,608);
    c4->Divide(3,1);
    
    for (int i=0;i<3;++i)
    {
        c4->cd(i+1);
        c4->cd(i+1)->Range(-5.887851,-1.930603,70.65421,17.37543);
        c4->cd(i+1)->SetFillColor(0);
        c4->cd(i+1)->SetBorderMode(0);
        c4->cd(i+1)->SetBorderSize(2);
        c4->cd(i+1)->SetGridx();
        c4->cd(i+1)->SetGridy();
        c4->cd(i+1)->SetLeftMargin(0.07692308);
        c4->cd(i+1)->SetTopMargin(0.07124352);
        c4->cd(i+1)->SetFrameBorderMode(0);
        c4->cd(i+1)->SetFrameBorderMode(0);
        eff_d->GetXaxis()->SetTitle("Tracking particle d_{0}(in cm)");
        eff_d->GetYaxis()->SetTitle("Efficiency");
        eff_d->Draw();
        
        if (i==0)
        {
            d_stubs_eff->SetMarkerStyle(20);
            d_stubs_eff->Draw("Psame");
        }
        
        if (i==1)
        {
            d_tower_eff->SetMarkerStyle(20);
            d_tower_eff->Draw("Psame");
        }
        
        if (i==2)
        {
            d_AM_eff->SetMarkerStyle(20);
            d_AM_eff->Draw("Psame");
            d_inclusive_eff->SetMarkerStyle(24);
            d_inclusive_eff->SetMarkerColor(4);
            d_inclusive_eff->Draw("Psame");
        }
    }
    
    c4->Update();
   
    
    
}



void acceptance_study(std::string filename, int nevt)
{
    // First get the data
    // and link the file
    
    TChain *FullI = new TChain("FullInfo");
    
    FullI->Add(filename.c_str());
    
    //
    //1/ FullInfo TREE content:
    //
    // https://github.com/sviret/HL_LHC/blob/master/Utils/AM_ana/track_eff.h
    //
    
    
    int n_part;                        // The total number of particles inducing at least one stub in the event
    
    std::vector<int>     *part_pdg;       // PDG id of the particles
    std::vector<int>     *part_nstubs;    // How many stubs are induced by the particle in the tracker?
    std::vector<int>     *part_nhits;     // How many different layers/disks are hit by the particle?
    std::vector<int>     *part_nhits_fe;  // How many different layers/disks are hit by the particle and pass the FE cuts?
    
    std::vector<float>   *part_pt;     // pt of the particles
    std::vector<float>   *part_rho;    // rho0 of the particles
    std::vector<float>   *part_d0;     // d0 of the particles
    std::vector<float>   *part_z0;     // z0 of the particles
    std::vector<float>   *part_eta;    // eta of the particles
    std::vector<float>   *part_phi;    // phi of the particles
    std::vector<float>   *part_dist;
   
    part_pdg=0;
    part_nstubs=0;
    part_nhits=0;
    part_nhits_fe=0;
    
    part_pt=0;
    part_rho=0;
    part_d0=0;
    part_z0=0;
    part_eta=0;
    part_phi=0;
    part_dist=0;
    
    
    TEfficiency *myEff = new TEfficiency();
    
    FullI->SetBranchAddress("n_part",       &n_part);
    FullI->SetBranchAddress("part_pdg",     &part_pdg);
    FullI->SetBranchAddress("part_nstubs",  &part_nstubs);
    FullI->SetBranchAddress("part_nhits",   &part_nhits);
    FullI->SetBranchAddress("part_nhits_fe",&part_nhits_fe);

    FullI->SetBranchAddress("part_pt",      &part_pt);
    FullI->SetBranchAddress("part_rho",     &part_rho);
    FullI->SetBranchAddress("part_d0",      &part_d0);
    FullI->SetBranchAddress("part_z0",      &part_z0);
    FullI->SetBranchAddress("part_eta",     &part_eta);
    FullI->SetBranchAddress("part_phi",     &part_phi);
    FullI->SetBranchAddress("part_dist",    &part_dist);

    int n_entries = FullI->GetEntries();

    
    const int nbins = 24;
    
    bool insec,inpatt,intc,intrack;
    
    float min_pt  = 0;
    float min_z0  = -15.;
    float max_z0  = 15.;

    float max_pt  = 10;
    float max_eta = 2.4;
    
    float bin_pt  = (max_pt-min_pt)/nbins;
    float bin_eta = max_eta/nbins;
    
    float sig_val[5][20][nbins];
    float count[5][20][nbins];
    
    float eff_val[nbins][nbins][4][3][4];
    
    float eff_vals[4][3][4];

    float abscissa[5][2][nbins];
    
    for (int i=0;i<nbins;++i)
    {
      for (int j=0;j<nbins;++j)
      {
	for (int k=0;k<4;++k)
	{
	  for (int l=0;l<3;++l)
	  {
	    for (int m=0;m<4;++m)
	    {
	      eff_vals[k][l][m]      = 0.;
	      eff_val[i][j][k][l][m] = 0.;
	    }
	  }
	}
      }
    }
    
    // Definition of eff_val values
    
    // eff_val[ptbin][etabin][type][nhits][j]

    // type = 0 all primaries (rho<1cm)
    // type = 1 only mu (rho<1cm)
    // type = 2 only e  (rho<1cm)
    // type = 3 only pi (rho<1cm)
    //
    // nhits = 0 <-> 4
    // nhits = 1 <-> 5
    // nhits = 2 <-> 6
    //
    // Then
    //
    // eff_val[][][][][0] : all the particles in the acceptance with at least 1 stub
    // eff_val[][][][][1] : all the particles in [0] with at least nhits stubs
    // eff_val[][][][][2] : all the particles in [1] with stubs in at least nhits layer/disks
    // eff_val[][][][][3] : all the particles in [2] with stubs passing FE in at least nhits layer/disks
    //
    //
    
    
    // Loop over events
    
    int bin[5];
    
    float npart_0=0;
    float npart_1=0;
    float npart_2=0;
    float npart_3=0;
    float npart_4=0;
    
    for (int j=0;j<std::min(n_entries,nevt);++j)
    {
        if (j%1000==0) cout << j << endl;
        
        int n_comb_tot = 0;
        
        FullI->GetEntry(j); // Get the tree info
        
        // First of all count the number of patterns in the sector
        
        int npatt=0;
        int npatt_full=0;
        int ptbin;
        int inacc;
        int sectype;

	int etabin;
	int pdid,type,nstubs;

	int ns,nh,nhf;

        for (int k=0;k<n_part;++k) // Loop over all the particles
        {
	  if (part_rho->at(k)>1) continue;
	  if (part_nstubs->at(k)<1) continue;
	  if (part_pt->at(k)<min_pt) continue;
	  if (part_pt->at(k)>max_pt) continue;
	  if (part_z0->at(k)>max_z0) continue;
	  if (part_z0->at(k)<min_z0) continue;
	  if (std::abs(part_eta->at(k))>max_eta) continue;

	  ptbin  = std::min(int((part_pt->at(k)-min_pt)/bin_pt),nbins-1);
	  etabin = std::min(int((std::abs(part_eta->at(k)))/bin_eta),nbins-1);
	  pdid   = std::abs(part_pdg->at(k));
	  ns     = std::min(part_nstubs->at(k),6);
	  nh     = std::min(part_nhits->at(k),6);
	  nhf    = std::min(part_nhits_fe->at(k),6);
	  
	  type=0;
	  if (pdid==13)  type=1;
	  if (pdid==11)  type=2;
	  if (pdid==211) type=3;
	  
	  for (int l=0;l<3;++l)
	  {
	    ++eff_val[ptbin][etabin][0][l][0];
	    if (type!=0) ++eff_val[ptbin][etabin][type][l][0];

	    if (part_pt->at(k)<3) continue;
	    ++eff_vals[0][l][0];
	    if (type!=0) ++eff_vals[type][l][0];
	  }
	  
	  if (ns<4) continue;
            
	  for (int l=0;l<ns-3;++l)
	  {
	    ++eff_val[ptbin][etabin][0][l][1];
	    if (type!=0) ++eff_val[ptbin][etabin][type][l][1];

	    if (part_pt->at(k)<3) continue;
	    ++eff_vals[0][l][1];
	    if (type!=0) ++eff_vals[type][l][1];
	  }
	  
	  if (nh<4) continue;            

	  for (int l=0;l<nh-3;++l)
	  {
	    ++eff_val[ptbin][etabin][0][l][2];
	    if (type!=0) ++eff_val[ptbin][etabin][type][l][2];

	    if (part_pt->at(k)<3) continue;
	    ++eff_vals[0][l][2];
	    if (type!=0) ++eff_vals[type][l][2];
	  }
	  
	  if (nhf<4) continue;            
	  
	  for (int l=0;l<nhf-3;++l)
	  {
	    ++eff_val[ptbin][etabin][0][l][3];
	    if (type!=0) ++eff_val[ptbin][etabin][type][l][3];

	    if (part_pt->at(k)<3) continue;
	    ++eff_vals[0][l][3];
	    if (type!=0) ++eff_vals[type][l][3];
	  }
        } // End of loop over particles
    } // End of loop on events


  std::vector<TH2F*> eff_all; 
  std::vector<TH2F*> eff_mu; 
  std::vector<TH2F*> eff_ele; 
  std::vector<TH2F*> eff_pi; 

  for (int i=0;i<3;++i)
  {
    TH2F *eff_all_0  = new TH2F("eff_a0","eff_a0",nbins,0,max_eta,nbins,min_pt,max_pt);
    TH2F *eff_all_1  = new TH2F("eff_a1","eff_a1",nbins,0,max_eta,nbins,min_pt,max_pt);
    TH2F *eff_all_2  = new TH2F("eff_a2","eff_a2",nbins,0,max_eta,nbins,min_pt,max_pt);
    TH2F *eff_all_3  = new TH2F("eff_a3","eff_a3",nbins,0,max_eta,nbins,min_pt,max_pt);

    eff_all.push_back(eff_all_0);
    eff_mu.push_back(eff_all_1);
    eff_ele.push_back(eff_all_2);
    eff_pi.push_back(eff_all_3);
  }
    

  for (int j=0;j<nbins;++j)
  {
    for (int k=0;k<nbins;++k)
    {
      for (int kk=0;kk<4;++kk)
      {
	for (int l=0;l<3;++l)
	{
	  if (eff_val[j][k][kk][l][0]==0) continue;
	  
	  for (int m=1;m<4;++m) eff_val[j][k][kk][l][m] /= eff_val[j][k][kk][l][0];
	  
	  if (kk==0) eff_all[l]->Fill((k+0.5)*bin_eta,min_pt+(j+0.5)*bin_pt,eff_val[j][k][kk][l][2]);
	  if (kk==1) eff_mu[l]->Fill((k+0.5)*bin_eta,min_pt+(j+0.5)*bin_pt,eff_val[j][k][kk][l][2]);
	  if (kk==2) eff_ele[l]->Fill((k+0.5)*bin_eta,min_pt+(j+0.5)*bin_pt,eff_val[j][k][kk][l][2]);
	  if (kk==3) eff_pi[l]->Fill((k+0.5)*bin_eta,min_pt+(j+0.5)*bin_pt,eff_val[j][k][kk][l][2]);
	}
      }
    }
  }

    /*
    
    cout << eff_vals[0][0] << " particles pass the kinematic cuts and have at least one stub" << endl;
    cout << eff_vals[0][1] << " of them have at least " << nhits << " stubs" << endl;
    cout << eff_vals[0][2] << " of them have stubs in at least " << nhits << " layer disks" << endl;
    cout << eff_vals[0][3] << " of them have stubs passing FE cuts in at least " << nhits << " layer disks" << endl;
    */


  char buffer[80];
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetNumberContours(90);
  gStyle->SetPalette(kBlueYellow);
  float low,up;
  
    
  TCanvas *c1 = new TCanvas("c1","All tracks efficiencies",201,77,870,808);
  c1->Divide(1,3);
    
  for (int i=0;i<3;++i)
  {
    c1->cd(i+1);
    c1->cd(i+1)->SetFillColor(0);
    c1->cd(i+1)->SetGridx();
    c1->cd(i+1)->SetGridy();
    eff_all[i]->GetXaxis()->SetLabelSize(0.05);
    eff_all[i]->GetXaxis()->SetTitleSize(0.05);
    eff_all[i]->GetYaxis()->SetLabelSize(0.05);
    eff_all[i]->GetYaxis()->SetTitleSize(0.05);
    eff_all[i]->GetYaxis()->SetTitleOffset(0.5);
    eff_all[i]->SetAxisRange(0., 1.0,"Z");
    eff_all[i]->GetYaxis()->SetTitle("Particle p_{T} (GeV)");
    eff_all[i]->GetXaxis()->SetTitle("Particle #eta");
    eff_all[i]->Draw("cont4z");
      
    low = 100*(myEff->ClopperPearson(eff_vals[0][i][0],eff_vals[0][i][2],0.68,0)-eff_vals[0][i][2]/eff_vals[0][i][0]);
    up  = 100*(myEff->ClopperPearson(eff_vals[0][i][0],eff_vals[0][i][2],0.68,1)-eff_vals[0][i][2]/eff_vals[0][i][0]);
    
    TPaveText *eta1 = new TPaveText(0.25,0.25,0.5,0.5,"br");
    sprintf(buffer,"#epsilon_{p_{T}>3}^{N_{h}>%d}=%.2f_{%.2f}^{+%.2f}%%",i+4,100*eff_vals[0][i][2]/eff_vals[0][i][0],low,up);
    eta1->AddText(buffer);
    eta1->SetTextFont(102);
    eta1->SetTextSize(0.055);
    eta1->SetFillColor(0);
    eta1->Draw(); 
  }
  
  c1->Update();

  TCanvas *c2 = new TCanvas("c2","Muons efficiencies",201,77,870,808);
  c2->Divide(1,3);
    
  for (int i=0;i<3;++i)
  {
    c2->cd(i+1);
    c2->cd(i+1)->SetFillColor(0);
    c2->cd(i+1)->SetGridx();
    c2->cd(i+1)->SetGridy();
    eff_mu[i]->GetXaxis()->SetLabelSize(0.05);
    eff_mu[i]->GetXaxis()->SetTitleSize(0.05);
    eff_mu[i]->GetYaxis()->SetLabelSize(0.05);
    eff_mu[i]->GetYaxis()->SetTitleSize(0.05);
    eff_mu[i]->GetYaxis()->SetTitleOffset(0.5);
    eff_mu[i]->SetAxisRange(0., 1.0,"Z");
    eff_mu[i]->GetYaxis()->SetTitle("Particle p_{T} (GeV)");
    eff_mu[i]->GetXaxis()->SetTitle("Particle #eta");
    eff_mu[i]->Draw("cont4z");
    
    low = 100*(myEff->ClopperPearson(eff_vals[1][i][0],eff_vals[1][i][2],0.68,0)-eff_vals[1][i][2]/eff_vals[1][i][0]);
    up  = 100*(myEff->ClopperPearson(eff_vals[1][i][0],eff_vals[1][i][2],0.68,1)-eff_vals[1][i][2]/eff_vals[1][i][0]);
    
    TPaveText *eta1 = new TPaveText(0.25,0.25,0.5,0.5,"br");
    sprintf(buffer,"#epsilon_{p_{T}>3}^{N_{h}>%d}=%.2f_{%.2f}^{+%.2f}%%",i+4,100*eff_vals[1][i][2]/eff_vals[1][i][0],low,up);
    eta1->AddText(buffer);
    eta1->SetTextFont(102);
    eta1->SetTextSize(0.055);
    eta1->SetFillColor(0);
    eta1->Draw(); 
  }
  
  c2->Update();


  TCanvas *c3 = new TCanvas("c3","Electrons efficiencies",201,77,870,808);
  c3->Divide(1,3);
    
  for (int i=0;i<3;++i)
  {
    c3->cd(i+1);
    c3->cd(i+1)->SetFillColor(0);
    c3->cd(i+1)->SetGridx();
    c3->cd(i+1)->SetGridy();
    eff_ele[i]->GetXaxis()->SetLabelSize(0.05);
    eff_ele[i]->GetXaxis()->SetTitleSize(0.05);
    eff_ele[i]->GetYaxis()->SetLabelSize(0.05);
    eff_ele[i]->GetYaxis()->SetTitleSize(0.05);
    eff_ele[i]->GetYaxis()->SetTitleOffset(0.5);
    eff_ele[i]->SetAxisRange(0., 1.0,"Z");
    eff_ele[i]->GetYaxis()->SetTitle("Particle p_{T} (GeV)");
    eff_ele[i]->GetXaxis()->SetTitle("Particle #eta");
    eff_ele[i]->Draw("cont4z");
    
    low = 100*(myEff->ClopperPearson(eff_vals[2][i][0],eff_vals[2][i][2],0.68,0)-eff_vals[2][i][2]/eff_vals[2][i][0]);
    up  = 100*(myEff->ClopperPearson(eff_vals[2][i][0],eff_vals[2][i][2],0.68,1)-eff_vals[2][i][2]/eff_vals[2][i][0]);
    
    TPaveText *eta1 = new TPaveText(0.25,0.25,0.5,0.5,"br");
    sprintf(buffer,"#epsilon_{p_{T}>3}^{N_{h}>%d}=%.2f_{%.2f}^{+%.2f}%%",i+4,100*eff_vals[2][i][2]/eff_vals[2][i][0],low,up);
    eta1->AddText(buffer);
    eta1->SetTextFont(102);
    eta1->SetTextSize(0.055);
    eta1->SetFillColor(0);
    eta1->Draw(); 
  }
  
  c3->Update();
    

  TCanvas *c4 = new TCanvas("c4","Pions efficiencies",201,77,870,808);
  c4->Divide(1,3);
    
  for (int i=0;i<3;++i)
  {
    c4->cd(i+1);
    c4->cd(i+1)->SetFillColor(0);
    c4->cd(i+1)->SetGridx();
    c4->cd(i+1)->SetGridy();
    eff_pi[i]->GetXaxis()->SetLabelSize(0.05);
    eff_pi[i]->GetXaxis()->SetTitleSize(0.05);
    eff_pi[i]->GetYaxis()->SetLabelSize(0.05);
    eff_pi[i]->GetYaxis()->SetTitleSize(0.05);
    eff_pi[i]->GetYaxis()->SetTitleOffset(0.5);
    eff_pi[i]->SetAxisRange(0., 1.0,"Z");
    eff_pi[i]->GetYaxis()->SetTitle("Particle p_{T} (GeV)");
    eff_pi[i]->GetXaxis()->SetTitle("Particle #eta");
    eff_pi[i]->Draw("cont4z");
    
    low = 100*(myEff->ClopperPearson(eff_vals[3][i][0],eff_vals[3][i][2],0.68,0)-eff_vals[3][i][2]/eff_vals[3][i][0]);
    up  = 100*(myEff->ClopperPearson(eff_vals[3][i][0],eff_vals[3][i][2],0.68,1)-eff_vals[3][i][2]/eff_vals[3][i][0]);
    
    TPaveText *eta1 = new TPaveText(0.25,0.25,0.5,0.5,"br");
    sprintf(buffer,"#epsilon_{p_{T}>3}^{N_{h}>%d}=%.2f_{%.2f}^{+%.2f}%%",i+4,100*eff_vals[3][i][2]/eff_vals[3][i][0],low,up);
    eta1->AddText(buffer);
    eta1->SetTextFont(102);
    eta1->SetTextSize(0.055);
    eta1->SetFillColor(0);
    eta1->Draw(); 
  }
  
  c4->Update();
}
