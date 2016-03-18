#include <vector>
#include <list>
#ifdef __MAKECINT__
#pragma link C++ class vector<vector<int> >+;
#pragma link C++ class vector<vector<float> >+;
#endif

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

/*
  ROOT macros for pattern reco efficiencies visualization

  Use:

  root[1]-> .L MakePR_EffPlots.C

  Then you can choose between the following methods 
  knowing that filename denotes the extracted ROOT file 
  containing the PR results (created by the 
  sectorMaker program).

  Solution 1 do the test for one sector
  in order to avoid efficiency due to the overlap area we consider only the tracks
  have at least nh hits ONLY in the sector considered (so the track hitting the sector center)

  root[2]-> do_pattern_effs(std::string filename, int sec_num, int nh)

  where:

  filename: name of the root file created using the AM ana tool (see 6.2.2. of the tutorial)
  sec_num : number of the sector you consider (defined here (http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCData))
  nh      : minimal number of distinct layer/hits with stubs induced by the track in order to enter the efficiency calc (denominator)

  /////////////////////////////////////////////////////////////////////////////

  Solution 2 do the test for all the sectors
  root[2]-> do_full_effs(std::string filename, int nh)

  where:

  filename: name of the root file created using the AM ana tool (see 6.2.2. of the tutorial)
  nh      : minimal number of distinct layer/hits with stubs induced by the track in order to enter the efficiency calc (denominator)


  Information about this macro may be found in the following page:

  http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto620 (Part 6.2.2)

  Author: viret@in2p3_dot_fr
  Date: 26/06/2014
  Maj. update: 18/03/2016
*/

//
// Method plotting pattern efficiencies  
// for a given sector
//
// =>filename is the name of the root file containing the results
// =>nh is the minimum number of layer/disks containing primary stubs 
// from the primary particle 
// =>sec_num is the sector number


void do_pattern_effs(std::string filename, int sec_num, int nh, int pdgid=-1, float ptcut=2., float ptmax=100., float d0cut=1.,float etamax=2.4)
{
  // First get the data
  // by merging all the available files

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
 
    
  int n_part;                        // The total number of particles inducing at least one stub in the event
    
  std::vector<int>     *part_pdg;    // PDG id of the particles
  std::vector< std::vector<int> >     *part_sec;
  std::vector<int>     *part_nsec;   // In how many trigger towers this particle hit more than 4 different layers/disks?
  std::vector<int>     *part_nhits;  // How many different layers/disks are hit by the particle?
  std::vector<int>     *part_npatt;  // How many patterns contains more than 4 stubs of the particle (in 4 different layers/disks)?
    
  std::vector<float>   *part_pt;     // pt of the particles
  std::vector<float>   *part_d0;     // d0 of the particles
  std::vector<float>   *part_z0;     // z0 of the particles
  std::vector<float>   *part_eta;    // eta of the particles
  std::vector<float>   *part_phi;    // phi of the particles

  TFile *file    = TFile::Open(filename.c_str());
  TTree *newtree = (TTree*)file->Get("FullInfo");
    
  newtree->SetBranchAddress("n_part",       &n_part);
  newtree->SetBranchAddress("part_pdg",     &part_pdg);
  newtree->SetBranchAddress("part_nsec",    &part_nsec);
  newtree->SetBranchAddress("part_nhits",   &part_nhits);
  newtree->SetBranchAddress("part_npatt",   &part_npatt);
  newtree->SetBranchAddress("part_pt",      &part_pt);
  newtree->SetBranchAddress("part_d0",      &part_d0);
  newtree->SetBranchAddress("part_z0",      &part_z0);
  newtree->SetBranchAddress("part_eta",     &part_eta);
  newtree->SetBranchAddress("part_phi",     &part_phi);
  newtree->SetBranchAddress("part_sec",     &part_sec);

  char buffer[80];
    
  const int nbin_eta = 200;
  const int nbin_phi = 200;

  TH2F *tracks     = new TH2F("tracks","tracks",nbin_eta,-2.4,2.4,nbin_phi,-3.15.,3.15);
  TH2F *eff_pat    = new TH2F("eff_p","eff_p",nbin_eta,-2.4,2.4,nbin_phi,-3.15.,3.15);
  TH2F *pt_eff     = new TH2F("eff_pt","eff_pt",100,0.,100.,100,0.,1.1);
  TH2F *pt_eff_z   = new TH2F("eff_pt_z","eff_pt_z",50,0.,10.,100,0.,1.1);

  float eff_stub[nbin_phi][nbin_eta];
  float eff_sector[nbin_phi][nbin_eta];
  float entries_stub[nbin_phi][nbin_eta];
  float eff_pat_m[nbin_phi][nbin_eta];
  float eff_sec_m[nbin_phi][nbin_eta];
  float entries_stub[nbin_phi][nbin_eta];

  for (int i=0;i<nbin_phi;++i)
  {
    for (int j=0;j<nbin_eta;++j)
    {
      eff_stub[i][j]     = 0.;
      eff_sector[i][j]   = 0.;
      entries_stub[i][j] = 0.;
      eff_pat_m[i][j]    = 0.;
      eff_sec_m[i][j]    = 0.;
    }
  }

  float tot_in_sec    = 0.;
  float tot_match     = 0.;

  float pt, eta, phi;
    
  int i_bin;
  int j_bin;

  float PI=4.*atan(1.);

  float pt_in_sec[100];
  float pt_in_pat[100];

  float pt_in_sec_z[50];
  float pt_in_pat_z[50];
  float good_patt_z[50];
  float tot_patt_z[50];

  for (int i=0;i<100;++i)
  {
    pt_in_sec[i]=0.;
    pt_in_pat[i]=0.;
  }
  
  for (int i=0;i<50;++i)
  {
    pt_in_sec_z[i]=0.;
    pt_in_pat_z[i]=0.;
  }

  int pt_bin, pt_bin_z;
  int n_sec=0,in_sec=0;

  int n_entries = newtree->GetEntries();
  
  // Loop over all the interesting particles found in the event

  for (int i=0;i<n_entries;++i)
  {
    newtree->GetEntry(i);
    
    if (i%10000==0) 
      cout << "Processed  event " << i << "/" << n_entries << endl;

    if (n_part==0) continue;

    for (int j=0;j<n_part;++j)
    {
      pt    = part_pt->at(j);
      eta   = part_eta->at(j);
      phi   = part_phi->at(j);

      if (abs(part_pdg->at(j))!=pdgid && pdgid!=-1) continue;
      if (fabs(pt)>ptmax) continue;
      if (fabs(eta)>etamax) continue;
      if (part_nhits->at(j)<nh) continue;
      if (fabs(pt)<ptcut) continue;
      if (fabs(part_d0->at(j))>d0cut) continue;
    
      pt_bin=static_cast<int>(pt);
      pt_bin_z=static_cast<int>(5*pt);

      i_bin = static_cast<int>(nbin_phi*(phi+PI)/(2*PI));
      j_bin = static_cast<int>(nbin_eta*(eta+2.2)/4.4);

      // Check that the track is in the sector acceptance 
      if (part_nsec->at(j)<1) continue;
    
      n_sec=0;
      in_sec=0;

      for (unsigned int jj=0;jj<part_sec->at(j).size();++jj)
      {
          if (part_sec->at(j).at(jj)!=sec_num)	++n_sec;
          if (part_sec->at(j).at(jj)==sec_num)	++in_sec;
      }
    
      // If n_sec is different from 0, it means that the 
      // track is in the acceptance of other trigger towers
      // 
      // This could bias the efficiency calculation so we don't account for it
      //

      if (in_sec==0) continue;
      if (n_sec!=0) continue;

      ++pt_in_sec[pt_bin];
      
      eff_sec_m[i_bin][j_bin]+=1.;
      
      ++tot_in_sec;
        
      if (pt<10) ++pt_in_sec_z[pt_bin_z];
        
      if (part_npatt->at(j)!=0) // A pattern was matched for this track
      {
          ++pt_in_pat[pt_bin];
      
          eff_pat_m[i_bin][j_bin]+=1.;
          tot_match+=1;
      
          if (pt<10)
          {
              ++pt_in_pat_z[pt_bin_z];
          }
      }
    }	
  }

  cout << tot_in_sec << " tracks of sector " << sec_num << " have more than " << nh 
       << " hits in the sector and no more than " << nh << " hits in any other sector" << endl;
  cout << tot_match << " were matched in a pattern..." << endl;
  cout << " Efficiency is " << 100*tot_match/tot_in_sec << "%" << endl;
  
  for (int i=0;i<nbin_phi;++i)
  {
    for (int j=0;j<nbin_eta;++j)   
    {
      tracks->Fill(-2.2+(j+0.5)*4.4/nbin_eta,-PI+(i+0.5)*2*PI/nbin_phi,eff_sec_m[i][j]);

      if (eff_sec_m[i][j]!=0.)
      {
          eff_pat_m[i][j] /= eff_sec_m[i][j];
          eff_pat->Fill(-2.2+(j+0.5)*4.4/nbin_eta,-PI+(i+0.5)*2*PI/nbin_phi,eff_pat_m[i][j]);
      }
    }
  }

  for (int i=0;i<100;++i)
  {
    if (pt_in_sec[i]!=0.)
    {
      pt_in_pat[i]/=pt_in_sec[i];
      pt_eff->Fill(i+0.5,pt_in_pat[i]);
    }
  }

  for (int i=0;i<50;++i)
  {
    if (pt_in_sec_z[i]!=0.)
    {
      pt_in_pat_z[i]/=pt_in_sec_z[i];
      pt_eff_z->Fill(0.2*(i+0.5),pt_in_pat_z[i]);
    }
  }

  char buffer[80];
  
  c1 = new TCanvas("c1","Overall efficiencies",6,102,1526,921);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c1->Range(-1.833856,-0.1286626,21.1442,1.157964);
  c1->Divide(2,2);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetLeftMargin(0.08);
  c1->SetRightMargin(0.05);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderMode(0);
  
  c1->cd(1);  
  eff_pat->GetXaxis()->SetTitle("#eta");
  eff_pat->GetXaxis()->SetLabelFont(42);
  eff_pat->GetXaxis()->SetLabelSize(0.03);
  eff_pat->GetXaxis()->SetTitleSize(0.035);
  eff_pat->GetXaxis()->SetTitleFont(42);
  eff_pat->GetYaxis()->SetTitle("#phi (in rad)");
  eff_pat->GetYaxis()->SetLabelFont(42);
  eff_pat->GetYaxis()->SetLabelSize(0.03);
  eff_pat->GetYaxis()->SetTitleSize(0.035);
  eff_pat->GetYaxis()->SetTitleFont(42);
  eff_pat->SetContour(99);
  eff_pat->Draw("cont4z");

  TPaveText *pte = new TPaveText(-0.57,0.61,0.58,0.68,"br");
  pte->SetFillColor(0);
  pte->SetTextFont(42);
  pte->SetTextSize(0.0349944);
  sprintf (buffer, "Pattern recognition efficiency in the sector");  
  TText *text = pte->AddText(buffer);
  pte->Draw();
  

  c1->cd(2); 
  c1->cd(2)->SetGridx();
  c1->cd(2)->SetGridy();
  pt_eff->SetMinimum(-0.1);
  pt_eff->SetMaximum(1.1);
  pt_eff->SetMarkerStyle(3);
  pt_eff->GetXaxis()->SetTitle("Particle p_{T} (in GeV/c)");
  pt_eff->GetXaxis()->SetLabelFont(42);
  pt_eff->GetXaxis()->SetLabelSize(0.035);
  pt_eff->GetXaxis()->SetTitleSize(0.035);
  pt_eff->GetXaxis()->SetTitleOffset(1.19);
  pt_eff->GetXaxis()->SetTitleFont(42);
  pt_eff->GetYaxis()->SetTitle("Proportion of patterns");
  pt_eff->GetYaxis()->SetLabelFont(42);
  pt_eff->GetYaxis()->SetLabelSize(0.035);
  pt_eff->GetYaxis()->SetTitleSize(0.035);
  pt_eff->GetYaxis()->SetTitleFont(42);
  pt_eff->SetMarkerStyle(20);
  pt_eff->Draw("P");

  leg = new TLegend(0.6,0.5,0.85,0.65);
  leg->AddEntry(pt_eff,"Pattern reco efficiency","p");
  leg->Draw();

  c1->cd(4); 
  c1->cd(4)->SetGridx();
  c1->cd(4)->SetGridy();
  pt_eff_z->SetMinimum(-0.1);
  pt_eff_z->SetMaximum(1.1);
  pt_eff_z->SetMarkerStyle(3);
  pt_eff_z->GetXaxis()->SetTitle("Particle p_{T} (in GeV/c)");
  pt_eff_z->GetXaxis()->SetLabelFont(42);
  pt_eff_z->GetXaxis()->SetLabelSize(0.035);
  pt_eff_z->GetXaxis()->SetTitleSize(0.035);
  pt_eff_z->GetXaxis()->SetTitleOffset(1.19);
  pt_eff_z->GetXaxis()->SetTitleFont(42);
  pt_eff_z->GetYaxis()->SetTitle("Proportion of patterns");
  pt_eff_z->GetYaxis()->SetLabelFont(42);
  pt_eff_z->GetYaxis()->SetLabelSize(0.035);
  pt_eff_z->GetYaxis()->SetTitleSize(0.035);
  pt_eff_z->GetYaxis()->SetTitleFont(42);
  pt_eff_z->SetMarkerStyle(20);
  pt_eff_z->Draw("P");

  leg = new TLegend(0.6,0.5,0.85,0.65);
  leg->AddEntry(pt_eff_z,"Pattern reco efficiency","p");
  leg->Draw();

  c1->cd(3);  
  tracks->GetXaxis()->SetTitle("#eta");
  tracks->GetXaxis()->SetLabelFont(42);
  tracks->GetXaxis()->SetLabelSize(0.03);
  tracks->GetXaxis()->SetTitleSize(0.035);
  tracks->GetXaxis()->SetTitleFont(42);
  tracks->GetYaxis()->SetTitle("#phi (in rad)");
  tracks->GetYaxis()->SetLabelFont(42);
  tracks->GetYaxis()->SetLabelSize(0.03);
  tracks->GetYaxis()->SetTitleSize(0.035);
  tracks->GetYaxis()->SetTitleFont(42);
  tracks->SetContour(99);
  tracks->Draw("cont4z");

  TPaveText *pte = new TPaveText(-0.57,0.61,0.58,0.68,"br");
  pte->SetFillColor(0);
  pte->SetTextFont(42);
  pte->SetTextSize(0.0349944);
  sprintf (buffer, "Tracks collected in the sector");  
  TText *text = pte->AddText(buffer);
  pte->Draw();


  c1->Modified();
  c1->Update();


  sprintf (buffer, "Pattern_eff_sector_%d_%d_hits_track.png",sec_num,nh);

  c1->Print(buffer);
  

}



//
// Method plotting pattern efficiencies  
// for the full tracker
//
// =>filename is the name of the root file containing the results
// =>nh is the minimum number of layer/disks containing primary stubs 
// from the primary particle 
// =>plot_fake is an option that puts the fake rate vs pt info over
// the plot. This is relevant only for PGun events so could be
// turned off otherwise 
//
// Each entry corresponds to a given primary particle of the 
// input file. So there could be many entries with the same evt ID



void do_full_effs(std::string filename, int nh, int pdgid=-1, float ptcut=2., float d0cut=1., bool plot_fake=false)
{
  // First get the data
  // by merging all the available files

  //  gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    
    
    int n_part;                        // The total number of particles inducing at least one stub in the event
    
    std::vector<int>     *part_pdg;    // PDG id of the particles
    std::vector< std::vector<int> >     *part_sec;
    std::vector<int>     *part_nsec;   // In how many trigger towers this particle hit more than 4 different layers/disks?
    std::vector<int>     *part_nhits;  // How many different layers/disks are hit by the particle?
    std::vector<int>     *part_npatt;  // How many patterns contains more than 4 stubs of the particle (in 4 different layers/disks)?
    std::vector<int>     *part_ntc;    // How many TCs contains more than 4 stubs of the particle (in 4 different layers/disks)?
    
    
    std::vector<float>   *part_pt;     // pt of the particles
    std::vector<float>   *part_rho;    // rho0 of the particles
    std::vector<float>   *part_d0;     // d0 of the particles
    std::vector<float>   *part_z0;     // z0 of the particles
    std::vector<float>   *part_eta;    // eta of the particles
    std::vector<float>   *part_phi;    // phi of the particles
    std::vector<float>   *part_dist;   // isolation wrt other primaries with pT>3GeV
    
    int   evt;        // Event number (for PU event, where there is more than 1 primary)
  int   nsec;       // The number of sectors containing at least 5 stubs of the prim. track
  float pt;         // The pT of the prim. track
  float d0;         // The d0 of the prim. track
  float z0;         // The d0 of the prim. track
  float z0_f;         // The d0 of the prim. track
  float eta;        // The eta of the prim. track
  float phi;        // The phi of the prim. track
  float pt_f;       // The pT of the prim. track
  float eta_f;      // The eta of the prim. track
  float phi_f;      // The phi of the prim. track
  int   pdg;        // The pdg id of the prim. track
  int   npatt;      // The number of patterns containing at least 4 stubs of the prim. track
  int   ntotpatt;   // The total number of tracks 
  int   ntrck;      // The number of tracks containing at least 4 stubs of the prim. track
  int   ntottrck;   // The total number of tracks 
  int   mult[500];  // The total number of stubs per sector 
  int   nhits;      // The total number of layers/disks hit by the prim track

    
    
    TFile *file    = TFile::Open(filename.c_str());
    TTree *newtree = (TTree*)file->Get("FullInfo");
    
    newtree->SetBranchAddress("n_part",       &n_part);
    newtree->SetBranchAddress("part_pdg",     &part_pdg);
    newtree->SetBranchAddress("part_nsec",    &part_nsec);
    newtree->SetBranchAddress("part_nhits",   &part_nhits);
    newtree->SetBranchAddress("part_npatt",   &part_npatt);
    newtree->SetBranchAddress("part_ntc",     &part_ntc);
    newtree->SetBranchAddress("part_pt",      &part_pt);
    newtree->SetBranchAddress("part_rho",     &part_rho);
    newtree->SetBranchAddress("part_d0",      &part_d0);
    newtree->SetBranchAddress("part_z0",      &part_z0);
    newtree->SetBranchAddress("part_eta",     &part_eta);
    newtree->SetBranchAddress("part_phi",     &part_phi);
    newtree->SetBranchAddress("part_sec",     &part_sec);
    newtree->SetBranchAddress("part_dist",    &part_dist);
    
  int n_entries = newtree->GetEntries();

  const int nbin_eta = 50;
  const int nbin_phi = 50;
  const int nbin_pt  = 50;
  const int nbin_z0  = 50;


  TH2F *tracks     = new TH2F("tracks","tracks",nbin_eta,-2.5,2.5,nbin_phi,-3.15.,3.15);
  TH2F *eff_sec    = new TH2F("eff_s","eff_s",nbin_eta,-2.5,2.5,nbin_phi,-3.15.,3.15);
  TH2F *eff_pat    = new TH2F("eff_p","eff_p",nbin_eta,-2.5,2.5,nbin_phi,-3.15.,3.15);
  TH2F *eff_trk    = new TH2F("eff_t","eff_t",nbin_eta,-2.5,2.5,nbin_phi,-3.15.,3.15);
  TH2F *eff_tot    = new TH2F("eff_o","eff_o",nbin_eta,-2.5,2.5,nbin_phi,-3.15.,3.15);


  TH2F *pt_patt_eff   = new TH2F("eff_patt_pt","eff_patt_pt",nbin_pt,0.,100.,200,0.,1.1);
  TH2F *pt_trck_eff   = new TH2F("eff_trck_pt","eff_trck_pt",nbin_pt,0.,100.,200,0.,1.1);
  TH2F *pt_tot_eff    = new TH2F("eff_tot_pt","eff_tot_pt",nbin_pt,0.,100.,200,0.,1.1);

  TH2F *eta_patt_eff  = new TH2F("eff_patt_eta","eff_patt_eta",nbin_eta,-2.5,2.5,200,0.,1.1);
  TH2F *eta_trck_eff  = new TH2F("eff_trck_eta","eff_trck_eta",nbin_eta,-2.5,2.5,200,0.,1.1);
  TH2F *eta_tot_eff   = new TH2F("eff_tot_eta","eff_tot_eta",nbin_eta,-2.5,2.5,200,0.,1.1);

  TH2F *z0_patt_eff  = new TH2F("eff_patt_z0","eff_patt_z0",nbin_z0,-20.,20.,200,0.,1.1);
  TH2F *z0_trck_eff  = new TH2F("eff_trck_z0","eff_trck_z0",nbin_z0,-20.,20.,200,0.,1.1);
  TH2F *z0_tot_eff   = new TH2F("eff_tot_z0","eff_tot_z0",nbin_z0,-20.,20.,200,0.,1.1);

  TH2F *pt_patt_eff_z = new TH2F("eff_patt_pt_z","eff_patt_pt_z_z",nbin_pt,0.,10.,200,0.,1.1);
  TH2F *pt_trck_eff_z = new TH2F("eff_trck_pt_z","eff_trck_pt_z",nbin_pt,0.,10.,200,0.,1.1);
  TH2F *pt_tot_eff_z  = new TH2F("eff_tot_pt_z","eff_tot_pt_z",nbin_pt,0.,10.,200,0.,1.1);


  TH2F *pt_resol   = new TH2F("pt_resol","pt_resol",nbin_pt,0.,100.,200,0.,1.);
  TH2F *eta_resol  = new TH2F("eta_resol","eta_resol",nbin_pt,0.,100.,200,0.,1.);
  TH2F *phi_resol  = new TH2F("phi_resol","phi_resol",nbin_pt,0.,100.,200,0.,1.);
  TH2F *z0_resol   = new TH2F("z0_resol","z0_resol",nbin_pt,0.,100.,200,0.,10.);

  TH1F *pt_d   = new TH1F("pt_d","pt_d",100,-0.15,0.15);
  TH1F *eta_d  = new TH1F("eta_d","eta_d",100,-0.25,0.25);
  TH1F *phi_d  = new TH1F("phi_d","phi_d",100,-0.05,0.05);
  TH1F *z0_d   = new TH1F("z0_d","z0_d",100,-2.,2.);

  // [0] : all the tracks in the bin
  // [1] : all the tracks with at least nhits in one sector
  // [2] : all the tracks matched in one pattern
  // [3] : all the tracks fitted

  float pt_eff[nbin_pt][4];
  float pt_eff_z[nbin_pt][4];
  float eta_eff[nbin_eta][4];
  float z0_eff[nbin_z0][4];
  float e_p_eff[nbin_phi][nbin_eta][4];

  float pt_res[nbin_pt][4];
  float phi_res[nbin_pt][4];
  float eta_res[nbin_pt][4];
  float z0_res[nbin_pt][4];

    for (int k=0;k<4;++k)
    {
        for (int i=0;i<nbin_phi;++i)
        {
            for (int j=0;j<nbin_eta;++j)
            {
                e_p_eff[i][j][k] = 0.;
            }
        }

        for (int i=0;i<nbin_eta;++i) eta_eff[i][k]  = 0.;
        for (int i=0;i<nbin_pt;++i)  pt_eff[i][k]   = 0.;
        for (int i=0;i<nbin_pt;++i)  pt_eff_z[i][k] = 0.;
        for (int i=0;i<nbin_z0;++i)  z0_eff[i][k]   = 0.;

        for (int i=0;i<nbin_pt;++i)  pt_res[i][k]   = 0.;
        for (int i=0;i<nbin_pt;++i) eta_res[i][k]  = 0.;
        for (int i=0;i<nbin_pt;++i) phi_res[i][k]  = 0.;
        for (int i=0;i<nbin_pt;++i)  z0_res[i][k]   = 0.;

    }

    int tot_prims     = 0;
    int tot_in_sec    = 0;
    int tot_match     = 0;
    int tot_patts     = 0;
    int tot_goodpatts = 0;
    int tot_match_t   = 0;
    int tot_trcks     = 0;
    int tot_goodtrcks = 0;

    int prev_evt = -1;

    int i_bin;
    int j_bin;
    
    float PI=4.*atan(1.);

    int eta_bin,phi_bin,z0_bin,pt_bin, pt_bin_z;
    int n_sec=0;
  
    // First loop to get the sector eta/phi acceptance

    for (int i=0;i<n_entries;++i)
    //for (int i=0;i<1000;++i)
    {

        if (i%10000==0)
            cout << "Processed " << i << "/" << n_entries << endl;

        newtree->GetEntry(i);

        for (int j=0;j<n_part;++j)
        {
        pdg   = part_pdg->at(j);
        pt    = part_pt->at(j);
        eta   = part_eta->at(j);
        phi   = part_phi->at(j);
        nhits = part_nhits->at(j);
        d0    = part_d0->at(j);
        z0    = part_z0->at(j);
        
        if (pdg!=pdgid && pdgid!=-1) continue;
        if (fabs(pt)>100) continue;
        if (fabs(eta)>2.5) continue;
        if (nhits<nh) continue;
        if (fabs(pt)<ptcut) continue;
        if (fabs(d0)>d0cut) continue;
        if (fabs(z0)>20) continue;

        pt_bin  = static_cast<int>(nbin_pt*pt/100);
        pt_bin_z= static_cast<int>(nbin_pt*pt/10);
        z0_bin  = static_cast<int>(nbin_z0*(z0+20.)/40.);
        eta_bin = static_cast<int>(nbin_eta*(eta+2.5)/5.);
        phi_bin = static_cast<int>(nbin_phi*(phi+PI)/(2*PI));

        n_sec   = part_nsec->at(j);

        tot_prims++;
        ++pt_eff[pt_bin][0];
        if (pt_bin_z<nbin_pt) ++pt_eff_z[pt_bin_z][0];
        ++eta_eff[eta_bin][0];
        ++z0_eff[z0_bin][0];
        ++e_p_eff[phi_bin][eta_bin][0];

        if (n_sec==0) continue; // The track is not in the acceptance

        ++tot_in_sec;
        ++pt_eff[pt_bin][1];
        if (pt_bin_z<nbin_pt) ++pt_eff_z[pt_bin_z][1];
        ++eta_eff[eta_bin][1];
        ++z0_eff[z0_bin][1];
        ++e_p_eff[phi_bin][eta_bin][1];
  
        if (part_npatt->at(j)!=0) // This prim is in at least one pattern
        {
            ++tot_match;
            ++pt_eff[pt_bin][2];
            if (pt_bin_z<nbin_pt) ++pt_eff_z[pt_bin_z][2];
            ++eta_eff[eta_bin][2];
            ++z0_eff[z0_bin][2];
            ++e_p_eff[phi_bin][eta_bin][2];
        }

        if (part_ntc->at(j)!=0) // This prim is in at least one TC
        {
            ++tot_match_t;
            ++pt_eff[pt_bin][3];
            if (pt_bin_z<nbin_pt) ++pt_eff_z[pt_bin_z][3];
            ++eta_eff[eta_bin][3];
            ++z0_eff[z0_bin][3];
            ++e_p_eff[phi_bin][eta_bin][3];
        }
    }
    }
  
    cout << tot_in_sec << "/" << tot_prims << " tracks have more than " << nh << " hits in the full tracker" << endl;
    cout << tot_match << " were matched in a pattern..." << endl;
    cout << tot_match_t << " were matched in a TC..." << endl;

    for (int i=0;i<nbin_phi;++i)
    {
        for (int j=0;j<nbin_eta;++j)
        {
            if (e_p_eff[i][j][0]!=0.) eff_sec->Fill(-2.5+(j+0.5)*5./nbin_eta,-PI+(i+0.5)*2*PI/nbin_phi,e_p_eff[i][j][1]/e_p_eff[i][j][0]);
            if (e_p_eff[i][j][1]!=0.) eff_pat->Fill(-2.5+(j+0.5)*5./nbin_eta,-PI+(i+0.5)*2*PI/nbin_phi,e_p_eff[i][j][2]/e_p_eff[i][j][1]);
            if (e_p_eff[i][j][2]!=0.)
            {
                if (e_p_eff[i][j][3]/e_p_eff[i][j][2] > 1)
                {
                    eff_trk->Fill(-2.5+(j+0.5)*5./nbin_eta,-PI+(i+0.5)*2*PI/nbin_phi,1.);
                }
                else
                {
                    eff_trk->Fill(-2.5+(j+0.5)*5./nbin_eta,-PI+(i+0.5)*2*PI/nbin_phi,e_p_eff[i][j][3]/e_p_eff[i][j][2]);
                }
            }

            if (e_p_eff[i][j][1]!=0.) eff_tot->Fill(-2.5+(j+0.5)*5./nbin_eta,-PI+(i+0.5)*2*PI/nbin_phi,e_p_eff[i][j][3]/e_p_eff[i][j][1]);
        }
    }

    for (int i=0;i<nbin_pt;++i)
    {
    if (pt_eff[i][0]>10)
    {
      if (pt_eff[i][1]!=0.) pt_patt_eff->Fill((i+0.5)*100./nbin_pt,pt_eff[i][2]/pt_eff[i][1]);
      if (pt_eff[i][2]!=0.) pt_trck_eff->Fill((i+0.5)*100./nbin_pt,pt_eff[i][3]/pt_eff[i][2]);
      if (pt_eff[i][1]!=0.) pt_tot_eff->Fill((i+0.5)*100./nbin_pt,pt_eff[i][3]/pt_eff[i][1]);
    }

    if (pt_eff_z[i][0]>10)
    {
      if (pt_eff_z[i][1]!=0.) pt_patt_eff_z->Fill((i+0.5)*10./nbin_pt,pt_eff_z[i][2]/pt_eff_z[i][1]);
      if (pt_eff_z[i][2]!=0.) pt_trck_eff_z->Fill((i+0.5)*10./nbin_pt,pt_eff_z[i][3]/pt_eff_z[i][2]);
      if (pt_eff_z[i][1]!=0.) pt_tot_eff_z->Fill((i+0.5)*10./nbin_pt,pt_eff_z[i][3]/pt_eff_z[i][1]);
    }
  }

  for (int i=0;i<nbin_eta;++i)
  {
    if (eta_eff[i][0]>10)
    {
      if (eta_eff[i][1]!=0.) eta_patt_eff->Fill(-2.5+(i+0.5)*5./nbin_eta,eta_eff[i][2]/eta_eff[i][1]);
      if (eta_eff[i][2]!=0.) eta_trck_eff->Fill(-2.5+(i+0.5)*5./nbin_eta,eta_eff[i][3]/eta_eff[i][2]);
      if (eta_eff[i][1]!=0.) eta_tot_eff->Fill(-2.5+(i+0.5)*5./nbin_eta,eta_eff[i][3]/eta_eff[i][1]);
    }
  }

  for (int i=0;i<nbin_z0;++i)
  {
    if (z0_eff[i][0]>10)
    {
      if (z0_eff[i][1]!=0.) z0_patt_eff->Fill(-20.+(i+0.5)*40./nbin_z0,z0_eff[i][2]/z0_eff[i][1]);
      if (z0_eff[i][2]!=0.) z0_trck_eff->Fill(-20.+(i+0.5)*40./nbin_z0,z0_eff[i][3]/z0_eff[i][2]);
      if (z0_eff[i][1]!=0.) z0_tot_eff->Fill(-20.+(i+0.5)*40./nbin_z0,z0_eff[i][3]/z0_eff[i][1]);
    }
  }

  char buffer[80];


 
  c1 = new TCanvas("c1","Overall efficiencies",6,102,1526,921);
  gStyle->SetOptStat(0);  
  gStyle->SetOptTitle(0);
  c1->Range(-1.833856,-0.1286626,21.1442,1.157964);
  c1->Divide(2,2);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetLeftMargin(0.08);
  c1->SetRightMargin(0.05);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderMode(0);
  
  
  c1->cd(1);  
  eff_sec->GetXaxis()->SetTitle("#eta");
  eff_sec->GetXaxis()->SetLabelFont(42);
  eff_sec->GetXaxis()->SetLabelSize(0.03);
  eff_sec->GetXaxis()->SetTitleSize(0.035);
  eff_sec->GetXaxis()->SetTitleFont(42);
  eff_sec->GetYaxis()->SetTitle("#phi (in rad)");
  eff_sec->GetYaxis()->SetLabelFont(42);
  eff_sec->GetYaxis()->SetLabelSize(0.03);
  eff_sec->GetYaxis()->SetTitleSize(0.035);
  eff_sec->GetYaxis()->SetTitleFont(42);
  eff_sec->SetContour(99);
  eff_sec->Draw("cont4z");

  TPaveText *pte = new TPaveText(-0.57,0.61,0.58,0.68,"br");
  pte->SetFillColor(0);
  pte->SetTextFont(42);
  pte->SetTextSize(0.0349944);
  sprintf (buffer, "Trigger tower efficiency");  
  TText *text = pte->AddText(buffer);
  pte->Draw();

  c1->cd(2);  
  eff_pat->GetXaxis()->SetTitle("#eta");
  eff_pat->GetXaxis()->SetLabelFont(42);
  eff_pat->GetXaxis()->SetLabelSize(0.03);
  eff_pat->GetXaxis()->SetTitleSize(0.035);
  eff_pat->GetXaxis()->SetTitleFont(42);
  eff_pat->GetYaxis()->SetTitle("#phi (in rad)");
  eff_pat->GetYaxis()->SetLabelFont(42);
  eff_pat->GetYaxis()->SetLabelSize(0.03);
  eff_pat->GetYaxis()->SetTitleSize(0.035);
  eff_pat->GetYaxis()->SetTitleFont(42);
  eff_pat->SetContour(99);
  eff_pat->Draw("cont4z");

  TPaveText *pte = new TPaveText(-0.57,0.61,0.58,0.68,"br");
  pte->SetFillColor(0);
  pte->SetTextFont(42);
  pte->SetTextSize(0.0349944);
  sprintf (buffer, "Pattern recognition efficiency");  
  TText *text = pte->AddText(buffer);
  pte->Draw();
  
  
  c1->cd(3);  
  eff_trk->GetXaxis()->SetTitle("#eta");
  eff_trk->GetXaxis()->SetLabelFont(42);
  eff_trk->GetXaxis()->SetLabelSize(0.03);
  eff_trk->GetXaxis()->SetTitleSize(0.035);
  eff_trk->GetXaxis()->SetTitleFont(42);
  eff_trk->GetYaxis()->SetTitle("#phi (in rad)");
  eff_trk->GetYaxis()->SetLabelFont(42);
  eff_trk->GetYaxis()->SetLabelSize(0.03);
  eff_trk->GetYaxis()->SetTitleSize(0.035);
  eff_trk->GetYaxis()->SetTitleFont(42);
  eff_trk->SetContour(99);
  eff_trk->Draw("cont4z");

  TPaveText *pte = new TPaveText(-0.57,0.61,0.58,0.68,"br");
  pte->SetFillColor(0);
  pte->SetTextFont(42);
  pte->SetTextSize(0.0349944);
  sprintf (buffer, "Track fit efficiency");  
  TText *text = pte->AddText(buffer);
  pte->Draw();

  
  c1->cd(4);  
  eff_tot->GetXaxis()->SetTitle("#eta");
  eff_tot->GetXaxis()->SetLabelFont(42);
  eff_tot->GetXaxis()->SetLabelSize(0.03);
  eff_tot->GetXaxis()->SetTitleSize(0.035);
  eff_tot->GetXaxis()->SetTitleFont(42);
  eff_tot->GetYaxis()->SetTitle("#phi (in rad)");
  eff_tot->GetYaxis()->SetLabelFont(42);
  eff_tot->GetYaxis()->SetLabelSize(0.03);
  eff_tot->GetYaxis()->SetTitleSize(0.035);
  eff_tot->GetYaxis()->SetTitleFont(42);
  eff_tot->SetContour(99);
  eff_tot->Draw("cont4z");

  TPaveText *pte = new TPaveText(-0.57,0.61,0.58,0.68,"br");
  pte->SetFillColor(0);
  pte->SetTextFont(42);
  pte->SetTextSize(0.0349944);
  sprintf (buffer, "Overall L1 track efficiency");  
  TText *text = pte->AddText(buffer);
  pte->Draw();


  c1->Modified();
  c1->Update();


  sprintf (buffer, "L1_2D_summary_%d_hits_track.png",nh);

  c1->Print(buffer);
 
  sprintf (buffer, "L1_2D_summary_%d_hits_track.C",nh);

  c1->SaveSource(buffer);


  //  char buffer[80];
  
  c2 = new TCanvas("c2","Overall efficiencies",6,102,1526,921);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c2->Range(-1.833856,-0.1286626,21.1442,1.157964);
  c2->Divide(2,2);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetBorderSize(2);
  c2->SetLeftMargin(0.08);
  c2->SetRightMargin(0.05);
  c2->SetFrameBorderMode(0);
  c2->SetFrameBorderMode(0);


  c2->cd(1); 
  c2->cd(1)->SetGridx();
  c2->cd(1)->SetGridy();
  pt_tot_eff->SetMinimum(-0.1);
  pt_tot_eff->SetMaximum(1.1);
  pt_tot_eff->SetMarkerStyle(3);
  pt_tot_eff->GetXaxis()->SetTitle("Particle p_{T} (in GeV/c)");
  pt_tot_eff->GetXaxis()->SetLabelFont(42);
  pt_tot_eff->GetXaxis()->SetLabelSize(0.035);
  pt_tot_eff->GetXaxis()->SetTitleSize(0.035);
  pt_tot_eff->GetXaxis()->SetTitleOffset(1.19);
  pt_tot_eff->GetXaxis()->SetTitleFont(42);
  pt_tot_eff->GetYaxis()->SetTitle("Proportion of tracks reconstructed");
  pt_tot_eff->GetYaxis()->SetLabelFont(42);
  pt_tot_eff->GetYaxis()->SetLabelSize(0.035);
  pt_tot_eff->GetYaxis()->SetTitleSize(0.035);
  pt_tot_eff->GetYaxis()->SetTitleFont(42);
  pt_tot_eff->GetYaxis()->SetRangeUser(-0.1,1.1);
  pt_tot_eff->SetMarkerStyle(20);
  pt_tot_eff->Draw();
  pt_patt_eff->SetMarkerStyle(24);
  pt_patt_eff->Draw("same");
  pt_trck_eff->SetMarkerStyle(3);
  //pt_trck_eff->Draw("same");


  TLatex Tl;
  Tl.SetTextSize(0.05);
  Tl.DrawLatex(30., 0.25, "CMS Preliminary Simulation");  
  c1->Modified();

  leg = new TLegend(0.7,0.85,0.95,1.);
  leg->AddEntry(pt_patt_eff,"AM pattern reco efficiency","p");
  //leg->AddEntry(pt_trck_eff,"Hough track fit efficiency","p");
  leg->AddEntry(pt_tot_eff,"Overall L1 track reco efficiency","p");
  leg->Draw();


  c2->cd(2); 
  c2->cd(2)->SetGridx();
  c2->cd(2)->SetGridy();
  pt_tot_eff_z->SetMinimum(-0.1);
  pt_tot_eff_z->SetMaximum(1.1);
  pt_tot_eff_z->SetMarkerStyle(3);
  pt_tot_eff_z->GetXaxis()->SetTitle("Particle p_{T} (in GeV/c)");
  pt_tot_eff_z->GetXaxis()->SetLabelFont(42);
  pt_tot_eff_z->GetXaxis()->SetLabelSize(0.035);
  pt_tot_eff_z->GetXaxis()->SetTitleSize(0.035);
  pt_tot_eff_z->GetXaxis()->SetTitleOffset(1.19);
  pt_tot_eff_z->GetXaxis()->SetTitleFont(42);
  pt_tot_eff_z->GetYaxis()->SetTitle("Proportion of tracks reconstructed");
  pt_tot_eff_z->GetYaxis()->SetLabelFont(42);
  pt_tot_eff_z->GetYaxis()->SetLabelSize(0.035);
  pt_tot_eff_z->GetYaxis()->SetTitleSize(0.035);
  pt_tot_eff_z->GetYaxis()->SetTitleFont(42);
  pt_tot_eff_z->GetYaxis()->SetRangeUser(-0.1,1.1);
  pt_tot_eff_z->SetMarkerStyle(20);
  pt_tot_eff_z->Draw();
  pt_patt_eff_z->SetMarkerStyle(24);
  pt_patt_eff_z->Draw("same");


  TLatex Tl;
  Tl.SetTextSize(0.05);
  Tl.DrawLatex(3., 0.25, "CMS Preliminary Simulation");  
  c1->Modified();

  leg = new TLegend(0.7,0.85,0.95,1.);
  leg->AddEntry(pt_patt_eff_z,"AM pattern reco efficiency","p");
  leg->AddEntry(pt_tot_eff_z,"Overall L1 track reco efficiency","p");
  leg->Draw();


  c2->cd(3); 
  c2->cd(3)->SetGridx();
  c2->cd(3)->SetGridy();
  eta_tot_eff->SetMinimum(-0.1);
  eta_tot_eff->SetMaximum(1.1);
  eta_tot_eff->SetMarkerStyle(3);
  eta_tot_eff->GetXaxis()->SetTitle("Particle #eta");
  eta_tot_eff->GetXaxis()->SetLabelFont(42);
  eta_tot_eff->GetXaxis()->SetLabelSize(0.035);
  eta_tot_eff->GetXaxis()->SetTitleSize(0.035);
  eta_tot_eff->GetXaxis()->SetTitleOffset(1.19);
  eta_tot_eff->GetXaxis()->SetTitleFont(42);
  eta_tot_eff->GetYaxis()->SetTitle("Proportion of tracks reconstructed");
  eta_tot_eff->GetYaxis()->SetLabelFont(42);
  eta_tot_eff->GetYaxis()->SetLabelSize(0.035);
  eta_tot_eff->GetYaxis()->SetTitleSize(0.035);
  eta_tot_eff->GetYaxis()->SetTitleFont(42);
  eta_tot_eff->GetYaxis()->SetRangeUser(-0.1,1.1);
  eta_tot_eff->SetMarkerStyle(20);
  eta_tot_eff->Draw();
  eta_patt_eff->SetMarkerStyle(24);
  eta_patt_eff->Draw("same");
  eta_trck_eff->SetMarkerStyle(3);
  //eta_trck_eff->Draw("same");


  TLatex Tl;
  Tl.SetTextSize(0.05);
  Tl.DrawLatex(-1., 0.25, "CMS Preliminary Simulation");  
  c1->Modified();

  leg = new TLegend(0.7,0.85,0.95,1.);
  leg->AddEntry(eta_patt_eff,"AM pattern reco efficiency","p");
  //leg->AddEntry(eta_trck_eff,"Hough track fit efficiency","p");
  leg->AddEntry(eta_tot_eff,"Overall L1 track reco efficiency","p");
  leg->Draw();


  c2->cd(4); 
  c2->cd(4)->SetGridx();
  c2->cd(4)->SetGridy();
  z0_tot_eff->SetMinimum(-0.1);
  z0_tot_eff->SetMaximum(1.1);
  z0_tot_eff->SetMarkerStyle(3);
  z0_tot_eff->GetXaxis()->SetTitle("Particle z0 (in cm)");
  z0_tot_eff->GetXaxis()->SetLabelFont(42);
  z0_tot_eff->GetXaxis()->SetLabelSize(0.035);
  z0_tot_eff->GetXaxis()->SetTitleSize(0.035);
  z0_tot_eff->GetXaxis()->SetTitleOffset(1.19);
  z0_tot_eff->GetXaxis()->SetTitleFont(42);
  z0_tot_eff->GetYaxis()->SetTitle("Proportion of tracks reconstructed");
  z0_tot_eff->GetYaxis()->SetLabelFont(42);
  z0_tot_eff->GetYaxis()->SetLabelSize(0.035);
  z0_tot_eff->GetYaxis()->SetTitleSize(0.035);
  z0_tot_eff->GetYaxis()->SetTitleFont(42);
  z0_tot_eff->GetYaxis()->SetRangeUser(-0.1,1.1);
  z0_tot_eff->SetMarkerStyle(20);
  z0_tot_eff->Draw();
  z0_patt_eff->SetMarkerStyle(24);
  z0_patt_eff->Draw("same");
  //z0_trck_eff->SetMarkerStyle(3);
  //z0_trck_eff->Draw("same");


  TLatex Tl;
  Tl.SetTextSize(0.05);
  Tl.DrawLatex(-7., 0.25, "CMS Preliminary Simulation");  
  c1->Modified();

  leg = new TLegend(0.7,0.85,0.95,1.);
  leg->AddEntry(z0_patt_eff,"AM pattern reco efficiency","p");
  //leg->AddEntry(z0_trck_eff,"Hough track fit efficiency","p");
  leg->AddEntry(z0_tot_eff,"Overall L1 track reco efficiency","p");
  leg->Draw();


  c2->Modified();
  c2->Update();


  sprintf (buffer, "Track_eff_%d_hits_track.png",nh);

  c2->Print(buffer);
 
  sprintf (buffer, "Track_eff_%d_hits_track.C",nh);

  c2->SaveSource(buffer);
}

