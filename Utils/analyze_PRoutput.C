#include <vector>
#ifdef __MAKECINT__
#pragma link C++ class vector<vector<int> >+;
#pragma link C++ class vector<vector<float> >+;
#endif

/*
  Small ROOT macro showing how to analyze the PR output.

  Use:

  root[1]-> .L PR_output.C
  root[2]-> am_result(filename,evtnum) 

  where filename is the name of the extracted ROOT file and evtnum the number of entry to analyze

  Information about this macro may be found in the following page:

  http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto (STEP V)

  Author: viret@in2p3_dot_fr
  Date: 02/10/2013

*/

void am_result(std::string filename, int evtnum)
{
  // First get the data
  // by merging all the available files

  TChain *L1TT            = new TChain("FullInfo");  

  L1TT->Add(filename.c_str());

//
//1/ FullInfo TREE content:
//
// https://github.com/sviret/HL_LHC/blob/master/Extractors/RecoExtractor/test/SectorMaker/sector_test.h
//


  int evt;             // Event number (for PU event, where there is more than 1 primary)
  int n_stub_total;    // The total number of stubs in the event
  int n_stub;          // The total number of stubs contained in matched patterns in the event

  std::vector<float>   *stub_x;      // x coordinates of ALL the stubs
  std::vector<float>   *stub_y;      // y coordinates of ALL the stubs
  std::vector<float>   *stub_z;      // z coordinates of ALL the stubs
  std::vector<int>     *stub_layer;  // layer number of ALL the stubs
  std::vector<int>     *stub_ladder; // ladder number of ALL the stubs
  std::vector<int>     *stub_module; // module number of ALL the stubs
  std::vector<int>     *stub_tp;     // tp index of ALL the stubs (in part_*** vectors of this tree!!!!)
  std::vector<int>     *stub_inpatt; // is the stub in a pattern (1) of not (0)?

  int n_part;                        // The total number of particles inducing at least one stub in the event

  std::vector<int>     *part_pdg;    // PDG id of the particles
  std::vector<int>     *part_nsec;   // In how many trigger towers this particle hit more than 4 different layers/disks?
  std::vector<int>     *part_nhits;  // How many different layers/disks are hit by the particle?
  std::vector<int>     *part_npatt;  // How many patterns contains more than 4 stubs of the particle (in 4 different layers/disks)?
  std::vector<float>   *part_pt;     // pt of the particles
  std::vector<float>   *part_rho;    // rho0 of the particles
  std::vector<float>   *part_z0;     // z0 of the particles
  std::vector<float>   *part_eta;    // eta of the particles 
  std::vector<float>   *part_phi;    // phi of the particles

  int n_patt;                        // The total number of patterns matched in the event

  // Sector id of all the patterns
  std::vector<int>                  *patt_sec;   

  // tp index of ALL the particles contained in the pattern (in part_*** vectors of this tree!!!!)
  std::vector< std::vector<int> >   *patt_parts; 

  // index of ALL the stubs contained in the pattern (in stub_*** vectors of this tree!!!!) 
  std::vector< std::vector<int> >   *patt_stubs; 

  L1TT->SetBranchAddress("evt",          &evt); 

  L1TT->SetBranchAddress("n_stub_total", &n_stub_total); 
  L1TT->SetBranchAddress("n_stub_inpat", &n_stub); 
  L1TT->SetBranchAddress("stub_x",       &stub_x); 
  L1TT->SetBranchAddress("stub_y",       &stub_y); 
  L1TT->SetBranchAddress("stub_z",       &stub_z); 
  L1TT->SetBranchAddress("stub_layer",   &stub_layer); 
  L1TT->SetBranchAddress("stub_ladder",  &stub_ladder);
  L1TT->SetBranchAddress("stub_module",  &stub_module);
  L1TT->SetBranchAddress("stub_tp",      &stub_tp);
  L1TT->SetBranchAddress("stub_inpatt",  &stub_inpatt);

  L1TT->SetBranchAddress("n_part",       &n_part); 
  L1TT->SetBranchAddress("part_pdg",     &part_pdg); 
  L1TT->SetBranchAddress("part_nsec",    &part_nsec); 
  L1TT->SetBranchAddress("part_nhits",   &part_nhits); 
  L1TT->SetBranchAddress("part_npatt",   &part_npatt); 
  L1TT->SetBranchAddress("part_pt",      &part_pt); 
  L1TT->SetBranchAddress("part_rho",     &part_rho);
  L1TT->SetBranchAddress("part_z0",      &part_z0);
  L1TT->SetBranchAddress("part_eta",     &part_eta);  
  L1TT->SetBranchAddress("part_phi",     &part_phi); 

  L1TT->SetBranchAddress("n_patt",       &n_patt); 
  L1TT->SetBranchAddress("patt_sec",     &patt_sec); 
  L1TT->SetBranchAddress("patt_parts",   &patt_parts); 
  L1TT->SetBranchAddress("patt_stubs",   &patt_stubs); 

  int n_entries = L1TT->GetEntries();

  if (evtnum>=n_entries) evtnum = n_entries-1;
  if (evtnum<0) evtnum = 0;


  L1TT->GetEntry(evtnum);

  cout <<endl;
  cout << "Analyzing event " << evt <<endl;
  cout << "where " << n_stub_total << " stub(s) were produced" <<endl;
  cout << n_part << " particle(s) have induced at least one stub in the tracker" <<endl;
  cout << "      " << n_stub << " stub(s) are contained in the " << n_patt << " pattern(s) matched in this event" <<endl;
  cout << endl;
  cout << "Now looping over the patterns... " <<endl;
  cout << endl;

  // Loop over patterns

  int idx_s;
  int idx_p;

  for (int k=0;k<n_patt;++k)
  {
    cout << "-------------------------------------------------"  <<endl;
    cout << "Pattern " << k+1 << " properties:"  <<endl;
    cout << "=> Sector id : " << patt_sec->at(k) <<endl;
    cout << "=> Number of stubs : " << patt_stubs->at(k).size() <<endl;
    cout << "=> Number of particles w/more than four stubs in the pattern : " << patt_parts->at(k).size() <<endl;

    if (patt_parts->at(k).size()==0)
    {
      cout << "!! FAKE PATTERN containing the following stubs: " <<endl;

      for (int kk=0;kk<patt_stubs->at(k).size();++kk)
      {
	idx_s = patt_stubs->at(k).at(kk);
	idx_p = stub_tp->at(idx_s);

	cout << " Stub " << kk+1 << endl;  
	cout << " X/Y/Z (in cm)       : " << stub_x->at(idx_s) 
	     << "/" << stub_y->at(idx_s) 
	     << "/" << stub_z->at(idx_s) << endl;

	if (stub_tp->at(idx_s)>=0)  // The cluster is matched
	{
	  cout << "    Matched with PART   : " << idx_p << endl;
	  cout << "    PTgen     : " << part_pt->at(idx_p) 
	       << endl;
	  
	  cout << "    PART origin (R/Z) : " << part_rho->at(idx_p) << " / " 
	       << part_z0->at(idx_p) << endl; 
	  cout << "    PART pdg code       : " << part_pdg->at(idx_p) << endl; 
	}
	else
	{
	  cout << "Unmatched" << endl;
	}
      }
    }
    else
    {
      for (int jj=0;jj<patt_parts->at(k).size();++jj)
      {
	cout << "Particle: " << jj+1 <<endl;
      }

      cout << "This patterns contains the following stubs: " <<endl;

      for (int kk=0;kk<patt_stubs->at(k).size();++kk)
      {
	idx_s = patt_stubs->at(k).at(kk);
	idx_p = stub_tp->at(idx_s);

	cout << " Stub " << kk+1 << endl;  
	cout << " X/Y/Z (in cm)       : " << stub_x->at(idx_s) 
	     << "/" << stub_y->at(idx_s) 
	     << "/" << stub_z->at(idx_s) << endl;

	if (stub_tp->at(idx_s)>=0)  // The cluster is matched
	{
	  cout << "    Matched with PART   : " << idx_p << endl;
	  cout << "    PTgen     : " << part_pt->at(idx_p) 
	       << endl;
	  
	  cout << "    PART origin (R/Z) : " << part_rho->at(idx_p) << " / " 
	       << part_z0->at(idx_p) << endl; 
	  cout << "    PART pdg code       : " << part_pdg->at(idx_p) << endl; 
	}
	else
	{
	  cout << "Unmatched" << endl;
	}
      }
    }
  }
}


void fake_ana(std::string filename, float ptmin, float d0max)
{
  // First get the data
  // by merging all the available files

  TChain *L1TT            = new TChain("FullInfo");  

  L1TT->Add(filename.c_str());

//
//1/ FullInfo TREE content:
//
// https://github.com/sviret/HL_LHC/blob/master/Extractors/RecoExtractor/test/SectorMaker/sector_test.h
//


  int evt;             // Event number (for PU event, where there is more than 1 primary)
  int n_stub_total;    // The total number of stubs in the event
  int n_stub;          // The total number of stubs contained in matched patterns in the event
  int n_stub_t;        // The total number of stubs contained in fitted tracks in the event

  std::vector<float>   *stub_x;      // x coordinates of ALL the stubs
  std::vector<float>   *stub_y;      // y coordinates of ALL the stubs
  std::vector<float>   *stub_z;      // z coordinates of ALL the stubs
  std::vector<int>     *stub_layer;  // layer number of ALL the stubs
  std::vector<int>     *stub_ladder; // ladder number of ALL the stubs
  std::vector<int>     *stub_module; // module number of ALL the stubs
  std::vector<int>     *stub_seg;    // segment number of ALL the stubs
  std::vector<float>   *stub_strip;  // strip number of ALL the stubs
  std::vector<int>     *stub_tp;     // tp index of ALL the stubs (in part_*** vectors of this tree!!!!)
  std::vector<int>     *stub_inpatt; // is the stub in a pattern (1) of not (0)?
  std::vector<int>     *stub_intrk;  // is the stub in a track (1) of not (0)?



  int n_part;                        // The total number of particles inducing at least one stub in the event

  std::vector<int>     *part_pdg;    // PDG id of the particles
  std::vector<int>     *part_nsec;   // In how many trigger towers this particle hit more than 4 different layers/disks?
  std::vector<int>     *part_nhits;  // How many different layers/disks are hit by the particle?
  std::vector<int>     *part_npatt;  // How many patterns contains more than 4 stubs of the particle (in 4 different layers/disks)?
  std::vector<int>     *part_ntrk;   // How many tracks contains more than 4 stubs of the particle (in 4 different layers/disks)?

  std::vector<float>   *part_pt;     // pt of the particles
  std::vector<float>   *part_rho;    // rho0 of the particles
  std::vector<float>   *part_z0;     // z0 of the particles
  std::vector<float>   *part_eta;    // eta of the particles 
  std::vector<float>   *part_phi;    // phi of the particles


  int n_patt;                        // The total number of patterns matched in the event

  // Sector id of all the patterns
  std::vector<int>                  *patt_sec;   

  // tp index of ALL the particles contained in the pattern (in part_*** vectors of this tree!!!!)
  std::vector< std::vector<int> >   *patt_parts; 

  // index of ALL the stubs contained in the pattern (in stub_*** vectors of this tree!!!!) 
  std::vector< std::vector<int> >   *patt_stubs; 


  int n_track;                        // The total number of L1 tracks matched in the event
  std::vector<int>                  *trk_sec;    // Sector id of all the tracks
  std::vector< std::vector<int> >   *trk_parts;  // tp index of ALL the particles contained in the pattern (in part_*** vectors of this tree!!!!)
  std::vector< std::vector<int> >   *trk_stubs;  // index of ALL the stubs contained in the track (in stub_*** vectors of this tree!!!!) 
  std::vector<float>                *trk_pt;     // pt of the track
  std::vector<float>                *trk_eta;    // eta of the track
  std::vector<float>                *trk_z;      // z0 of the track
  std::vector<float>                *trk_phi;    // phi of the track


  std::vector<int>    stub_used;  


  L1TT->SetBranchAddress("evt",          &evt); 

  L1TT->SetBranchAddress("n_stub_total", &n_stub_total); 
  L1TT->SetBranchAddress("n_stub_inpat", &n_stub); 
  L1TT->SetBranchAddress("n_stub_intrk", &n_stub_t); 
  L1TT->SetBranchAddress("stub_x",       &stub_x); 
  L1TT->SetBranchAddress("stub_y",       &stub_y); 
  L1TT->SetBranchAddress("stub_z",       &stub_z); 
  L1TT->SetBranchAddress("stub_layer",   &stub_layer); 
  L1TT->SetBranchAddress("stub_ladder",  &stub_ladder);
  L1TT->SetBranchAddress("stub_module",  &stub_module);
  L1TT->SetBranchAddress("stub_segment", &stub_seg);
  L1TT->SetBranchAddress("stub_strip",   &stub_strip);
  L1TT->SetBranchAddress("stub_tp",      &stub_tp);
  L1TT->SetBranchAddress("stub_inpatt",  &stub_inpatt);
  L1TT->SetBranchAddress("stub_intrack", &stub_intrk);

  L1TT->SetBranchAddress("n_part",       &n_part); 
  L1TT->SetBranchAddress("part_pdg",     &part_pdg); 
  L1TT->SetBranchAddress("part_nsec",    &part_nsec); 
  L1TT->SetBranchAddress("part_nhits",   &part_nhits); 
  L1TT->SetBranchAddress("part_npatt",   &part_npatt); 
  L1TT->SetBranchAddress("part_ntrk",    &part_ntrk); 
  L1TT->SetBranchAddress("part_pt",      &part_pt); 
  L1TT->SetBranchAddress("part_rho",     &part_rho);
  L1TT->SetBranchAddress("part_z0",      &part_z0);
  L1TT->SetBranchAddress("part_eta",     &part_eta);  
  L1TT->SetBranchAddress("part_phi",     &part_phi); 

  L1TT->SetBranchAddress("n_patt",       &n_patt); 
  L1TT->SetBranchAddress("patt_sec",     &patt_sec); 
  L1TT->SetBranchAddress("patt_parts",   &patt_parts); 
  L1TT->SetBranchAddress("patt_stubs",   &patt_stubs); 

  L1TT->SetBranchAddress("n_tc",      &n_track);
  L1TT->SetBranchAddress("tc_sec",      &trk_sec);
  L1TT->SetBranchAddress("tc_parts",    &trk_parts);
  L1TT->SetBranchAddress("tc_stubs",    &trk_stubs);
  L1TT->SetBranchAddress("tc_pt",       &trk_pt);
  L1TT->SetBranchAddress("tc_phi",      &trk_phi);
  L1TT->SetBranchAddress("tc_z",        &trk_z);
  L1TT->SetBranchAddress("tc_eta",      &trk_eta);

  int n_entries = L1TT->GetEntries();
  
  // Loop over patterns

  int idx_s;
  int idx_p;

  float sec_stubs[48];
  float sec_stubs_2[48];
  float sec_patt[48];
  float sec_trk[48];

  float sec_stubs_fake[48];
  float sec_stubs_fake_2[48];
  float sec_patt_fake[48];
  float sec_trk_fake[48];

  float sec_stubs_fake_e[48];
  float sec_stubs_fake_2_e[48];
  float sec_patt_fake_e[48];
  float sec_trk_fake_e[48];


  float eta_stubs[6];
  float eta_stubs_2[6];
  float eta_patt[6];
  float eta_trk[6];

  float eta_stubs_fake[6];
  float eta_stubs_fake_2[6];
  float eta_patt_fake[6];
  float eta_trk_fake[6];

  float eta_stubs_fake_e[6];
  float eta_stubs_fake_2_e[6];
  float eta_patt_fake_e[6];
  float eta_trk_fake_e[6];



  for (int k=0;k<48;++k)
  {
    sec_stubs[k]          = 0.;
    sec_stubs_fake[k]     = 0.;
    sec_stubs_2[k]        = 0.;
    sec_stubs_fake_2[k]   = 0.;
    sec_patt[k]           = 0.;
    sec_patt_fake[k]      = 0.;
    sec_trk[k]            = 0.;
    sec_trk_fake[k]       = 0.;
    sec_stubs_fake_e[k]   = 0.;
    sec_stubs_fake_2_e[k] = 0.;
    sec_patt_fake_e[k]    = 0.;
    sec_trk_fake_e[k]     = 0.;
  }

  for (int k=0;k<6;++k)
  {
    eta_stubs[k]          = 0.;
    eta_stubs_fake[k]     = 0.;
    eta_stubs_2[k]        = 0.;
    eta_stubs_fake_2[k]   = 0.;
    eta_patt[k]           = 0.;
    eta_patt_fake[k]      = 0.;
    eta_trk[k]            = 0.;
    eta_trk_fake[k]       = 0.;
    eta_stubs_fake_e[k]   = 0.;
    eta_stubs_fake_2_e[k] = 0.;
    eta_patt_fake_e[k]    = 0.;
    eta_trk_fake_e[k]     = 0.;
  }




  int ngp;
  int range;


//  for (int j=0;j<n_entries ;++j)
  for (int j=0;j<150 ;++j)
  {
    L1TT->GetEntry(j);

    stub_used.clear();  
    for (int k=0;k<n_stub_total;++k) stub_used.push_back(0); 

    for (int k=0;k<n_patt;++k)
    {
      range = static_cast<int>(patt_sec->at(k)/8); 

      ngp=0;
      sec_patt[patt_sec->at(k)] += 1.;
      eta_patt[range]           += 1.;

      if (patt_parts->at(k).size()==0) 
      {
          sec_patt_fake[patt_sec->at(k)] += 1.;
          eta_patt_fake[range]           += 1.;
      }
      else
      {
          for (int kk=0;kk<patt_parts->at(k).size();++kk)
          {
              idx_p = patt_parts->at(k).at(kk);

              if (part_pt->at(idx_p)>ptmin && part_rho->at(idx_p)<d0max) ++ngp;
          }

          if (ngp==0)
          {
              sec_patt_fake[patt_sec->at(k)] += 1.;
              eta_patt_fake[range]           += 1.;
          }
      }


      for (int kk=0;kk<patt_stubs->at(k).size();++kk)
      {
          idx_s = patt_stubs->at(k).at(kk);
          idx_p = stub_tp->at(idx_s);

          if (stub_used.at(idx_s)!=0) continue;
          stub_used.at(idx_s)=1;

          sec_stubs[patt_sec->at(k)] += 1.;
          eta_stubs[range]           += 1.;

          if (idx_p<0)
          {
              sec_stubs_fake[patt_sec->at(k)] += 1.;
              eta_stubs_fake[range]           += 1.;
          }
          else
          {
              if (part_pt->at(idx_p)<ptmin || part_rho->at(idx_p)>d0max)
              {
                  sec_stubs_fake[patt_sec->at(k)] += 1.;
                  eta_stubs_fake[range]           += 1.;
              }
          }
        }
    }
    
    stub_used.clear();  
    for (int k=0;k<n_stub_total;++k) stub_used.push_back(0); 

    for (int k=0;k<n_track;++k)
    {
      if (trk_pt->at(k)<0.9*ptmin) continue;

      range = static_cast<int>(trk_sec->at(k)/8); 
      ngp=0;
      sec_trk[trk_sec->at(k)] += 1.;
      eta_trk[range]          += 1.;

      if (trk_parts->at(k).size()==0)
      {
	sec_trk_fake[trk_sec->at(k)] += 1.;
	eta_trk_fake[range]          += 1.;
      }
      else
      {
	for (int kk=0;kk<trk_parts->at(k).size();++kk)
	{
	  idx_p = trk_parts->at(k).at(kk);
	  
	  if (part_pt->at(idx_p)>ptmin && part_rho->at(idx_p)<d0max) ++ngp; 
	}

	if (ngp==0)
	{
	  sec_trk_fake[trk_sec->at(k)] += 1.;
	  eta_trk_fake[range]          += 1.;
	}
      }



      for (int kk=0;kk<trk_stubs->at(k).size();++kk)
      {
	idx_s = trk_stubs->at(k).at(kk);
	idx_p = stub_tp->at(idx_s);

	if (stub_used.at(idx_s)!=0) continue;
	stub_used.at(idx_s)=1;	

	sec_stubs_2[trk_sec->at(k)] += 1.;
	eta_stubs_2[range]          += 1.;

	if (idx_p<0)
	{
	  sec_stubs_fake_2[trk_sec->at(k)] += 1.;
	  eta_stubs_fake_2[range]          += 1.;
	}
	else
	{
	  if (part_pt->at(idx_p)<ptmin || part_rho->at(idx_p)>d0max)
	  {
	    sec_stubs_fake_2[trk_sec->at(k)] += 1.;
	    eta_stubs_fake_2[range]          += 1.;
	  }
	}
      }
    }


  }

  float patt_rate[6];
  float trk_rate[6];

  float patt_rate_e[6];
  float trk_rate_e[6];

  for (int k=0;k<48;++k)
  {
    sec_stubs_fake_e[k]   = sqrt((sec_stubs_fake[k]/sec_stubs[k]*(1-sec_stubs_fake[k]/sec_stubs[k]))/sec_stubs[k]);
    sec_stubs_fake_2_e[k] = sqrt((sec_stubs_fake_2[k]/sec_stubs_2[k]*(1-sec_stubs_fake_2[k]/sec_stubs_2[k]))/sec_stubs_2[k]); 
    sec_patt_fake_e[k]    = sqrt((sec_patt_fake[k]/sec_patt[k]*(1-sec_patt_fake[k]/sec_patt[k]))/sec_patt[k]);
    sec_trk_fake_e[k]     = sqrt((sec_trk_fake[k]/sec_trk[k]*(1-sec_trk_fake[k]/sec_trk[k]))/sec_trk[k]);
  }

  for (int k=0;k<6;++k)
  {
    eta_stubs_fake_e[k]   = sqrt((eta_stubs_fake[k]/eta_stubs[k]*(1-eta_stubs_fake[k]/eta_stubs[k]))/eta_stubs[k]);
    eta_stubs_fake_2_e[k] = sqrt((eta_stubs_fake_2[k]/eta_stubs_2[k]*(1-eta_stubs_fake_2[k]/eta_stubs_2[k]))/eta_stubs_2[k]); 
    eta_patt_fake_e[k]    = sqrt((eta_patt_fake[k]/eta_patt[k]*(1-eta_patt_fake[k]/eta_patt[k]))/eta_patt[k]);
    eta_trk_fake_e[k]     = sqrt((eta_trk_fake[k]/eta_trk[k]*(1-eta_trk_fake[k]/eta_trk[k]))/eta_trk[k]);

    eta_stubs_fake[k]     = eta_stubs_fake[k]/eta_stubs[k];
    eta_stubs_fake_2[k]   = eta_stubs_fake_2[k]/eta_stubs_2[k];
    eta_patt_fake[k]      = eta_patt_fake[k]/eta_patt[k];
    eta_trk_fake[k]       = eta_trk_fake[k]/eta_trk[k];


    patt_rate[k]   = eta_patt[k]/(8*n_entries);
    trk_rate[k]    = eta_trk[k]/(8*n_entries);
    patt_rate_e[k] = sqrt(patt_rate[k]);
    trk_rate_e[k]  = sqrt(trk_rate[k]);

    //    cout << k << " / " << patt_rate[k] << " / " << patt_rate_e[k] << endl;
    //    cout << k << " / " << trk_rate[k] << " / " << trk_rate_e[k] << endl;
  }


  const int nbin_eta = 6;
  const int nbin_phi = 8;

  TH2F *eff_s_p    = new TH2F("eff_s_p","eff_s_p",nbin_eta,0,6,nbin_phi,0.,8);
  TH2F *eff_s_t    = new TH2F("eff_s_t","eff_s_t",nbin_eta,0,6,nbin_phi,0.,8);
  TH2F *eff_p      = new TH2F("eff_p","eff_p",nbin_eta,0,6,nbin_phi,0.,8);
  TH2F *eff_t      = new TH2F("eff_t","eff_t",nbin_eta,0,6,nbin_phi,0.,8);

  for (int k=0;k<6;++k)
  {
    for (int l=0;l<8;++l)
    {
      eff_s_p->Fill(k+0.5,l+0.5,sec_stubs_fake[8*k+l]/sec_stubs[8*k+l]);
      eff_s_t->Fill(k+0.5,l+0.5,sec_stubs_fake_2[8*k+l]/sec_stubs_2[8*k+l]);
      eff_p->Fill(k+0.5,l+0.5,sec_patt_fake[8*k+l]/sec_patt[8*k+l]);
      eff_t->Fill(k+0.5,l+0.5,sec_trk_fake[8*k+l]/sec_trk[8*k+l]);
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
  eff_s_p->GetXaxis()->SetTitle("#eta sector number");
  eff_s_p->GetXaxis()->SetLabelFont(42);
  eff_s_p->GetXaxis()->SetLabelSize(0.03);
  eff_s_p->GetXaxis()->SetTitleSize(0.035);
  eff_s_p->GetXaxis()->SetTitleFont(42);
  eff_s_p->GetYaxis()->SetTitle("#phi sector number");
  eff_s_p->GetYaxis()->SetLabelFont(42);
  eff_s_p->GetYaxis()->SetLabelSize(0.03);
  eff_s_p->GetYaxis()->SetTitleSize(0.035);
  eff_s_p->GetYaxis()->SetTitleFont(42);
  eff_s_p->SetContour(99);
  eff_s_p->Draw("colz");
  TPaveText *pte = new TPaveText(1.7,8.4,4.5,9.0,"br");
  pte->SetFillColor(0);
  pte->SetTextFont(42);
  pte->SetTextSize(0.0349944);
  sprintf (buffer, "Fake stub proportion (PR step)");  
  TText *text = pte->AddText(buffer);
  pte->Draw();
  
  c1->cd(2);  
  eff_s_t->GetXaxis()->SetTitle("#eta");
  eff_s_t->GetXaxis()->SetLabelFont(42);
  eff_s_t->GetXaxis()->SetLabelSize(0.03);
  eff_s_t->GetXaxis()->SetTitleSize(0.035);
  eff_s_t->GetXaxis()->SetTitleFont(42);
  eff_s_t->GetYaxis()->SetTitle("#phi (in rad)");
  eff_s_t->GetYaxis()->SetLabelFont(42);
  eff_s_t->GetYaxis()->SetLabelSize(0.03);
  eff_s_t->GetYaxis()->SetTitleSize(0.035);
  eff_s_t->GetYaxis()->SetTitleFont(42);
  eff_s_t->SetContour(99);
  eff_s_t->Draw("colz");

  TPaveText *pte = new TPaveText(1.7,8.4,4.5,9.0,"br");
  pte->SetFillColor(0);
  pte->SetTextFont(42);
  pte->SetTextSize(0.0349944);
  sprintf (buffer, "Fake stub proportion (FIT step)");  
  TText *text = pte->AddText(buffer);
  pte->Draw();
  
  
  c1->cd(3);  
  eff_p->GetXaxis()->SetTitle("#eta");
  eff_p->GetXaxis()->SetLabelFont(42);
  eff_p->GetXaxis()->SetLabelSize(0.03);
  eff_p->GetXaxis()->SetTitleSize(0.035);
  eff_p->GetXaxis()->SetTitleFont(42);
  eff_p->GetYaxis()->SetTitle("#phi (in rad)");
  eff_p->GetYaxis()->SetLabelFont(42);
  eff_p->GetYaxis()->SetLabelSize(0.03);
  eff_p->GetYaxis()->SetTitleSize(0.035);
  eff_p->GetYaxis()->SetTitleFont(42);
  eff_p->SetContour(99);
  eff_p->Draw("colz");

  TPaveText *pte = new TPaveText(1.7,8.4,4.5,9.0,"br");
  pte->SetFillColor(0);
  pte->SetTextFont(42);
  pte->SetTextSize(0.0349944);
  sprintf (buffer, "Fake pattern proportion");  
  TText *text = pte->AddText(buffer);
  pte->Draw();

  
  c1->cd(4);  
  eff_t->GetXaxis()->SetTitle("#eta");
  eff_t->GetXaxis()->SetLabelFont(42);
  eff_t->GetXaxis()->SetLabelSize(0.03);
  eff_t->GetXaxis()->SetTitleSize(0.035);
  eff_t->GetXaxis()->SetTitleFont(42);
  eff_t->GetYaxis()->SetTitle("#phi (in rad)");
  eff_t->GetYaxis()->SetLabelFont(42);
  eff_t->GetYaxis()->SetLabelSize(0.03);
  eff_t->GetYaxis()->SetTitleSize(0.035);
  eff_t->GetYaxis()->SetTitleFont(42);
  eff_t->SetContour(99);
  eff_t->Draw("colz");

  TPaveText *pte = new TPaveText(1.7,8.4,4.5,9.0,"br");
  pte->SetFillColor(0);
  pte->SetTextFont(42);
  pte->SetTextSize(0.0349944);
  sprintf (buffer, "Fake track proportion");  
  TText *text = pte->AddText(buffer);
  pte->Draw();
  

  c1->Modified();
  c1->Update();
  /*

  sprintf (buffer, "L1_2D_summary_%d_hits_track.png",nh);

  c1->Print(buffer);
 
  sprintf (buffer, "L1_2D_summary_%d_hits_track.C",nh);

  c1->SaveSource(buffer);
  */

  float eta[6];
  float eta_e[6];



  for (int k=0;k<6;++k)
  {
    eta[k]                = k+0.5;
    eta_e[k]              = 0.;
  }

  TH2F *cadre    = new TH2F("cadre","cadre",6,0,6,200,0.,1.);


  c2 = new TCanvas("c2","Fake stub proportions",200,10,700,500);
  c2->SetGrid();
  c2->GetFrame()->SetBorderSize(12);

  st_patt = new TGraphErrors(6,eta,eta_stubs_fake,eta_e,eta_stubs_fake_e);
  st_patt->SetMarkerStyle(20);

  st_trk = new TGraphErrors(6,eta,eta_stubs_fake_2,eta_e,eta_stubs_fake_2_e);
  st_trk->SetMarkerStyle(24);

  cadre->Draw();
  st_patt->Draw("LPsame");
  st_trk->Draw("LPsame");

  leg = new TLegend(0.7,0.85,1.,1.);
  leg->AddEntry(st_patt,"Fake stub proportion after PR","p");
  leg->AddEntry(st_trk,"Fake stub proportion after FIT","p");
  leg->Draw();
  
  c2->Modified();
  c2->Update();




  c3 = new TCanvas("c3","Fake proportions",200,10,700,500);
  c3->SetGrid();
  c3->GetFrame()->SetBorderSize(12);

  f_patt = new TGraphErrors(6,eta,eta_patt_fake,eta_e,eta_patt_fake_e);
  f_patt->SetMarkerStyle(20);

  f_trk = new TGraphErrors(6,eta,eta_trk_fake,eta_e,eta_trk_fake_e);
  f_trk->SetMarkerStyle(24);

  cadre->Draw();
  f_patt->Draw("LPsame");
  f_trk->Draw("LPsame");

  leg = new TLegend(0.7,0.85,1.,1.);
  leg->AddEntry(f_patt,"Fake pattern proportion","p");
  leg->AddEntry(f_trk,"Fake track proportion","p");
  leg->Draw();
  
  c3->Modified();
  c3->Update();


  TH2F *cadre2    = new TH2F("cadre2","cadre2",6,0,6,200,0.,1000.);
  TH2F *cadre3    = new TH2F("cadre3","cadre3",6,0,6,200,0.,20.);


  c4 = new TCanvas("c4","Rates",200,10,700,500);
  c4->SetGrid();
  c4->GetFrame()->SetBorderSize(12);
  c4->SetLogy();

  r_patt = new TGraphErrors(6,eta,patt_rate,eta_e,patt_rate_e);
  r_patt->SetMarkerStyle(20);

  cadre2->Draw();
  cadre2->GetXaxis()->SetLabelFont(42);
  cadre2->GetXaxis()->SetLabelSize(0.025);
  cadre2->GetXaxis()->SetTitleSize(0.03);
  cadre2->GetXaxis()->SetTitleFont(42);
  cadre2->GetYaxis()->SetLabelFont(42);
  cadre2->GetYaxis()->SetLabelOffset(0.004);
  cadre2->GetYaxis()->SetLabelSize(0.025);
  cadre2->GetYaxis()->SetTitleSize(0.03);
  cadre2->GetYaxis()->SetTitle("Average pattern rate per sector (in patterns/sec/BX)"); 
  cadre2->GetXaxis()->SetTitle("Sector #eta range"); 
  r_patt->Draw("LPsame");
  
  c4->Modified();
  c4->Update();

  c5 = new TCanvas("c5","Rates",200,10,700,500);
  c5->SetGrid();
  c5->GetFrame()->SetBorderSize(12);
  c5->SetLogy();

  r_trk = new TGraphErrors(6,eta,trk_rate,eta_e,trk_rate_e);
  r_trk->SetMarkerStyle(24);

  cadre3->Draw();
  cadre3->GetXaxis()->SetLabelFont(42);
  cadre3->GetXaxis()->SetLabelSize(0.025);
  cadre3->GetXaxis()->SetTitleSize(0.03);
  cadre3->GetXaxis()->SetTitleFont(42);
  cadre3->GetYaxis()->SetLabelFont(42);
  cadre3->GetYaxis()->SetLabelOffset(0.004);
  cadre3->GetYaxis()->SetLabelSize(0.025);
  cadre3->GetYaxis()->SetTitleSize(0.03);
  cadre3->GetYaxis()->SetTitle("Average track rate per sector (in tracks/sec/BX)"); 
  cadre3->GetXaxis()->SetTitle("Sector #eta range"); 
  r_patt->Draw("LPsame");
  r_trk->Draw("LPsame");
  
  c5->Modified();
  c5->Update();
}
