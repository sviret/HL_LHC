// Class for the stub windows determination
// For more info, look at the header file

#include "windows.h"

// Main constructor

windows::windows(std::string file_r, std::string file_e, float pminm, float pmaxm, float pmine, float pmaxe, float prop, float CIC_lim)
{
  m_pmin=pminm;
  m_pmax=pmaxm;
  m_prop=prop;
  m_lim=CIC_lim;
    
  windows::initTuple(file_r,file_e); // Open the input files (PU0 PGUN and PUXXX data sequence
  windows::reset();    
  windows::initVars(); 
  windows::get_rates();              // Determine the ave. stub rate in function of stub window, for the 108 regions
  windows::get_losses();             // Determine the good stub losses in function of stub window, for the 108 regions
    
  windows::get_effs(0);              // Determine the muon stub efficiency in function of stub window, for the 108 regions

  m_pmin=pmine;
  m_pmax=pmaxe;
  windows::get_effs(1);              // Determine the electron stub efficiency in function of stub window, for the 108 regions
  windows::get_result(1,0);          // Determine the muon-based tuning
  windows::get_result(1,1);          // Determine the electron-based tuning
  windows::print_result(0);          // Print the TIGHT SW tuning
  windows::print_result(1);          // Print the LOOSE SW tuning
    
    m_ratetree->Fill();
    m_outfile->Write();
    delete L1TT;
    
    delete m_outfile;}


//////////////////////////////////////////////
//
// 
//////////////////////////////////////////////
void windows::get_rates()
{
  // Initialize some params
  
  double fact   = 1./static_cast<float>(8*n_r); // Normalization factor

  double n_mod_tilt[6]      = {31,35,39,24,24,24};
  double n_lad_tilt[6]      = {18,26,36,48,60,78};
  double n_mod_endcap[15]   = {20,24,24,28,32,32,36,40,40,44,52,60,64,72,76};

  // Then loop over events
  
  int rank;

  for (int j=0;j<n_r;++j)
  {
    Losses->GetEntry(j);
    
    //if (j%1000==0) cout << j << endl;

    for (int k=0;k<6;++k)
    {
      for (int l=0;l<40;++l)
      {
          rank=0;
          if (k<=2)
          {
              if (l<12) rank=12-l;
              if (l>=n_mod_tilt[k]-12) rank=l+13-n_mod_tilt[k];
          }

          for (int m=0;m<16;++m) // Loop over all the SW indexes (max is 8*2)
          {
              if (k>2)
              {
                  barrel_w[k][rank][m][0] += (fact/(n_mod_tilt[k]*n_lad_tilt[k]))*m_rate_b[k][l][m][2][0];
              }
              else
              {
                  if (rank==0)
                  {
                      barrel_w[k][rank][m][0] += (fact/((n_mod_tilt[k]-24)*n_lad_tilt[k]))*m_rate_b[k][l][m][2][0];
                  }
                  else
                  {
                      barrel_w[k][rank][m][0] += fact/(2*n_lad_tilt[k])*m_rate_b[k][l][m][2][0];
                  }
              }
          }
      }
    }

    for (int k=0;k<5;++k)
    {
      for (int l=0;l<15;++l)
      {
          for (int m=0;m<16;++m) // Loop over all the particles
          {
              if (k<2)
              {
                  disk_w[k][l][m][0] += fact/(2*n_mod_endcap[l])*m_rate_d[k][l][m][2][0];
              }
              else
              {
                  if (l>=12) continue;
                  disk_w[k][l][m][0] += fact/(2*n_mod_endcap[l+3])*m_rate_d[k][l][m][2][0];
              }
          }
        }
    }
  } // End of loop over events

    // Finally we normalize the rates
    

    for (int i=0;i<16;++i)
    {
        for (int j=0;j<6;++j)
        {
            for (int l=0;l<13;++l) barrel_w[j][l][i][0] /= barrel_w[j][l][15][0];
        }
        
        for (int j=0;j<5;++j)
        {
            for (int l=0;l<15;++l) disk_w[j][l][i][0] /= disk_w[j][l][15][0];
        }
    }

    
    
  //  windows::print_result();
}

void windows::get_losses()
{
  // Initialize some params
  
  double fact   = 100./static_cast<float>(n_r); // Normalization factor

  double n_mod_tilt[6]      = {31,35,39,24,24,24};
  // Then loop over events
  
  int rank;

  for (int j=0;j<n_r;++j)
  {
    Losses->GetEntry(j);

    for (int k=0;k<6;++k)
    {
      for (int l=0;l<40;++l)
      {
          rank=0;
          if (k<=2)
          {
              if (l<12) rank=12-l;
              if (l>=n_mod_tilt[k]-12) rank=l+13-n_mod_tilt[k];
          }

          for (int m=0;m<16;++m) // Loop over all the particles
          {
              if (k>2)
              {
                  barrel_w[k][0][m][1] += (fact/n_mod_tilt[k])*m_ovflow_b[k][l][m][2][1]/m_ovflow_b[k][l][m][0][1];
                  barrel_w[k][0][m][3] += (fact/n_mod_tilt[k])*m_loss_b[k][l][m][2][1];
              }
              else
              {
                  if (rank==0)
                  {
                      barrel_w[k][rank][m][1] += (fact/(n_mod_tilt[k]-24))*m_ovflow_b[k][l][m][2][1]/m_ovflow_b[k][l][m][0][1];
                      barrel_w[k][rank][m][3] += (fact/(n_mod_tilt[k]-24))*m_loss_b[k][l][m][2][1];
                  }
                  else
                  {
                      barrel_w[k][rank][m][1] += fact*m_ovflow_b[k][l][m][2][1]/m_ovflow_b[k][l][m][0][1];
                      barrel_w[k][rank][m][3] += fact*m_loss_b[k][l][m][2][1];

                  }
              }
          }
      }
    }

    for (int k=0;k<5;++k)
    {
        for (int l=0;l<15;++l)
        {
            for (int m=0;m<16;++m) // Loop over all the particles
            {
                if (m_ovflow_d[k][l][m][0][1]==0) continue;

                disk_w[k][l][m][1] += fact*m_ovflow_d[k][l][m][2][1]/m_ovflow_d[k][l][m][0][1];
                disk_w[k][l][m][3] += fact*m_loss_d[k][l][m][2][1];
            }
        }
    }
  } // End of loop over events

  delete Losses;
}

void windows::get_effs(int ptype)
{
  // Initialize some params
  
  // Then loop over events
  
  int stub_per_lay[10][20];

    (ptype==0)
    ? m_ptype = 13
    : m_ptype = 11;
    
  for (int j=0;j<n_e;++j)
  {
    L1TT->GetEntry(j);
     
    for (int j=0;j<10;++j)
    { 
      for (int k=0;k<20;++k)
      {
	stub_per_lay[j][k]=0;
      }
    } 
    // Loop over stubs
    
    if (m_stub==0) continue;

    for (int k=0;k<m_stub;++k)
    {
      if (std::abs(m_stub_pdg->at(k))!=m_ptype) continue;
      if (m_stub_tp->at(k)>=10 || m_stub_tp->at(k)<0) continue;

      float pt_GEN = sqrt(m_stub_pxGEN->at(k)*m_stub_pxGEN->at(k)+m_stub_pyGEN->at(k)*m_stub_pyGEN->at(k));            
      if (pt_GEN<m_pmin || pt_GEN>m_pmax) continue;

      int modid = m_stub_layer->at(k)-5;

      ++stub_per_lay[m_stub_tp->at(k)][modid];
    }
    
    int ladder=0,layer=0;
    bool isba;

    for (int k=0;k<m_stub;++k)
    {
      if (std::abs(m_stub_pdg->at(k))!=m_ptype) continue;
      if (m_stub_tp->at(k)>=10 || m_stub_tp->at(k)<0) continue;

      layer  = m_stub_layer->at(k)-5;
      isba=false;

      // Here we add a cut to remove multi-cluster stubs. Electron-tuning is sufficient to account for them, and they would bias
      // the tight tuning definition for muons

      if (stub_per_lay[m_stub_tp->at(k)][layer]>=3) continue;
          
      float pt_GEN = sqrt(m_stub_pxGEN->at(k)*m_stub_pxGEN->at(k)+m_stub_pyGEN->at(k)*m_stub_pyGEN->at(k));            
      if (pt_GEN<m_pmin || pt_GEN>m_pmax) continue;

      int idx   = static_cast<int>(std::abs(m_stub_deltas->at(k)*2));
      if (idx>=16) idx=15;
      ladder = 0;
 
      if (m_stub_type->at(k)==1) // Tilted -
      {
          ladder = 12-m_stub_module->at(k);
      }
            
      if (m_stub_type->at(k)==2) // Tilted +
      {
          if (layer==0) ladder = m_stub_module->at(k)-18;
          if (layer==1) ladder = m_stub_module->at(k)-22;
          if (layer==2) ladder = m_stub_module->at(k)-26;
      }

      if (layer<=5)
      {
          isba=true;
          for (int l=0;l<idx+1;++l) ++barrel_w[layer][ladder][l][2+2*ptype];
      }
            
      if (isba) continue;

      // Now look at disks

      layer  = (layer-6)%7;
      ladder = m_stub_ladder->at(k);

      for (int l=0;l<idx+1;++l) ++disk_w[layer][ladder][l][2+2*ptype];
                                    
    } // End for (int k=0;k<m_stub;++k)
  }   // End of loop over events

  for (int j=0;j<6;++j)
  {
    for (int l=0;l<13;++l)
    {
      int ntot = barrel_w[j][l][0][2+2*ptype];

      if (ntot==0) continue;

      for (int i=0;i<16;++i)
      {
          barrel_w[j][l][i][2+2*ptype] = (ntot-barrel_w[j][l][i][2+2*ptype])/ntot;
      }
    }
  }      

  for (int j=0;j<5;++j)
  {
    for (int l=0;l<15;++l) 
    {
      int ntot = disk_w[j][l][0][2+2*ptype];

      if (ntot==0) continue;

      for (int i=0;i<16;++i)
      {
          disk_w[j][l][i][2+2*ptype] = (ntot-disk_w[j][l][i][2+2*ptype])/ntot;
      }
    }
  }  
}


void windows::initVars()
{
  for (int k=0;k<5;++k)
  {
    for (int i=0;i<16;++i)
    {
      for (int j=0;j<6;++j)
      {
          for (int l=0;l<13;++l) barrel_w[j][l][i][k] = 0;
      }
      
      for (int j=0;j<5;++j)
      {
          for (int l=0;l<15;++l) disk_w[j][l][i][k] = 0;
      }
    }
  }

    for (int k=0;k<2;++k)
    {
        for (int j=0;j<6;++j)
        {
            for (int l=0;l<13;++l) barrel_tune[j][l][k] = 0;
        }
            
        for (int j=0;j<5;++j)
        {
            for (int l=0;l<15;++l) disk_tune[j][l][k] = 0;
        }
    }
    

    for (int i=0;i<16;++i)
    {
      for (int j=0;j<2;++j)
      {
	for (int k=0;k<6;++k)
	{
	  for (int l=0;l<13;++l)
	  {
	    loss_b[k][l][i][j] = 0.;
	  }
	}

	for (int k=0;k<5;++k)
	{
	  for (int l=0;l<15;++l)
	  {
	    loss_d[k][l][i][j] = 0.;
	  }
	}
      }
    }
}

void windows::get_result(int ltype, int ptype)
{
    float cut,loss,eff,raeff,raeff_p,rate;
    int i,it,it1,it2;
    
    for(it=0;it<6;it++)
    {
        for(i=0 ; i< 16 ; i++)
        {
            loss = barrel_w[it][0][i][ltype];
            eff  = barrel_w[it][0][i][2+2*ptype];
            rate = barrel_w[it][0][i][0];
            raeff= eff*sqrt(1-rate);
            raeff_p = raeff;
            
            if (i>0) raeff_p = barrel_w[it][0][i-1][2+2*ptype]*sqrt(1-barrel_w[it][0][i-1][0]);
            
            if (loss>m_lim)
            {
                i--;
                break;
            }
            
            if (eff>=m_prop && raeff<raeff_p)
            {
                if (barrel_w[it][0][i-1][2+2*ptype]>=m_prop)
                {
                    i--;
                    break;
                }
            }
            
        }
        
        cut=float(i)/2.;
        barrel_tune[it][0][ptype] = cut;
    }
    
    // Special case barrel tilted
 
    
    for(it1=0;it1<3;it1++)
    {
        for(it2=1;it2<13;it2++)
        {
            for(i=0 ; i< 16 ; i++)
            {
                loss = barrel_w[it1][it2][i][ltype];
                eff  = barrel_w[it1][it2][i][2+2*ptype];
                rate = barrel_w[it1][it2][i][0];
                raeff= eff*sqrt(1-rate);
                raeff_p = raeff;
                
                if (i>0) raeff_p = barrel_w[it1][it2][i-1][2+2*ptype]*sqrt(1-barrel_w[it1][it2][i-1][0]);
                
                if (loss>m_lim)
                {
                    i--;
                    break;
                }
                
                if (eff>=m_prop && raeff<raeff_p)
                {
                    if (barrel_w[it1][it2][i-1][2+2*ptype]>=m_prop)
                    {
                        i--;
                        break;
                    }
                }
            }
            
            cut=float(i)/2;
            barrel_tune[it1][it2][ptype] = cut;
        }
    }

    // Then finally the endcap cuts
    
    for(it1=0;it1<5;it1++)
    {
        for(it2=0;it2<15;it2++)
        {
            if (it1>1 && it2>11) continue;
            
            for(i=0 ; i< 16 ; i++)
            {
                loss = disk_w[it1][it2][i][ltype];
                eff  = disk_w[it1][it2][i][2+2*ptype];
                rate = disk_w[it1][it2][i][0];
                raeff= eff*sqrt(1-rate);
                raeff_p = raeff;
                
                if (i>0) raeff_p = disk_w[it1][it2][i-1][2+2*ptype]*sqrt(1-disk_w[it1][it2][i-1][0]);
                
                if (loss>m_lim)
                {
                    i--;
                    break;
                }
                
                if (eff>=m_prop && raeff<raeff_p)
                {
                    if (disk_w[it1][it2][i-1][2+2*ptype]>=m_prop)
                    {
                        i--;
                        break;
                    }
                }
            }
            
            cut=float(i)/2;
            disk_tune[it1][it2][ptype] = cut;
        }
    }//end loop on endcaps
}

void windows::print_result(int type)
{
    cout << endl;

    
    (type==0) // Muon-based tuning
    ? cout << "# --> Tight SW tuning (Muon-based)"<< endl
    : cout << "# --> Loose SW tuning (Muon/Electron-based)"<< endl;

    cout << "# --> Good stub losses should be below " << m_lim << "%"<< endl;
      
    cout << "process.TTStubAlgorithm_official_Phase2TrackerDigi_.BarrelCut = cms.vdouble( 0, " ;

    float cut_mu,cut_ele,cut;
    int it,it1,it2;
  
    for(it=0;it<6;it++)
    {
        cut_mu  = barrel_tune[it][0][0];
        cut_ele = barrel_tune[it][0][1];
      
        (it<1)
        ? cut = cut_mu
        : cut = std::max(cut_mu,type*cut_ele);
      
        if (cut>7) cut=7;
      
        cout << cut;
		
        if (it<5)
        {
            cout << ", ";
        }else{
            cout << ") " << endl;
        }
    }

    // Special case barrel tilted
        
    cout << "process.TTStubAlgorithm_official_Phase2TrackerDigi_.TiltedBarrelCutSet = cms.VPSet(" << endl;
    cout << "cms.PSet( TiltedCut = cms.vdouble( 0 ) ),"<< endl;
        
    for(it1=0;it1<3;it1++)
    {
        cout << "cms.PSet( TiltedCut = cms.vdouble( 0, ";
        
        for(it2=1;it2<13;it2++)
        {
            cut_mu  = barrel_tune[it1][it2][0];
            cut_ele = barrel_tune[it1][it2][1];
        
            (it1<1)
            ? cut = cut_mu
            : cut = std::max(cut_mu,type*cut_ele);
        
            if (cut>7) cut=7;
            cout << cut;
        
            if (it2<12)
            {
                cout << ", ";
            }else{
                cout << ") )," << endl;
            }
        }
    }
    cout <<")" << endl;

    // Then finally the endcap cuts
  
    cout << "process.TTStubAlgorithm_official_Phase2TrackerDigi_.EndcapCutSet = cms.VPSet(" << endl;
    cout << "cms.PSet( EndcapCut = cms.vdouble( 0 ) ),"<< endl;

    for(it1=0;it1<5;it1++)
    {
        cout << "cms.PSet( EndcapCut = cms.vdouble( 0, ";
    
        for(it2=0;it2<15;it2++)
        {
            if (it1>1 && it2>11) continue;
        
            cut_mu  = disk_tune[it1][it2][0];
            cut_ele = disk_tune[it1][it2][1];
        
            (it1<2 && it2<3)
            ? cut = cut_mu
            : cut = std::max(cut_mu,type*cut_ele);
        
            if (cut>7) cut=7;
            cout << cut;
                            
            if ((it2<14 && it1<=1) || (it2<11 && it1>1))
            {
                cout << ", ";
            }else{
                cout << ") )," << endl;
            }
        }
    }//end loop on endcaps
    cout <<")" << endl;

    cout << endl; // Then we print the final eff/losses couples
}


void windows::reset()
{
  m_stub=0;
  m_stub_layer->clear();
  m_stub_ladder->clear();
  m_stub_tp->clear();
  m_stub_pxGEN->clear();
  m_stub_pyGEN->clear();
  m_stub_etaGEN->clear();
  m_stub_module->clear();
  m_stub_type->clear();
  m_stub_deltas->clear();
  m_stub_pdg->clear();
}


void windows::initTuple(std::string in_r,std::string in_e)
{

  L1TT   = new TChain("TkStubs");
  Losses = new TChain("Trigger_SUM");

  L1TT->Add(in_e.c_str());
  Losses->Add(in_r.c_str());
    
  n_e = L1TT->GetEntries();
  n_r = Losses->GetEntries();  
  
  m_stub_layer   = new  std::vector<int>;
  m_stub_ladder  = new  std::vector<int>;
  m_stub_module  = new  std::vector<int>;
  m_stub_type    = new  std::vector<int>;
  m_stub_tp      = new  std::vector<int>;
  m_stub_pdg     = new  std::vector<int>;
  m_stub_deltas  = new  std::vector<float>;
  m_stub_pxGEN   = new  std::vector<float>;    
  m_stub_pyGEN   = new  std::vector<float>;    
  m_stub_etaGEN  = new  std::vector<float>;    
  
  //  L1TT->SetBranchStatus("*",0);
  L1TT->SetBranchAddress("L1TkSTUB_n",       &m_stub);
  L1TT->SetBranchAddress("L1TkSTUB_layer",   &m_stub_layer);
  L1TT->SetBranchAddress("L1TkSTUB_module",  &m_stub_module);
  L1TT->SetBranchAddress("L1TkSTUB_ladder",  &m_stub_ladder);
  L1TT->SetBranchAddress("L1TkSTUB_type",    &m_stub_type);
  L1TT->SetBranchAddress("L1TkSTUB_tp",      &m_stub_tp);
  L1TT->SetBranchAddress("L1TkSTUB_deltas",  &m_stub_deltas);
  L1TT->SetBranchAddress("L1TkSTUB_pxGEN",   &m_stub_pxGEN);
  L1TT->SetBranchAddress("L1TkSTUB_pyGEN",   &m_stub_pyGEN);
  L1TT->SetBranchAddress("L1TkSTUB_etaGEN",  &m_stub_etaGEN);
  L1TT->SetBranchAddress("L1TkSTUB_pdgID",   &m_stub_pdg);
  
  //  Losses->SetBranchStatus("*",0);
  Losses->SetBranchAddress("BAR_OVFLOW_I",       &m_ovflow_b);
  Losses->SetBranchAddress("DIS_OVFLOW_I",       &m_ovflow_d);
  Losses->SetBranchAddress("BAR_MULT",       &m_rate_b);
  Losses->SetBranchAddress("DIS_MULT",       &m_rate_d);
  Losses->SetBranchAddress("BAR_LOSS",       &m_loss_b);
  Losses->SetBranchAddress("DIS_LOSS",       &m_loss_d);
    
  m_outfile  = new TFile("stub_windows.root","recreate");
  m_ratetree = new TTree("Windows","Stub Windows info");
  
  m_ratetree->Branch("barrel_summary",  &barrel_w,  "barrel_w[6][13][16][5]/F");
  m_ratetree->Branch("disk_summary",    &disk_w,    "disk_w[5][15][16][5]/F");
}
