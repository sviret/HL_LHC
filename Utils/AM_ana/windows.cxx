// Class for the stub windows determination
// For more info, look at the header file

#include "windows.h"

// Main constructor

windows::windows(std::string filename, float pmin, float pmax, float prop, int ptype)
{
  m_pmin=pmin;
  m_pmax=pmax;
  m_prop=prop;
  m_ptype=ptype;

  windows::initTuple(filename);
  windows::initVars();
  windows::get_windows();
}


//////////////////////////////////////////////
//
// 
//////////////////////////////////////////////

void windows::get_windows()
{
  // Initialize some params
  
  int n_entries = L1TT->GetEntries();

  // Then loop over events
  
  int stub_per_lay_0[20];
  int stub_per_lay_1[20];

  for (int j=0;j<n_entries;++j)
  {
    if (j%10000==0) cout << j <<endl;
    
    L1TT->GetEntry(j,1);

    for (int k=0;k<20;++k)
    {
      stub_per_lay_0[k]=0;
      stub_per_lay_1[k]=0;
    }
        
    // Loop over stubs
        
    for (int k=0;k<m_stub;++k)
    {
      if (std::abs(m_stub_pdg.at(k))!=m_ptype) continue;

      if (m_stub_tp.at(k)!=0 && m_stub_tp.at(k)!=1) continue;
      if (m_stub_tp.at(k)==0) ++stub_per_lay_0[m_stub_layer.at(k)-5];
      if (m_stub_tp.at(k)==1) ++stub_per_lay_1[m_stub_layer.at(k)-5];
    }

    int ladder=0,layer=0;
    bool isba;

    for (int k=0;k<m_stub;++k)
    {
      if (std::abs(m_stub_pdg.at(k))!=m_ptype) continue;

      layer  = m_stub_layer.at(k)-5;
      isba=false;
      if (m_stub_tp.at(k)!=0 && m_stub_tp.at(k)!=1) continue;
      if (m_stub_tp.at(k)==0 && stub_per_lay_0[layer]!=1) continue;
      if (m_stub_tp.at(k)==1 && stub_per_lay_1[layer]!=1) continue;
          
      float pt_GEN = sqrt(m_stub_pxGEN.at(k)*m_stub_pxGEN.at(k)+m_stub_pyGEN.at(k)*m_stub_pyGEN.at(k));            
      if (pt_GEN<m_pmin || pt_GEN>m_pmax) continue;

      int idx   = static_cast<int>(std::abs(m_stub_deltas.at(k)*2));
      if (idx>=40) idx=39;
      ladder = 0;
 
      if (m_stub_type.at(k)==1) // Tilted -
      {
	ladder = 12-m_stub_module.at(k);
      }
            
      if (m_stub_type.at(k)==2) // Tilted +
      {
	if (layer==0) ladder = m_stub_module.at(k)-19;
	if (layer==1) ladder = m_stub_module.at(k)-23;
	if (layer==2) ladder = m_stub_module.at(k)-27;
      }

      if (layer<=5)
      {
	isba=true;
	for (int l=0;l<idx+1;++l) ++barrel_w[layer][ladder][l];	   
      }
            
      if (isba) continue;

      // Now look at disks

      layer  = (layer-6)%7;
      ladder = m_stub_ladder.at(k);

      for (int l=0;l<idx+1;++l) ++disk_w[layer][ladder][l];
                                    
    } // End for (int k=0;k<m_stub;++k)
  }   // End of loop over events
  
  windows::print_result();
}


void windows::initVars()
{
  for (int i=0;i<6;++i)
  {
    for (int j=0;j<13;++j)
    {
      for (int k=0;k<40;++k) barrel_w[i][j][k] = 0;
    }
  }

  for (int i=0;i<5;++i)
  {
    for (int j=0;j<15;++j)
    {
      for (int k=0;k<40;++k) disk_w[i][j][k] = 0;
    }
  }
}

void windows::print_result()
{
  cout << "# "<< m_pmin << "GeV/c Threshold" << endl;
  cout << "# --> Keep " << static_cast<int>(100*m_prop) << "% of the stubs between " << m_pmin << " and " << m_pmax << " GeV/c" << endl;
  cout << "process.TTStubAlgorithm_official_Phase2TrackerDigi_.BarrelCut = cms.vdouble( 0, " ;
    
  float cut;
  int i,it,it1,it2;
  
  for(it=0;it<6;it++)
  {
    for(i=0 ; i< 40 ; i++)
    {
      int ntot = barrel_w[it][0][0];
      int nup  = ntot-barrel_w[it][0][i];

      if(nup >= m_prop*ntot) break; // Limit is passed
    }

    cut=float(i)/2.;                
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
      for(i=0 ; i< 40 ; i++)
      {
	int ntot = barrel_w[it1][it2][0];
	int nup  = ntot-barrel_w[it1][it2][i];

	if(nup >= m_prop*ntot) break; // Limit is passed
      }
	    
      cut=float(i)/2;
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
      for(i=0 ; i< 40 ; i++)
      {
	int ntot = disk_w[it1][it2][0];
	int nup  = ntot-disk_w[it1][it2][i];
            
	if(nup >= m_prop*ntot) break; // Limit is passed
      }

      cut=float(i)/2;
      cout << cut;
                            
      if (it2<14)
      {
	cout << ", ";
      }else{
	cout << ") )," << endl;
      }
    }
  }//end loop on endcaps
  cout <<")" << endl;    
}


void windows::reset()
{
}


void windows::initTuple(std::string in)
{

    L1TT  = new TChain("TkStubs");

    L1TT->Add(in.c_str());

    pm_stub_layer=&m_stub_layer;
    pm_stub_ladder=&m_stub_ladder;
    pm_stub_module=&m_stub_module;
    pm_stub_deltas=&m_stub_deltas;
    pm_stub_pxGEN=&m_stub_pxGEN;
    pm_stub_pyGEN=&m_stub_pyGEN;
    pm_stub_etaGEN=&m_stub_etaGEN;
    pm_stub_tp=&m_stub_tp;
    pm_stub_pdg=&m_stub_pdg;
    pm_stub_type=&m_stub_type;

    L1TT->SetBranchStatus("*",0);

    L1TT->SetBranchAddress("L1TkSTUB_n",       &m_stub); 
    L1TT->SetBranchAddress("L1TkSTUB_tp",      &pm_stub_tp);
    L1TT->SetBranchAddress("L1TkSTUB_pxGEN",   &pm_stub_pxGEN);
    L1TT->SetBranchAddress("L1TkSTUB_pyGEN",   &pm_stub_pyGEN);
    L1TT->SetBranchAddress("L1TkSTUB_etaGEN",  &pm_stub_etaGEN);
    L1TT->SetBranchAddress("L1TkSTUB_layer",   &pm_stub_layer);
    L1TT->SetBranchAddress("L1TkSTUB_module",  &pm_stub_module);
    L1TT->SetBranchAddress("L1TkSTUB_ladder",  &pm_stub_ladder);
    L1TT->SetBranchAddress("L1TkSTUB_type",    &pm_stub_type);
    L1TT->SetBranchAddress("L1TkSTUB_deltas",  &pm_stub_deltas);
    L1TT->SetBranchAddress("L1TkSTUB_pdgID",   &pm_stub_pdg);
}
