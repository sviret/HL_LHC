// Class for the efficiencies determination
// For more info, look at the header file

#include "efficiencies.h"

// Main constructor

efficiencies::efficiencies(std::string filename, std::string outfile, int ptype)
{
    m_type = ptype;

    efficiencies::initTuple(filename,outfile);
    efficiencies::reset();
    efficiencies::initVars();
    efficiencies::get_efficiencies();
}


//////////////////////////////////////////////
//
// 
//////////////////////////////////////////////

void efficiencies::get_efficiencies()
{
    // Initialize some params
 
    double PTGEN=0;
    double D0GEN=0;

    int pix_i;
    int clus_i;
    int stub_i;
  
    int i_lay;
    int i_mod;
    int i_seg;
    float i_eta;

    int evtID;
    std::vector<int> stID;
    std::vector<int> pix_evtID;
    std::vector<int> pix_stID;

    int hit_on_lay[30];
    int clus_on_lay[30];
    int stub_on_lay[30];

    int n_entries = L1TT->GetEntries();
    bool inTP;

    // Then loop over events

    bool m_dbg =false;
  
    std::vector<int> pix_tp;
    std::vector<int> pix_cl;
    std::vector<int> pix_st;



    for (int j=0;j<n_entries;++j)
    //for (int j=0;j<1000;++j)
    {
        efficiencies::reset();

        PTGEN=0;
        D0GEN=0;

        pix_i=-1;
        clus_i=-1;
        stub_i=-1;

        if (j%10000==0) cout << j <<endl;

        // Principle is the following:
        // We loop over all digis,
        // and look if those digis have induced CLUSTER/STUBS
        //
        // Efficiencies are then simply N_object/N_digis

        MC->GetEntry(j,1); // Get the MC info
        Pix->GetEntry(j,1);
        L1TT->GetEntry(j,1);

        // This code is for cluster/stub efficiency calculation
        //
	
	if (m_pclus==0) continue; // No digi, skip
        
	pix_tp.clear();
	pix_cl.clear();
	pix_st.clear();

        for (int k=0;k<m_pclus;++k) // Loop over TPs
        {
	  pix_tp.push_back(-1);
	  pix_cl.push_back(-1);
	  pix_st.push_back(-1);
	}

        //for (int k=0;k<1;++k) // Loop over TPs
        for (int k=0;k<m_part_n;++k) // Loop over TPs
        {
	  // cout << m_part_pdg->at(k) << endl;

	  if (abs(m_part_pdg->at(k)) != m_type) continue;

            // Thie TP corresponds to what we are looking for
            // we continue

            
            //D0GEN = fabs(m_part_y->at(k)*cos(atan2(m_part_py->at(k),m_part_px->at(k)))-m_part_x->at(k)*sin(atan2(m_part_py->at(k),m_part_px->at(k))));
            D0GEN = sqrt(m_part_x->at(k)*m_part_x->at(k)+m_part_y->at(k)*m_part_y->at(k));
         

            PTGEN = sqrt(m_part_px->at(k)*m_part_px->at(k)+m_part_py->at(k)*m_part_py->at(k));

            if (PTGEN < 0.5) continue; // Not an interesting
	    if (D0GEN > 1.) continue; // Not a primary, skip

            evtID = m_part_evtId->at(k);
            stID  = m_part_stId->at(k);

	    //	    if (stID.size()!=1) continue;

            for (int i=0;i<30;++i)
            {
                hit_on_lay[i]=0;
                stub_on_lay[i]=0;
                clus_on_lay[i]=0;
            }

            if (m_dbg)
            {
                cout << endl;
                cout << "Dealing with TP " << k << " of " << m_part_n << endl;
                cout << "TP PDG code is " << m_part_pdg->at(k) << endl;
                cout << "This event contains " << m_pclus << " pixel digis" << endl;
            }

            i_eta = m_part_eta->at(k);

            for (int l=0;l<m_pclus;++l) // Loop over TP digis
            {
        
                if (m_pixclus_layer->at(l)<5) continue; // Don't care about pixel hits

                pix_evtID = m_pixclus_evtID->at(l);
                inTP=false;

                for (unsigned int ll=0;ll<pix_evtID.size();++ll)
                {
                    if (inTP) break;

                    if (m_dbg)
                        cout << l << "/" << evtID << "/" << pix_evtID.at(ll) << endl;

		    if (pix_evtID.at(ll) == evtID) inTP=true;
                }

                if (!inTP) continue;

                if (m_dbg) cout << "In the same evt ID" << endl;
	
                inTP=false;
                pix_stID = m_pixclus_simhitID->at(l);
                
                if (m_dbg) cout << pix_stID.size() << "/" << stID.size() << endl;

                for (unsigned int ll=0;ll<pix_stID.size();++ll)
                {
                    if (inTP) break;

                    for (unsigned int lll=0;lll<stID.size();++lll)
                    {
                        if (m_dbg) cout << pix_stID.at(ll) << "/" << stID.at(lll) << endl;
                        if (pix_stID.at(ll) == stID.at(lll)) inTP=true;
                    }
                }
                
                if (!inTP) continue;


		pix_tp.at(l) = k;


                // This pixel digi comes from the TP
	
                i_lay = m_pixclus_layer->at(l);
                i_mod = m_pixclus_module->at(l);

                ++hit_on_lay[i_lay-5];
		
		if (m_pixclus_type->at(l)<3 && i_lay<=10)
		{
		  //	  cout << "/" << 20+i_lay-5 << endl;
		  ++hit_on_lay[20+i_lay-5];		  
		}

		if (i_lay>10 && m_pixclus_ladder->at(l)<9)
		{
		  //		  cout << "///" << 25+(i_lay-11)%7 << endl;
		  ++hit_on_lay[25+(i_lay-11)%7];
		}
		

                if (m_dbg)
                {
                    cout << endl;
                    cout << "   DIGI    " << endl;
                    cout << "   Coords " << i_lay
                    << " / " << m_pixclus_ladder->at(l)
                    << " / " << i_mod 
		    << " / " << i_eta << endl;
                }

                if (PTGEN<20) entries_pt[i_lay-5][static_cast<int>(5*PTGEN)]+=1;
                if (fabs(i_eta)<2.5 && PTGEN>20) entries_eta[i_lay-5][static_cast<int>(25+10*i_eta)]+=1;

                pix_i  = l;  // The index of the induced DIGI
                clus_i = -1; // The index of the induced CLUSTER
                stub_i = -1; // The index of the induced STUB

                i_seg = m_pixclus_column->at(pix_i);
                
                if (PTGEN<20) digi_pt[i_lay-5][static_cast<int>(5*PTGEN)]+=1;
                if (fabs(i_eta)<2.5 && PTGEN>20) digi_eta[i_lay-5][static_cast<int>(25+10*i_eta)]+=1;
	
                // We have the digis, now look at the clusters
                // Then the official one
                for (int i=0;i<m_tkclus;++i)
                {
                    if (clus_i!=-1) break;
                    if (m_tkclus_nstrips->at(i)>4) continue;
                    if (m_tkclus_layer->at(i)!=i_lay) continue;
                    if (m_tkclus_module->at(i)!=i_mod) continue;
                    if (m_tkclus_seg->at(i)!=i_seg) continue;
                    if (m_tkclus_ladder->at(i)!=m_pixclus_ladder->at(l)) continue;
                    if (m_tkclus_bottom->at(i)!=m_pixclus_bottom->at(l)) continue;
                    if (m_tkclus_tp->at(i)!=k) continue;
                    
                    if (m_dbg)
                        cout << m_tkclus_strip->at(i) << "/" << m_pixclus_row->at(pix_i) << "/" << m_tkclus_nstrips->at(i) << endl;
                    
                    if (m_dbg)
                        cout << m_tkclus_seg->at(i) << "/" << i_seg << endl;
                    
                    if (fabs(m_tkclus_strip->at(i)-m_pixclus_row->at(pix_i))>=m_tkclus_nstrips->at(i)) continue;
                    clus_i = i;
                }

                if (clus_i!=-1)
                {
//                    cout << k << "/" << m_tkclus_tp->at(clus_i) << endl;

		  pix_cl.at(l)=1;

		  ++clus_on_lay[i_lay-5];
    	  
                    if (m_dbg)
                    {
                        cout << endl;
                        cout << "   CLUS FOUND !!   " << endl;
                    }
              
                    if (PTGEN<20) clus_pt[i_lay-5][static_cast<int>(5*PTGEN)]+=1;
                    if (fabs(i_eta)<2.5 && PTGEN>20) clus_eta[i_lay-5][static_cast<int>(25+10*i_eta)]+=1;

                    // We have the clusters, now look at the stubs
                    for (int i=0;i<m_tkstub;++i)
                    {
                        if (stub_i!=-1) break;
                        if (m_tkstub_layer->at(i)!=i_lay) continue;
                        
                        //if (m_tkstub_tp->at(i)!=k) continue;
               
                        if (m_tkstub_clust1->at(i)==clus_i) stub_i = i;
                        if (m_tkstub_clust2->at(i)==clus_i) stub_i = i;
                  
//                        stub_i = i;
                        
//                        if (m_tkstub_tp->at(i)!=k)
//                            cout << k << "/" << m_tkstub_tp->at(i) << " / " << m_tkclus_tp->at(m_tkstub_clust1->at(i)) << " / " << m_tkclus_tp->at(m_tkstub_clust2->at(i)) << " / " << m_tkclus_tp->at(clus_i) << " / " << m_tkstub_clust2->at(i)  <<" / " << m_tkstub_clust1->at(i) << " / " << clus_i  <<endl;
                    }

                    if (stub_i!=-1)
                    {
                        if (m_dbg)
                        {
                            cout << endl;
                            cout << "   STUB FOUND !!   " << endl;
                        }
                  
                        ++stub_on_lay[i_lay-5];
    

			
			if (m_tkstub_type->at(stub_i)<3 && i_lay<=10)
			{

			  ++stub_on_lay[20+i_lay-5];
			}

			if (m_tkstub_layer->at(stub_i)>10 && m_tkstub_ladder->at(stub_i)<9)
			{

			  ++stub_on_lay[25+(i_lay-11)%7];
			}
			

			pix_st.at(l)=1;

                        if (PTGEN<20) stub_pt[i_lay-5][static_cast<int>(5*PTGEN)]+=1;
                        if (fabs(i_eta)<2.5 && PTGEN>20) stub_eta[i_lay-5][static_cast<int>(25+10*i_eta)]+=1;
                    }
                }
		else
		{
		  if (m_dbg)
		    {
		      cout << endl;
		      cout << "   NOCLUS !!   " << endl;
		    }
                  
		}
		
            } // End of loop over SimHits

            for (int i=0;i<30;++i)
            {
                if (hit_on_lay[i]!=0)
                {
                    if (PTGEN<20) entries_pt_lay[i][static_cast<int>(5*PTGEN)]+=1;
                    if (fabs(i_eta)<2.5 && PTGEN>20) entries_eta_lay[i][static_cast<int>(25+10*i_eta)]+=1;
                }

		if (clus_on_lay[i]!=0)
                {
                    if (PTGEN<20) clus_pt_lay[i][static_cast<int>(5*PTGEN)]+=1;
                    if (fabs(i_eta)<2.5 && PTGEN>20) clus_eta_lay[i][static_cast<int>(25+10*i_eta)]+=1;
                }

                if (stub_on_lay[i]!=0)
                {
                    if (PTGEN<20) stub_pt_lay[i][static_cast<int>(5*PTGEN)]+=1;
                    if (fabs(i_eta)<2.5 && PTGEN>20) stub_eta_lay[i][static_cast<int>(25+10*i_eta)]+=1;
                }
            }
        } // End of loop over TPs

	/*
	cout << endl;
	cout << "This event contains " << m_pclus << " pixel digis" << endl;
        for (int k=0;k<m_pclus;++k) // Loop over TPs
        {
	  // if (m_pixclus_layer->at(k)!=10) continue;

	  cout << k 
	       << " / "  << m_pixclus_layer->at(k)
	       << " / "  << m_pixclus_ladder->at(k)
	       << " / "  << m_pixclus_module->at(k)
	       << " / "  << m_pixclus_bottom->at(k)
	       << " / "  << m_pixclus_row->at(k) 
	       << " / "  << m_pixclus_column->at(k)
	       << " --> "  << pix_tp.at(k)
	       << " / "  <<   pix_cl.at(k)
	       << " / "  <<   pix_st.at(k) << endl;
	}

	cout << "Clusters" << endl;
	for (int k=0;k<m_tkclus;++k)
	{
	  //	  if (m_tkclus_layer->at(k)!=10) continue;

	  cout << k 
	       << " / "  << m_tkclus_layer->at(k)
	       << " / "  << m_tkclus_ladder->at(k)
	       << " / "  << m_tkclus_module->at(k)
	       << " / "  << m_tkclus_bottom->at(k)
	       << " / "  << m_tkclus_strip->at(k) 
	       << " / "  << m_tkclus_nstrips->at(k)
	       << " / "  << m_tkclus_seg->at(k)
	       << " --> "  << m_tkclus_tp->at(k) << endl;
	}	


	cout << "Unmatched digis" << endl;
        for (int k=0;k<m_pclus;++k) // Loop over TPs
        {
	  if (pix_tp.at(k)==-1 || pix_cl.at(k)==1) continue;

	  //	  if (m_pixclus_layer->at(k)!=10) continue;

	  cout << k 
	       << " / "  << m_pixclus_layer->at(k)
	       << " / "  << m_pixclus_ladder->at(k)
	       << " / "  << m_pixclus_module->at(k)
	       << " / "  << m_pixclus_bottom->at(k)
	       << " / "  << m_pixclus_row->at(k) 
	       << " / "  << m_pixclus_column->at(k)
	       << " --> "  << pix_tp.at(k)
	       << " / "  <<   pix_cl.at(k)
	       << " / "  <<   pix_st.at(k) << endl;
	}

	cout << "Unmatched clusters" << endl;
	for (int k=0;k<m_tkclus;++k)
	{
	  if (m_tkclus_tp->at(k)!=-1) continue;

	  //	  if (m_tkclus_layer->at(k)!=10) continue;

	  cout << k 
	       << " / "  << m_tkclus_layer->at(k)
	       << " / "  << m_tkclus_ladder->at(k)
	       << " / "  << m_tkclus_module->at(k)
	       << " / "  << m_tkclus_bottom->at(k)
	       << " / "  << m_tkclus_strip->at(k) 
	       << " / "  << m_tkclus_nstrips->at(k)
	       << " / "  << m_tkclus_seg->at(k)
	       << " --> "  << m_tkclus_tp->at(k)
	       << " --> "  << m_tkclus_pdg->at(k) << endl;
	}
	*/

    } // End of loop over events

    // Finally get the efficiencies

    for (int i=0;i<30;++i)
    {
        for (int j=0;j<100;++j)
        {
            if (entries_pt[i][j]>20.)
            {
                digi_pt[i][j] = digi_pt[i][j]/entries_pt[i][j];
                clus_pt[i][j] = clus_pt[i][j]/entries_pt[i][j];
                stub_pt[i][j] = stub_pt[i][j]/entries_pt[i][j];
            }

            if (entries_pt_lay[i][j]>0.)
            {
                stub_pt_lay[i][j] = stub_pt_lay[i][j]/entries_pt_lay[i][j];
		clus_pt_lay[i][j] = clus_pt_lay[i][j]/entries_pt_lay[i][j];
            }
        }

        for (int j=0;j<50;++j)
        {
            if (entries_eta[i][j]>20.)
            {
                digi_eta[i][j] = digi_eta[i][j]/entries_eta[i][j];
                clus_eta[i][j] = clus_eta[i][j]/entries_eta[i][j];
                stub_eta[i][j] = stub_eta[i][j]/entries_eta[i][j];
            }
            
            if (entries_eta_lay[i][j]>0.)
            {
                stub_eta_lay[i][j] = stub_eta_lay[i][j]/entries_eta_lay[i][j];
                clus_eta_lay[i][j] = clus_eta_lay[i][j]/entries_eta_lay[i][j];
            }
        }
    }

    m_tree->Fill();
    m_outfile->Write();
    m_outfile->Close();
}


void efficiencies::initVars()
{

  for (int j=0;j<100;++j) pt_val[j]=0.2*j;
  for (int j=0;j<50;++j)  eta_val[j]=-2.5+0.1*j;

  for (int i=0;i<30;++i) 
  {
      for (int j=0;j<100;++j)
      {
          entries_pt[i][j]=0.;
          entries_pt_lay[i][j]=0.;

          digi_pt[i][j]=0.;
          clus_pt[i][j]=0.;
          stub_pt[i][j]=0.;
          clus_pt_lay[i][j]=0.;
          stub_pt_lay[i][j]=0.;
      }

      for (int j=0;j<50;++j)
      {
          entries_eta[i][j]=0.;
          entries_eta_lay[i][j]=0.;

          digi_eta[i][j]=0.;
          clus_eta[i][j]=0.;
          stub_eta[i][j]=0.;
          clus_eta_lay[i][j]=0.;
          stub_eta_lay[i][j]=0.;
      }
  }
}

void efficiencies::reset()
{
  m_pixclus_row->clear();
  m_pixclus_column->clear();
  m_pixclus_layer->clear();
  m_pixclus_module->clear();
  m_pixclus_bottom->clear();
  m_pixclus_ladder->clear();
  m_pixclus_x->clear();
  m_pixclus_y->clear();
  m_pixclus_type->clear();
  m_pixclus_simhitID->clear();
  m_pixclus_evtID->clear();

  m_part_hits->clear();
  m_part_pdg->clear();
  m_part_px->clear();
  m_part_py->clear();
  m_part_eta->clear();
  m_part_x->clear();
  m_part_y->clear();
  m_part_evtId->clear();
  m_part_stId->clear();
   
  m_tkclus_nstrips->clear();
  m_tkclus_layer->clear(); 
  m_tkclus_module->clear();								       
  m_tkclus_bottom->clear();								       
  m_tkclus_ladder->clear();				      
  m_tkclus_seg->clear();   
  m_tkclus_strip->clear();
  m_tkstub_layer->clear();
  m_tkstub_type->clear();
  m_tkstub_ladder->clear();
  m_tkstub_clust1->clear();
  m_tkstub_clust2->clear();
  m_tkstub_tp->clear();
  m_tkstub_type->clear();
  m_tkclus_tp->clear();
  m_tkclus_pdg->clear();
}


void efficiencies::initTuple(std::string in,std::string out)
{

  L1TT  = new TChain("TkStubs"); 
  Pix   = new TChain("Pixels");
  MC    = new TChain("MC");

  L1TT->Add(in.c_str());
  Pix->Add(in.c_str());
  MC->Add(in.c_str()); 

  m_pixclus_row      = new std::vector<int>;      
  m_pixclus_column   = new std::vector<int>;      
  m_pixclus_layer    = new std::vector<int>;    
  m_pixclus_module   = new std::vector<int>; 
  m_pixclus_ladder   = new std::vector<int>;    
  m_pixclus_type     = new std::vector<int>;  
  m_pixclus_bottom   = new std::vector<int>;  
  m_pixclus_x        = new std::vector<float>;  
  m_pixclus_y        = new std::vector<float>; 
  m_pixclus_simhitID = new std::vector< std::vector<int> >;
  m_pixclus_evtID    = new std::vector< std::vector<int> >; 

  m_part_hits   = new std::vector<int>; 
  m_part_pdg    = new std::vector<int>; 
  m_part_px     = new std::vector<float>; 
  m_part_py     = new std::vector<float>; 
  m_part_eta    = new std::vector<float>;   
  m_part_x      = new std::vector<float>;   
  m_part_y      = new std::vector<float>;   
  m_part_evtId  = new std::vector<int>;
  m_part_stId   = new std::vector< std::vector<int> >;  

  m_tkclus_layer   = new  std::vector<int>;
  m_tkclus_module  = new  std::vector<int>;  
  m_tkclus_bottom  = new  std::vector<int>;  
  m_tkclus_ladder  = new  std::vector<int>;  
  m_tkclus_seg     = new  std::vector<int>; 
  m_tkclus_strip   = new  std::vector<float>; 
  m_tkclus_nstrips = new  std::vector<int>; 
  m_tkclus_pdg     = new  std::vector<int>; 
  
  m_tkstub_layer   = new  std::vector<int>; 
  m_tkstub_type    = new  std::vector<int>; 
  m_tkstub_ladder  = new  std::vector<int>; 
  m_tkstub_clust1  = new  std::vector<int>;  
  m_tkstub_clust2  = new  std::vector<int>;  
  m_tkstub_tp      = new  std::vector<int>; 
  m_tkclus_tp      = new  std::vector<int>;
    

  Pix->SetBranchStatus("*",0);
  MC->SetBranchStatus("*",0);
  L1TT->SetBranchStatus("*",0);

  Pix->SetBranchStatus("PIX_n",1);
  Pix->SetBranchStatus("PIX_layer",1);   
  Pix->SetBranchStatus("PIX_module",1);  
  Pix->SetBranchStatus("PIX_bottom",1);  
  Pix->SetBranchStatus("PIX_ladder",1);  
  Pix->SetBranchStatus("PIX_x",1);       
  Pix->SetBranchStatus("PIX_y",1);       
  Pix->SetBranchStatus("PIX_row",1);     
  Pix->SetBranchStatus("PIX_column",1);  
  Pix->SetBranchStatus("PIX_type",1);  
  Pix->SetBranchStatus("PIX_simhitID",1);     
  Pix->SetBranchStatus("PIX_evtID",1); 

  MC->SetBranchStatus("subpart_n",1);    
  MC->SetBranchStatus("subpart_x",1);    
  MC->SetBranchStatus("subpart_y",1);    
  MC->SetBranchStatus("subpart_px",1);   
  MC->SetBranchStatus("subpart_py",1);   
  MC->SetBranchStatus("subpart_eta",1);
  MC->SetBranchStatus("subpart_pdgId",1);  
  MC->SetBranchStatus("subpart_evtId",1); 
  MC->SetBranchStatus("subpart_stId",1); 

  L1TT->SetBranchStatus("L1TkSTUB_n",1);      
  L1TT->SetBranchStatus("L1TkSTUB_layer",1);  
  L1TT->SetBranchStatus("L1TkSTUB_type",1);  
  L1TT->SetBranchStatus("L1TkSTUB_ladder",1); 
  L1TT->SetBranchStatus("L1TkSTUB_tp",1);     
  L1TT->SetBranchStatus("L1TkSTUB_clust1",1); 
  L1TT->SetBranchStatus("L1TkSTUB_clust2",1); 

  L1TT->SetBranchStatus("L1TkCLUS_n",1);      
  L1TT->SetBranchStatus("L1TkCLUS_layer",1);  
  L1TT->SetBranchStatus("L1TkCLUS_module",1); 
  L1TT->SetBranchStatus("L1TkCLUS_bottom",1);
  L1TT->SetBranchStatus("L1TkCLUS_ladder",1); 
  L1TT->SetBranchStatus("L1TkCLUS_seg",1);
  L1TT->SetBranchStatus("L1TkCLUS_tp",1);
  L1TT->SetBranchStatus("L1TkCLUS_pdgID",1);
  L1TT->SetBranchStatus("L1TkCLUS_strip",1);  
  L1TT->SetBranchStatus("L1TkCLUS_nstrip",1); 
  
  Pix->SetBranchAddress("PIX_n",         &m_pclus);
  Pix->SetBranchAddress("PIX_layer",     &m_pixclus_layer);
  Pix->SetBranchAddress("PIX_module",    &m_pixclus_module);
  Pix->SetBranchAddress("PIX_bottom",    &m_pixclus_bottom);
  Pix->SetBranchAddress("PIX_ladder",    &m_pixclus_ladder);
  Pix->SetBranchAddress("PIX_x",         &m_pixclus_x);
  Pix->SetBranchAddress("PIX_y",         &m_pixclus_y);
  Pix->SetBranchAddress("PIX_row",       &m_pixclus_row);
  Pix->SetBranchAddress("PIX_type",      &m_pixclus_type);
  Pix->SetBranchAddress("PIX_column",    &m_pixclus_column);  
  Pix->SetBranchAddress("PIX_simhitID",  &m_pixclus_simhitID);
  Pix->SetBranchAddress("PIX_evtID",     &m_pixclus_evtID);

  MC->SetBranchAddress("subpart_n",        &m_part_n);    
  MC->SetBranchAddress("subpart_x",        &m_part_x);    
  MC->SetBranchAddress("subpart_y",        &m_part_y);    
  MC->SetBranchAddress("subpart_px",       &m_part_px);   
  MC->SetBranchAddress("subpart_py",       &m_part_py);   
  MC->SetBranchAddress("subpart_eta",      &m_part_eta);  
  MC->SetBranchAddress("subpart_pdgId",    &m_part_pdg);  
  MC->SetBranchAddress("subpart_stId",    &m_part_stId);
  MC->SetBranchAddress("subpart_evtId",    &m_part_evtId);

  L1TT->SetBranchAddress("L1TkSTUB_n",         &m_tkstub);    
  L1TT->SetBranchAddress("L1TkSTUB_layer",     &m_tkstub_layer);
  L1TT->SetBranchAddress("L1TkSTUB_ladder",    &m_tkstub_ladder);
  L1TT->SetBranchAddress("L1TkSTUB_type",      &m_tkstub_type);
  L1TT->SetBranchAddress("L1TkSTUB_tp",        &m_tkstub_tp);
  L1TT->SetBranchAddress("L1TkSTUB_clust1",    &m_tkstub_clust1);
  L1TT->SetBranchAddress("L1TkSTUB_clust2",    &m_tkstub_clust2);


  L1TT->SetBranchAddress("L1TkCLUS_n",         &m_tkclus);
  L1TT->SetBranchAddress("L1TkCLUS_layer",     &m_tkclus_layer);
  L1TT->SetBranchAddress("L1TkCLUS_module",    &m_tkclus_module);
  L1TT->SetBranchAddress("L1TkCLUS_bottom",    &m_tkclus_bottom);
  L1TT->SetBranchAddress("L1TkCLUS_ladder",    &m_tkclus_ladder);
  L1TT->SetBranchAddress("L1TkCLUS_seg",       &m_tkclus_seg);
  L1TT->SetBranchAddress("L1TkCLUS_tp",        &m_tkclus_tp);
  L1TT->SetBranchAddress("L1TkCLUS_strip",     &m_tkclus_strip);
  L1TT->SetBranchAddress("L1TkCLUS_nstrip",    &m_tkclus_nstrips);
  L1TT->SetBranchAddress("L1TkCLUS_pdgID",     &m_tkclus_pdg);

  m_outfile  = new TFile(out.c_str(),"recreate");
  m_tree     = new TTree("Efficiencies","Efficiencies info");

  m_tree->Branch("pt_val",&pt_val,"pt_val[100]");
  m_tree->Branch("eta_val",&eta_val,"eta_val[50]");

  m_tree->Branch("digi_pt",&digi_pt,"digi_pt[30][100]");
  m_tree->Branch("digi_eta",&digi_eta,"digi_eta[30][50]");
  m_tree->Branch("clus_pt",&clus_pt,"clus_pt[30][100]");
  m_tree->Branch("clus_eta",&clus_eta,"clus_eta[30][50]");
  m_tree->Branch("clus_pt_lay",&clus_pt_lay,"clus_pt_lay[30][100]");
  m_tree->Branch("clus_eta_lay",&clus_eta_lay,"clus_eta_lay[30][50]");
  m_tree->Branch("stub_pt",&stub_pt,"stub_pt[30][100]");
  m_tree->Branch("stub_eta",&stub_eta,"stub_eta[30][50]");
  m_tree->Branch("stub_pt_lay",&stub_pt_lay,"stub_pt_lay[30][100]");
  m_tree->Branch("stub_eta_lay",&stub_eta_lay,"stub_eta_lay[30][50]");
}
