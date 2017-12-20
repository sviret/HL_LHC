// Class for the L1 track reco efficiency tests
// For more info, look at the header file

#include "track_eff.h"

track_eff::track_eff(std::string filename, std::string secfilename, 
		     std::string outfile, int nevt, float pt_min, 
                     float d0_min, bool dbg) :
stub_x(new std::vector<float>),
stub_y(new std::vector<float>),
stub_z(new std::vector<float>),
stub_layer(new std::vector<int>),
stub_ladder(new std::vector<int>),
stub_module(new std::vector<int>),
stub_seg(new std::vector<int>),
stub_strip(new std::vector<float>),
stub_did(new std::vector<int>),
stub_sw(new std::vector<float>),
stub_sw_c(new std::vector<float>),
stub_sw_truth(new std::vector<float>),
stub_cw1(new std::vector<int>),
stub_cw2(new std::vector<int>),
stub_cn1(new std::vector<int>),
stub_cn2(new std::vector<int>),
stub_tp(new std::vector<int>),
stub_fe(new std::vector<int>),
stub_ptGEN(new std::vector<float>),
stub_d0GEN(new std::vector<float>),
stub_r0GEN(new std::vector<float>),
stub_etaGEN(new std::vector<float>),
stub_z0GEN(new std::vector<float>),
stub_insec(new std::vector<int>),
stub_inpatt(new std::vector<int>),
stub_intc(new std::vector<int>),
stub_intrk(new std::vector<int>),
stub_sec(new std::vector< std::vector<int> >),
part_pdg(new std::vector<int>),
part_nsec(new std::vector<int>),
part_nhits_fe(new std::vector<int>),
part_nstubs(new std::vector<int>),
part_nhits(new std::vector<int>),
part_npatt(new std::vector<int>),
part_ntc(new std::vector<int>),
part_ntrk(new std::vector<int>),
part_pt(new std::vector<float>),
part_rho(new std::vector<float>),
part_d0(new std::vector<float>),
part_z0(new std::vector<float>),
part_eta(new std::vector<float>),
part_phi(new std::vector<float>),
part_dist(new std::vector<float>),
part_sec(new std::vector< std::vector<int> >),
patt_sec(new std::vector<int>),
patt_id(new std::vector<int>),
patt_miss(new std::vector<int>),
patt_insec(new std::vector<int>),
patt_insec_full(new std::vector<int>),
patt_comb(new std::vector<int>),
patt_gcomb(new std::vector<int>),
patt_parts(new std::vector< std::vector<int> >),
patt_stubs(new std::vector< std::vector<int> >),
patt_rank(new std::vector<int>),
patt_rank_full(new std::vector<int>),
tc_sec(new std::vector<int>),
tc_id(new std::vector<int>),
tc_parts(new std::vector< std::vector<int> >),
tc_stubs(new std::vector< std::vector<int> >),
tc_pt(new std::vector<float>),
tc_eta(new std::vector<float>),
tc_z(new std::vector<float>),
tc_phi(new std::vector<float>),
tc_pt_t(new std::vector<float>),
tc_eta_t(new std::vector<float>),
tc_z_t(new std::vector<float>),
tc_phi_t(new std::vector<float>),
tc_PDG_t(new std::vector<int>),
trk_sec(new std::vector<int>),
trk_parts(new std::vector< std::vector<int> >),
trk_stubs(new std::vector< std::vector<int> >),
trk_pt(new std::vector<float>),
trk_eta(new std::vector<float>),
trk_z(new std::vector<float>),
trk_phi(new std::vector<float>),
trk_pt_t(new std::vector<float>),
trk_eta_t(new std::vector<float>),
trk_z_t(new std::vector<float>),
trk_phi_t(new std::vector<float>),
trk_PDG_t(new std::vector<int>),
dtc_id(new std::vector<int>),
dtc_mult(new std::vector<int>),
dtc_mult_FE(new std::vector<int>)
{  
  m_dbg    = dbg;
  m_d0_min = d0_min;
  m_pt_min = pt_min;
  has_patt = true;

  m_tilted = true;
    
  int n_tilted_rings[6];
  int n_flat_rings[6];
  
  for (int i=0; i < 6; ++i) n_tilted_rings[i]=0;
  for (int i=0; i < 6; ++i) n_flat_rings[i]=0;
  
  if (m_tilted)
  {
    n_tilted_rings[0]=12;
    n_tilted_rings[1]=12;
    n_tilted_rings[2]=12;
    n_flat_rings[0]=7;
    n_flat_rings[1]=11;
    n_flat_rings[2]=15;
  }
    
  for (int i=0; i < 6; ++i)
  {
    for (int j=0; j < 3; ++j)
    {
      limits[i][j]=0;
      
      if (n_tilted_rings[i]==0) continue;
      
      limits[i][j]=(j%2)*n_flat_rings[i]+(j>0)*n_tilted_rings[i];
    }
  }  
    
  cout << "PT A" << endl;
  track_eff::initTuple(filename,outfile);
  cout << "PT B" << endl;  

  if (has_patt)
  {
    if (!track_eff::convert_towers(secfilename)) return; 
  }
  else
  {
    if (!track_eff::convert_cabling(secfilename)) return; 
  }

  cout << "PT C" << endl;
  track_eff::do_test(nevt); // Launch the test loop over n events
}


/////////////////////////////////////////////////////////////////////////////////
//
// ==> track_eff::do_test(int nevt)
//
// Main method, where the efficiency calculations are made
//
/////////////////////////////////////////////////////////////////////////////////

void track_eff::do_test(int nevt)
{
  int id,dtcid;

  if (!has_patt) m_sec_mult=1; // We don't deal with trigger towers then

  const int m_nsec = m_sec_mult; // How many sectors are in the file (from convert_tower)

  //  cout << nevt << " / " << m_L1TT->GetEntries() << endl;

  int ndat = std::min(nevt,static_cast<int>(m_L1TT->GetEntries())); // How many events will we test

  cout << "Starting a test loop over " << ndat << " events..." << endl;
  if (has_patt) cout << "... using " << m_nsec << " trigger sectors..." << endl;

  int is_sec_there[m_nsec];
  int ladder,module,layer;
  int n_per_lay[20];
  int n_per_lay_patt[20];
  int nhits_p;

  int patt_psec[48];
  int patt_psec_full[48];
  int n_comb;
  float sval,lay;
  int n_inlayer[20];
  
  m_primaries.clear();
  
  std::vector<int> the_sectors;
  std::vector<int> prim_stubs;
  
  std::vector<int> stubs;
  std::vector<int> parts;
  std::vector<int> msec;
  std::vector<int> psec;
  
  std::vector<std::vector <float> > lay_content;
  std::vector<float> stub_list;
  
  std::vector<int> dtc_list;
  std::vector<int> dtc_list_FE;
  
  // Loop over the events
    
  for (int i=0;i<ndat;++i)
  {
    for (int j=0;j<48;++j) patt_psec[j]=0;
    for (int j=0;j<48;++j) patt_psec_full[j]=0;

    track_eff::reset();

    dtc_list.clear();
    dtc_list_FE.clear();
	
    for (int j=0;j<2000;++j) dtc_list.push_back(0);	
    for (int j=0;j<2000;++j) dtc_list_FE.push_back(0);	

    nb_patterns =0;
    nb_tcs =0;
    nb_tracks =0;
    
    m_L1TT->GetEntry(i);
    m_PIX->GetEntry(i);
    if (has_patt) m_PATT->GetEntry(i);

    if (i%1000==0)
      cout << "Processed " << i << "/" << ndat << endl;

    if (m_stub == 0) continue; // No stubs, don't go further

    m_npu = m_PU;

    evt         = i;
    n_stub_total=m_stub;
    n_clus_total=m_clus;
    n_pix_total =m_pix;
    n_patt      =nb_patterns;
    n_tc        =nb_tcs;
    n_track     =nb_tracks;
    
    parts.clear();
    msec.clear();
    psec.clear();

    // Here we fill the pattern and tracks containers
    // the info is already there

    if (nb_patterns>0)
    {
      for(int kk=0;kk<nb_patterns;kk++)
      {
	patt_psec[m_pattsecid.at(kk)]+=1;
	if (m_pattmiss.at(kk)==0) patt_psec_full[m_pattsecid.at(kk)]+=1;
	patt_rank->push_back(int(patt_psec[m_pattsecid.at(kk)]));
	(m_pattmiss.at(kk)==0)
	  ? patt_rank_full->push_back(patt_psec_full[m_pattsecid.at(kk)])
	  : patt_rank_full->push_back(-1);
	patt_sec->push_back(m_pattsecid.at(kk));
	patt_id->push_back(m_pattid.at(kk));
	patt_miss->push_back(m_pattmiss.at(kk));
	patt_parts->push_back(parts);
	patt_stubs->push_back(m_pattlinks.at(kk));
	
	n_comb=1;
	lay_content.clear();
        
	for (int kb=0;kb<20;++kb)
	{
	  stub_list.clear();
	  lay_content.push_back(stub_list);
	}
        
	for (unsigned int kb=0;kb<patt_stubs->at(kk).size();++kb) // Loop over stubs in the pattern
	{
	  sval = m_stub_strip[patt_stubs->at(kk).at(kb)];
	  lay  = m_stub_layer[patt_stubs->at(kk).at(kb)]-5;                
	  lay_content.at(lay).push_back(sval);
	}
                
	// We now know how many stubs per layer are in the pattern
        
	for (int kb=0;kb<20;++kb)
	{
	  n_inlayer[kb] = lay_content.at(kb).size();
	  if (n_inlayer[kb]==0) continue;
	  n_comb = n_comb*n_inlayer[kb];
	}

	n_patt_comb+=n_comb;
	patt_comb->push_back(n_comb);
	patt_gcomb->push_back(0);
      }

      for(int kk=0;kk<nb_patterns;kk++)
      {
	patt_insec->push_back(patt_psec[m_pattsecid.at(kk)]);
	patt_insec_full->push_back(patt_psec_full[m_pattsecid.at(kk)]);
      }
    }
      
    // Look at the track candidates then
    if (nb_tcs>0)
    {
      for(int kk=0;kk<nb_tcs;kk++)
      {
	tc_sec->push_back(m_tcsecid.at(kk));
	tc_id->push_back(m_tcid.at(kk));
	tc_pt->push_back(m_tcpt.at(kk));
	tc_eta->push_back(m_tceta.at(kk));
	tc_phi->push_back(m_tcphi.at(kk));
	tc_z->push_back(m_tcz.at(kk));
	tc_parts->push_back(parts);
	tc_stubs->push_back(m_tclinks.at(kk));
	tc_pt_t->push_back(0);
	tc_eta_t->push_back(0);
	tc_phi_t->push_back(0);
	tc_z_t->push_back(0);
	tc_PDG_t->push_back(0);
      }
    }

    // And finally at the L1Tracks
    if (nb_tracks>0)
    {
      for(int kk=0;kk<nb_tracks;kk++)
      {
	trk_sec->push_back(m_trksecid.at(kk));
	trk_pt->push_back(m_trkpt.at(kk));
	trk_eta->push_back(m_trketa.at(kk));
	trk_phi->push_back(m_trkphi.at(kk));
	trk_z->push_back(m_trkz.at(kk));
	trk_parts->push_back(parts);
	trk_stubs->push_back(m_trklinks.at(kk));
	trk_pt_t->push_back(0);
	trk_eta_t->push_back(0);
	trk_phi_t->push_back(0);
	trk_z_t->push_back(0);
	trk_PDG_t->push_back(0);
      }
    }

    //
    // First, we make a loop on the stubs in order to check
    // how many interesting tracks we have in the event
    //
    // By interesting we mean a track with pt>m_pt_min GeV/c with an IP<m_d0_min cm
    //
    // The info is then stored in the m_primaries vector
    //
    // IMPORTANT NOTE: we base the efficiency calculation on interesting tracks
    // which have induced at least one stub. Interesting tracks which are not
    // inducing any stub are not accounted for.
    //

    bool already_there = 0;
        
    float eta_s,r_s,thick,pitch;
    float tilt,theta0,pt;
    int nchip,ncic,ncn1,ncn2;
    bool isPS;
    
    for (int j=0;j<m_stub;++j)
    {
      msec.clear();
      psec.clear();
      msec.push_back(-1);
      
      layer  = m_stub_layer[j];
      ladder = m_stub_ladder[j];
      module = m_stub_module[j];
      id     = 10000*layer+100*ladder+module; // Get the module ID

      nchip=0;
      ncic=0;
      ncn1=0;
      ncn2=0;
      isPS = false;

      if (layer==5 && ladder>12 && ladder<19) n_stub_total_c++;

      if (layer<=7 || (layer>=11 && ladder<=9)) isPS = true;
      if (layer>=11 && (layer-11)%7>1 && ladder>=6) isPS =false;


      for (int jj=0;jj<m_stub;++jj)
      {
	if (layer   != m_stub_layer[jj]) continue;
	if (ladder  != m_stub_ladder[jj]) continue;
	if (module  != m_stub_module[jj]) continue;
	if (!isPS && (m_stub_segment[j] != m_stub_segment[jj])) continue;
	if (isPS && (m_stub_segment[j]/16 != m_stub_segment[jj])/16) continue;
	++ncic;
	
	if (m_stub_chip[j] != m_stub_chip[jj]) continue;
	++nchip;	   
      }
      /*
      for (int jj=0;jj<m_clus;++jj)
      {	    
	//if (m_stub_tp[j] != m_clus_tp[jj]) continue;
	if (layer   != m_clus_layer[jj]) continue;
	if (ladder  != m_clus_ladder[jj]) continue;
	if (module  != m_clus_module[jj]) continue;
	if (!isPS && (m_clus_segment[jj] != m_stub_segment[j])) continue;
	if (isPS && m_clus_bottom[jj] && (m_stub_segment[j]/16 != m_clus_segment[jj])/16) continue;
	if (isPS && !m_clus_bottom[jj] && (m_stub_segment[j]/16 != m_clus_segment[jj])) continue;
	
	//cout << m_stub_chip[j] << " / " << int(m_clus_strip[jj]/120) << endl;

	if (isPS && m_clus_bottom[jj] && (m_stub_chip[j] != int(m_clus_strip[jj]/120))) continue;
	if (isPS && !m_clus_bottom[jj] && (m_stub_chip[j] != int(m_clus_strip[jj]/127))) continue;
	if (!isPS && (m_stub_chip[j] != int(m_clus_strip[jj]/127))) continue;

	(m_clus_bottom[jj]) 
	  ? ncn1+=m_clus_cw[jj]
	  : ncn2+=m_clus_cw[jj];
      }
      */
      if (has_patt)
      {
	if (m_modules.at(id).size()>0)
	{
	  for (unsigned int kk=0;kk<m_modules.at(id).size();++kk) // In which sector the module is
	  {
	    msec.push_back(m_modules.at(id).at(kk));

	    if (std::abs(m_stub_sw[j])>7) continue;

	    ++dtc_list.at(m_modules.at(id).at(kk));
	  }
	}
      }          
      else
      {
	if ( m_modules_to_DTC.at(id).size()>0)
	{
	  dtcid  = m_modules_to_DTC.at(id).at(0);
	  
	  ++dtc_list.at(dtcid);
	  
	  if (m_stub_c[j]<100) ++dtc_list_FE.at(dtcid);
	}
      }
  
      eta_s  =-log(tan(atan2(sqrt(m_stub_y[j]*m_stub_y[j]+m_stub_x[j]*m_stub_x[j]),m_stub_z[j])/2));
      r_s    = sqrt(m_stub_x[j]*m_stub_x[j]+m_stub_y[j]*m_stub_y[j]);
 
      thick  = 0.18;
      pitch  = 0.09;
 
      tilt   = 0;
      if (m_stub_layer[j]<11) tilt=2.*atan(1.);
 
      if (m_stub_layer[j]<8 || (m_stub_layer[j]>10 && m_stub_ladder[j]<9))
      {
	thick = 0.16;
	pitch = 0.1;
	
	if (m_stub_layer[j]>=10) thick=0.4;
      }

      if (m_stub_layer[j]==5) thick=0.26;
      if (m_stub_layer[j]>=11 && m_stub_ladder[j]==9)  thick=0.4;
      if (m_stub_layer[j]==13 && m_stub_ladder[j]==10) thick=0.4;
      if (m_stub_layer[j]==20 && m_stub_ladder[j]==10) thick=0.4;
      if (m_stub_layer[j]==14 && m_stub_ladder[j]==10) thick=0.4;
      if (m_stub_layer[j]==21 && m_stub_ladder[j]==10) thick=0.4;
      if (m_stub_layer[j]==15 && m_stub_ladder[j]==10) thick=0.4;
      if (m_stub_layer[j]==22 && m_stub_ladder[j]==10) thick=0.4;
      if (m_stub_layer[j]==15 && m_stub_ladder[j]==11) thick=0.4;
      if (m_stub_layer[j]==22 && m_stub_ladder[j]==11) thick=0.4;
      
      theta0=2*atan(exp(-eta_s));

      stub_x->push_back(m_stub_x[j]);
      stub_y->push_back(m_stub_y[j]);
      stub_z->push_back(m_stub_z[j]);
      stub_layer->push_back(m_stub_layer[j]);
      stub_ladder->push_back(m_stub_ladder[j]);
      stub_module->push_back(m_stub_module[j]);
      stub_seg->push_back(m_stub_segment[j]);
      stub_strip->push_back(m_stub_strip[j]);
      stub_chip->push_back(1000*ncic+nchip);
      stub_sw->push_back(m_stub_sw[j]);
      //stub_sw_c->push_back(m_stub_sw_c[j]);
      stub_sw_c->push_back(0);
      stub_sw_truth->push_back(0);
      stub_cw1->push_back(-1);
      stub_cw2->push_back(-1);
      stub_cn1->push_back(ncn1);
      stub_cn2->push_back(ncn2);
      stub_tp->push_back(-1);
      stub_fe->push_back(0);
      /*
      if (ncn1+ncn2>12 && layer>7 && layer<11) stub_fe->at(j) = 1;
      if (ncn1+ncn2<2 && layer>7 && layer<11) stub_fe->at(j) = 1;
      
      if (ncn1+ncn2<1 && layer<=7) stub_fe->at(j) = 1;
      if (ncn1+ncn2>20 && layer<=7) stub_fe->at(j) = 1;

      if (ncn1+ncn2<1 && layer>10 && (layer-11)%7>1 && ladder<6) stub_fe->at(j) = 1;
      if (ncn1+ncn2>15 && layer>10 && (layer-11)%7>1 && ladder<6) stub_fe->at(j) = 1;
    
    if (ncn1+ncn2<1 && layer>10 && (layer-11)%7>1 && ladder>=6) stub_fe->at(j) = 1;
    if (ncn1+ncn2>10 && layer>10 && (layer-11)%7>1 && ladder>=6) stub_fe->at(j) = 1;
  
  if (ncn1+ncn2<1 && layer>10 && (layer-11)%7<=1 && ladder<9) stub_fe->at(j) = 1;
  if (ncn1+ncn2>18 && layer>10 && (layer-11)%7<=1 && ladder<9) stub_fe->at(j) = 1;
if (ncn1+ncn2<1 && layer>10 && (layer-11)%7<=1 && ladder>=9) stub_fe->at(j) = 1;
if (ncn1+ncn2>8 && layer>10 && (layer-11)%7<=1 && ladder>=9) stub_fe->at(j) = 1;
      */

      if (isPS && nchip>3) stub_fe->at(j) = 1;
      if (!isPS && nchip>2) stub_fe->at(j) = 1;

      //     if (m_stub_c[j]>300) stub_fe->at(j) = 1;
      //      if (m_stub_c[j]<300 && m_stub_c[j]>50) stub_fe->at(j) = 2;

      stub_ptGEN->push_back(-1);
      stub_d0GEN->push_back(-1);
      stub_r0GEN->push_back(-1);
      stub_etaGEN->push_back(-1);
      stub_z0GEN->push_back(-1);
      stub_insec->push_back(msec.size()-1);
      stub_inpatt->push_back(0);
      stub_intrk->push_back(0);
      stub_intc->push_back(0);
      stub_sec->push_back(msec);

      if (m_stub_tp[j]<0) continue; // Bad stub (unmatched
      
      theta0=2*atan(exp(-m_stub_etaGEN[j]));

      stub_ptGEN->at(j)   = abs(m_stub_pdg[j])/m_stub_pdg[j]*sqrt(m_stub_pxGEN[j]*m_stub_pxGEN[j]+m_stub_pyGEN[j]*m_stub_pyGEN[j]);
      stub_d0GEN->at(j)   = m_stub_Y0[j]*cos(atan2(m_stub_pyGEN[j],m_stub_pxGEN[j]))-m_stub_X0[j]*sin(atan2(m_stub_pyGEN[j],m_stub_pxGEN[j]));
      stub_r0GEN->at(j)   = sqrt(m_stub_X0[j]*m_stub_X0[j]+m_stub_Y0[j]*m_stub_Y0[j]);
      stub_etaGEN->at(j)  = m_stub_etaGEN[j];
      stub_z0GEN->at(j)   = m_stub_Z0[j];
      stub_sw_truth->at(j)= 1/stub_ptGEN->at(j)*0.057*sin(theta0)/cos(theta0-tilt)*r_s*thick/pitch;
      
      if (m_clus_bottom[m_stub_clust1[j]])
      {
	stub_cw1->at(j)    = m_clus_cw[m_stub_clust1[j]];
	stub_cw2->at(j)    = m_clus_cw[m_stub_clust2[j]];
      }
      else
      {
	stub_cw2->at(j)    = m_clus_cw[m_stub_clust1[j]];
	stub_cw1->at(j)    = m_clus_cw[m_stub_clust2[j]];
      }

      //stub_cw1->at(j)     = m_clus_cw[m_stub_clust1[j]];
      //stub_cw2->at(j)     = m_clus_cw[m_stub_clust2[j]];


      // Basic particle selection (pt and d0 cuts)
      if (sqrt(m_stub_pxGEN[j]*m_stub_pxGEN[j]+m_stub_pyGEN[j]*m_stub_pyGEN[j])<m_pt_min) continue;      
      if (fabs(stub_r0GEN->at(j))>m_d0_min) continue;

      already_there = false;



      for (unsigned int k=0;k<m_primaries.size();++k) // Check if it's already been found
      {
	if (m_primaries.at(k).at(0)!=m_stub_tp[j]) continue;
	
	stub_tp->at(j)=k;
	already_there = true;
	m_primaries.at(k).push_back(j); // If yes, link the stub to the corresponding part
	break;
      }

      if (already_there == true) continue;

      // Here we have a new primary, we create a new entry

      ++n_part;
      part_pdg->push_back(m_stub_pdg[j]);
      part_nsec->push_back(0);
      part_sec->push_back(psec);
      part_nhits->push_back(0);
      part_nstubs->push_back(0);
      part_nhits_fe->push_back(0);
      part_npatt->push_back(0);
      part_ntc->push_back(0);
      part_ntrk->push_back(0);
      part_pt->push_back(sqrt(m_stub_pxGEN[j]*m_stub_pxGEN[j]+m_stub_pyGEN[j]*m_stub_pyGEN[j]));
      part_rho->push_back(sqrt(m_stub_X0[j]*m_stub_X0[j]+m_stub_Y0[j]*m_stub_Y0[j]));
      part_d0->push_back(m_stub_Y0[j]*cos(atan2(m_stub_pyGEN[j],m_stub_pxGEN[j]))-m_stub_X0[j]*sin(atan2  (m_stub_pyGEN[j],m_stub_pxGEN[j])));
      part_z0->push_back(m_stub_Z0[j]);
      part_eta->push_back(m_stub_etaGEN[j]);
      part_phi->push_back(atan2(m_stub_pyGEN[j],m_stub_pxGEN[j]));
      part_dist->push_back(-1);
      
      prim_stubs.clear();
      prim_stubs.push_back(m_stub_tp[j]); // The TP id
      prim_stubs.push_back(j); // The first stub ID
      m_primaries.push_back(prim_stubs);
      
      stub_tp->at(j)=m_primaries.size()-1; // Link the stub to the local tree now
    } // End of loop over stubs

    for (int k=0;k<2000;++k)
    {
      if (dtc_list.at(k)==0) continue;

      dtc_id->push_back(k);
      dtc_mult->push_back(dtc_list.at(k));
      dtc_mult_FE->push_back(dtc_list_FE.at(k));
    }
    
    if (m_primaries.size()>300)
      std::cout << "Found " << m_primaries.size()
		<< " interesting particles inducing at least one stub in the event" << i << std::endl;

    // Look at isolation
    
    float dist_min, dist, deta, dphi;
    
    for (int k=0;k<n_part;++k)
    {
      dist_min = 1000;
      if (part_pt->at(k)  < 3) continue;
      if (fabs(part_rho->at(k)) > 1) continue;

      for (int l=0;l<n_part;++l)
      {
	if (k==l) continue;
	if (part_pt->at(l)  < 3) continue;
	if (part_rho->at(l) > 1) continue;
	
	deta = part_eta->at(k)-part_eta->at(l);
	dphi = part_phi->at(k)-part_phi->at(l);
	
	if (part_phi->at(k)*part_phi->at(l)<0 && fabs(dphi)>4*atan(1))
	{
	  if (part_phi->at(k)<0) dphi = part_phi->at(k) + 8*atan(1) -part_phi->at(l) ;
	  if (part_phi->at(l)<0) dphi = part_phi->at(k) - 8*atan(1) -part_phi->at(l) ;
	}

	dist = sqrt(deta*deta+dphi*dphi);
	if (dist < dist_min) dist_min = dist;
      }
      
      part_dist->at(k) = dist_min;
    }

    //
    // Then, in the second loop, we test all the interesting particles
    // which are inducing at least 4 stubs in the detector
    // in order to check if they have matched a pattern or a track
    //
    
    int idx = 0;
    int ns=0;
    int ns_f=0;
    int discrep;
      
    for (unsigned int k=0;k<m_primaries.size();++k)
    {
      if (m_primaries.at(k).size()<1) continue; // Less then 4 stubs, give up this one
                                                // as this is the lowest possible threshold

      for (int j=0;j<20;++j) n_per_lay[j]=0;
      for (int j=0;j<20;++j) n_per_lay_patt[j]=0;
      for (int j=0;j<m_nsec;++j) is_sec_there[j]=0;
            
      int nhits=0;
      int nhits_fe=0;
      nsec=0;
      eta=0.;
      pt=0.;
      phi=0.;
      d0=0.;
      z0=0.;
      int npatt=-1;
      int ntc=-1;
      int ntrack=-1;
      eta_f=0.;
      pt_f=0.;
      phi_f=0.;
      z0_f=0.;
      nstubs=0;
      nstubs_f=0;
      
      // Here we make a loop over prim stubs to compute the value of NHits
      //
      // NHits is the number of distinct layers/disks in which the track has let stubs
      // this corresponds to the number of active layers of the track, in AM langage
      //
                  
      for (unsigned int j=1;j<m_primaries.at(k).size();++j)
      {	
	idx    = m_primaries.at(k).at(j); // stub index

	// First of all we compute the ID of the stub's module
	layer  = m_stub_layer[idx];
	ladder = m_stub_ladder[idx];
	module = m_stub_module[idx];
	id     = 10000*layer+100*ladder+module; // Get the module ID
	                
	if (has_patt)
	{
	  if (m_modules.at(id).size()>0)
	  {
	    // Just to be precise, sectors are numbered starting from 1 in convert()
	    
	    for (unsigned int kk=0;kk<m_modules.at(id).size();++kk) // In which sector the module is
	    {
	      ++is_sec_there[m_modules.at(id).at(kk)];
	    }
	  }
	}
                
	if (stub_fe->at(idx)==0) ++n_per_lay_patt[layer-5];

	++n_per_lay[layer-5];
	pt  = sqrt(m_stub_pxGEN[idx]*m_stub_pxGEN[idx]+m_stub_pyGEN[idx]*m_stub_pyGEN[idx]);
	eta = m_stub_etaGEN[idx];
	phi = atan2(m_stub_pyGEN[idx],m_stub_pxGEN[idx]);
	d0  = m_stub_Y0[j]*cos(atan2(m_stub_pyGEN[j],m_stub_pxGEN[j]))-m_stub_X0[j]*sin(atan2(m_stub_pyGEN[j],m_stub_pxGEN[j]));
	z0  = m_stub_Z0[idx];
	pdg = m_stub_pdg[idx];
      }
            
      // For the sector efficiency we proceed as follow
      // We check if the track has let at least 4 stubs in
      // 4 different layer/disk (nhits>=4). If yes, we compute
      // the number of sectors containing at least 4 of those
      // stubs (nsec)
      
      // First we get the number of different layers/disks hit by the track (nhits)
      for (int kk=0;kk<20;++kk)
      {
	if (n_per_lay[kk]!=0) ++nhits;
	if (n_per_lay_patt[kk]!=0) ++nhits_fe;
      }
      
      // Then we get the number of sectors containing more than 4 primary hits
      for (int kk=0;kk<m_nsec;++kk)
      {
	if (is_sec_there[kk]>=4) part_sec->at(k).push_back(kk);
	if (is_sec_there[kk]>=4) ++nsec;
      }

      part_nsec->at(k)  = nsec;
      part_nhits->at(k) = nhits;
      part_nstubs->at(k)= static_cast<int>(m_primaries.at(k).size());
      part_nhits_fe->at(k) = nhits_fe;	  
      
      // Finally we do the pattern/TC/track loops

      npatt    = 0;           // The patterns containing at least 4 prim hits
      ntc      = 0;           // The tcs containing at least 4 prim hits
      ntrack   = 0;           // The tracks containing at least 4 prim hits
      
      if (nsec==0)
      {
	part_npatt->at(k) = npatt;
	part_ntc->at(k)   = ntc;
	part_ntrk->at(k)  = ntrack;
	
	continue;
      }

      if (nb_patterns>0)
      {
	for(int kk=0;kk<nb_patterns;kk++)
	{
	  for (int j=0;j<20;++j) n_per_lay_patt[j]=0;
	  for (int j=0;j<20;++j) n_per_lay[j]=0;
	  
	  // First we count the number of prim stubs in the pattern
	  
	  if (m_pattlinks.at(kk).size()==0) continue;
	  
	  for(unsigned int m=0;m<m_pattlinks.at(kk).size();m++) // Loop over pattern stubs
	  {
	    n_per_lay[m_stub_layer[m_pattlinks.at(kk).at(m)]-5]++;
	    
	    if (m_stub_tp[m_pattlinks.at(kk).at(m)]==m_primaries.at(k).at(0)) // Stub belongs to the primary
	      n_per_lay_patt[m_stub_layer[m_pattlinks.at(kk).at(m)]-5]++;

	    if (stub_inpatt->at(m_pattlinks.at(kk).at(m))==0)
	    {
	      stub_inpatt->at(m_pattlinks.at(kk).at(m))=1;
	      ++n_stub_patt;
	    }
	  }
	  
	  nhits_p=0;
	  discrep=0;
	  
	  for (int jk=0;jk<20;++jk)
	  {
	    if (n_per_lay_patt[jk]!=0) ++nhits_p;
	    if (n_per_lay_patt[jk]==0 && n_per_lay[jk]!=0) ++discrep;
	  }
	  
	  if (nhits_p>=4)
	  {
	    ++npatt; // More than 4, the pattern is good
	    patt_parts->at(kk).push_back(k);
	  }

	  if (discrep==0)
	  {
	    patt_gcomb->at(kk) = patt_gcomb->at(kk)+1;
	  }
	}
      } // End of loop over patterns
      
      if (nb_tcs>0)
      {
	for(int kk=0;kk<nb_tcs;kk++)
	{
	  ns      = 0;
	  ns_f    = 0;
	  
	  // First we count the number of prim stubs in the track
	  
	  if (m_tclinks.at(kk).size()==0) continue;
	  
	  for(unsigned int m=0;m<m_tclinks.at(kk).size();m++) // Loop over track stubs
	  {
	    ++ns;
	    if (m_stub_tp[m_tclinks.at(kk).at(m)]==m_primaries.at(k).at(0)) ++ns_f;
	    
	    if (stub_intc->at(m_tclinks.at(kk).at(m))==0)
	    {
	      stub_intc->at(m_tclinks.at(kk).at(m))=1;
	      ++n_stub_tc;
	    }
	  }

	  if ((ns-ns_f)<=1)
	  {
	    tc_parts->at(kk).push_back(k);
	    tc_pt_t->at(kk)  = pt;
	    tc_eta_t->at(kk) = eta;
	    tc_phi_t->at(kk) = phi;
	    tc_z_t->at(kk)   = z0;
	    tc_PDG_t->at(kk) = pdg;
	    ++ntc;
	  }
	}
      } // End of loop over tracks
      
      if (nb_tracks>0)
      {
	for(int kk=0;kk<nb_tracks;kk++)
	{
	  ns      = 0;
	  ns_f    = 0;
	  
	  // First we count the number of prim stubs in the track
	
	  if (m_trklinks.at(kk).size()==0) continue;
                    
	  for(unsigned int m=0;m<m_trklinks.at(kk).size();m++) // Loop over track stubs
	  {
	    ++ns;
	    if (m_stub_tp[m_trklinks.at(kk).at(m)]==m_primaries.at(k).at(0)) ++ns_f;

	    if (stub_intrk->at(m_trklinks.at(kk).at(m))==0)
	    {
	      stub_intrk->at(m_trklinks.at(kk).at(m))=1;
	      ++n_stub_trk;
	    }
	  }

	  if ((ns-ns_f)<=1)
	  {
	    trk_parts->at(kk).push_back(k);
	    trk_pt_t->at(kk)  = pt;
	    trk_eta_t->at(kk) = eta;
	    trk_phi_t->at(kk) = phi;
	    trk_z_t->at(kk)   = z0;
	    trk_PDG_t->at(kk) = pdg;
	    ++ntrack;
	  }
	}
      } // End of loop over tracks

      part_npatt->at(k) = npatt;
      part_ntc->at(k)   = ntc;
      part_ntrk->at(k)  = ntrack;
      
    } // End of loop over primaries
    
    m_finaltree->Fill();
  }

  m_outfile->Write();
  m_outfile->Close();
}

/////////////////////////////////////////////////////////////////////////////////
//
// ==> track_eff::initTuple(std::string test,std::string out)
//
// This method opens the input  roottuple
//
/////////////////////////////////////////////////////////////////////////////////


void track_eff::initTuple(std::string test,std::string out)
{
  m_PIX    = new TChain("Pixels");    // Tree containing the Digi info
  m_L1TT   = new TChain("TkStubs");   // Tree containing the Stub info 
  m_PATT   = new TChain("L1tracks");  // Tree containing the L1Track reco info

  // Input data file 

  std::size_t found = test.find(".root");

  // Case 1, it's a root file
  if (found!=std::string::npos)
  {
    m_L1TT->Add(test.c_str());
    m_PIX->Add(test.c_str());
    if (has_patt) m_PATT->Add(test.c_str());
  }
  else // This is a list provided into a text file
  {
    std::string STRING;
    std::ifstream in(test.c_str());
    if (!in) 
    {
      std::cout << "Please provide a valid data filename list" << std::endl; 
      return;
    }    
  
    while (!in.eof()) 
    {
      getline(in,STRING);

      found = STRING.find(".root");
      if (found!=std::string::npos)
      {
	m_L1TT->Add(STRING.c_str());   
	m_PIX->Add(STRING.c_str());   
	if (has_patt) m_PATT->Add(STRING.c_str());   
      }
    }

    in.close();
  }
  
  pm_clus_layer=&m_clus_layer;
  pm_clus_ladder=&m_clus_ladder;
  pm_clus_module=&m_clus_module;
  pm_clus_segment=&m_clus_segment;
  pm_clus_bottom=&m_clus_bottom;
  pm_clus_tp=&m_clus_tp;

  pm_stub_layer=&m_stub_layer;
  pm_stub_ladder=&m_stub_ladder;
  pm_stub_module=&m_stub_module;
  pm_stub_segment=&m_stub_segment;
  pm_stub_strip=&m_stub_strip;
  pm_stub_chip=&m_stub_chip;
  pm_stub_sw=&m_stub_sw;
  pm_stub_sw_c=&m_stub_sw_c;
  pm_stub_did=&m_stub_did;
  pm_stub_pxGEN=&m_stub_pxGEN;
  pm_stub_pyGEN=&m_stub_pyGEN;
  pm_stub_etaGEN=&m_stub_etaGEN;
  pm_stub_tp=&m_stub_tp;
  pm_stub_pdg=&m_stub_pdg;
  pm_stub_X0=&m_stub_X0;
  pm_stub_Y0=&m_stub_Y0;
  pm_stub_Z0=&m_stub_Z0;
  pm_stub_x=&m_stub_x;
  pm_stub_y=&m_stub_y;
  pm_stub_z=&m_stub_z;
  pm_stub_c=&m_stub_c;
  pm_clus_x=&m_clus_x;
  pm_clus_y=&m_clus_y;
  pm_clus_z=&m_clus_z;
  pm_stub_clust2=&m_stub_clust2;
  pm_stub_clust1=&m_stub_clust1;
  pm_clus_cw=&m_clus_cw;
  pm_clus_strip=&m_clus_strip;

  pm_pattlinks=&m_pattlinks;
  pm_pattsecid=&m_pattsecid;
  pm_pattmiss=&m_pattmiss;
  pm_pattsmult=&m_pattsmult;
  pm_pattid=&m_pattid;
    
  pm_tclinks=&m_tclinks;
  pm_tcsecid=&m_tcsecid;
  pm_tcpt=&m_tcpt;
  pm_tceta=&m_tceta;
  pm_tcphi=&m_tcphi;
  pm_tcz=&m_tcz;    
  pm_tcid=&m_tcid;
    
  pm_trklinks=&m_trklinks;
  pm_trksecid=&m_trksecid;
  pm_trkpt=&m_trkpt;
  pm_trketa=&m_trketa;
  pm_trkphi=&m_trkphi;
  pm_trkz=&m_trkz;

  m_PIX->SetBranchAddress("PIX_nPU",            &m_PU); 
  m_PIX->SetBranchAddress("PIX_n",              &m_pix); 


  m_L1TT->SetBranchAddress("L1Tkevt",            &m_evtid);
  m_L1TT->SetBranchAddress("L1TkCLUS_n",         &m_clus); 
  m_L1TT->SetBranchAddress("L1TkCLUS_tp",        &pm_clus_tp);
  m_L1TT->SetBranchAddress("L1TkCLUS_layer",     &pm_clus_layer);
  m_L1TT->SetBranchAddress("L1TkCLUS_ladder",    &pm_clus_ladder);
  m_L1TT->SetBranchAddress("L1TkCLUS_module",    &pm_clus_module);
  m_L1TT->SetBranchAddress("L1TkCLUS_seg",       &pm_clus_segment);
  m_L1TT->SetBranchAddress("L1TkCLUS_bottom",    &pm_clus_bottom);
  m_L1TT->SetBranchAddress("L1TkCLUS_strip",     &pm_clus_strip);
  //  m_L1TT->SetBranchAddress("L1TkCLUS_n",         &m_clus); 
  m_L1TT->SetBranchAddress("L1TkSTUB_n",         &m_stub);
  m_L1TT->SetBranchAddress("L1TkSTUB_layer",     &pm_stub_layer);
  m_L1TT->SetBranchAddress("L1TkSTUB_ladder",    &pm_stub_ladder);
  m_L1TT->SetBranchAddress("L1TkSTUB_module",    &pm_stub_module);
  m_L1TT->SetBranchAddress("L1TkSTUB_seg",       &pm_stub_segment);
  m_L1TT->SetBranchAddress("L1TkSTUB_strip",     &pm_stub_strip);
  m_L1TT->SetBranchAddress("L1TkSTUB_detid",     &pm_stub_did);
  m_L1TT->SetBranchAddress("L1TkSTUB_chip",      &pm_stub_chip);
  m_L1TT->SetBranchAddress("L1TkSTUB_cor",       &pm_stub_c);
  m_L1TT->SetBranchAddress("L1TkSTUB_deltas",    &pm_stub_sw);
  m_L1TT->SetBranchAddress("L1TkSTUB_deltasf",   &pm_stub_sw_c);
  m_L1TT->SetBranchAddress("L1TkSTUB_pxGEN",     &pm_stub_pxGEN);
  m_L1TT->SetBranchAddress("L1TkSTUB_pyGEN",     &pm_stub_pyGEN);
  m_L1TT->SetBranchAddress("L1TkSTUB_X0",        &pm_stub_X0);
  m_L1TT->SetBranchAddress("L1TkSTUB_Y0",        &pm_stub_Y0);
  m_L1TT->SetBranchAddress("L1TkSTUB_Z0",        &pm_stub_Z0);
  m_L1TT->SetBranchAddress("L1TkSTUB_x",         &pm_stub_x);
  m_L1TT->SetBranchAddress("L1TkSTUB_y",         &pm_stub_y);
  m_L1TT->SetBranchAddress("L1TkSTUB_z",         &pm_stub_z);
  m_L1TT->SetBranchAddress("L1TkSTUB_etaGEN",    &pm_stub_etaGEN);
  m_L1TT->SetBranchAddress("L1TkSTUB_tp",        &pm_stub_tp);
  m_L1TT->SetBranchAddress("L1TkSTUB_pdgID",     &pm_stub_pdg);
  m_L1TT->SetBranchAddress("L1TkSTUB_clust1",    &pm_stub_clust1);
  m_L1TT->SetBranchAddress("L1TkSTUB_clust2",    &pm_stub_clust2);
  //m_L1TT->SetBranchAddress("L1TkCLUS_x",         &pm_clus_x);
  //m_L1TT->SetBranchAddress("L1TkCLUS_y",         &pm_clus_y);
  //m_L1TT->SetBranchAddress("L1TkCLUS_z",         &pm_clus_z);
  m_L1TT->SetBranchAddress("L1TkCLUS_nstrip",    &pm_clus_cw);

  if (has_patt)
  {
    m_PATT->SetBranchAddress("L1PATT_n",           &nb_patterns);
    m_PATT->SetBranchAddress("L1PATT_links",       &pm_pattlinks);
    m_PATT->SetBranchAddress("L1PATT_secid",       &pm_pattsecid);
    m_PATT->SetBranchAddress("L1PATT_nmiss",       &pm_pattmiss);
    m_PATT->SetBranchAddress("L1PATT_pattid",      &pm_pattid);
      
    m_PATT->SetBranchAddress("L1TC_n",             &nb_tcs);
    m_PATT->SetBranchAddress("L1TC_links",         &pm_tclinks);
    m_PATT->SetBranchAddress("L1TC_secid",         &pm_tcsecid);
    m_PATT->SetBranchAddress("L1TC_pt",            &pm_tcpt);
    m_PATT->SetBranchAddress("L1TC_eta",           &pm_tceta);
    m_PATT->SetBranchAddress("L1TC_phi",           &pm_tcphi);
    m_PATT->SetBranchAddress("L1TC_z",             &pm_tcz);
    m_PATT->SetBranchAddress("L1TC_pattid",        &pm_tcid);
      
    m_PATT->SetBranchAddress("L1TRK_n",            &nb_tracks);
    m_PATT->SetBranchAddress("L1TRK_links",        &pm_trklinks);
    m_PATT->SetBranchAddress("L1TRK_secid",        &pm_trksecid);
    m_PATT->SetBranchAddress("L1TRK_pt",           &pm_trkpt);
    m_PATT->SetBranchAddress("L1TRK_eta",          &pm_trketa);
    m_PATT->SetBranchAddress("L1TRK_phi",          &pm_trkphi);
    m_PATT->SetBranchAddress("L1TRK_z",            &pm_trkz);
  }

  // Output file definition (see the header)

  m_outfile = new TFile(out.c_str(),"recreate");
 
  m_finaltree = new TTree("FullInfo","");

  m_finaltree->Branch("evt",          &evt); 
  m_finaltree->Branch("npu",          &m_npu); 

  if (m_dbg)
  {
  m_finaltree->Branch("n_pix_total", &n_pix_total); 
  m_finaltree->Branch("n_clus_total", &n_clus_total); 
  m_finaltree->Branch("n_clus_total_c", &n_clus_total_c); 
  m_finaltree->Branch("n_stub_total", &n_stub_total); 
  m_finaltree->Branch("n_stub_total_c", &n_stub_total_c); 
  m_finaltree->Branch("n_stub_inpat", &n_stub_patt); 
  m_finaltree->Branch("n_stub_intc",  &n_stub_tc); 
  m_finaltree->Branch("n_stub_intrk", &n_stub_trk); 
  m_finaltree->Branch("stub_x",       &stub_x); 
  m_finaltree->Branch("stub_y",       &stub_y); 
  m_finaltree->Branch("stub_z",       &stub_z);
  m_finaltree->Branch("stub_ldid",    &stub_did);
  m_finaltree->Branch("stub_layer",   &stub_layer); 
  m_finaltree->Branch("stub_ladder",  &stub_ladder);
  m_finaltree->Branch("stub_module",  &stub_module);
  m_finaltree->Branch("stub_segment", &stub_seg);
  m_finaltree->Branch("stub_chip",    &stub_chip);
  m_finaltree->Branch("stub_strip",   &stub_strip);
  m_finaltree->Branch("stub_sw",      &stub_sw);
  m_finaltree->Branch("stub_sw_off",  &stub_sw_c);
  m_finaltree->Branch("stub_sw_truth",&stub_sw_truth);
  m_finaltree->Branch("stub_cw1",     &stub_cw1);
  m_finaltree->Branch("stub_cw2",     &stub_cw2);
  m_finaltree->Branch("stub_cn1",     &stub_cn1);
  m_finaltree->Branch("stub_cn2",     &stub_cn2);
  m_finaltree->Branch("stub_tp",      &stub_tp);
  m_finaltree->Branch("stub_fe",      &stub_fe);
  m_finaltree->Branch("stub_ptGEN",   &stub_ptGEN);
  m_finaltree->Branch("stub_d0GEN",   &stub_d0GEN);
  m_finaltree->Branch("stub_r0GEN",   &stub_r0GEN);
  m_finaltree->Branch("stub_etaGEN",  &stub_etaGEN);
  m_finaltree->Branch("stub_z0GEN",   &stub_z0GEN);
  m_finaltree->Branch("stub_sec",     &stub_sec);
  m_finaltree->Branch("stub_insec",   &stub_insec);
  m_finaltree->Branch("stub_inpatt",  &stub_inpatt);
  m_finaltree->Branch("stub_intc",    &stub_intc);
  m_finaltree->Branch("stub_intrack", &stub_intrk);
  }

  m_finaltree->Branch("dtc_id",       &dtc_id); 
  m_finaltree->Branch("dtc_mult",     &dtc_mult); 
  m_finaltree->Branch("dtc_mult_FE",  &dtc_mult_FE); 

  m_finaltree->Branch("n_part",       &n_part); 
  m_finaltree->Branch("part_pdg",     &part_pdg); 
  m_finaltree->Branch("part_nsec",    &part_nsec); 
  m_finaltree->Branch("part_sec",     &part_sec);
  m_finaltree->Branch("part_nhits",   &part_nhits); 
  m_finaltree->Branch("part_nstubs",  &part_nstubs); 
  m_finaltree->Branch("part_nhits_fe",&part_nhits_fe); 
  m_finaltree->Branch("part_npatt",   &part_npatt); 
  m_finaltree->Branch("part_ntc",     &part_ntc); 
  m_finaltree->Branch("part_ntrk",    &part_ntrk); 
  m_finaltree->Branch("part_pt",      &part_pt); 
  m_finaltree->Branch("part_rho",     &part_rho);
  m_finaltree->Branch("part_d0",      &part_d0);
  m_finaltree->Branch("part_z0",      &part_z0);
  m_finaltree->Branch("part_eta",     &part_eta); 
  m_finaltree->Branch("part_phi",     &part_phi); 
  m_finaltree->Branch("part_dist",    &part_dist);  

  if (m_dbg)
  {
  m_finaltree->Branch("n_patt",       &n_patt); 
  m_finaltree->Branch("n_combs",      &n_patt_comb); 
  m_finaltree->Branch("patt_sec",     &patt_sec);
  m_finaltree->Branch("patt_id",      &patt_id);
  m_finaltree->Branch("patt_miss",    &patt_miss);
  m_finaltree->Branch("patt_parts",   &patt_parts); 
  m_finaltree->Branch("patt_stubs",   &patt_stubs);
  m_finaltree->Branch("patt_comb",    &patt_comb); 
  m_finaltree->Branch("patt_gcomb",   &patt_gcomb); 
  m_finaltree->Branch("patt_insec",   &patt_insec); 
  m_finaltree->Branch("patt_insec_full",  &patt_insec_full);
  m_finaltree->Branch("patt_rank",        &patt_rank);
  m_finaltree->Branch("patt_rank_full",   &patt_rank_full);
       
  m_finaltree->Branch("n_tc",         &n_tc);
  m_finaltree->Branch("tc_sec",       &tc_sec);
  m_finaltree->Branch("tc_id",        &tc_id);
  m_finaltree->Branch("tc_parts",     &tc_parts); 
  m_finaltree->Branch("tc_stubs",     &tc_stubs); 
  m_finaltree->Branch("tc_pt",        &tc_pt); 
  m_finaltree->Branch("tc_phi",       &tc_phi); 
  m_finaltree->Branch("tc_z",         &tc_z); 
  m_finaltree->Branch("tc_eta",       &tc_eta); 
  m_finaltree->Branch("tc_pt_t",      &tc_pt_t); 
  m_finaltree->Branch("tc_phi_t",     &tc_phi_t); 
  m_finaltree->Branch("tc_z_t",       &tc_z_t); 
  m_finaltree->Branch("tc_eta_t",     &tc_eta_t);
  m_finaltree->Branch("tc_PDG_t",     &tc_PDG_t);
  }

  m_finaltree->Branch("n_track",      &n_track);
  m_finaltree->Branch("trk_sec",      &trk_sec); 
  m_finaltree->Branch("trk_parts",    &trk_parts); 
  m_finaltree->Branch("trk_stubs",    &trk_stubs); 
  m_finaltree->Branch("trk_pt",       &trk_pt); 
  m_finaltree->Branch("trk_phi",      &trk_phi); 
  m_finaltree->Branch("trk_z",        &trk_z); 
  m_finaltree->Branch("trk_eta",      &trk_eta); 
  m_finaltree->Branch("trk_pt_t",     &trk_pt_t); 
  m_finaltree->Branch("trk_phi_t",    &trk_phi_t); 
  m_finaltree->Branch("trk_z_t",      &trk_z_t); 
  m_finaltree->Branch("trk_eta_t",    &trk_eta_t);
  m_finaltree->Branch("trk_PDG_t",    &trk_PDG_t);
}


/////////////////////////////////////////////////////////////////////////////////
//
// ==> track_eff::convert(std::string sectorfilename) 
//
// Here we retrieve info from the TKLayout CSV file containing the sector definition
//
// This file contains, for each sector, the ids of the modules contained in the sector
//
// The role of this method is to create the opposite, ie a vector containing, for every module the list of sectors belonging to it
//
/////////////////////////////////////////////////////////////////////////////////

bool track_eff::convert_cabling(std::string cablingfilename)
{
  int modid,lay,lad,mod,disk,type,dtcid;
    
  m_sec_mult = 0;
    
  std::vector<int> module;
    
  m_modules_to_DTC.clear();
  m_DTC_to_modules.clear();
    
  for (unsigned int i=0;i<230000;++i)
  {
    module.clear();
    m_modules_to_DTC.push_back(module);
  }

  for (unsigned int i=0;i<2000;++i)
  {
    module.clear();
    m_DTC_to_modules.push_back(module);
  }
    
  std::string STRING;
  std::ifstream in(cablingfilename.c_str());
  if (!in)
  {
    std::cout << "Please provide a valid csv sector filename" << std::endl;
    return false;
  }
    
  int npar = 0;

  while (!in.eof())
  {
    ++m_sec_mult;
    getline(in,STRING);
    std::istringstream ss(STRING);
    npar = 0;
    
    while (ss)
    {
      std::string s;
      if (!getline( ss, s, ',' )) break;
      
      ++npar;
      if (npar==1) modid = atoi(s.c_str());
      if (npar==2) dtcid = atoi(s.c_str());
    }
            
    std::bitset<32> detid = modid; // Le detid
            
    int rmodid; // Le modid que l'on utilise
            
    if (detid[25]) // barrel;
    {
      lay  = 8*detid[23]+4*detid[22]+2*detid[21]+detid[20]+4;
      type = 2*detid[19]+detid[18];
      
      if (type==3) // Pas tilté
      {
	lad  = 128*detid[17]+64*detid[16]+32*detid[15]+16*detid[14]+
	  8*detid[13]+4*detid[12]+2*detid[11]+detid[10]-1;
	mod  = 128*detid[9]+64*detid[8]+32*detid[7]+16*detid[6]+
	  8*detid[5]+4*detid[4]+2*detid[3]+detid[2]-1+limits[lay-5][type-1];
      }
      else // tilté
      {
	mod  = 128*detid[17]+64*detid[16]+32*detid[15]+16*detid[14]+
	  8*detid[13]+4*detid[12]+2*detid[11]+detid[10]-1+limits[lay-5][type-1];
	lad  = 128*detid[9]+64*detid[8]+32*detid[7]+16*detid[6]+
	  8*detid[5]+4*detid[4]+2*detid[3]+detid[2]-1;
      }
    }
    else // endcap
    {
      disk  = 8*detid[21]+4*detid[20]+2*detid[19]+detid[18];
      lay   = 10+disk+abs(2-(2*detid[24]+detid[23]))*7;
      lad   = 32*detid[17]+16*detid[16]+8*detid[15]+4*detid[14]+2*detid[13]+detid[12]-1;          
      mod   = 128*detid[9]+64*detid[8]+32*detid[7]+16*detid[6]+
	8*detid[5]+4*detid[4]+2*detid[3]+detid[2]-1;
    }
              
    rmodid = 10000*lay+100*lad+mod;
    
    module = m_modules_to_DTC.at(rmodid);
    module.push_back(dtcid);
    
    m_modules_to_DTC.at(rmodid) = module;
   
    module = m_DTC_to_modules.at(dtcid);
    module.push_back(rmodid);
    
    m_DTC_to_modules.at(dtcid) = module;
  }
    
  int n_mods=0, n_DTCs=0;
  
  for (unsigned int i=0;i<230000;++i)
  {
    if (m_modules_to_DTC.at(i).size()!=0) ++n_mods;
  }
  
  for (unsigned int i=0;i<2000;++i)
  {
    if (m_DTC_to_modules.at(i).size()!=0) ++n_DTCs;
  }
  
  std::cout << "Found " << n_mods << " modules despatched into " << n_DTCs << " DTCs" << endl;
    
  in.close();
        
  return true;
}

bool track_eff::convert_towers(std::string towerfilename)
{
  int modid,lay,lad,mod,disk,type;
    
  m_sec_mult = 0;
  
  std::vector<int> module;
  
  m_modules.clear();
    
  for (unsigned int i=0;i<230000;++i)
  {
    module.clear();
    m_modules.push_back(module);
  }
    
  std::string STRING;
  std::ifstream in(towerfilename.c_str());
  if (!in)
  {
    std::cout << "Please provide a valid csv sector filename" << std::endl;
    return false;
  }
  
  int npar = 0;
  
  while (!in.eof())
  {
    ++m_sec_mult;
    getline(in,STRING);
    
    if (m_sec_mult<1) continue;
        
    std::istringstream ss(STRING);
    npar = 0;
        
    while (ss)
    {
      std::string s;
      if (!getline( ss, s, ',' )) break;
      
      ++npar;
      if (npar<=2) continue;
            
      modid = atoi(s.c_str());
            
      std::bitset<32> detid = modid; // Le detid
            
      int rmodid; // Le modid que l'on utilise
            
      if (detid[25]) // barrel;
      {
	lay  = 8*detid[23]+4*detid[22]+2*detid[21]+detid[20]+4;
	type = 2*detid[19]+detid[18];
        
	if (type==3) // Pas tilté
	{
	  lad  = 128*detid[17]+64*detid[16]+32*detid[15]+16*detid[14]+
	    8*detid[13]+4*detid[12]+2*detid[11]+detid[10]-1;
	  mod  = 128*detid[9]+64*detid[8]+32*detid[7]+16*detid[6]+
	    8*detid[5]+4*detid[4]+2*detid[3]+detid[2]-1+limits[lay-5][type-1];
	}
	else // tilté
	{
	  mod  = 128*detid[17]+64*detid[16]+32*detid[15]+16*detid[14]+
	    8*detid[13]+4*detid[12]+2*detid[11]+detid[10]-1+limits[lay-5][type-1];
	  lad  = 128*detid[9]+64*detid[8]+32*detid[7]+16*detid[6]+
	    8*detid[5]+4*detid[4]+2*detid[3]+detid[2]-1;
	}
      }
      else // endcap
      {
	disk  = 8*detid[21]+4*detid[20]+2*detid[19]+detid[18];
	lay   = 10+disk+abs(2-(2*detid[24]+detid[23]))*7;
	lad   = 32*detid[17]+16*detid[16]+8*detid[15]+4*detid[14]+2*detid[13]+detid[12]-1;
	
	mod  = 128*detid[9]+64*detid[8]+32*detid[7]+16*detid[6]+
	  8*detid[5]+4*detid[4]+2*detid[3]+detid[2]-1;
      }            
      
      rmodid = 10000*lay+100*lad+mod;
            
      module = m_modules.at(rmodid);
      module.push_back(m_sec_mult-1);
      
      m_modules.at(rmodid) = module;
    }
  }
    
  in.close();
  
  m_sec_mult -= 1;
  
  return true;
}


void track_eff::reset() 
{
  m_primaries.clear(); 

  m_npu=0;

  n_pix_total=0; 
  n_clus_total=0; 
  n_stub_total=0; 
  n_stub_total_c=0; 
  n_clus_total_c=0; 
  n_stub_patt=0; 
  n_stub_tc=0;
  n_stub_trk=0;
 
  stub_x->clear();
  stub_y->clear();
  stub_z->clear();
  stub_layer->clear();   
  stub_ladder->clear();  
  stub_module->clear();
  stub_seg->clear();
  stub_did->clear();
  stub_chip->clear();    
  stub_strip->clear();    
  stub_sw->clear();
  stub_sw_c->clear();
  stub_sw_truth->clear();
  stub_cw1->clear();
  stub_cw2->clear();
  stub_cn1->clear();
  stub_cn2->clear();
  stub_tp->clear();
  stub_fe->clear();

  stub_insec->clear(); 
  stub_inpatt->clear(); 
  stub_intc->clear(); 
  stub_intrk->clear();  
  stub_sec->clear();  

  stub_ptGEN->clear();
  stub_d0GEN->clear();  
  stub_r0GEN->clear();  
  stub_etaGEN->clear();   
  stub_z0GEN->clear();

  n_part=0; 
  part_pdg->clear();   
  part_nsec->clear();   
  part_sec->clear();   
  part_dist->clear();   
  part_nhits->clear();   
  part_nhits_fe->clear();   
  part_nstubs->clear();   
  part_npatt->clear();   
  part_ntc->clear();  
  part_ntrk->clear();   
  part_pt->clear();   
  part_rho->clear();  
  part_d0->clear();  
  part_z0->clear();  
  part_eta->clear();    
  part_phi->clear();   

  n_patt=0; 
  n_patt_comb=0;
  patt_sec->clear();
  patt_id->clear();
  patt_miss->clear();   
  patt_parts->clear();   
  patt_stubs->clear(); 
  patt_comb->clear(); 
  patt_gcomb->clear(); 
  patt_insec->clear(); 
  patt_insec_full->clear();
  patt_rank->clear();
  patt_rank_full->clear();

  n_tc=0;
  tc_parts->clear();   
  tc_stubs->clear();
  tc_sec->clear();
  tc_id->clear();
  tc_pt->clear();   
  tc_phi->clear();   
  tc_z->clear();   
  tc_eta->clear(); 
  tc_pt_t->clear();   
  tc_phi_t->clear();   
  tc_z_t->clear();   
  tc_eta_t->clear();    
  tc_PDG_t->clear();

  n_track=0;
  trk_parts->clear();   
  trk_stubs->clear();
  trk_sec->clear();   
  trk_pt->clear();   
  trk_phi->clear();   
  trk_z->clear();   
  trk_eta->clear();   
  trk_pt_t->clear();   
  trk_phi_t->clear();   
  trk_z_t->clear();   
  trk_eta_t->clear();    
  trk_PDG_t->clear();

  dtc_id->clear();
  dtc_mult->clear();
  dtc_mult_FE->clear();
}
