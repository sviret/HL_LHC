#ifndef TRACK_EFF_H
#define TRACK_EFF_H

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <cmath>
#include <bitset>
#include <stdio.h>  
#include <stdlib.h> 

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include <fstream>
#include <string>
#include <sstream> 

///////////////////////////////////
//
//
// Base class for L1 track reconstruction efficiency using the official chain
//
// The role of this class is to check that all particles are contained in at least one 
// trigger sector, to count the stub multiplicities per sector, and finally to 
// compute the pattern bank and track fit efficiencies
//
// This class is using the official object and is starting from rootuple created
// in step 6 of the tutorial
//
// filename    : the name and directory of the input ROOT file containing the data to analyze
// secfilename : the name and directory of the input csv file  containing the sectors definition
// outfile     : the name of the output ROOT file containing the efficiency results 
// ptmin       : the minimal pt of the tracks to be tested
// d0max       : the maximal |d0| of the tracks to be tested
//           
// nevt        : the number of particles to test
// dbg         : full output (1) or only partial one (0)
//
// Info about the code:
//
//  http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto620 (STEP 6)
//
//  Author: viret@in2p3_dot_fr
//  Date: 25/02/2014
//  Maj. rev: 08/03/2016
//
///////////////////////////////////

using namespace std;

class track_eff
{
 public:

  track_eff(std::string filename, std::string secfilename, 
	    std::string outfile, int nevt, float pt_min, 
	    float d0_min, bool dbg);

  void   do_test(int nevt);    
  void   initTuple(std::string test,std::string out);
  bool   convert_towers(std::string towerfilename);
  bool   convert_cabling(std::string cablingfilename);
  void   reset();
    
 private:

  bool m_dbg;
  bool has_patt;
  float m_pt_min;
  float m_d0_min;
  bool m_tilted;
  int limits[6][3];
  
  TFile  *m_outfile;
  TChain *m_L1TT;
  TTree  *m_finaltree;
  TChain *m_PATT;
  TChain *m_PIX;

  // List of stubs from primaries (one vector of stub ids per prim track) 
  // for a given event 

  std::vector< std::vector<int> >   m_primaries;

  // Coding conventions for barrel and endcap module IDs
  
  // We define a barrel ID and an endcap ID as follows:

  // B_id = (layer-1)*10000 + ladder*100 + module (0 to 57523)

  // E_id = (disk-1)*10000 + ladder*100 + module (0 to 131377)

  // Disks 0 to 6  (towards positive Z)
  // Disks 7 to 13 (towards negative Z)

  int m_sec_mult;


  std::vector< std::vector<int> >   m_modules;         // Table providing for a module the list of trigger towers it belongs to
  std::vector< std::vector<int> >   m_modules_to_DTC;  // Table providing for a module the list of DTC it belongs to (should have only one entry)
  std::vector< std::vector<int> >   m_DTC_to_modules;  // Table providing for a DTC the list of modules belonging to it

  int m_evtid;

  int m_npu,m_PU;

  int m_pix;
  int m_clus;
  int m_stub;
  std::vector<int>   m_stub_layer;
  std::vector<int>   m_stub_ladder;
  std::vector<int>   m_stub_module;
  std::vector<int>   m_stub_segment;
  std::vector<int>   m_stub_did;
  std::vector<float> m_stub_strip;
  std::vector<float> m_stub_sw;
  std::vector<float> m_stub_sw_c;
  std::vector<int>   m_stub_tp;
  std::vector<int>   m_stub_pdg;
  std::vector<int>   m_stub_chip;
  std::vector<float> m_stub_pxGEN;
  std::vector<float> m_stub_pyGEN;
  std::vector<float> m_stub_etaGEN;
  std::vector<float> m_stub_X0;
  std::vector<float> m_stub_Y0;
  std::vector<float> m_stub_Z0;
  std::vector<float> m_stub_x;
  std::vector<float> m_stub_y;
  std::vector<float> m_stub_z;
  std::vector<float> m_clus_x;
  std::vector<float> m_clus_y;
  std::vector<float> m_clus_z;
  std::vector<float> m_stub_c;
  std::vector<int>   m_stub_clust1;
  std::vector<int>   m_stub_clust2;
  std::vector<int>   m_clus_cw;
  std::vector<float>   m_clus_strip;

  std::vector<int>   m_clus_layer;
  std::vector<int>   m_clus_ladder;
  std::vector<int>   m_clus_module;
  std::vector<int>   m_clus_segment;
  std::vector<int>   m_clus_bottom;
  std::vector<int>   m_clus_tp;

  std::vector<int>   *pm_stub_layer;
  std::vector<int>   *pm_stub_ladder;
  std::vector<int>   *pm_stub_module;
  std::vector<int>   *pm_stub_segment;
  std::vector<int>   *pm_stub_did;
  std::vector<int>   *pm_stub_chip;
  std::vector<float> *pm_stub_strip;
  std::vector<float> *pm_stub_sw;
  std::vector<float> *pm_stub_sw_c;
  std::vector<int>   *pm_stub_tp;
  std::vector<int>   *pm_stub_pdg;
  std::vector<float> *pm_stub_pxGEN;
  std::vector<float> *pm_stub_pyGEN;
  std::vector<float> *pm_stub_etaGEN;
  std::vector<float> *pm_stub_X0;
  std::vector<float> *pm_stub_Y0;
  std::vector<float> *pm_stub_Z0;
  std::vector<float> *pm_stub_x;
  std::vector<float> *pm_stub_y;
  std::vector<float> *pm_stub_z;
  std::vector<float> *pm_clus_x;
  std::vector<float> *pm_clus_y;
  std::vector<float> *pm_clus_z;
  std::vector<float> *pm_stub_c;
  std::vector<int>   *pm_stub_clust1;
  std::vector<int>   *pm_stub_clust2;
  std::vector<int>   *pm_clus_cw;
  std::vector<float>   *pm_clus_strip;

  std::vector<int>   *pm_clus_layer;
  std::vector<int>   *pm_clus_ladder;
  std::vector<int>   *pm_clus_module;
  std::vector<int>   *pm_clus_segment;
  std::vector<int>   *pm_clus_bottom;
  std::vector<int>   *pm_clus_tp;

  // The output tree has one entry per event

  int   evt;        // Event number (for PU event, where there is more than 1 primary)
  int   pdg;        // The pdg ID of the particle
  int   nsec;       // The number of sectors containing at least 5 stubs of the prim. track
  float eta;        // The eta of the prim. track
  float phi;        // The phi of the prim. track
  float pt_f;       // The fitted pT of the prim. track
  float eta_f;      // The fitted eta of the prim. track
  float phi_f;      // The fitted phi of the prim. track
  float d0;         // The origin radius
  float z0;         // The origin z
  float z0_f;       // The origin z

  int nstubs;
  int nstubs_f;

  int n_pix_total;
  int n_clus_total;
  int n_stub_total_c;
  int n_clus_total_c;
  int n_stub_total;   // The total number of stubs in the event
  int n_stub_patt;    // The total number of stubs contained in matched patterns in the event
  int n_stub_tc;      // The total number of stubs contained in tcs in the event
  int n_stub_trk;     // The total number of stubs contained in reco tracks in the event

  std::vector<float>   *stub_x;      // x coordinates of ALL the stubs |
  std::vector<float>   *stub_y;      // y coordinates of ALL the stubs |-> Bottom cluster
  std::vector<float>   *stub_z;      // z coordinates of ALL the stubs |
  std::vector<int>     *stub_layer;  // layer number of ALL the stubs
  std::vector<int>     *stub_ladder; // ladder number of ALL the stubs
  std::vector<int>     *stub_module; // module number of ALL the stubs
  std::vector<int>     *stub_seg;    // segment number of ALL the stubs
  std::vector<float>   *stub_strip;  // strip number of ALL the stubs
  std::vector<int>     *stub_did;    // detid number of ALL the stubs
  std::vector<int>     *stub_chip;
  std::vector<float>   *stub_sw;     // stub width at FE level of ALL the stubs
  std::vector<float>   *stub_sw_c;   // stub width at offline level of ALL the stubs
  std::vector<float>   *stub_sw_truth;// stub width at offline level of ALL the stubs
  std::vector<int>     *stub_cw1;
  std::vector<int>     *stub_cw2;
  std::vector<int>     *stub_cn1;
  std::vector<int>     *stub_cn2;
  std::vector<int>     *stub_tp;     // tp index of ALL the stubs (in part_*** vectors of this tree!!!!)
  std::vector<int>     *stub_fe;     // 0 passing CIC, 1 passing CBC/MPA, 2 not passing CBC/MPA 
  std::vector<float>   *stub_ptGEN;  // PT gen of the particle inducing the stub (-1 if bad)
  std::vector<float>   *stub_d0GEN;  // d0 gen of the particle inducing the stub (-1 if bad)
  std::vector<float>   *stub_r0GEN;  // r0 gen of the particle inducing the stub (-1 if bad)
  std::vector<float>   *stub_etaGEN; // eta gen of the particle inducing the stub (-1 if bad)
  std::vector<float>   *stub_z0GEN;  // z0 gen of the particle inducing the stub (-1 if bad)
  std::vector<int>     *stub_insec;  // is the stub in trigger towers (1,2,...) or not (0)?
  std::vector<int>     *stub_inpatt; // is the stub in a pattern (1) or not (0)?
  std::vector<int>     *stub_intc;   // is the stub in a tc (1) or not (0)?
  std::vector<int>     *stub_intrk;  // is the stub in a track (1) or not (0)?
  std::vector< std::vector<int> >     *stub_sec;    // the list of sectors in which the stub is

  int n_part;                        // The total number of particles inducing at least one stub in the event

  std::vector<int>     *part_pdg;    // PDG id of the particles
  std::vector<int>     *part_nsec;   // In how many trigger towers this particle hit more than 4 different layers/disks?
  std::vector<int>     *part_nhits_fe;  // 
  std::vector<int>     *part_nstubs;  // 
  std::vector<int>     *part_nhits;  // How many different layers/disks are hit by the particle?
  std::vector<int>     *part_npatt;  // How many patterns contains more than 4 stubs of the particle (in 4 different layers/disks)?
  std::vector<int>     *part_ntc;    // How many tcs of size N contains more than N-1 stubs of the particle?
  std::vector<int>     *part_ntrk;   // How many tracks of size N contains more than N-1 stubs of the particle?
  std::vector<float>   *part_pt;     // pt of the particles
  std::vector<float>   *part_rho;    // rho0 of the particles
  std::vector<float>   *part_d0;     // IP of the particles
  std::vector<float>   *part_z0;     // z0 of the particles
  std::vector<float>   *part_eta;    // eta of the particles 
  std::vector<float>   *part_phi;    // phi of the particles
  std::vector<float>   *part_dist;   // Isolation (min r/phi dist of the particle wrt another high pt (>3 GeV/c) primary)
  std::vector< std::vector<int> >     *part_sec;    // the list of sectors in which the particle is

  int n_patt;                        // The total number of patterns matched in the event
  int n_patt_comb;                   // The total number of combinations from the matched patterns
  std::vector<int>                  *patt_sec;        // Sector id of all the patterns
  std::vector<int>                  *patt_id;
  std::vector<int>                  *patt_miss;       // Number of superstrip missing on the pattern (max 1)
  std::vector<int>                  *patt_insec;      // Total number of matched roads in the tower containing the pattern
  std::vector<int>                  *patt_insec_full; // Total number of complete matched roads (no missing strip)
  std::vector<int>                  *patt_comb;       // Number of combinations (if running PCA without TCB)
  std::vector<int>                  *patt_gcomb;      //
  std::vector< std::vector<int> >   *patt_parts; // tp index of ALL the particles contained in the pattern (in part_*** vectors of this tree!!!!)
  std::vector< std::vector<int> >   *patt_stubs; // index of ALL the stubs contained in the pattern (in stub_*** vectors of this tree!!!!)
  std::vector<int>                  *patt_rank;
  std::vector<int>                  *patt_rank_full;


  int n_tc;                        // The total number of TCs in the event
  std::vector<int>                  *tc_sec;    // Sector id of all the TCs
  std::vector<int>                  *tc_id;
  std::vector< std::vector<int> >   *tc_parts;  // tp index of ALL the particles contained in the TC (in part_*** vectors of this tree!!!!)
  std::vector< std::vector<int> >   *tc_stubs;  // index of ALL the stubs contained in the TC (in stub_*** vectors of this tree!!!!)
  std::vector<float>                *tc_pt;     // pt of the tc
  std::vector<float>                *tc_eta;    // eta of the tc
  std::vector<float>                *tc_z;      // z0 of the tc
  std::vector<float>                *tc_phi;    // phi of the tc
  std::vector<float>                *tc_pt_t;     // pt of the tc (TRUTH)
  std::vector<float>                *tc_eta_t;    // eta of the tc (TRUTH)
  std::vector<float>                *tc_z_t;      // z0 of the tc (TRUTH)
  std::vector<float>                *tc_phi_t;    // phi of the tc (TRUTH)
  std::vector<int>                  *tc_PDG_t;    // PDG id of the tc (TRUTH) (0 if unmatched) 

  int n_track;                        // The total number of L1 tracks matched in the event
  std::vector<int>                  *trk_sec;    // Sector id of all the tracks
  std::vector< std::vector<int> >   *trk_parts;  // tp index of ALL the particles contained in the pattern (in part_*** vectors of this tree!!!!)
  std::vector< std::vector<int> >   *trk_stubs;  // index of ALL the stubs contained in the track (in stub_*** vectors of this tree!!!!) 
  std::vector<float>                *trk_pt;     // pt of the track
  std::vector<float>                *trk_eta;    // eta of the track
  std::vector<float>                *trk_z;      // z0 of the track
  std::vector<float>                *trk_phi;    // phi of the track
  std::vector<float>                *trk_pt_t;     // pt of the track
  std::vector<float>                *trk_eta_t;    // eta of the track
  std::vector<float>                *trk_z_t;      // z0 of the track
  std::vector<float>                *trk_phi_t;    // phi of the track 
  std::vector<int>                  *trk_PDG_t;    // PDG id of the track 

  // Informations contained in the pattern file

  int nb_patterns; // Number of patterns in event evtID
  int nb_tcs; // Number of patterns in event evtID

//  int event_id;    // Event ID

  std::vector<int>                *dtc_id;
  std::vector<int>                *dtc_mult;
  std::vector<int>                *dtc_mult_FE;


  // Info id dbg==false (CMSSW fast mode)

  std::vector< std::vector<int> > m_links; // Links to the stub ids of the pattern 
  std::vector<int>                m_secid; // Link to the sector ids of the pattern 

//  std::vector< std::vector<int> > *pm_links;
//  std::vector<int>                *pm_secid;

  // Informations contained in the pattern tree of official file

  int nb_tracks; // Number of L1 tracks in event

  std::vector< std::vector<int> > m_pattlinks; // Links to the stub ids of the pattern 
  std::vector<int>                m_pattsecid; // Link to the sector ids of the pattern
  std::vector<int>                m_pattid;
  std::vector<int>                m_pattmiss;
  std::vector< std::vector<int> > m_pattsmult;
  std::vector< std::vector<int> > m_tclinks;  // Links to the stub ids of the track 
  std::vector<int>                m_tcsecid;  // Link to the sector ids of the track
  std::vector<int>                m_tcid;
  std::vector<float>              m_tcpt;
  std::vector<float>              m_tceta; 
  std::vector<float>              m_tcz; 
  std::vector<float>              m_tcphi; 
  std::vector< std::vector<int> > m_trklinks;  // Links to the stub ids of the track 
  std::vector<int>                m_trksecid;  // Link to the sector ids of the track
  std::vector<float>              m_trkpt; 
  std::vector<float>              m_trketa; 
  std::vector<float>              m_trkz; 
  std::vector<float>              m_trkphi; 

  std::vector< std::vector<int> > *pm_pattlinks; 
  std::vector<int>                *pm_pattsecid;
  std::vector<int>                *pm_pattid;
  std::vector<int>                *pm_pattmiss;
  std::vector< std::vector<int> > *pm_pattsmult;
  std::vector< std::vector<int> > *pm_trklinks;
  std::vector<int>                *pm_trksecid; 
  std::vector<int>                *pm_tcid;
  std::vector<float>              *pm_trkpt;
  std::vector<float>              *pm_trketa; 
  std::vector<float>              *pm_trkz; 
  std::vector<float>              *pm_trkphi; 
  std::vector< std::vector<int> > *pm_tclinks;
  std::vector<int>                *pm_tcsecid; 
  std::vector<float>              *pm_tcpt; 
  std::vector<float>              *pm_tceta; 
  std::vector<float>              *pm_tcz; 
  std::vector<float>              *pm_tcphi;

};

#endif

