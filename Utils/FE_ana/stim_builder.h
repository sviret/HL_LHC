#ifndef STIM_BUILDER_H
#define STIM_BUILDER_H


#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <set>
#include <iomanip>   
#include <random>  

#include <stdio.h>    
#include <stdlib.h>     /* srand, rand */
#include <time.h>    

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include <fstream>
#include <string>
#include <sstream> 
#include <bitset>         // std::bitset

using namespace std;

///////////////////////////////////
//
//
// Base class for pattern generation for the concentrator chip
//
// Input infos are :
//
// filename : the name and directory of the input ROOT file containing the STUB information
// outfile  : the name of the output ROOT file containing the pattern info (a text file with binary info is 
//            also produced) 
// npatt    : the number of patterns you want to generate 
//
// Info about the code:
//
//  ...TBD...
// 
//
//  Author: viret@in2p3_dot_fr
//  Date: 12/11/2013
//
///////////////////////////////////



class stim_builder
{
 public:

  stim_builder(std::string filenameRAW, std::string filenameTRG, std::string outfile,
               int npatt, int layer, int ladder, int module, int tower, std::string sector,
               int L1prop, int TRGsize);

  void  initVars();
  void  initTuple(std::string inRAW,std::string inTRG,std::string out);
  void  get_stores(int nevts,bool conc);
  bool  convert(std::string sectorfilename);

  void  fill_TRG_block(std::vector<std::vector<int> > stubs, bool ps, bool conc, int BXid);

  void  fill_CONC_TRG_header(int BXid, int MPA);
    
    void  ana_pix(int lay,int lad,int mod, std::vector<int> digits);
  //void  do_stub(int lay,int lad,int mod);
    void do_stub(std::vector<int> stubs);
    
    
    std::bitset<3> convert_bend_PS(int layer,int ring, float bend);
    std::bitset<4> convert_bend_2S(int layer,int ring, float bend);
    
 private:

  TChain *L1TT;     // The trees containing the input data
  TChain *PIX;      // The trees containing the input data
  TChain *MC;      // The trees containing the input data
    
  // Coding conventions for barrel and endcap module IDs
  
  // We define a barrel ID and an endcap ID as follows:

  // B_id = (layer-1)*10000 + ladder*100 + module (0 to 57523)

  // E_id = (disk-1)*10000 + ladder*100 + module (0 to 131377)

  // Disks 0 to 6  (towards positive Z)  (0 to  4 for 5D)
  // Disks 7 to 13 (towards negative Z)  (7 to 11 for 5D)

  // Here are the parameters needed from the data
  // Details on these might be found on
  //
  // https://github.com/sviret/HL_LHC/blob/master/Extractors/RecoExtractor/interface/L1TrackTrigger_analysis.h
  //

    int m_ntp;
    std::vector<float>  m_part_px;
    std::vector<float>  m_part_py;
    std::vector<float>  m_part_eta;
    std::vector<float>  m_part_x;
    std::vector<float>  m_part_y;
    std::vector<float>  m_part_z;
    
    std::vector<float>  *pm_part_px;
    std::vector<float>  *pm_part_py;
    std::vector<float>  *pm_part_eta;
    std::vector<float>  *pm_part_x;
    std::vector<float>  *pm_part_y;
    std::vector<float>  *pm_part_z;
    
    int m_pix;
    int m_npu;
    std::vector<int>    m_pix_layer;
    std::vector<int>    m_pix_ladder;
    std::vector<int>    m_pix_module;
    std::vector<int>    m_pix_row;
    std::vector<int>    m_pix_col;
    std::vector<int>    m_pix_ncol;
    std::vector<int>    m_pix_bot;
    std::vector<float>  m_pix_ch;
    std::vector<float>  m_pix_x;
    std::vector<float>  m_pix_y;
    std::vector<float>  m_pix_z;
    
    std::vector<int>    *pm_pix_layer;
    std::vector<int>    *pm_pix_ladder;
    std::vector<int>    *pm_pix_module;
    std::vector<int>    *pm_pix_row;
    std::vector<int>    *pm_pix_col;
    std::vector<int>    *pm_pix_ncol;
    std::vector<int>    *pm_pix_bot;
    std::vector<float>  *pm_pix_x;
    std::vector<float>  *pm_pix_y;
    std::vector<float>  *pm_pix_z;
    
    int m_clus;
    std::vector<int>   m_clus_layer;
    std::vector<int>   m_clus_ladder;
    std::vector<int>   m_clus_module;
    std::vector<int>   m_clus_nrows;
    std::vector<int>   m_clus_bot;
    std::vector<int>   m_clus_tp;
    std::vector<std::vector<int> >   m_clus_pix;
    std::vector<std::vector<int> >   m_clus_mult;
    std::vector<int>    m_clus_nseg;
    
    std::vector<int>   *pm_clus_layer;
    std::vector<int>   *pm_clus_ladder;
    std::vector<int>   *pm_clus_module;
    std::vector<int>   *pm_clus_nrows;
    std::vector<int>   *pm_clus_bot;
    std::vector<int>   *pm_clus_tp;
    std::vector<std::vector<int> >   *pm_clus_pix;
    std::vector<std::vector<int> >   *pm_clus_mult;
    std::vector<int>   *pm_clus_nseg;
    
    
    int m_stub;
    std::vector<int>    m_stub_layer;
    std::vector<int>    m_stub_ladder;
    std::vector<int>    m_stub_module;
    std::vector<float>  m_stub_pxGEN;  // px generated of stub i (in GeV/c)
    std::vector<float>  m_stub_pyGEN;  // py generated of stub i (in GeV/c)
    std::vector<float>  m_stub_etaGEN; // eta generated of stub i (in GeV/c)
    std::vector<float>  m_stub_X0;     // x origin of particle inducing stub i (in cm)
    std::vector<float>  m_stub_Y0;     // y origin of particle inducing stub i (in cm)
    std::vector<float>  m_stub_Z0;     // z origin of particle inducing stub i (in cm)
    std::vector<float>  m_stub_PHI0;   // phi origin of particle inducing stub i (in rad)
    std::vector<float>  m_stub_pt;
    std::vector<int>    m_stub_tp;
    std::vector<float>  m_stub_deltas;
    std::vector<float>  m_stub_strip;
    std::vector<float>  m_stub_z;
    std::vector<int>    m_stub_seg;
    std::vector<int>    m_stub_chip;
    std::vector<int>    m_stub_clust1;
    std::vector<int>    m_stub_clust2;
    
    
    std::vector<int>    *pm_stub_layer;
    std::vector<int>    *pm_stub_ladder;
    std::vector<int>    *pm_stub_module;
    std::vector<float>  *pm_stub_pxGEN;  // px generated of stub i (in GeV/c)
    std::vector<float>  *pm_stub_pyGEN;  // py generated of stub i (in GeV/c)
    std::vector<float>  *pm_stub_etaGEN; // eta generated of stub i (in GeV/c)
    std::vector<float>  *pm_stub_X0;     // x origin of particle inducing stub i (in cm)
    std::vector<float>  *pm_stub_Y0;     // y origin of particle inducing stub i (in cm)
    std::vector<float>  *pm_stub_Z0;     // z origin of particle inducing stub i (in cm)
    std::vector<float>  *pm_stub_pt;
    std::vector<int>    *pm_stub_tp;
    std::vector<float>  *pm_stub_deltas;
    std::vector<float>  *pm_stub_strip;
    std::vector<float>  *pm_stub_z;
    std::vector<int>    *pm_stub_seg;
    std::vector<int>    *pm_stub_chip;
    std::vector<int>    *pm_stub_clust1;
    std::vector<int>    *pm_stub_clust2;


  std::multimap< int, std::vector<int> > m_chip_trig;  //
  std::multimap< int, std::vector<int> > m_conc_trig;  //

  std::multimap< int, std::vector<int> > m_chip_raw;  //
  std::multimap< int, std::vector<int> > m_conc_raw;  //

  std::multimap< int, std::vector<int> >::const_iterator m_iter; 
  std::multimap< int, std::vector<int> >::const_iterator m_iter2; 

  std::vector<std::multimap< int, std::vector<int> > > m_data_trig;
  std::vector<std::multimap< int, std::vector<int> > > m_data_raw;


  int m_lay;  // the layer number (5/6/7 for the MPA, 8/9/10 for the CBC)
  int m_lad; 
  int m_mod; 

  std::vector< std::vector<int> >   m_modules; 

  std::vector< int >   m_chips; 
  std::vector< int >   m_concs;

  TFile *m_outfile;  // The output file
  TTree *m_tri_tree; // 
  TTree *m_raw_tree; //
  TTree *m_raw_summary; //
 
  std::vector< std::vector<int> >   trig_sequence; 

  double m_L1prop;

  int m_CICsize;
    int m_tower;
    int m_TRGsize;
  int m_PHYsize;

  int bend_bit_MPA;
  int bend_bit_CBC;

  int m_tri_bx;
  int m_tri_size;
  int m_tri_size_anders;
  int m_tri_chip;
  int m_tri_nstubs;
  int m_tri_nstubs_s;
  int m_tri_nstubs_g;
  int m_tri_nstubs_gs;

  int m_raw_bx;
  int m_raw_mbits;
  int m_raw_chip;
  int m_raw_size;
  int m_raw_FIFO_FULL;
  int m_tri_lay;
  int m_tri_lad;
  int m_tri_mod;
  int m_raw_lay;
  int m_raw_lad;
  int m_raw_mod;
  int m_raw_np;
  int m_raw_ns;

  std::vector<int>   *m_tri_data;
  std::vector<int>   *m_raw_data;

    std::vector<int>   m_evt_pix;
    std::vector<int>   m_evt_clu;
    std::vector<int>   m_evt_stu;
    std::vector<int>   m_evt_tp;

    ofstream FE_TRG_IN;
    ofstream FE_TRG_OUT;
    ofstream CIC_TRG_OUT;

};

#endif

