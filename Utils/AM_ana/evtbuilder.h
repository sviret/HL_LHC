#ifndef EVTBUILDER_H
#define EVTBUILDER_H


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



class evtbuilder
{
 public:

  evtbuilder(std::string filenameRAW, std::string filenameTRG, std::string outfile, 
	     int npatt, int rate, int layer, int ladder, int module, std::string sector, 
	     bool RAW, bool TRIG, int npblock, bool conc,float L1prop,bool dbg);

  void  initVars();
  void  initTuple(std::string inRAW,std::string inTRG,std::string out);
  void  get_stores(int nevts,bool conc);
  bool  convert(std::string sectorfilename); 


  void  fill_TRG_block(std::vector< std::vector<int> > stubs, bool spars, bool conc, int BXid,int CIClim);


  void  fill_RAW_block(std::vector<int> digis,bool spars,int BXid);
  void  fill_CONC_RAW_block(std::vector<int> digis,bool spars,int BXid);

  void  fill_RAW_header_CBC(int L1id);
  void  fill_RAW_header_MPA(int L1id);

  void  fill_CONC_TRG_header(int BXid, int MPA);
  void  fill_CONC_RAW_header(int L1id);


 private:

  TChain *L1TT;     // The trees containing the input data
  TChain *PIX;      // The trees containing the input data

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


  int m_pix;
  int m_npu;
  std::vector<int>    m_pix_layer;
  std::vector<int>    m_pix_ladder;
  std::vector<int>    m_pix_module;
  std::vector<int>    m_pix_row;
  std::vector<int>    m_pix_col;
  std::vector<int>    m_pix_bot;
  std::vector<int>    m_pix_nrow;
  std::vector<int>    m_pix_ncol;
  std::vector<float>  m_pix_z;

  std::vector<int>    *pm_pix_layer;
  std::vector<int>    *pm_pix_ladder;
  std::vector<int>    *pm_pix_module;
  std::vector<int>    *pm_pix_row;
  std::vector<int>    *pm_pix_col;
  std::vector<int>    *pm_pix_bot;
  std::vector<int>    *pm_pix_nrow;
  std::vector<int>    *pm_pix_ncol;
  std::vector<float>  *pm_pix_z;

  int m_stub;
  std::vector<int>    m_stub_layer;
  std::vector<int>    m_stub_ladder;
  std::vector<int>    m_stub_module;
  std::vector<float>  m_stub_pt;
  std::vector<float>  m_stub_deltas;
  std::vector<float>  m_stub_strip;
  std::vector<int>    m_stub_seg;
  std::vector<int>    m_stub_chip;
  std::vector<float>  m_stub_X0;
  std::vector<float>  m_stub_Y0;
  std::vector<float>  m_stub_pxGEN;
  std::vector<float>  m_stub_pyGEN;
  std::vector<float>  m_stub_z;

  std::vector<int>    *pm_stub_layer;
  std::vector<int>    *pm_stub_ladder;
  std::vector<int>    *pm_stub_module;
  std::vector<float>  *pm_stub_pt;
  std::vector<float>  *pm_stub_deltas;
  std::vector<float>  *pm_stub_strip;
  std::vector<int>    *pm_stub_seg;
  std::vector<int>    *pm_stub_chip;
  std::vector<float>  *pm_stub_X0;
  std::vector<float>  *pm_stub_Y0;
  std::vector<float>  *pm_stub_pxGEN;
  std::vector<float>  *pm_stub_pyGEN;
  std::vector<float>  *pm_stub_z;

  std::multimap< int, std::vector<int> > m_chip_trig;  //

  std::multimap< int, std::vector<int> > m_chip_trig1;  //
  std::multimap< int, std::vector<int> > m_chip_trig2;  //

  std::multimap< int, std::vector<int> > m_chip_raw;  //
  std::multimap< int, std::vector<int> >::const_iterator m_iter; 
  std::multimap< int, std::vector<int> >::const_iterator m_iter2; 

  std::multimap< int, std::vector<int> > m_raw_FIFO;  //

  std::vector<std::multimap< int, std::vector<int> > > m_data_trig;
  std::vector<std::multimap< int, std::vector<int> > > m_data_raw;

  std::multimap< int, std::vector<int> > m_chip_FIFOs;

  std::vector<int>  *m_raw_chip_fifo;
  std::vector<int>  *m_raw_chip_bx;

  float m_raw_chip_slope;
  float m_raw_chip_slope_err;
  float m_raw_chip_int;
  float m_raw_chip_int_err;

  bool m_dbg;

  int m_rate; // the input L1 rate, in kHz
  int m_lay;  // the layer number (5/6/7 for the MPA, 8/9/10 for the CBC)
  int m_lad;
  int m_mod;

  std::vector< std::vector<int> >   m_modules; 

  std::vector< int >   m_chips; 
  std::vector< int >   m_concs;

  TFile *m_outfile;  // The output file
  TTree *m_tri_tree; // 
  TTree *m_tri_summary; // 
  TTree *m_raw_tree; //
  TTree *m_raw_summary; //
 
  std::vector< std::vector<int> >   trig_sequence; 

  double m_L1prop;

  int m_CICsize;

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

  int m_tri_nstubs_f[40][3][2];

  float m_tri_b_nstubs[6][40][40][3][2];
  float m_los_b_nstubs[6][40][40][3][2];
  //int   m_lp_b_nstubs[6][40][40][3][2];
  float   m_ovflow_b_nstubs[6][40][40][3][2];

  float m_tri_d_nstubs[5][15][40][3][2];
  float m_los_d_nstubs[5][15][40][3][2];
  //int   m_lp_d_nstubs[5][15][40][3][2];
  float   m_ovflow_d_nstubs[5][15][40][3][2];


  float m_sw[40];

  int m_raw_bx;
  int m_raw_chip;
  int m_raw_size;
  int m_raw_mbits;
  int m_raw_FIFO_FULL;
  int m_raw_FIFO_SIZE;
  int m_npblock;
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

  bool m_write_raw;
  bool m_write_trg;
  bool m_write_out;

};

#endif

