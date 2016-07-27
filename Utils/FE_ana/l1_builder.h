#ifndef L1_BUILDER_H
#define L1_BUILDER_H


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
//  Date: 26/08/2015
//
///////////////////////////////////



class l1_builder
{
 public:

    l1_builder(std::string filenameRAW, std::string outfile, std::string sectorfile,
               int npatt, int rate, int layer, int ladder, int module, int npblock, int error);
    void  initVars();
    void  initTuple(std::string inRAW,std::string out);
    void  get_stores(int nevts);

    void  fill_RAW_block(std::vector<int> digis,bool spars,int BXid);
    void  fill_RAW_header_CBC(int L1id);
    void  fill_RAW_header_MPA(int L1id);
    
    void  fill_CONC_RAW_block(std::vector<int> digis,bool MPA,int BXid,bool unsparsified);
    void  fill_CONC_RAW_header(int L1id);

    bool  convert(std::string sectorfilename);
    
    int   badbit();

 private:

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

    bool unsp;

    int m_error_FE;
    float m_error_prop;

    int m_error_rdm;
    int m_error_hdr;
    int m_error_siz;

    int m_pix;
    int m_npu;
    std::vector<int>    m_pix_layer;
    std::vector<int>    m_pix_ladder;
    std::vector<int>    m_pix_module;
    std::vector<int>    m_pix_row;
    std::vector<int>    m_pix_col;

    std::vector<int>    *pm_pix_layer;
    std::vector<int>    *pm_pix_ladder;
    std::vector<int>    *pm_pix_module;
    std::vector<int>    *pm_pix_row;
    std::vector<int>    *pm_pix_col;
    
    std::multimap< int, std::vector<int> > m_conc_raw;  //
    std::multimap< int, std::vector<int> > m_chip_raw;  //
    std::multimap< int, std::vector<int> >::const_iterator m_iter;
    std::multimap< int, std::vector<int> >::const_iterator m_iter2;
    std::multimap< int, std::vector<int> >::const_iterator m_iter3;

    std::vector<std::multimap< int, std::vector<int> > > m_data_raw;

    int m_lay;  // the layer number (5/6/7 for the MPA, 8/9/10 for the CBC)
    int m_lad;
    int m_mod;

    int m_MPA_L1_delay;
    int m_CBC_L1_delay;
    int m_CIC_L1_delay;
    
    int m_MPA_FIFO_depth;
    int m_CIC_FIFO_depth;
    
    int m_MPA_FIFO_size;
    int m_CIC_FIFO_size;
    
    std::multimap< int, std::vector<int> > m_raw_FIFO;
    std::multimap< int, std::vector<int> > m_chip_FIFOs;
    std::multimap< int, std::vector<int> > m_words;
    
    std::multimap< int, std::vector<int> > m_c_raw_FIFO;
    std::multimap< int, std::vector<int> > m_c_chip_FIFOs;
    std::multimap< int, std::vector<int> > m_c_words;
    
    std::vector<int>  *m_raw_chip_fifo;
    std::vector<int>  *m_raw_chip_bx;
    
    std::vector<int>  *m_raw_conc_fifo;
    std::vector<int>  *m_raw_conc_bx;
    
    float m_raw_chip_slope;
    float m_raw_chip_slope_err;
    float m_raw_chip_int;
    float m_raw_chip_int_err;

    int m_rate; // the input L1 rate, in kHz

    std::vector< std::vector<int> >   m_modules;

    std::vector< int >   m_chips;
    std::vector< int >   m_concs;

    TFile *m_outfile;  // The output file
    TTree *m_raw_tree; //
    TTree *m_raw_summary; //
 
    int m_raw_cic;
    int m_raw_bx;
    int m_raw_chip;
    int m_raw_size;
    int m_raw_mbits;
    int m_raw_FIFO_FULL;
    int m_raw_FIFO_SIZE;
    int m_npblock;
    int m_raw_lay;
    int m_raw_lad;
    int m_raw_mod;
    int m_raw_np;
    int m_raw_ns;
    
    std::vector<int>   *m_raw_data;
    
    ofstream FE_L1_IN;
    ofstream FE_L1_OUT;
    ofstream FE_L1_OUT_E;
    ofstream CIC_L1_OUT;
};

#endif

