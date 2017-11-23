#ifndef WINDOWS_H
#define WINDOWS_H


#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <cmath>



#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include <fstream>
#include <string>

using namespace std;

///////////////////////////////////
//
//
// Base class for stub efficiencies calculation (in stubs/module/BX)
//
// This code computes the efficiencies per module or the efficiencies per layer
// both for official and private stub producer
//
// For infos about the efficiency definition, have a look at the following presentation:
//
// https://indico.cern.ch/getFile.py/access?contribId=0&resId=0&materialId=slides&confId=263068
//
// Input infos are :
//
// filename : the name and directory of the input ROOT file containing the STUB information
// outfile  : the name of the output ROOT file containing the rates 
// ptype    : the pdg ID of the particle type you want to test 
//
// Info about the code:
//
//  http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto (STEP III)
//
//  Author: viret@in2p3_dot_fr
//  Date: 23/05/2013
//
///////////////////////////////////



class windows
{
 public:

  windows(std::string file_r, std::string file_e, float pmin, float pmax, float prop, float CIC_lim, int ptype);

  void  get_rates();  // The main method
  void  get_losses();  // The main method
  void  get_effs();  // The main method
  void  initVars();
  void  reset();
  void  initTuple(std::string in_r,std::string in_e);
  void  print_result();
  
 private:
  
  float m_pmin;
  float m_pmax;
  float m_prop;
  int   m_ptype;
  float m_lim;

  TChain *L1TT;      // The trees
  TChain *Losses;      // The trees
  

  TFile *m_outfile;  // The output file
  TTree *m_ratetree; // The tree containing the rate information

  // Here are the parameters needed from the data
  // Details on these might be found on
  //
  // https://github.com/sviret/HL_LHC/blob/master/Extractors/RecoExtractor/interface/L1TrackTrigger_analysis.h
  //
  
  int n_e, n_r;

  int m_stub;

  std::vector<float>  *m_stub_pxGEN;
  std::vector<float>  *m_stub_pyGEN;
  std::vector<float>  *m_stub_etaGEN;
  std::vector<int>    *m_stub_layer;
  std::vector<int>    *m_stub_module;
  std::vector<int>    *m_stub_ladder;
  std::vector<int>    *m_stub_type;
  std::vector<float>  *m_stub_deltas;
  std::vector<int>    *m_stub_tp;
  std::vector<int>    *m_stub_pdg;

  float barrel_w[6][13][16][3];
  float disk_w[5][15][16][3];   
  
  float   m_ovflow_b[6][40][40][3][2];
  float   m_ovflow_d[5][15][40][3][2];
     
  float   m_rate_b[6][40][40][3][2];
  float   m_rate_d[5][15][40][3][2];

  float   loss_b[6][13][40][2];
  float   loss_d[5][15][40][2];
};

#endif

