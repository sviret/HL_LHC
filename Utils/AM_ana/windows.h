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
// Base class for stub window tuning definition 
//
// This code computes the efficiencies per module or the efficiencies per layer
// both for official and private stub producer
//
// For infos about the SW definition, have a look at the following document:
//
// https://www.dropbox.com/s/sc1ifhfbcv5pts3/Stub%20windows%20tuning-%20a%20tutorial.pdf?dl=0
//
//
//  Author: s.viret@ipnl_dot_in2p3_dot_fr
//  Date: 20/12/2017
//
///////////////////////////////////



class windows
{
 public:

  windows(std::string file_r, std::string file_e, float pminm, float pmaxm, float pmine, float pmaxe, float prop, float CIC_lim);

  void  get_rates();  // The main method
  void  get_losses();  // The main method
  void  get_effs(int ptype);  // The main method
  void  initVars();
  void  reset();
  void  initTuple(std::string in_r,std::string in_e);
  void  get_result(int ltype, int ptype);
  void  print_result(int type);
    
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
  // https://github.com/sviret/HL_LHC/blob/Tools_for_10_0_0/Extractors/RecoExtractor/interface/StubExtractor.h
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

  float barrel_w[6][13][16][5];
  float disk_w[5][15][16][5];
  
  float   m_ovflow_b[6][40][40][3][2];
  float   m_ovflow_d[5][15][40][3][2];
     
  float   m_rate_b[6][40][40][3][2];
  float   m_rate_d[5][15][40][3][2];
 
  float   m_loss_b[6][40][40][3][2];
  float   m_loss_d[5][15][40][3][2];
    
  float   loss_b[6][13][40][2];
  float   loss_d[5][15][40][2];
    
  float barrel_tune[6][13][2];
  float disk_tune[5][15][2];
    
};

#endif

