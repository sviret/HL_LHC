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
// Base class for stub windows calculation
//
// This code computes the stubs windows, based on an input file created without 
// windows (and preferably PU0 and particle gun)
//
//
// Input infos are :
//
// filename : the name and directory of the input ROOT file containing the STUB information
// pmin     : minimum pT of the particle to be considered
// pmax     : maximum pT of the particle to be considered
// prop     : proportion of stubs you want to keep between pmin and pmax (bet. 0 and 1)
// ptype    : the pdg ID of the particle type you want to test 
//
// Info about the code:
//
//  TBD
//
//  Author: s.viret@ipnl_dot_in2p3_dot_fr
//
///////////////////////////////////



class windows
{
    public:

    windows(std::string filename, float pmin, float pmax, float prop, int ptype);
    
    ~windows(){}

    void  get_windows();  // The main method
    void  initVars();
    void  reset();
    void  initTuple(std::string in);
    void  print_result();
    
    private:

    float m_pmin;
    float m_pmax;
    float m_prop;
    int   m_ptype;

    TChain *L1TT;      // The trees
  
    // Here are the parameters needed from the data
    // Details on these might be found on
    //
    // https://github.com/sviret/HL_LHC/blob/master/Extractors/RecoExtractor/interface/L1TrackTrigger_analysis.h
    //

    int m_stub;
    
    std::vector<float>  m_stub_pxGEN;
    std::vector<float>  m_stub_pyGEN;
    std::vector<float>  m_stub_etaGEN;
    std::vector<int>    m_stub_layer;
    std::vector<int>    m_stub_module;
    std::vector<int>    m_stub_ladder;
    std::vector<int>    m_stub_type;
    std::vector<float>  m_stub_deltas;
    std::vector<int>    m_stub_tp;
    std::vector<int>    m_stub_pdg;
    
    std::vector<float>  *pm_stub_pxGEN;
    std::vector<float>  *pm_stub_pyGEN;
    std::vector<float>  *pm_stub_etaGEN;
    std::vector<int>    *pm_stub_layer;
    std::vector<int>    *pm_stub_module;
    std::vector<int>    *pm_stub_ladder;
    std::vector<int>    *pm_stub_type;
    std::vector<float>  *pm_stub_deltas;
    std::vector<int>    *pm_stub_tp;
    std::vector<int>    *pm_stub_pdg;

    float barrel_w[6][13][40];
    float disk_w[5][15][40];    
};

#endif

