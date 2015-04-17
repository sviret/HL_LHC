#ifndef FILTER_H
#define FILTER_H

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <cmath>
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
// Base class for tracker coordinates filtering 
//
// filename    : the name and directory of the input ROOT file to filter
// secfilename : the name and directory of the csv file containing the sectors definition
// outfile     : the name of the output ROOT file containing the filtered events
//           
// secid       : the sector number
//
// Info about the sector definition:
//
//  http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCData
//
///////////////////////////////////

using namespace std;

class filter
{
 public:

  filter(std::string filename, std::string secfilename, 
	 std::string outfile, int secid);

  void   do_filter(int secid);    
  void   initTuple(std::string test,std::string out);
  bool   convert(std::string sectorfilename); 
    
 private:

  TFile  *m_outfile;
  TChain *m_L1TT;

  TTree  *m_efftree;

  std::vector< std::vector<int> >   m_modules; 

  int m_sec_mult;

  int    m_layer;
  int    m_module;
  int    m_ladder;
  int    m_row;
  int    m_column;
  float  m_x;
  float  m_y;
  float  m_z;
};

#endif

