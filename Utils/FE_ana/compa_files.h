#ifndef COMPA_FILES_H
#define COMPA_FILES_H


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

#include <fstream>
#include <string>
#include <sstream> 
#include <bitset>         // std::bitset

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cmath>
#include <vector>

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



class compa_files
{
 public:

  compa_files(std::string verilog, std::string simu, std::string outfile, int type, int MPA);

  void compa_TRG();
  void compa_L1();

  void do_ana(std::vector<int> sec_ver,std::vector<int> sec_sim, int n);

 private:

  std::string m_out; 
  std::string m_insim; 
  std::string m_indat; 

  std::ofstream m_outbinary; // txt file containing the output sequences
  std::ifstream fver;
  std::ifstream fsim;

  int m_isMPA;

};

#endif

