#include <iostream>
#include <fstream>
#include <iomanip>

// Internal includes

#include "patterngen.h"
#include "evtbuilder.h"
#include "stim_builder.h"
#include "l1_builder.h"
#include "compa_files.h"
#include "jobparams.h"
#include "TROOT.h"

using namespace std;

///////////////////////////////////
//
//
// Base code for the FE data analysis tool
//
// This code is extensively documented on the following page:
//
// http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.FETuto
//
//
//  Author: viret@in2p3_dot_fr
//  Date       : 24/03/2013
//  Last rev.  : 25/08/2014
//
///////////////////////////////////

int main(int argc, char** argv) {

  // Necessary lines to make branches containing vectors
  gROOT->ProcessLine(".L Loader.C+");

  // read jobParams
  jobparams params(argc,argv);

  // Depending on the option chosen, process the information


  // Option 1: generate a serie of patterns for the concentrator in a text file
  if (params.option()=="pattgen_CONC")
  {
    patterngen* my_pgen = new patterngen(params.inputfile(),params.outfile(),params.nevt(),params.rate(),0,0,0);
    delete my_pgen;
  }
  

  // Option 2: generate a serie of events for the MPA chip in a text file
  if (params.option()=="pattgen_MPA")
  {
    patterngen* my_pgen = new patterngen(params.inputfile(),params.outfile(),params.nevt(),0,params.lay(),params.lad(),params.mod());
    delete my_pgen;
  }

  // Option 3: generate a serie of events into a root file
  

  if (params.option()=="evtbuild_TRG_CONC")
  {
    evtbuilder* my_pgen = new evtbuilder(params.inputfile(),params.inputfileTRG(),params.outfile(),params.nevt(),params.rate(),-1,params.testfile(),false,true,params.l1size(),params.mpabend(),params.cbcbend(),true,params.prop(),params.cicsize());
    delete my_pgen;
  }

  if (params.option()=="evtbuild_L1_CONC")
  {
    evtbuilder* my_pgen = new evtbuilder(params.inputfile(),params.inputfileTRG(),params.outfile(),params.nevt(),params.rate(),-1,params.testfile(),true,false,params.l1size(),params.mpabend(),params.cbcbend(),true,0,8);
    delete my_pgen;
  }


  if (params.option()=="evtbuild_TRG_FE")
  {
    evtbuilder* my_pgen = new evtbuilder(params.inputfile(),params.inputfileTRG(),params.outfile(),params.nevt(),params.rate(),-1,params.testfile(),false,true,params.l1size(),params.mpabend(),params.cbcbend(),false,params.prop(),8);
    delete my_pgen;
  }

  if (params.option()=="evtbuild_L1_FE")
  {
    evtbuilder* my_pgen = new evtbuilder(params.inputfile(),params.inputfileTRG(),params.outfile(),params.nevt(),params.rate(),-1,params.testfile(),true,false,params.l1size(),params.mpabend(),params.cbcbend(),false,0,8);
    delete my_pgen;
  }

  if (params.option()=="build_stim_TRG")
  {
    stim_builder* my_stim = new stim_builder(params.inputfile(),params.inputfileTRG(),params.outfile(),params.nevt(),params.lay(),params.lad(),params.mod(),-1,params.testfile(),params.prop(),params.cicsize());
    delete my_stim;
  }
    
  if (params.option()=="build_stim_RAW")
  {
    l1_builder* my_stim = new l1_builder(params.inputfile(),params.outfile(),params.testfile(),params.nevt(),params.rate(),params.lay(),params.lad(),params.mod(),params.cicsize(),params.l1size());
    delete my_stim;
  }
    
  if (params.option()=="compa_out_TRG")
  {
    compa_files* my_comp = new compa_files(params.inputfile(),params.inputfileTRG(),params.outfile(),0,0);
    delete my_comp;
  }

  if (params.option()=="compa_out_L1_MPA")
  {
    compa_files* my_comp = new compa_files(params.inputfile(),params.inputfileTRG(),params.outfile(),1,1);
    delete my_comp;
  }

  if (params.option()=="compa_out_L1_CBC")
  {
    compa_files* my_comp = new compa_files(params.inputfile(),params.inputfileTRG(),params.outfile(),1,0);
    delete my_comp;
  }
  return 0;
}
