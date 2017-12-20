#include <iostream>
#include <fstream>
#include <iomanip>

// Internal includes

#include "rates.h"
//#include "TC_win.h"
#include "track_eff.h"
#include "evtbuilder.h"
#include "efficiencies.h"
#include "windows.h"
#include "jobparams.h"
#include "TROOT.h"
//#include "latency.h"

using namespace std;

///////////////////////////////////
//
//
// Base code for the AM analysis tool
//
// This code is extensively documented on the AM-tracking tutorial page:
//
// http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto620
//
//
//  Author: viret@in2p3_dot_fr
//  Date       : 23/05/2013
//  Maj. update: 13/11/2013
//
///////////////////////////////////

int main(int argc, char** argv) {

  // Necessary lines to make branches containing vectors
  gROOT->ProcessLine(".L Loader.C+");

  // read jobParams
  jobparams params(argc,argv);

  // Depending on the option chosen, process the information


  // Option 1: just do the rate calculation
  // Documented in part 3.2.2 of the mentioned tutorial page
  if (params.option()=="rates")
  {
    rates* my_rates = new rates(params.inputfile(),params.outfile(),params.dbg());
    delete my_rates;    
  }
  
  // Option 2: basic stub efficiencies calculation
  // Documented in part 3.2.3 of the official tutorial page
  if (params.option()=="stub_eff")
  {
      efficiencies* my_effs = new efficiencies(params.inputfile(),params.outfile(),params.type());
      delete my_effs;
  }


  // Option 3: pattern reco efficiency from official stubs 
  // Documented in part 6.2.2 of the official tutorial page
  if (params.option()=="L1track_eff")
  {
    track_eff* my_test = new track_eff(params.testfile(),params.inputfile(),
				       params.outfile(),params.nevt(),params.ptmin(),
				       params.qmax(),params.dbg());

    delete my_test;
  }

  //    // Option 4: pattern reco efficiency from official stubs
  //    // Documented in part 6.2.2 of the official tutorial page
  //    if (params.option()=="latency")
  //    {
  //      latency* my_test = new latency(params.testfile(),params.inputfile(),
     //                                    params.outfile(),params.nevt(),12,10,15,10000,1);
  //                                    params.outfile(),params.nevt(),12,10,15,10000,params.ptmin(),1);
  //     
  //    delete my_test;
  //}
    
    
    // Option 5: stub windows calculation
    //
  
    if (params.option()=="stub_win")
    {
      windows* my_effs = new windows(params.inputfile(),params.outfile(),
                                     params.sptminmu(),params.sptmaxmu(),params.sptminele(),params.sptmaxele(),
                                     params.qmax(),params.maxlosses());
      delete my_effs;
    }
    


//  if (params.option()=="get_windows")
//  {
//    TC_win* my_test = new TC_win(params.testfile(),params.outfile(),params.nevt(),params.ptmin(),
//				 params.qmax(),params.dbg());
//
//    delete my_test;
 // }
    // Option 6: generate a serie of events into a root file
    
    if (params.option()=="losses_TRG_CONC") // TRG block for CIC
    {
        std::cout << "Evaluating the trigger data losses for the whole tracker: " << std::endl;
        
        evtbuilder* my_pgen = new evtbuilder(params.inputfile(),params.inputfileQCD(),params.outfile(),
                                             params.nevt(),100,-1,-1,-1,
                                             params.testfile(),false,true,0,
                                             true,params.qmax(),0);
        delete my_pgen;
    }

  return 0;
}
