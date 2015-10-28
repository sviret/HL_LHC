// Class for the data filtering
// For more info, look at the header file

#include "filter.h"

filter::filter(std::string filename, std::string secfilename, 
	       std::string outfile, int secid)
{  
  filter::initTuple(filename,outfile);

  if (!filter::convert(secfilename)) return; // Don't go further if there is no sector file

  filter::do_filter(secid); // Launch the filter loop
}

/////////////////////////////////////////////////////////////////////////////////
//
// ==> filter::do_filter(int secid)
//
// Main method, where the filtering is made
//
/////////////////////////////////////////////////////////////////////////////////

void filter::do_filter(int secid)
{
  const int m_nsec = m_sec_mult; // How many sectors are in the file

  int ndat = m_L1TT->GetEntries(); // How many events will we test

  cout << "Starting a filtering loop over " << ndat << " coordinates..." << endl;
  cout << "... using " << m_nsec << " trigger sectors..." << endl;
  cout << "Looking coords in sector " << secid << "..." << endl;

  bool is_there;
  int modid;
  int layer;

  int ladder;
  int module;

  // Loop over the events
 
  for (int i=0;i<ndat;++i)
  {    
    m_L1TT->GetEntry(i);

    if (i%1000000==0) 
      cout << "Processed " << i << "/" << ndat << endl;
 
    layer  = m_layer; 
    ladder = m_ladder-1; 
    module = m_module; 
    //
    if(layer<=7)
      module=module/2;
      
    modid  = layer*10000+ladder*100+module; 

    if (m_modules.at(modid).size()<=1) continue;

    is_there=false;

    for (unsigned int kk=1;kk<m_modules.at(modid).size();++kk) // In which sector the module is
    {
      if (m_modules.at(modid).at(kk)==secid) is_there=true;  
    }

    if (!is_there) continue;

    m_efftree->Fill(); // If yes fill the skimmed tree

  }

  m_outfile->Write();
  m_outfile->Close();
}


/////////////////////////////////////////////////////////////////////////////////
//
// ==> filter::initTuple(std::string test,std::string out)
//
// This method opens and creates the differe rootuples involved
//
/////////////////////////////////////////////////////////////////////////////////


void filter::initTuple(std::string test,std::string out)
{
  m_L1TT   = new TChain("StripCoords"); 

  // Input data file 

  std::size_t found = test.find(".root");

  // Case 1, it's a root file
  if (found!=std::string::npos)
  {
    m_L1TT->Add(test.c_str());
  }
  else // This is a list provided into a text file
  {
    std::string STRING;
    std::ifstream in(test.c_str());
    if (!in) 
    {
      std::cout << "Please provide a valid data filename list" << std::endl; 
      return;
    }    
  
    while (!in.eof()) 
    {
      getline(in,STRING);

      found = STRING.find(".root");
      if (found!=std::string::npos) m_L1TT->Add(STRING.c_str());   
    }

    in.close();
  }

  m_L1TT->SetBranchAddress("layer",        &m_layer);
  m_L1TT->SetBranchAddress("module",       &m_module);
  m_L1TT->SetBranchAddress("ladder",       &m_ladder);
  m_L1TT->SetBranchAddress("strip",        &m_row);
  m_L1TT->SetBranchAddress("segment",      &m_column);
  m_L1TT->SetBranchAddress("x",            &m_x);
  m_L1TT->SetBranchAddress("y",            &m_y);
  m_L1TT->SetBranchAddress("z",            &m_z);


  // Output file definition (see the header)

  m_outfile = new TFile(out.c_str(),"recreate");

  m_efftree = new TTree("StripCoords","");

  m_efftree->Branch("layer",     &m_layer);
  m_efftree->Branch("module",    &m_module);
  m_efftree->Branch("ladder",    &m_ladder);
  m_efftree->Branch("strip",     &m_row);
  m_efftree->Branch("segment",   &m_column);
  m_efftree->Branch("x",         &m_x);
  m_efftree->Branch("y",         &m_y);
  m_efftree->Branch("z",         &m_z);
}



/////////////////////////////////////////////////////////////////////////////////
//
// ==> filter::convert(std::string sectorfilename) 
//
// Here we retrieve info from the TKLayout CSV file containing the sector definition
//
// This file contains, for each sector, the ids of the modules contained in the sector
//
// The role of this method is to create the opposite, ie a vector containing, for every module the list of sectors belonging to it
//
/////////////////////////////////////////////////////////////////////////////////

bool filter::convert(std::string sectorfilename) 
{
  int n_rods[6] = {16,24,34,48,62,76};

  int modid,lay,lad,mod;

  m_sec_mult = 0;

  std::vector<int> module;

  m_modules.clear();

  for (int i=0;i<230000;++i)
  {
    module.clear();
    module.push_back(-1);
    m_modules.push_back(module);
  }

  std::string STRING;
  std::ifstream in(sectorfilename.c_str());
  if (!in) 
  {
    std::cout << "Please provide a valid csv sector filename" << std::endl; 
    return false;
  }    
  
  int npar = 0;

  while (!in.eof()) 
  {
    ++m_sec_mult;
    getline(in,STRING);

    if (m_sec_mult<2) continue;

    std::istringstream ss(STRING);
    npar = 0;
    while (ss)
    {
      std::string s;
      if (!getline( ss, s, ',' )) break;

      ++npar;
      if (npar<=2) continue;

      modid = atoi(s.c_str());
      
      lay   = int(modid/10000); 
      modid-= 10000*lay;
      lad   = int(modid/100); 
      modid-= 100*lad;
      mod   = modid; 

      ///////
      // This hack is temporary and is due to a numbering problem in the TkLayout tool
      if (lay<=10) lad = (lad+n_rods[lay-5]/4)%(n_rods[lay-5]);
     // if (lay<=7)  mod = mod/2;
      ///////

      modid = 10000*lay+100*lad+mod;

      m_modules.at(modid).push_back(m_sec_mult-2);
    }
  }

  in.close();

  m_sec_mult -= 1;

  return true;
}

