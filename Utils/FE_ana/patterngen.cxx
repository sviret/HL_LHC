// Base class for simple pattern generation
// 
// Output is stored in binary text files corresponding to the front-end data format 

// Meaning of the parameters
//
//
// filename : the input data files, we use official stubs from extracted ROOT files
//            see how to produce that in part 3.1 of this page:
//            http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto620
// 
// outfile  : name of the output text file containing binary info
// npatt    : how many BXs should be produced (will be fitted correctly afterwards)
// rate     : input L1A rate (for the L1 raw block writing


#include "patterngen.h"

// Main constructor

patterngen::patterngen(std::string filename, std::string outfile,
		       int npatt, int rate, int lay, int lad, int mod)
{
  m_rate = rate;

  patterngen::initVars();
  
  if (npatt == 0) //
  {
    if (rate!=0)
    {
      patterngen::initTuple(filename,outfile,0);

      std::cout << "...Reading output file " << filename 
		<< " and putting info in " << outfile << std::endl;
      
      patterngen::get_CONC_output(filename);
    }
  }
  else // Create patterns or events
  {
    patterngen::initTuple(filename,outfile,1);

    if (rate!=0)  // Concentrator
    {
      std::cout << "...Creating " << npatt 
		<< " patterns with an L1 rate of " 
		<< m_rate << " kHz" << std::endl;

      patterngen::get_CONC_input(npatt);
    }
    else // MPA
    {      
      std::cout << "...Creating " << npatt 
		<< " events for the MPA " << std::endl;

      patterngen::get_MPA_input(npatt,lay,lad,mod);
    }
  }
}


//////////////////////////////////////////////
//
// Macro creating the input data stream to serve as input 
// for the MPA verilog model
//
//////////////////////////////////////////////

void patterngen::get_MPA_input(int nevt,int lay, int lad, int mod)
{
  int ladder,module,strip;
  int seg;

  ladder= -1;
  module= -1; 

  int idx=-1;
  int its=-1;
  int itp=-1;

  int modid;
  std::vector<int> m_digi_list;
  std::vector<int> m_mod_list;

  float ptGEN;
  float d0GEN;
  bool isGOOD;

  for (int j=0;j<nevt;++j)
    //for (int j=0;j<1;++j)
  { 
    L1TT->GetEntry(j); 
    PIX->GetEntry(j);
    MC->GetEntry(j); 

    m_pix_idx.clear();
    m_mod_list.clear();

    std::cout << std::endl;
    std::cout << "##################################################" << std::endl;
    std::cout << "#" << std::endl;
    std::cout << "# Event  " << j << std::endl;
    std::cout << "# NPU  = " << m_npu << std::endl;


    // First we loop over the pixels

    for (int i=0;i<m_pix;++i)
    {
      if (m_pix_layer[i]!=lay) continue; // Only interested into one layer
      
      ladder= m_pix_ladder[i];
      if (ladder!=lad && lad != -1) continue;

      module= static_cast<int>((m_pix_module[i]-1)/2); 
      if (module!=mod && mod != -1) continue;


      modid = ladder*100 + module;

      m_iter = m_pix_idx.find(modid);

      m_digi_list.clear();

      //      cout << lay << " / " << ladder << " / " << module << " / " 
      //	   << m_pix_row[i] << " / " << m_pix_col[i] << " / " << m_pix_ch[i] << endl;

      if (m_iter == m_pix_idx.end()) // New mod, add it
      {
	m_digi_list.push_back(i);
	m_pix_idx.insert(std::make_pair(modid,m_digi_list));
	m_mod_list.push_back(modid);
      }
      else
      {
	m_digi_list = m_iter->second;
	m_digi_list.push_back(i);
	m_pix_idx.erase(m_iter->first);
	m_pix_idx.insert(std::make_pair(modid,m_digi_list));
      }
    }

    // We now have a list of the modules containing pixels, along 
    // with the pixels contained into them
    //
    // m_pix_idx contains, for each touched module in the layer
    // a pair composed of the module ID number and a vector of int (pix indices)
    //
    // m_mod_list is just the list of modules

    if (m_mod_list.size()==0) continue;
    
    for (unsigned int i=0;i<m_mod_list.size();++i) // Loop over the list of touched modules
    {
      // We first got the entry in m_pix_idx

      m_iter = m_pix_idx.find(m_mod_list.at(i));

      m_digi_list = m_iter->second; // The pixel digi list
      modid       = m_iter->first;  // The module ID

      ladder      = modid/100;
      module      = 2*(modid-100*ladder)+1; 

      std::cout << std::endl;
      std::cout << "__________________________________________________" << std::endl;
      std::cout << "Layer/Disk  " << lay << std::endl;	 
      std::cout << "Ladder/Ring " << ladder << std::endl;	 
      std::cout << "Module      " << module << std::endl;	 
     
      // Link the stub/cluster info
      patterngen::ana_pix(lay,ladder,modid-100*ladder,m_digi_list); 

      // Start the printing loop
      for (unsigned int k=0;k<m_digi_list.size();++k) // Loop over all digis
      {       
	its=-1;
	itp=-1;
	idx   = m_digi_list.at(k); // Digi index in the original rootuple
	seg   = m_pix_col[idx]; 
	strip = m_pix_row[idx]; 
       
	// Got it or not in the flagged content from ana_pix
	for (unsigned int kk=0;kk<m_evt_pix.size();++kk)
	{
	  if (idx!=m_evt_pix[kk]) continue;
	  its = m_evt_stu[kk];
	  itp = m_evt_tp[kk];

	  //	  if (its!=-1)
	  //	    cout << kk << " / " << its << " / " << itp << endl;
	}

	isGOOD=false;

	// Look if the pix is induced by a primary particle
	if (itp!=-1)
	{
	  ptGEN = sqrt(m_part_px[itp]*m_part_px[itp]+m_part_py[itp]*m_part_py[itp]);
	  d0GEN = sqrt(m_part_x[itp]*m_part_x[itp]+m_part_y[itp]*m_part_y[itp]); 

	  if (ptGEN>2 && d0GEN<0.2) isGOOD=true;
	}

	// Start to print
	(m_pix_module[idx]%2==1)
	  ? std::cout << "pixeldigi " << std::setw(5) << strip << " " << std::setw(3) << seg 
	  : std::cout << "stripdigi " << std::setw(5) << strip << " " << std::setw(3) << seg; 
	
	if (itp==-1) // Unmatched digi
	{
	  std::cout << " // 0   0.00   0.00 // " 
		    << std::fixed << std::setprecision(2)
		    << std::setw(6) << atan2(m_pix_y[idx],m_pix_x[idx]) << " " 
		    << std::setw(6) << sqrt(m_pix_x[idx]*m_pix_x[idx]+m_pix_y[idx]*m_pix_y[idx]) << " " 
		    << std::setw(6) << m_pix_z[idx]
		    << std::fixed << std::setprecision(0)
		    << std::endl;
	}
	else // Matched digi
	{
	  if (its==-1)  std::cout << " // 0 "; // No stub
	  if (its!=-1 && isGOOD)  std::cout << " // 2 "; // In a stub from a good part 
	  if (its!=-1 && !isGOOD) std::cout << " // 1 "; // In a stub from a not good part
	
	  std:: cout << std::fixed << std::setprecision(2)
		     << std::setw(6) << ptGEN << " " 
		     << std::setw(6) << d0GEN << " // "  
		     << std::setw(6) << atan2(m_pix_y[idx],m_pix_x[idx]) << " " 
		     << std::setw(6) << sqrt(m_pix_x[idx]*m_pix_x[idx]+m_pix_y[idx]*m_pix_y[idx]) << " " 
		     << std::setw(6) << m_pix_z[idx]
		     << std::fixed << std::setprecision(0)
		     << std::endl;
	}
      }

      // Finally print the stub info
      do_stub(lay,ladder,(module-1)/2);
    }    
  }
}




//
// This method creates, for one module, a list of the pixels used 
// to make a cluster and a stub
//
// Three tables are built. One for the digi flagged, one for the corresp.
// clusters, and one for the stubs when applicable

void patterngen::ana_pix(int lay,int lad,int mod, std::vector<int> digits)
{
  m_evt_pix.clear();
  m_evt_clu.clear();
  m_evt_stu.clear();
  m_evt_tp.clear();

  int c1 = -1;
  int c2 = -1;

  std::vector<int> list_pix;
  std::vector<int> list_pix_coords;

  int row,col,idx;
  bool found;

  // Loop over clusters in the studied module
  for (int i=0;i<m_clus;++i)
  {      
    if (m_clus_layer[i]!=lay) continue;
    if (m_clus_ladder[i]!=lad) continue;
    if (int(m_clus_module[i]-1)/2!=mod) continue;

    list_pix.clear();
    list_pix_coords.clear();
    list_pix_coords = m_clus_pix[i]; // List of pixels making the cluster

    // Try to match them with the list of digits effectively recorded in the module
    for (unsigned int j=0;j<list_pix_coords.size()/2;++j)
    {
      row=list_pix_coords.at(2*j);
      col=list_pix_coords.at(2*j+1);

      found=false;

      for (unsigned int k=0;k<digits.size();++k)
      {
	if (found) continue;

	idx = digits.at(k);

	if (m_pix_module[idx]!=m_clus_module[i]) continue;
	if (m_pix_row[idx]!=row) continue;
	if (m_pix_col[idx]!=col) continue;

	found=true;

	m_evt_pix.push_back(idx);
	m_evt_clu.push_back(i);
	m_evt_stu.push_back(-1);
	m_evt_tp.push_back(m_clus_tp[i]);
      }
    }
  }

  // Then we look if we found a stub with this cluster
  for (int i=0;i<m_stub;++i)
  {      
    if (m_stub_layer[i] !=lay) continue;
    if (m_stub_ladder[i]!=lad-1) continue;
    if (m_stub_module[i]!=mod) continue;

    // First of all we compute the ID of the stub's module

    c1       = m_stub_clust1[i];
    c2       = m_stub_clust2[i];

    for (unsigned int j=0;j<m_evt_clu.size();++j)
    {
      if (m_evt_clu.at(j)==c1) m_evt_stu.at(j)=i;
      if (m_evt_clu.at(j)==c2) m_evt_stu.at(j)=i;
    }
  }
}

void patterngen::do_stub(int lay,int lad,int mod)
{
  for (int i=0;i<m_stub;++i)
  {  
    if (m_stub_layer[i]!=lay) continue;
    if (m_stub_ladder[i]!=lad-1) continue;
    if (m_stub_module[i]!=mod) continue;

    // First of all we compute the ID of the stub's module

    std::cout << "stub: " 
	      << std::setw(6) << std::fixed << std::setprecision(1) << m_stub_strip[i] << " " 
	      << std::setw(5) << m_stub_deltas[i]     
	      << " " << std::setprecision(2) 
	      << std::setw(7) <<  m_stub_z[i] << std::endl;
  }
}

//////////////////////////////////////////////
//
// Macro creating the input stream to test 
// concentrator's Verilog model 
//
// A more complex evolution is in evt_builder.cxx
//
//////////////////////////////////////////////

void patterngen::get_CONC_input(int npatt)
{
  // Initialize some params
 
  int B_id; // The detector module IDs (defined in the header)

  int layer,ladder,module,strip,chip;

  int seg;

  int n_entries = L1TT->GetEntries();

  // Then loop over events

  std::vector<int> m_stub_list;
  std::vector<int> m_digi_list;

  std::vector<std::vector <int> > evt_mix;

  std::vector<int> evt_signature;
  std::vector<int> evt_tmp;

  evt_mix.clear();
  m_data_trig.clear();
  m_data_raw.clear();

  int n_seq    = npatt;

  int pt = 0;

  //
  // We make a first loop over all entries available
  // in order to build data stores, using the digis (L1raw block)
  // and the stubs (trigger block)
  // 

  cout << "Entering loop 1, producing the big data stores..." << endl;

  for (int j=0;j<n_entries;++j)
  {    
    if (j%100==0)
      cout << j << endl;

    m_chip_trig.clear();
    m_chip_raw.clear();

    L1TT->GetEntry(j); 
    PIX->GetEntry(j); 

    // First loop over the digis

    for (int i=0;i<m_pix;++i)
    {
      layer = m_pix_layer[i]; 
      ladder= m_pix_ladder[i]-1;

      if (layer<8 || (layer>10 && ladder<9)) continue; // Select only PS modules

      module= static_cast<int>((m_pix_module[i]-1)/2); 
      seg   = m_pix_col[i]; 
      strip = m_pix_row[i]; 
      chip  = static_cast<int>(strip/127)+seg*8;
      strip = strip%128+((m_pix_module[i]-1)%2)*128;

      B_id = (layer-5)*10000 + ladder*100 + module;
      B_id = 100*B_id+chip;

      // Look if this chip has already been touched in this event
      m_iter = m_chip_raw.find(B_id);

      m_digi_list.clear();

      if (m_iter == m_chip_raw.end()) // New chip, add it
      {
	m_digi_list.push_back(strip);
	m_chip_raw.insert(std::make_pair(B_id,m_digi_list));
      }
      else // Otherwise complement
      {
	m_digi_list = m_iter->second;
	m_digi_list.push_back(strip);
	m_chip_raw.erase(m_iter->first);
	m_chip_raw.insert(std::make_pair(B_id,m_digi_list));
      }
    } // End of digi collection


    // Start the stub loop

    for (int i=0;i<m_stub;++i)
    {  
      // First of all we compute the ID of the stub's module

      layer = m_stub_layer[i]; 
      ladder= m_stub_ladder[i]; 

      if (layer<8 || (layer>10 && ladder<9)) continue; // Select only PS modules

      module= m_stub_module[i];
      seg   = m_stub_seg[i]; 
      chip  = m_stub_chip[i]+seg*8;; 
      strip = m_stub_strip[i]+1; // Code between 1 and 127 to avoid 00000000 position
      pt    = static_cast<int>(2*m_stub_deltas[i]);
         
      B_id = (layer-5)*10000 + ladder*100 + module;
      B_id = 100*B_id+chip;

      // Look if this chip has already been touched in this event
      m_iter = m_chip_trig.find(B_id);

      m_stub_list.clear();

      if (m_iter == m_chip_trig.end()) // New chip, add it
      {
	m_stub_list.push_back(strip);
	m_stub_list.push_back(pt);
	m_chip_trig.insert(std::make_pair(B_id,m_stub_list));
      }
      else // Otherwise complement
      {
	m_stub_list = m_iter->second;
	m_stub_list.push_back(strip);
	m_stub_list.push_back(pt);
	m_chip_trig.erase(m_iter->first);
	m_chip_trig.insert(std::make_pair(B_id,m_stub_list));
      }
    }  // End of digi collection

    // Finally we add both collection to the stores
    
    m_data_trig.push_back(m_chip_trig);
    m_data_raw.push_back(m_chip_raw);
  }

  // Now we have two collections one containing the raw data, and one the trigger data for all event.

  // We now create a set of 8 BX trains

  // Trigger data is sent at every BX
  // Raw data is sent every n BX, where n depends on the L1 rate

  // If the L1 rate is given in kHz, we have n=40000/rate
  //
  // For the moment we produce it only for a concentrator chip


  int mod_num  = 31020;
  int seg_side = 1;
  int idx      = 100*mod_num+8*seg_side;
  int delta    = static_cast<int>(40000/m_rate);

  int n_data = m_data_trig.size();

  std::cout << "We will trig an L1 A every " << delta << " BXs" << std::endl;

  int BX_ID    = 0;
  int Block_ID = 0;

  srand (time(NULL));

  int trg_evnum;
  int raw_evnum;
  int raw_rank;


  int raw_data[8][266]; // The L1 data blocks of the 8 chips

  for (int i=0;i<12;++i) // Header
  { 
    for (int j=0;j<8;++j) raw_data[j][i]=1;
  }

  for (int i=12;i<266;++i) // Data
  { 
    for (int j=0;j<8;++j) raw_data[j][i]=0;
  }

  std::multimap< int, std::vector<int> > trg_evt;
  std::multimap< int, std::vector<int> > raw_evt;

  // Create the sequences, and write them on the fly 
  // in the output file
  
  for (int i=0;i<n_seq;++i)
  { 
    m_outbinary << "//Event_train " << Block_ID << "\n";

    for (int j=0;j<8;++j) // Loop over events
    { 
      patterngen::initVars();

      trg_evnum = rand()%n_data;
      trg_evt   = m_data_trig.at(trg_evnum);

      if (BX_ID%delta==0) // Time to send a L1
      {
	m_raw_bx=BX_ID;

	for (int k=0;k<12;++k)
	{ 
	  for (int kk=0;kk<8;++kk) m_raw_chp[kk][k]=1;
	}

	for (int k=12;k<266;++k)
	{ 
	  for (int kk=0;kk<8;++kk) raw_data[kk][k]=0;
	}

	raw_evnum= rand()%n_data;
	raw_rank=0;
	raw_evt = m_data_raw.at(raw_evnum);

	for (int k=0;k<8;++k)
	{ 
	  m_digi_list.clear();

	  m_iter = raw_evt.find(idx+k);

	  if (m_iter != raw_evt.end())  
	  {
	    m_digi_list = m_iter->second;

	    for (unsigned int kk=0;kk<m_digi_list.size();++kk)
	    {
	      m_raw_chp[k][12+m_digi_list[kk]]=1;
	      raw_data[k][12+m_digi_list[kk]]=1;
	    }
	  }
	}

	m_raw_tree->Fill();
      }

      // First write the trigger word

      m_outbinary << "--CBC_TRIG_content_for_event " << std::bitset<3>(j) << " "  << BX_ID << " " << trg_evnum << " " << "\n";
      
      m_tri_bx=BX_ID;

      for (int k=0;k<8;++k)
      {
	m_iter = trg_evt.find(idx+k);
	
	if (m_iter != trg_evt.end())  
	{
	  m_stub_list.clear();
	  m_outbinary << "Chip " << std::bitset<3>(k) << "  " ;

	  m_stub_list = m_iter->second;

	  for (unsigned int kk=0;kk<m_stub_list.size()/2;++kk)
	  {
	    if (kk>=3) continue; // No sorting for the moment
	    
	    std::bitset<8> pos = m_stub_list.at(2*kk);
	    std::bitset<5> off = abs(2*m_stub_list.at(2*kk+1));

	    m_outbinary << pos << off << " " ;

	    //    cout << pos << " # " << off << endl ;

	    for (unsigned int ik=13*kk;ik<13*(kk+1);++ik)
	    {
	      if (ik%13<8) m_tri_chp[k][ik] = pos[7-ik%13]; // 
	      if (ik%13>=8) m_tri_chp[k][ik] = off[12-ik%13]; // 
	    }
	  }

	  m_outbinary << "\n";      
	}
      } // End of trigger data loop on chip

      m_tri_tree->Fill();


      m_outbinary << "**CBC_RAW_content_for_event_sent_at_BX " << BX_ID-raw_rank << " " << raw_rank << " " << raw_evnum << "\n";

      // Then the L1 block, when applicable

      if (raw_rank<=33)
      {
	for (int k=0;k<8;++k)
	{
	  m_outbinary << "RawChip " << std::bitset<3>(k) << "  " ;
	  
	  for (int kk=8*raw_rank;kk<8*(raw_rank+1);++kk)
	  {
	    if (kk<266) m_outbinary << raw_data[k][kk] ;
	  }
	  
	  m_outbinary << "\n";      
	}
      }
  
      ++raw_rank;
      ++BX_ID;

      m_outbinary << "\n";  
   } // End of loop over events

    ++Block_ID;
    m_outbinary  << "#########################################" << "\n";

  }
  

  m_outbinary.close();

  m_outfile->Write();
  
  delete L1TT;
  delete PIX;
  delete m_outfile;
}

//////////////////////////////////////////////////////////////////////
//
// This method reads and decode the concentrator output
// then prepare an output file to compare it to the input data 
// file initially sent
//
// First decode the trigger block, then decode the L1 raw block, either sparsified or not 
//
//
//////////////////////////////////////////////////////////////////////

void patterngen::get_CONC_output(std::string in)
{

  // Open the input data file

  std::string STRING,L1BX,L1FRAG;

  std::ifstream in2(in.c_str());
  if (!in2)
  {
    std::cout << "Please provide a valid data filename list" << std::endl;
    return;
  }

  int sparsified = 0;

  std::size_t found; 

  found = in.find("unsparsified");
  if (found!=std::string::npos) 
  {
    std::cout << "Dealing with unsparsified raw data" << std::endl;
  }
  else
  {
    sparsified = 1;
    std::cout << "Dealing with sparsified raw data" << std::endl;
  }

  patterngen::initVars();


  int ntr=0;
  int bx = 0;
  char *evt_num  = new char[3];
  char *cbc_word = new char[43];
  char *cbc_raw  = new char[11];

  int chip_num  = 0;
  int chip_prev = 0;
  int n_clus    = 0;
  int BX_prev   = -1;
  int BX_L1     = 0;
  int frag_L1   = 0;
  int bit_count = 0;
  int nbits     = 0; 

  int new_L1    = 0;
  int new_trg   = 0;

  std::vector<int> spars_raw_tmp;
  
  int pos;
  int wdth;


  // Read the concentrator file line per line

  while (!in2.eof())
  {
    getline(in2,STRING);
    
    // New concentrator block
    found = STRING.find("//Event_train");

    if (found!=std::string::npos) 
    {	
      if (new_trg==1) // The previous event was empty
      {
	m_tri_tree->Fill();
	new_trg = 0;
      }

      STRING = STRING.substr(found+13);
      STRING.erase(std::remove(STRING.begin(), STRING.end(), ' '), STRING.end());
      ntr = std::stoi(STRING);
      bx  = 8*ntr; 
    }


    // Dealing with a Trigger fragment
    found = STRING.find("--CBC_TRIG_content_for_event");

    if (found!=std::string::npos) 
    {	
      if (new_trg==1) // The previous event was empty
      {
	m_tri_tree->Fill();
	new_trg = 0;
      }

      new_trg = 1;

      for (int i=0;i<8;++i) // New event, reset the trigger part
      {
	for (int j=0;j<40;++j)  m_tri_chp[i][j]  = 0.;
      }

      STRING = STRING.substr(found+28);
      STRING.erase(std::remove(STRING.begin(), STRING.end(), ' '), STRING.end());
      std::strcpy(evt_num, STRING.c_str());
   
      m_tri_bx = bx + 4*int(evt_num[0]-'0')+2*int(evt_num[1]-'0')+int(evt_num[2]-'0'); 
    }

    // Chip trigger word within an event (40 bits from CBC)
    found = STRING.find("Chip");

    if (found!=std::string::npos && new_trg==1) 
    {
      for (int j=0;j<43;++j)  cbc_word[j] = '0';
      STRING = STRING.substr(found+4);
      STRING.erase(std::remove(STRING.begin(), STRING.end(), ' '), STRING.end());
      std::strcpy(cbc_word, STRING.c_str());

      chip_num = 4*int(cbc_word[0]-'0')+2*int(cbc_word[1]-'0')+int(cbc_word[2]-'0');

      for (int j=0;j<40;++j)
      {
	(int(cbc_word[j+3]-'0')<0) 
	  ? m_tri_chp[chip_num][j] = 0 
	  : m_tri_chp[chip_num][j] = int(cbc_word[j+3]-'0');
      }
      
      if (chip_num==8)
      {
	m_tri_tree->Fill();
	new_trg = 0;
      }
    }
  

    // Dealing with a L1Raw block
    found = STRING.find("**CBC_RAW_content_for_event_sent_at_BX");

    if (found!=std::string::npos) 
    {

      // First handle the unsparsified case
      if (sparsified==0) 
      {	
	if (new_L1==0) // Start of a new L1 event 
	{
	  new_L1    = 1;
	  bit_count = 0;

	  for (int i=0;i<8;++i)
	  {
	    for (int j=0;j<266;++j)  m_raw_chp[i][j] = 0;
	  }
	}

	L1BX   = STRING.substr(found+39,14);
	L1FRAG = STRING.substr(found+53);
	L1BX.erase(std::remove(L1BX.begin(), L1BX.end(), ' '), L1BX.end());
	L1FRAG.erase(std::remove(L1FRAG.begin(), L1FRAG.end(), ' '), L1FRAG.end());
      
	BX_L1   = std::stoi(L1BX);
	frag_L1 = std::stoi(L1FRAG);
	
	// Recovering the raw word
	getline(in2,STRING);
	
	found = STRING.find("RawChip");
	if (found!=std::string::npos) 
	{
	  STRING = STRING.substr(found+8);
	  STRING.erase(std::remove(STRING.begin(), STRING.end(), ' '), STRING.end());
	  std::strcpy(cbc_raw, STRING.c_str());
	  
	  chip_num = 4*int(cbc_raw[0]-'0')+2*int(cbc_raw[1]-'0')+int(cbc_raw[2]-'0');
	  
	  nbits=8; // Number of bits to acount for in the word
	  
	  // Special cases (EOF)
	  if (frag_L1==34  && chip_num==0) nbits=2; 
	  if (frag_L1==34  && chip_num==1) nbits=6; 
	  if (frag_L1==67  && chip_num==1) nbits=4; 
	  if (frag_L1==67  && chip_num==2) nbits=4;
	  if (frag_L1==100 && chip_num==2) nbits=6; 
	  if (frag_L1==100 && chip_num==3) nbits=2; 
	  if (frag_L1==167 && chip_num==4) nbits=2; 
	  if (frag_L1==167 && chip_num==5) nbits=6; 
	  if (frag_L1==200 && chip_num==5) nbits=4; 
	  if (frag_L1==200 && chip_num==6) nbits=4;
	  if (frag_L1==233 && chip_num==6) nbits=6; 
	  if (frag_L1==233 && chip_num==7) nbits=2;  
	  
	  if (chip_num==7 && frag_L1==266)
	  {
	    std::cout << "Dealing with fragment " << frag_L1%33 
		      << " of chip " << chip_num << " so extracting " 
		      <<  nbits << " bits (" << bit_count << " to " 
		      << bit_count+nbits << ")" << std::endl;    
	  }

	  for (int j=0;j<nbits;++j)
	  {
	    m_raw_chp[chip_num][(j+bit_count)%266] = int(cbc_raw[j+3]-'0');
	    
	    if (chip_num==7 && frag_L1==266)
	    {
	      std::cout << j << " /  " << (j+bit_count)%266 << " /  " <<  int(cbc_raw[j+3]-'0') << " //  " ;
	    }
	  }

	  if (chip_num==7 && frag_L1==266)
	  {
	    std::cout  << std::endl;    
	  }

	  bit_count+=nbits;
	  
	  if (frag_L1==266) //EOF
	  {
	    new_L1    = 0;
	    m_raw_tree->Fill();
	  }
	}
	else
	{
	  cout << " Problem !! There should be a word here!!!" << endl;
	}
    
	m_raw_bits = 266*8;
	m_raw_bx   = BX_L1;
      }
      else // Sparsified
      {

	L1BX   = STRING.substr(found+39,14);
	L1FRAG = STRING.substr(found+53);
	L1BX.erase(std::remove(L1BX.begin(), L1BX.end(), ' '), L1BX.end());
	L1FRAG.erase(std::remove(L1FRAG.begin(), L1FRAG.end(), ' '), L1FRAG.end());
	
	BX_L1   = std::stoi(L1BX);
	frag_L1 = std::stoi(L1FRAG);
	
	//	cout << BX_L1 << " / " << frag_L1 << endl;

	if (BX_L1!=BX_prev)  // Start of a new L1 event 
	{	  
	  m_raw_bx   = BX_prev;
	  BX_prev    = BX_L1;
	  new_L1     = 1;
	  m_raw_bits = bit_count;

	  cout << bit_count << " bits were collected " << endl;
	  cout << " now analyzing them " << endl;
	  
	  chip_prev = -1;
	  
	  for (int i=0;i<8;++i)
	  {
	    for (int j=0;j<266;++j)  m_raw_chp[i][j] = 0;
	  }
	  
	  for (int i=0;i<bit_count;++i)
	  {
	    cout << spars_raw_tmp.at(2*i+1) << "/";
	  }
	  cout << endl;
	  
	  for (int i=0;i<bit_count;++i)
	  {
	    chip_num = spars_raw_tmp.at(2*i);
	    
	    if (chip_num!=chip_prev) // New chip
	    {
	      n_clus   = 0;

	      for (int j=i;j<i+12;++j) // Read header (12 bits)
		m_raw_chp[chip_num][j-i] = spars_raw_tmp.at(2*j+1);

	      for (int j=i+12;j<i+17;++j) // Number of clusters (5 bits)
	      {
		if (i+16-j != 0)
		{
		  n_clus += pow(2*spars_raw_tmp.at(2*j+1),i+16-j);
		}
		else if (spars_raw_tmp.at(2*j+1)!=0)
		{
		  n_clus += 1;
		}
	      }

	      i+=16; // Skip header + nclus block

	      //	      cout << n_clus << " clusters in chip " << chip_num << endl;
	    
	      // Expect 11 bits per stored cluster, 8 for position and 3 for width

	      if (n_clus!=0) 
	      {
		for (int j=0;j<n_clus;++j) 
		{
		  pos=0;
		  wdth=0;
		  
		  for (int k=i+1;k<=i+8;++k) // Read position (8 bits)
		  { 
		    if (i+8-k != 0)
		    {
		      pos += pow(2*spars_raw_tmp.at(2*k+1),i+8-k);
		    }
		    else if (spars_raw_tmp.at(2*k+1)!=0)
		    {
		      pos += 1;
		    }
		  }

		  for (int k=i+9;k<=i+11;++k) // Read width (3 bits)
		  { 
		    if (i+11-k != 0)
		    {
		      wdth += pow(2*spars_raw_tmp.at(2*k+1),i+11-k);
		    }
		    else if (spars_raw_tmp.at(2*k+1)!=0)
		    {
		      wdth += 1;
		    }
		  }
		  
		  i+=11;

		  cout << j << " / " << pos << " / " << wdth << endl;

		  // Pos is coded from 1

		  for (int k=pos+11;k<pos+11+wdth;++k) m_raw_chp[chip_num][k]=1;
		}
	      } // End of cluster decoding
	    }
	  }

	  bit_count=0;
	  spars_raw_tmp.clear();
	 
	  if (m_raw_bx!=-1) m_raw_tree->Fill();
	}

	// Recovering the raw word
	getline(in2,STRING);

	found = STRING.find("RawChip");

	if (found!=std::string::npos) 
        {
	  STRING = STRING.substr(found+8);
	  STRING.erase(std::remove(STRING.begin(), STRING.end(), ' '), STRING.end());
	  std::strcpy(cbc_raw, STRING.c_str());
	  
	  chip_num = 4*int(cbc_raw[0]-'0')+2*int(cbc_raw[1]-'0')+int(cbc_raw[2]-'0');
	  
	  for (int j=0;j<8;++j)
	  {
	    if (int(cbc_raw[j+3]-'0')<0) break;
	    
	    ++bit_count;
	    spars_raw_tmp.push_back(chip_num);
	    spars_raw_tmp.push_back(int(cbc_raw[j+3]-'0'));
	  }
	}
	else
	{
	  cout << " Problem !! There should be a word here!!!" << endl;
	}    
      }
    }
  }
 
  in2.close();

  m_outfile->Write();
  
  delete m_outfile;
}




/////////////////////////////////////////////////////////////
//
// Basic methods, initializations,...
//
/////////////////////////////////////////////////////////////

void patterngen::initVars()
{
  m_tri_bx=0;
  m_raw_bx=0;
  m_tri_chip=0;
  m_raw_chip=0;
  m_raw_bits=0;

  for (int i=0;i<8;++i)
  { 
    for (int j=0;j<40;++j)  m_tri_chp[i][j]  = 0.;
    for (int j=0;j<266;++j) m_raw_chp[i][j]  = 0;
  }

}

//
// ==> patterngen::initTuple(std::string in,std::string out,int type)
//
// Initialize the rootuple to read and to create, also initialize the text files if necessary 
// 

void patterngen::initTuple(std::string in,std::string out,int type)
{

  m_outfile  = new TFile(out.c_str(),"recreate");

  if (type >= 1)
  {
    m_tri_tree = new TTree("Trigger_bef","L1Trigger words");
    m_raw_tree = new TTree("Raw_bef","Raw data words");
  }
  else
  {
    m_tri_tree = new TTree("Trigger_aft","L1Trigger words");
    m_raw_tree = new TTree("Raw_aft","Raw data words");
  }

  m_tri_tree->Branch("TRI_BX",         &m_tri_bx,      "TRI_BX/I");
  m_tri_tree->Branch("TRI_CHIP",       &m_tri_chp,     "TRI_CHIP[8][40]/I");
  
  m_raw_tree->Branch("RAW_BX",         &m_raw_bx,      "RAW_BX/I");
  m_raw_tree->Branch("RAW_BITS",       &m_raw_bits,    "RAW_BITS/I");
  m_raw_tree->Branch("RAW_CHIP",       &m_raw_chp,     "RAW_CHIP[8][266]/I");

  if (type == 0) return;
 
  L1TT   = new TChain("TkStubs");
  PIX    = new TChain("Pixels"); 
  MC     = new TChain("MC");   
 
  // Input data file
  
  std::size_t found = in.find(".root");
  
  // Case 1, it's a root file
  if (found!=std::string::npos)
  {
    L1TT->Add(in.c_str());
    PIX->Add(in.c_str());
    MC->Add(in.c_str());
  }
  else // This is a list provided into a text file
  {
    std::string STRING;
    std::ifstream in2(in.c_str());
    if (!in2)
    {
      std::cout << "Please provide a valid data filename list" << std::endl;
      return;
    }
  
    while (!in2.eof())
    {
      getline(in2,STRING);
      
      found = STRING.find(".root");
      if (found!=std::string::npos) 
      {	
	L1TT->Add(STRING.c_str());
	PIX->Add(STRING.c_str());
	MC->Add(STRING.c_str());
      }
    }
    
    in2.close();
  }


  pm_part_px=&m_part_px;
  pm_part_py=&m_part_py;
  pm_part_eta=&m_part_eta;
  pm_part_x=&m_part_x;
  pm_part_y=&m_part_y;
  pm_part_z=&m_part_z;

  MC->SetBranchAddress("subpart_n",        &m_ntp);    
  MC->SetBranchAddress("subpart_x",        &pm_part_x);    
  MC->SetBranchAddress("subpart_y",        &pm_part_y);    
  MC->SetBranchAddress("subpart_z",        &pm_part_z);    
  MC->SetBranchAddress("subpart_px",       &pm_part_px);   
  MC->SetBranchAddress("subpart_py",       &pm_part_py);   
  MC->SetBranchAddress("subpart_eta",      &pm_part_eta);  
  
  pm_pix_layer=&m_pix_layer;
  pm_pix_ladder=&m_pix_ladder;
  pm_pix_module=&m_pix_module;
  pm_pix_row=&m_pix_row;
  pm_pix_col=&m_pix_col;
  pm_pix_x=&m_pix_x;
  pm_pix_y=&m_pix_y;
  pm_pix_z=&m_pix_z;
  pm_pix_ch=&m_pix_ch;

  PIX->SetBranchAddress("PIX_n",         &m_pix);
  PIX->SetBranchAddress("PIX_nPU",       &m_npu);
  PIX->SetBranchAddress("PIX_layer",     &pm_pix_layer);
  PIX->SetBranchAddress("PIX_ladder",    &pm_pix_ladder);
  PIX->SetBranchAddress("PIX_module",    &pm_pix_module);
  PIX->SetBranchAddress("PIX_row",       &pm_pix_row);
  PIX->SetBranchAddress("PIX_column",    &pm_pix_col);
  PIX->SetBranchAddress("PIX_charge",    &pm_pix_ch);
  PIX->SetBranchAddress("PIX_x",         &pm_pix_x);
  PIX->SetBranchAddress("PIX_y",         &pm_pix_y);
  PIX->SetBranchAddress("PIX_z",         &pm_pix_z);

  pm_stub_layer=&m_stub_layer;
  pm_stub_ladder=&m_stub_ladder;
  pm_stub_module=&m_stub_module;
  pm_stub_pt=&m_stub_pt;
  pm_stub_tp=&m_stub_tp;
  pm_stub_deltas=&m_stub_deltas;
  pm_stub_strip=&m_stub_strip;
  pm_stub_seg=&m_stub_seg;
  pm_stub_z=&m_stub_z;
  pm_clus_nseg=&m_clus_nseg;
  pm_clus_pix=&m_clus_pix;
  pm_clus_mult=&m_clus_mult;
  pm_stub_chip=&m_stub_chip;
  pm_stub_clust1=&m_stub_clust1;
  pm_stub_clust2=&m_stub_clust2;
  pm_stub_pxGEN=&m_stub_pxGEN;
  pm_stub_pyGEN=&m_stub_pyGEN;
  pm_stub_etaGEN=&m_stub_etaGEN;
  pm_stub_X0=&m_stub_X0;
  pm_stub_Y0=&m_stub_Y0;
  pm_stub_Z0=&m_stub_Z0;

  pm_clus_layer=&m_clus_layer;
  pm_clus_ladder=&m_clus_ladder;
  pm_clus_module=&m_clus_module;
  pm_clus_tp=&m_clus_tp;
 
  L1TT->SetBranchAddress("L1TkSTUB_n",         &m_stub);
  L1TT->SetBranchAddress("L1TkSTUB_layer",     &pm_stub_layer);
  L1TT->SetBranchAddress("L1TkSTUB_ladder",    &pm_stub_ladder);
  L1TT->SetBranchAddress("L1TkSTUB_module",    &pm_stub_module);
  L1TT->SetBranchAddress("L1TkSTUB_pt",        &pm_stub_pt);
  L1TT->SetBranchAddress("L1TkSTUB_z",         &pm_stub_z);
  L1TT->SetBranchAddress("L1TkSTUB_tp",        &pm_stub_tp);
  L1TT->SetBranchAddress("L1TkSTUB_deltas",    &pm_stub_deltas);
  L1TT->SetBranchAddress("L1TkSTUB_strip",     &pm_stub_strip);
  L1TT->SetBranchAddress("L1TkSTUB_seg",       &pm_stub_seg);
  L1TT->SetBranchAddress("L1TkSTUB_chip",      &pm_stub_chip);
  L1TT->SetBranchAddress("L1TkSTUB_clust1",    &pm_stub_clust1);
  L1TT->SetBranchAddress("L1TkSTUB_clust2",    &pm_stub_clust2);
  L1TT->SetBranchAddress("L1TkCLUS_PS",        &pm_clus_nseg);
  
  L1TT->SetBranchAddress("L1TkCLUS_n",         &m_clus);
  L1TT->SetBranchAddress("L1TkCLUS_PIX",       &pm_clus_pix);
  L1TT->SetBranchAddress("L1TkCLUS_MULT",      &pm_clus_mult);
  L1TT->SetBranchAddress("L1TkCLUS_tp",        &pm_clus_tp);
  L1TT->SetBranchAddress("L1TkCLUS_layer",     &pm_clus_layer);
  L1TT->SetBranchAddress("L1TkCLUS_ladder",    &pm_clus_ladder);
  L1TT->SetBranchAddress("L1TkCLUS_module",    &pm_clus_module);
  
  L1TT->SetBranchAddress("L1TkSTUB_pxGEN",     &pm_stub_pxGEN);
  L1TT->SetBranchAddress("L1TkSTUB_pyGEN",     &pm_stub_pyGEN);
  L1TT->SetBranchAddress("L1TkSTUB_etaGEN",    &pm_stub_etaGEN);
  L1TT->SetBranchAddress("L1TkSTUB_X0",        &pm_stub_X0);
  L1TT->SetBranchAddress("L1TkSTUB_Y0",        &pm_stub_Y0);
  L1TT->SetBranchAddress("L1TkSTUB_Z0",        &pm_stub_Z0);
 
  m_outbinary.open("concentrator_input.txt");
  
}


