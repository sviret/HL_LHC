// Base class for Phase II tracker FrontEnd tests stimulus builder
//
// Main methods for the trigger block construction
// (see stim_builder_L1 for the L1 block)
//
//
// Meaning of the parameters
//
//
// filenameRAW : the input data files for L1-like events
// filenameTRG : the input data files for minbias main events
// outfile     : name of the output root file containing some crosscheck datas
// npatt       : how many BXs should be produced (will be fitted correctly afterwards)
// layer       : restrict info to a given layer, or not (give -1 in this case)
// ladder      : restrict info to a given ladder, or not (give -1 in this case)
// module      : restrict info to a given module, or not (give -1 in this case)
// tower       : restrict info to a given trigger tower, or not (give -1 in this case)
// sector      : the input file defining the modules
// L1prop      : The proportion of rich events you want in the TRG flow (in %)
// TRGsize     : The size of the CIC trigger block (in bits)


#include "stim_builder.h"

stim_builder::stim_builder(std::string filenameRAW, std::string filenameTRG, std::string outfile, 
			   int npatt, int layer, int ladder, int module, int tower, std::string sector,
			   int L1prop, int TRGsize)
{
    m_lay        = layer;
    m_lad        = ladder;
    m_mod        = module;
    m_tri_data   = new  std::vector<int>;
    m_raw_data   = new  std::vector<int>;
    m_L1prop     = L1prop;
    m_CICsize    = 8;
    bend_bit_MPA = 3;
    bend_bit_CBC = 4;
    m_tower      = tower;
    m_TRGsize    = TRGsize;
    
    stim_builder::initVars();                                 // Initialize everything
    stim_builder::convert(sector);                            // Get the trigger tower module info
    stim_builder::initTuple(filenameRAW,filenameTRG,outfile); // Initialize ROOT stuff
    stim_builder::get_stores(npatt-npatt%8,0);                // Prepare the data store where to pickup info
}

/////////////////////////////////////////////////////////
//
//
// get_stores() is the core of the event builder
//
// This method is divided into 3 steps. In the first step we build a data store by picking events 
// in CMSSW samples, then we generate random sequences of events (L1 and TRIG). Finally, depending on the 
// type of data requested, we produce a root file with the output and a text file with the corresponding sequence
//
/////////////////////////////////////////////////////////

void stim_builder::get_stores(int nevts, bool conc)
{
 
    int B_id, B_id_conc; // The detector module IDs (defined in the header)

    int layer,ladder,module,strip,chip,seg,nseg;

    int n_entries_TRG = L1TT->GetEntries();

    bool  isPS        = false;
    
    int idx=-1;
    int its=-1;
    int itp=-1;

    // Then loop over events

    std::vector<int> m_stub_list;
    std::vector<int> m_digi_list;

    m_data_trig.clear();
    m_data_raw.clear();

    //
    // We make a first loop over all entries available
    // in order to build data stores, using the digis (L1raw block)
    // and the stubs (trigger block)
    // 
    // 2 data samples are defined for that, one with m_L1prop % physics events, 
    // for the trigger words, and one with pure physics events, for L1 data transmission 
    //
    // You can choose the proportion of complex events you want in the TRG block, in % 
    //

    int mixpar = 0;
    int store_size;

    store_size = n_entries_TRG;

    cout << "--> Entering loop 1, producing the big data stores for " 
	 << store_size << " events..." << endl;

    if (m_chips.size()==0 || m_concs.size()==0)
    {
      cout << "You haven't selected any module... " << endl;
      cout << "Try again... " << endl;
    }
    
    cout << "--> Considering " 
	 << m_chips.size() << " FE chips and..."
	 << m_concs.size() << " CIC chips..." << endl;


    for (int j=0;j<store_size;++j)
    {    
      if (j%10==0)
	cout << "Processing event " <<  j << "/" << store_size << endl;

      m_chip_trig.clear();
      m_chip_raw.clear();
      m_conc_trig.clear();
      m_conc_raw.clear();

      if (rand()%100>m_L1prop) // Mbias event in the trigger block
      {
	mixpar=rand()%(n_entries_TRG-m_PHYsize)+m_PHYsize;
	L1TT->GetEntry(mixpar);
	PIX->GetEntry(mixpar);
      }
      else // Phys event (stored in the s2nd half of the TRG sample)
      {
	mixpar=rand()%m_PHYsize; 
	L1TT->GetEntry(mixpar); 
	PIX->GetEntry(mixpar); 
      }
   
      // Initialize the map with dummy values for all referenced modules
      //
      // Please keep in mind that this map is 8 times larger in the FE output case
      // Maps are quite heavy, so don't produce a lot of FE events
      
      for (unsigned int i=0;i<m_chips.size();++i) // All the chips if one asks the FE output
      {
	m_digi_list.clear();
	m_digi_list.push_back(mixpar);
	m_chip_raw.insert(std::make_pair(m_chips.at(i),m_digi_list));
	m_chip_trig.insert(std::make_pair(m_chips.at(i),m_digi_list));
      }

      for (unsigned int i=0;i<m_concs.size();++i) // All the concentrators if one asks the CONC output (8 times less)
      {
	m_digi_list.clear();
	m_digi_list.push_back(mixpar);
	m_conc_raw.insert(std::make_pair(m_concs.at(i),m_digi_list));
	m_conc_trig.insert(std::make_pair(m_concs.at(i),m_digi_list));
      }

      // First loop over the digis (raw block)

      for (int i=0;i<m_pix;++i) // Loop over pixels
      {
        isPS  = false;
        layer = m_pix_layer[i];
	
        if (layer!=m_lay && m_lay!=-1) continue; // By default we loop over all layers (-1)
	
        ladder= m_pix_ladder[i]-1;
      
        if (ladder!=m_lad && m_lad!=-1) continue; // By default we loop over all ladders (-1)
        if (layer<8 || (layer>10 && ladder<9)) isPS = true;

        module= static_cast<int>(m_pix_module[i]-1);
        
        if (module!=m_mod && m_mod!=-1) continue; // By default we loop over all modules (-1)
        seg   = m_pix_col[i];
        strip = m_pix_row[i];
	nseg  = m_pix_ncol[i]; 

	if (isPS && nseg==32) // Pixel size
	{
	  chip  = static_cast<int>(strip/120)+(seg/16)*8;
	  strip = strip%120;
	}
	else if (isPS && m_pix_module[i]%2==0) // Ger the chip number for the PS-S
	{
	  chip  = static_cast<int>(strip/120)+seg*8;
	  strip = strip%120+120;
	}
	else // For the 2S
	{
	  chip  = static_cast<int>(strip/127)+seg*8;
	  strip = strip%127+(1-m_pix_bot[i])*127;
	}

	//   cout << m_pix_module[i] << " / " << strip << " / " << seg << " / " << chip << endl;       
	// cout << m_pix_module[i] << " // " << strip << endl;
	
        B_id = layer*1000000 + ladder*10000 + module*100 + chip; // Finally get the FE chip ID

        // Look if this chip has already been touched in this event
        m_iter  = m_chip_raw.find(B_id);
	
        if (m_iter == m_chip_raw.end()) // Unknown chip, this can happen because csv file is only considering chips involved in TRG
        {
            continue; // We don't consider them for the moment...
            m_digi_list.clear();
            m_digi_list.push_back(mixpar);
            m_chip_raw.insert(std::make_pair(B_id,m_digi_list));
            m_iter = m_chip_raw.find(B_id);
        }

       // cout << j << " PIX " << B_id << " / " << m_pix_module[i] << " / " << m_pix_row[i] << " / " << m_pix_col[i] << " / " << chip << endl;
        
        m_digi_list.clear();
        m_digi_list = m_iter->second;
        m_digi_list.push_back(i);

        m_chip_raw.erase(m_iter->first);
        m_chip_raw.insert(std::make_pair(B_id,m_digi_list));

      } // End of digi collection
 
      // return;
      
      for (int i=0;i<m_stub;++i) // Loop over stubs
      {  
	// First of all we compute the ID of the stub's module
        isPS     = false;
        layer    = m_stub_layer[i];
	
        if (layer!=m_lay && m_lay!=-1) continue; // By default we loop over all layers (-1)
	
        ladder   = m_stub_ladder[i]-1;
      
        if (ladder!=m_lad && m_lad!=-1) continue; // By default we loop over all ladders (-1)
        if (layer<8 || (layer>10 && ladder<9)) isPS = true;
      
        module   = m_stub_module[i]-1;
      
        if (module!=m_mod && m_mod!=-1) continue; // By default we loop over all modules (-1)
      
        seg      = m_stub_seg[i];
      
        (isPS)
        ? chip  =  m_stub_chip[i]+(seg/16)*8
        : chip  =  m_stub_chip[i]+seg*8;

        (isPS)
        ? strip = int(2*m_stub_strip[i])%240  // Bet 0 and 255
        : strip = int(2*m_stub_strip[i])%254+1; // Code between 1 and 254 to avoid 00000000 position

        B_id      = layer*1000000 + ladder*10000 + module*100 + chip;
        B_id_conc = B_id - chip%8; // Here we get the concentrator ID
	
        // Look if this chip has already been touched in this event
        m_iter  = m_chip_trig.find(B_id);
        m_iter2 = m_conc_trig.find(B_id_conc);
	
        if (m_iter == m_chip_trig.end()) // Unknown chip???
        {
            cout << "Wow, this should not happen, the chip ref " << B_id << " is unknown!!!" << endl;
            continue;
        }
        else // Otherwise complement
        {
         //   cout << j << " STUB " << B_id << " / " << m_stub_module[i] << " / " << m_stub_strip[i] << " / " << m_stub_seg[i] << " / " << chip << endl;
            
            m_stub_list.clear();
            m_stub_list = m_iter->second;
            m_stub_list.push_back(i);

            m_chip_trig.erase(m_iter->first);
            m_chip_trig.insert(std::make_pair(B_id,m_stub_list));
        }

        if (m_iter2 == m_conc_trig.end()) // Unknown chip???
        {
            cout << "Wow, this should not happen, the CIC ref " << B_id_conc << " is unknown!!!" << endl;
            continue;
        }
        else // Otherwise complement
        {
            m_stub_list.clear();
            m_stub_list = m_iter2->second;
            m_stub_list.push_back(i);

            m_conc_trig.erase(m_iter2->first);
            m_conc_trig.insert(std::make_pair(B_id_conc,m_stub_list));
        }
      }  // End of stub collection

      // Finally we add both collection to the stores
    
      m_data_trig.push_back(m_chip_trig);
      m_data_trig.push_back(m_conc_trig);
      m_data_raw.push_back(m_chip_raw);
      
    } // End of storage loop

    // Now we have two collections one containing the raw data,
    // and one the trigger data for all event at the concentrator level

    // Trigger data is sent at every BX

    cout << "--> Entering loop 2, producing the random Trig/L1 sequence..." << endl;

    int n_trig = m_data_raw.size();
    int trg_evnum;

    std::vector<int> trig_seq;

    trig_seq.clear();

    cout << "... Generate the sequence of " << nevts << " by picking up events in the store..." << std::endl;
    
    for (int i=0;i<nevts;++i)
    { 
      trg_evnum = rand()%n_trig;
      trig_seq.push_back(trg_evnum);
    }

    // Now we have a sequence of nevts, with L1A distributed randomly 
    // following the trigger rules
    
    // Next stage consists in creating a ROOT/txt file containing this sequence 
    // of events, with the correspondind data for each Concentrator chip
    
    cout << "--> Entering loop 3, producing the final root file..." << endl;

    // Write the trigger data 
    // Everything is written down on a 2BX basis for CBC/MPA   
    
    char buffer1[200];
    char buffer2[200];
    char buffer3[200];
    
    if (m_tower!=-1)
    {
      sprintf(buffer1, "trg_FE_IN_Tow%d_%d.txt",m_tower,nevts);
      sprintf(buffer2, "trg_FE_OUT_Tow%d_%d.txt",m_tower,nevts);
      sprintf(buffer3, "trg_CIC_OUT_Tow%d_%d.txt",m_tower,nevts);
    }
    else
    {
      if (m_lay!=-1 && m_lad!=-1 && m_mod!=-1)
      {
	sprintf(buffer1, "trg_FE_IN_L%dR%dM%d_%d.txt",m_lay,m_lad,m_mod,nevts);
	sprintf(buffer2, "trg_FE_OUT_L%dR%dM%d_%d.txt",m_lay,m_lad,m_mod,nevts);
	sprintf(buffer3, "trg_CIC_OUT_L%dR%dM%d_%d.txt",m_lay,m_lad,m_mod,nevts);
      }
      else if (m_lay!=-1 && m_lad!=-1)
      {
	sprintf(buffer1, "trg_FE_IN_L%dR%d_%d.txt",m_lay,m_lad,nevts);
	sprintf(buffer2, "trg_FE_OUT_L%dR%d_%d.txt",m_lay,m_lad,nevts);
	sprintf(buffer3, "trg_CIC_OUT_L%dR%d_%d.txt",m_lay,m_lad,nevts);
      }
      else if (m_lay!=-1 && m_mod!=-1)
      {
	sprintf(buffer1, "trg_FE_IN_L%dM%d_%d.txt",m_lay,m_mod,nevts);
	sprintf(buffer2, "trg_FE_OUT_L%dM%d_%d.txt",m_lay,m_mod,nevts);
	sprintf(buffer3, "trg_CIC_OUT_L%dM%d_%d.txt",m_lay,m_mod,nevts);
      }
      else if (m_lad!=-1 && m_mod!=-1)
      {
	sprintf(buffer1, "trg_FE_IN_R%dM%d_%d.txt",m_lad,m_mod,nevts);
	sprintf(buffer2, "trg_FE_OUT_R%dM%d_%d.txt",m_lad,m_mod,nevts);
	sprintf(buffer3, "trg_CIC_OUT_R%dM%d_%d.txt",m_lad,m_mod,nevts);
      }
      else if (m_lay!=-1)
      {
	sprintf(buffer1, "trg_FE_IN_L%dR%d_%d.txt",m_lay,m_lad,nevts);
	sprintf(buffer2, "trg_FE_OUT_L%dR%d_%d.txt",m_lay,m_lad,nevts);
	sprintf(buffer3, "trg_CIC_OUT_L%dR%d_%d.txt",m_lay,m_lad,nevts);
      }
      else if (m_mod!=-1)
      {
	sprintf(buffer1, "trg_FE_IN_L%dM%d_%d.txt",m_lay,m_mod,nevts);
	sprintf(buffer2, "trg_FE_OUT_L%dM%d_%d.txt",m_lay,m_mod,nevts);
	sprintf(buffer3, "trg_CIC_OUT_L%dM%d_%d.txt",m_lay,m_mod,nevts);
      }
      else if (m_lad!=-1)
      {
	sprintf(buffer1, "trg_FE_IN_R%dM%d_%d.txt",m_lad,m_mod,nevts);
	sprintf(buffer2, "trg_FE_OUT_R%dM%d_%d.txt",m_lad,m_mod,nevts);
	sprintf(buffer3, "trg_CIC_OUT_R%dM%d_%d.txt",m_lad,m_mod,nevts);
      }
      else
      {
	sprintf(buffer1, "trg_FE_IN_ALL_%d.txt",nevts);
	sprintf(buffer2, "trg_FE_OUT_ALL_%d.txt",nevts);
	sprintf(buffer3, "trg_CIC_OUT_ALL_%d.txt",nevts);
      }
    }
    
    FE_TRG_IN.open(buffer1);
    FE_TRG_OUT.open(buffer2);
    CIC_TRG_OUT.open(buffer3);

    FE_TRG_IN << "Digital input to the FE chip.\n";
    
    
    FE_TRG_OUT << "Digital output of the FE chip.\n";
    FE_TRG_OUT << "\n";
    FE_TRG_OUT << "Format defined in:\n";
    FE_TRG_OUT << "https://espace.cern.ch/Tracker-Upgrade/Electronics/CIC/Shared_Documents/Data_formats/CIC_IO_Formats_v2.pdf\n";
    
    CIC_TRG_OUT << "Digital output of the CIC chip.\n";
    CIC_TRG_OUT << "\n";
    CIC_TRG_OUT << "Format defined in:\n";
    CIC_TRG_OUT << "https://espace.cern.ch/Tracker-Upgrade/Electronics/CIC/Shared_Documents/Data_formats/CIC_IO_Formats_v2.pdf\n";
    
    float ptGEN;
    float d0GEN;
    bool isGOOD;
    
    int evtID;
    
    for (int i=0;i<nevts;++i) // Create the sequence
    {
      if (i%10==0) cout << i << endl;
      
      m_chip_raw  = m_data_raw.at(trig_seq.at(i)); // Pick up the event in the store
      
      evtID = m_chip_raw.begin()->second.at(0);
      
      PIX->GetEntry(evtID);
      MC->GetEntry(evtID);
      
      FE_TRG_IN << "\n";
      FE_TRG_IN << "##################################################\n";
      FE_TRG_IN << "#\n";
      FE_TRG_IN << "# Event  " << i << "\n";
      FE_TRG_IN << "# NPU  = " << m_npu << "\n";

      if (i%8==0)
      {
	CIC_TRG_OUT << "\n";
	CIC_TRG_OUT << "##################################################\n";
	CIC_TRG_OUT << "#\n";
	CIC_TRG_OUT << "# Events  " << i << " to " << i+7 << "\n";
        
	m_chip_trig   = m_data_trig.at(2*trig_seq.at(i)+1);
        
	for ( m_iter = m_chip_trig.begin(); m_iter != m_chip_trig.end();++m_iter )
        {
	  trig_sequence.clear();
          
	  isPS        = false;
	  m_tri_bx    = i;
	  m_tri_chip  = m_iter->first;
	  m_tri_lay   = m_tri_chip/1000000;
	  m_tri_lad   = (m_tri_chip-1000000*m_tri_lay)/10000;
	  m_tri_mod   = (m_tri_chip-1000000*m_tri_lay-10000*m_tri_lad)/100;
          
	  trig_sequence.push_back(m_iter->second);
	      
	  if (m_tri_lay<8 || (m_tri_lay>10 && m_tri_lad<9)) isPS = true;
          
	  for (int j=1;j<8;++j)
	  {
	    trig_sequence.push_back((m_data_trig.at(2*trig_seq.at(i+j)+1).find(m_tri_chip))->second);
	  }
              
	  // CIC_TRG_OUT << "\n";
	  // CIC_TRG_OUT << "__________________________________________________\n";
	  // CIC_TRG_OUT << "Layer/Disk  " << m_tri_lay << "\n";
	  // CIC_TRG_OUT << "Ladder/Ring " << m_tri_lad << "\n";
	  // CIC_TRG_OUT << "Module      " << m_tri_mod << "\n";
	  // CIC_TRG_OUT << "CIC Chip    " << (m_tri_chip%100)/8 << "\n";
          
	  CIC_TRG_OUT << m_tri_bx << " / " << m_tri_chip << " -> ";
	  stim_builder::fill_TRG_block(trig_sequence,isPS,1,i);
	}
      }

      if (i%2==0)
      {
	FE_TRG_OUT << "\n";
	FE_TRG_OUT << "##################################################\n";
	FE_TRG_OUT << "#\n";
	FE_TRG_OUT << "# Events  " << i << "/" << i+1 << "\n";
        
	m_chip_trig   = m_data_trig.at(2*trig_seq.at(i));
	
	for ( m_iter = m_chip_trig.begin(); m_iter != m_chip_trig.end();++m_iter )
        {
	  trig_sequence.clear();
          
	  isPS        = false;
	  m_tri_bx    = i;
	  m_tri_chip  = m_iter->first;
	  m_tri_lay   = m_tri_chip/1000000;
	  m_tri_lad   = (m_tri_chip-1000000*m_tri_lay)/10000;
	  m_tri_mod   = (m_tri_chip-1000000*m_tri_lay-10000*m_tri_lad)/100;
              
	  trig_sequence.push_back(m_iter->second);
          
	  if (m_tri_lay<8 || (m_tri_lay>10 && m_tri_lad<9)) isPS = true;
	      
	  if (isPS) // MPA case, data is sent over 2BXs
	  {
	    trig_sequence.push_back((m_data_trig.at(2*trig_seq.at(i+1)).find(m_tri_chip))->second);
            
	    // FE_TRG_OUT << "\n";
	    // FE_TRG_OUT << "__________________________________________________\n";
	    // FE_TRG_OUT << "Layer/Disk  " << m_tri_lay << "\n";
	    // FE_TRG_OUT << "Ladder/Ring " << m_tri_lad << "\n";
	    // FE_TRG_OUT << "Module      " << m_tri_mod << "\n";
	    // FE_TRG_OUT << "Chip        " << m_tri_chip%100 << "\n";
            
	    FE_TRG_OUT << m_tri_bx << " / " << m_tri_chip << " -> ";
	    stim_builder::fill_TRG_block(trig_sequence,isPS,0,i);
	  }
	  else // CBC case, one BX basis
	  {
	    FE_TRG_OUT << m_tri_bx << " / " << m_tri_chip << " -> ";
	    stim_builder::fill_TRG_block(trig_sequence,isPS,0,i);
	    trig_sequence.clear();
	    m_tri_bx    = i+1;
	    trig_sequence.push_back((m_data_trig.at(2*trig_seq.at(i+1)).find(m_tri_chip))->second);
	    FE_TRG_OUT << m_tri_bx << " / " << m_tri_chip << " -> ";
	    stim_builder::fill_TRG_block(trig_sequence,isPS,0,i+1);
	  }
	}
      }

      L1TT->GetEntry(evtID);
      for ( m_iter = m_chip_raw.begin(); m_iter != m_chip_raw.end();++m_iter )
      {
	m_raw_chip = m_iter->first;
	m_digi_list= m_iter->second;
	m_raw_lay  = m_raw_chip/1000000;
	m_raw_lad  = (m_raw_chip-1000000*m_raw_lay)/10000;
	m_raw_mod  = (m_raw_chip-1000000*m_raw_lay-10000*m_raw_lad)/100;
        
	if (m_raw_lay<8 || (m_raw_lay>10 && m_raw_lad<9)) isPS = true;
        
	FE_TRG_IN << "\n";
	FE_TRG_IN << "__________________________________________________\n";
	FE_TRG_IN << "Layer/Disk  " << m_raw_lay << "\n";
	FE_TRG_IN << "Ladder/Ring " << m_raw_lad << "\n";
	FE_TRG_IN << "Module      " << m_raw_mod << "\n";
	FE_TRG_IN << "Chip        " << m_raw_chip%100 << "\n";
       
	// Link the stub/cluster info
	
	if (m_digi_list.size()-1 == 0) continue;
          
	ana_pix(m_raw_lay,m_raw_lad,m_raw_mod,m_digi_list);
	
	// Start the printing loop
	for (unsigned int k=1;k<m_digi_list.size();++k) // Loop over all digis
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
	  if (isPS)
	  {
	    (m_pix_ncol[idx]==32)
	      ? FE_TRG_IN << "pixeldigi " << std::setw(5) << strip << " " << std::setw(3) << seg
	      : FE_TRG_IN << "stripdigi " << std::setw(5) << strip << " " << std::setw(3) << seg;
	  }
	  else
	  {
	    (m_pix_bot[idx]==1)
	      ? FE_TRG_IN << "stripdigi bot " << std::setw(5) << strip << " " << std::setw(3) << seg
	      : FE_TRG_IN << "stripdigi top " << std::setw(5) << strip << " " << std::setw(3) << seg;
	  }
              
	  if (itp==-1) // Unmatched digi
	  {
	    FE_TRG_IN << " // 0   0.00   0.00 // "
		      << std::fixed << std::setprecision(2)
		      << std::setw(6) << atan2(m_pix_y[idx],m_pix_x[idx]) << " "
		      << std::setw(6) << sqrt(m_pix_x[idx]*m_pix_x[idx]+m_pix_y[idx]*m_pix_y[idx]) << " "
		      << std::setw(6) << m_pix_z[idx]
		      << std::fixed << std::setprecision(0)
		      << "\n";
	  }
	  else // Matched digi
	  {
	    if (its==-1)  FE_TRG_IN << " // 0 "; // No stub
	    if (its!=-1 && isGOOD)  FE_TRG_IN << " // 2 "; // In a stub from a good part
	    if (its!=-1 && !isGOOD) FE_TRG_IN << " // 1 "; // In a stub from a not good part
            
	    FE_TRG_IN << std::fixed << std::setprecision(2)
		      << std::setw(6) << ptGEN << " " 
		      << std::setw(6) << d0GEN << " // "  
		      << std::setw(6) << atan2(m_pix_y[idx],m_pix_x[idx]) << " " 
		      << std::setw(6) << sqrt(m_pix_x[idx]*m_pix_x[idx]+m_pix_y[idx]*m_pix_y[idx]) << " " 
		      << std::setw(6) << m_pix_z[idx]
		      << std::fixed << std::setprecision(0)
		      << "\n";
	  }
	}
          
	m_chip_trig   = m_data_trig.at(2*trig_seq.at(i)); // The stubs
	m_digi_list = (m_chip_trig.find(m_raw_chip))->second;
          
	if (m_digi_list.size()-1 == 0) continue;

	do_stub(m_digi_list);
      }
    }
    
    cout << "End of event loop" << endl;
    FE_TRG_IN.close();
    FE_TRG_OUT.close();
    CIC_TRG_OUT.close();
    
    m_outfile->Write();
    
    delete L1TT;
    delete PIX;
    delete MC;
    delete m_outfile;
}




/////////////////////////////////////////////////////////////
//
// Basic methods, initializations,...
//
/////////////////////////////////////////////////////////////

void stim_builder::initVars()
{
  m_tri_bx=0;
  m_raw_bx=0;
  m_tri_chip=0;
  m_raw_chip=0;
  m_tri_size=0;
  m_tri_size_anders=0;
  m_raw_size=0;
  m_raw_mbits=0;
  m_raw_np=0;
  m_raw_ns=0;

  m_tri_lay=0;
  m_tri_lad=0;
  m_tri_mod=0;
  m_tri_nstubs=0;
  m_tri_nstubs_s=0;
  m_tri_nstubs_g=0;
  m_tri_nstubs_gs=0;
  m_raw_lay=0;
  m_raw_lad=0;
  m_raw_mod=0;
  m_raw_FIFO_FULL=0;

  m_tri_data->clear();
  m_raw_data->clear();

  m_raw_chip=0;
}

//
// This method creates, for one module, a list of the pixels used
// to make a cluster and a stub
//
// Three tables are built. One for the digi flagged, one for the corresp.
// clusters, and one for the stubs when applicable

void stim_builder::ana_pix(int lay,int lad,int mod, std::vector<int> digits)
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
        if (m_clus_ladder[i]!=lad+1) continue;
        if (int(m_clus_module[i])!=mod+1) continue;
        
        list_pix.clear();
        list_pix_coords.clear();
        list_pix_coords = m_clus_pix[i]; // List of pixels making the cluster
        
        // Try to match them with the list of digits effectively recorded in the module
        for (unsigned int j=0;j<list_pix_coords.size()/2;++j)
        {
            row=list_pix_coords.at(2*j);
            col=list_pix_coords.at(2*j+1);
            
            found=false;

            for (unsigned int k=1;k<digits.size();++k)
            {
                if (found) continue;
                
                idx = digits.at(k);
                
                if (m_pix_module[idx]!=m_clus_module[i]) continue;
                if (m_pix_row[idx]!=row)                 continue;
                if (m_pix_col[idx]!=col)                 continue;
                
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
        if (m_stub_ladder[i]!=lad+1) continue;
        if (m_stub_module[i]!=mod+1) continue;
        
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

void stim_builder::do_stub(std::vector<int> stubs)
{
    int idx;
    
    for (unsigned int i=1;i<stubs.size();++i)
    {
        idx = stubs.at(i);

        // First of all we compute the ID of the stub's module
        
        FE_TRG_IN << "stub: "
            << std::setw(6) << std::fixed << std::setprecision(1) << m_stub_strip[idx] << " "
            << std::setw(5) << m_stub_deltas[idx]
            << " " << std::setprecision(2)
            << std::setw(7) <<  m_stub_seg[idx] << "\n";
    }
}

void stim_builder::initTuple(std::string inRAW,std::string inTRG,std::string out)
{
    m_outfile    = new TFile(out.c_str(),"recreate");

    m_tri_tree    = new TTree("Trigger_FE","L1Trigger words after FE");
    m_raw_tree    = new TTree("Raw_FE","Raw data words after FE");
    m_raw_summary = new TTree("Raw_SUM","Raw data summary info");

    m_tri_tree->Branch("TRI_BX",         &m_tri_bx,      "TRI_BX/I");
    m_tri_tree->Branch("TRI_CHP",        &m_tri_chip,    "TRI_CHP/I");
    m_tri_tree->Branch("TRI_LAY",        &m_tri_lay,     "TRI_LAY/I");
    m_tri_tree->Branch("TRI_LAD",        &m_tri_lad,     "TRI_LAD/I");
    m_tri_tree->Branch("TRI_MOD",        &m_tri_mod,     "TRI_MOD/I");
    m_tri_tree->Branch("TRI_WORD",       &m_tri_data);
    m_tri_tree->Branch("TRI_SIZE",       &m_tri_size,    "TRI_SIZE/I");
    m_tri_tree->Branch("TRI_NSTUBS",     &m_tri_nstubs);
    m_tri_tree->Branch("TRI_NSTUBS_S",   &m_tri_nstubs_s);
    m_tri_tree->Branch("TRI_NSTUBS_G",   &m_tri_nstubs_g);
    m_tri_tree->Branch("TRI_NSTUBS_GS",  &m_tri_nstubs_gs);

    m_raw_tree->Branch("RAW_BX",         &m_raw_bx,      "RAW_BX/I");
    m_raw_tree->Branch("RAW_CHP",        &m_raw_chip,    "RAW_CHP/I");
    m_raw_tree->Branch("RAW_LAY",        &m_raw_lay,     "RAW_LAY/I");
    m_raw_tree->Branch("RAW_LAD",        &m_raw_lad,     "RAW_LAD/I");
    m_raw_tree->Branch("RAW_MOD",        &m_raw_mod,     "RAW_MOD/I");
    m_raw_tree->Branch("RAW_WORD",       &m_raw_data);
    m_raw_tree->Branch("RAW_SIZE",       &m_raw_size,    "RAW_SIZE/I");
    m_raw_tree->Branch("RAW_MBITS",      &m_raw_mbits,   "RAW_MBITS/I");

    m_raw_tree->Branch("RAW_NPCLUS",     &m_raw_np);
    m_raw_tree->Branch("RAW_NSCLUS",     &m_raw_ns);

    L1TT   = new TChain("TkStubs");
    PIX    = new TChain("Pixels");
    MC     = new TChain("MC");
    
    std::size_t found = inRAW.find(".root");
  
    // Case 1, it's a root file
    if (found!=std::string::npos)
    {
        PIX->Add(inRAW.c_str());
        L1TT->Add(inRAW.c_str());
        MC->Add(inRAW.c_str());
    }
    else // This is a list provided into a text file
    {
        std::string STRING;
        std::ifstream in2(inRAW.c_str());
        if (!in2)
        {
            std::cout << "Please provide a valid RAW data filename list" << std::endl;
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

    m_PHYsize = L1TT->GetEntries();

    found = inTRG.find(".root");
  
    // Case 1, it's a root file
    if (found!=std::string::npos)
    {
        PIX->Add(inTRG.c_str());
        L1TT->Add(inTRG.c_str());
        MC->Add(inTRG.c_str());
    }
    else // This is a list provided into a text file
    {
        std::string STRING;
        std::ifstream in2(inTRG.c_str());
        if (!in2)
        {
            std::cout << "Please provide a valid TRG data filename list" << std::endl;
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
    pm_pix_ncol=&m_pix_ncol;
    pm_pix_bot=&m_pix_bot;
    pm_pix_x=&m_pix_x;
    pm_pix_y=&m_pix_y;
    pm_pix_z=&m_pix_z;

    PIX->SetBranchAddress("PIX_n",         &m_pix);
    PIX->SetBranchAddress("PIX_nPU",       &m_npu);
    PIX->SetBranchAddress("PIX_layer",     &pm_pix_layer);
    PIX->SetBranchAddress("PIX_ladder",    &pm_pix_ladder);
    PIX->SetBranchAddress("PIX_module",    &pm_pix_module);
    PIX->SetBranchAddress("PIX_row",       &pm_pix_row);
    PIX->SetBranchAddress("PIX_column",    &pm_pix_col);
    PIX->SetBranchAddress("PIX_ncolumn",   &pm_pix_ncol);
    PIX->SetBranchAddress("PIX_bottom",    &pm_pix_bot);
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
    pm_clus_bot=&m_clus_bot;
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
    //    L1TT->SetBranchAddress("L1TkCLUS_MULT",      &pm_clus_mult);
    L1TT->SetBranchAddress("L1TkCLUS_tp",        &pm_clus_tp);
    L1TT->SetBranchAddress("L1TkCLUS_layer",     &pm_clus_layer);
    L1TT->SetBranchAddress("L1TkCLUS_ladder",    &pm_clus_ladder);
    L1TT->SetBranchAddress("L1TkCLUS_module",    &pm_clus_module);
    L1TT->SetBranchAddress("L1TkCLUS_bottom",    &pm_clus_bot);    

    L1TT->SetBranchAddress("L1TkSTUB_pxGEN",     &pm_stub_pxGEN);
    L1TT->SetBranchAddress("L1TkSTUB_pyGEN",     &pm_stub_pyGEN);
    L1TT->SetBranchAddress("L1TkSTUB_etaGEN",    &pm_stub_etaGEN);
    L1TT->SetBranchAddress("L1TkSTUB_X0",        &pm_stub_X0);
    L1TT->SetBranchAddress("L1TkSTUB_Y0",        &pm_stub_Y0);
    L1TT->SetBranchAddress("L1TkSTUB_Z0",        &pm_stub_Z0);

}





/////////////////////////////////////////////////////////////////////////////////
//
// ==> convert(std::string sectorfilename) 
//
// Here we retrieve info from the TKLayout CSV file containing the sector definition
//
// This file contains, for each sector, the ids of the modules and chips contained in the sector sec_num
//
// The role of this method is to create the opposite, ie a vector containing, for every module the list of sectors belonging to it
//
/////////////////////////////////////////////////////////////////////////////////

bool stim_builder::convert(std::string sectorfilename)
{
  std::vector<int> module;

  m_modules.clear();
  m_chips.clear();
  m_concs.clear();

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

  int m_sec_mult = 0;
  int npar = 0;

  int mlay, mlad, mmod;

  while (!in.eof()) 
  {
    ++m_sec_mult;
    getline(in,STRING);

    if (m_sec_mult<2) continue;
    if (m_sec_mult-2!=m_tower && m_tower!=-1) continue;

    std::istringstream ss(STRING);
    npar = 0;
    while (ss)
    {
      std::string s;
      if (!getline( ss, s, ',' )) break;

      ++npar;
      if (npar<=2) continue;

      mlay = int(atoi(s.c_str())/10000);
      mlad = int(atoi(s.c_str())-10000*mlay)/100;
      mmod = int(atoi(s.c_str())-10000*mlay-100*mlad);

      if (m_tower==-1)
      {
          if (mlay!=m_lay && m_lay!=-1) continue;
          if (mlad!=m_lad && m_lad!=-1) continue;
          if (mmod!=m_mod && m_mod!=-1) continue;
      }
        
      m_modules.at(atoi(s.c_str())).push_back(m_sec_mult-2);
    }
  }

  in.close();

  for (int i=0;i<230000;++i)
  {
    if (m_modules.at(i).size()<=1) continue;

    for (int j=0;j<16;++j) m_chips.push_back(100*i+j);
    m_concs.push_back(100*i);
    m_concs.push_back(100*i+8);
  }

  return true;
}

//
// List of method writing the data blocks, according to the format defined in
// the following document:
//
// https://espace.cern.ch/Tracker-Upgrade/Electronics/CIC/Shared%20Documents/Data%20formats/CIC_IO_Formats_v2.pdf
//

//
// Trigger block
//

// stubs contains the idx of the stubs belonging to the chip (in the tree)


void stim_builder::fill_TRG_block(std::vector<std::vector<int> > stubs, bool ps, bool conc, int BXid)
{
    m_tri_data->clear();

    // Then go for the data

    int nstubs[8]; // Number of stubs stored in the chips (need 8 for the CIC case)
    int ovflow[8]; // Chip in overflow or not (need 8 for the CIC case)
    int chip_id=-1;
    
    int ovflow_CIC = 0;
    
    int seqsize = static_cast<int>(stubs.size());
    std::vector<int> stublist;
    
    std::vector<float> stub_register;
    
    int nstubs_tot = 0;
    
    int lim_stub_FE  = 0;
    int lim_stub_CIC = 0;
    
    (ps)
    ? lim_stub_FE = 5
    : lim_stub_FE = 3;
    
    (ps)
    ? lim_stub_CIC = static_cast<int>((m_TRGsize-26)/21)
    : lim_stub_CIC = static_cast<int>((m_TRGsize-26)/18);
    
    if (conc)
    {
        stub_register.clear();
        
        (ps)
        ? stim_builder::fill_CONC_TRG_header(BXid%3564,1)
        : stim_builder::fill_CONC_TRG_header(BXid%3564,0);
    }
    else
    {
        // Start with synchro bit only for the PS

        if (ps) // For the PS one puts the number of stubs in the first BX
        {          // will be completed at the end
            m_tri_data->push_back(1); // Always start with synchro bit 1
            for (int j=0;j<3;++j) m_tri_data->push_back(0);
        }
    }
  
    int bit_bx = 3; // Number of bits necessary to precise the offset number in the CIC word
 
    int lay,lad;
    
    for (int j=0;j<8;++j) ovflow[j]=0;
    
    
    for (int i=0;i<seqsize;++i) // Loop over the block of events
    {
        // Here we initialize the stub per chip counters

        if (ps && i%2==0) // For the MPA, it's every 2BXs
        {
            for (int j=0;j<8;++j) nstubs[j]=0;
        }

        if (!ps) // For the CBC, it's every BX
        {
            for (int j=0;j<8;++j) nstubs[j]=0;
        }
        
        stublist = stubs.at(i); // Get the stubs indexes for the chip and evt i
        if (stublist.size()==1) continue;
        
        // Then loop over the stubs contained in the event

        int idx;
        int strip,chip,seg;
        int evtID = stublist.at(0);
        
        L1TT->GetEntry(evtID);
        
        for (unsigned int kk=1;kk<stublist.size();++kk)
        {
            idx    = stublist.at(kk);
            seg    = m_stub_seg[idx];
            lay    = m_stub_layer[idx];
            lad    = m_stub_ladder[idx];

            if (ps)
            {
                chip  =  m_stub_chip[idx];
                strip = int(2*m_stub_strip[idx])%240;  // Bet 0 and 255
            }
            else
            {
                chip  =  m_stub_chip[idx];
                strip = int(2*m_stub_strip[idx])%254+1;
            }
        
            // Here we apply the chip limits
 
            if (nstubs[chip]>=lim_stub_FE)
            {
                ovflow[chip]=1;
                continue; // No CHIP sorting for the moment
            }
            
            if (nstubs_tot>=lim_stub_CIC) ovflow_CIC=1;
            
            ++nstubs[chip];
            ++nstubs_tot;
            
            chip_id = chip;
            
            // For the CIC case we build a register first, in order to sort the stubs afterwards
            
            if (conc)
            {
	      //	      cout << i << " // " << chip << " // " << nstubs_tot << endl;
                stub_register.push_back(i);
                stub_register.push_back(chip);
                stub_register.push_back(strip);
                stub_register.push_back(seg%16);
                stub_register.push_back(m_stub_deltas[idx]);
                continue;
            }
            
            std::bitset<8> pos = strip;                     // Stub position
            std::bitset<4> col = seg%16;                    // Z position (for PS module)

            // Then the position
            for (int j=0;j<8;++j)
            {
                if (!conc && m_tri_data->size()==40 && ps) m_tri_data->push_back(0); // Second synchro bit
                (ps)
                ? m_tri_data->push_back(pos[7-j])
                : m_tri_data->push_back(pos[j]); // LSB first in the cbc
            }

            // The bend

            if (ps)
            {
                std::bitset<3> off = convert_bend_PS(lay,lad,m_stub_deltas[idx]);
                if (m_tri_data->size()==40 && !conc) m_tri_data->push_back(0); // Second synchro bit
                for (int j=0;j<bend_bit_MPA;++j) m_tri_data->push_back(off[2-j]);
            }
            else
            {
                std::bitset<4> off = convert_bend_2S(lay,lad,m_stub_deltas[idx]);
                for (int j=0;j<bend_bit_CBC;++j) m_tri_data->push_back(off[j]);
            }
            
            // Finally the column for MPA side
            if (ps)
            {
                for (int j=0;j<4;++j)
                {
                    if (!conc && m_tri_data->size()==40 && ps) m_tri_data->push_back(0); // Second synchro bit
                    m_tri_data->push_back(col[3-j]);
                }
            }
        }
        
        if (ps && i%2==0 && !conc) // Here we have the number of stubs in the first BX for the MPA
        {
            std::bitset<3> cnt= std::min(nstubs[chip_id],lim_stub_FE);
            for (int j=0;j<3;++j) m_tri_data->at(j+1)=cnt[2-j];
        }
    }
    
    // Here we perform the reordering
    //
    // This is a bit awkward
    //
    // Priority is given to the lowest bend, then in case of equal bends:
    //
    // Priority given to higher BXs, in groups of 2Bxs. Means that 6/7 will arrive
    // before 0/1
    //
    // Then in the same BXs group, priority is given to higher chip numbers, 7 will arrive before 4, for example
    //
    // Finally, in case of equal chip numbers, priority is given to higher strip number
    //
    // There is no sorting on the Z position
    
    
    int CBC_word[40];
    int idx=-1;
    float min_bend, min_off, min_add, min_chp,tp1, tp2, tp3, tp4, tp5;
    
    if (stub_register.size()!=0 && conc)
    {
        for (unsigned int i=0;i<stub_register.size()/5;++i)
        {
	  //	  cout << idx << " /÷/ " << i << endl;
            
	  if (int(i)>=lim_stub_CIC) continue;

	  idx      = i;
	  min_bend = fabs(stub_register.at(5*i+4));
	  min_off  = int(stub_register.at(5*i)/2);
	  min_chp  = stub_register.at(5*i+1);
	  min_add  = stub_register.at(5*i+2);
          
	  if (i!=stub_register.size()/5-1)
          {
	    for (unsigned int j=i+1;j<stub_register.size()/5;++j)
	    {
	      if (min_bend>fabs(stub_register.at(5*j+4))) // lowest bend
	      {
		min_bend = fabs(stub_register.at(5*j+4));
		min_off  = int(stub_register.at(5*j)/2);
		min_chp  = stub_register.at(5*j+1);
		min_add  = stub_register.at(5*j+2);
		idx=j;
	      }
	      else if (min_bend==fabs(stub_register.at(5*j+4)) && min_off<int(stub_register.at(5*j)/2)) // highest offset
	      {
		min_off  = int(stub_register.at(5*j)/2);
		min_chp  = stub_register.at(5*j+1);
		min_add  = stub_register.at(5*j+2);
		idx=j;
	      }
	      else if (min_bend==fabs(stub_register.at(5*j+4))
		       && min_off==int(stub_register.at(5*j)/2)
		       && min_chp<stub_register.at(5*j+1)
		       ) // highest offset
	      {
		min_chp  = stub_register.at(5*j+1);
		min_add  = stub_register.at(5*j+2);
		idx=j;
	      }
	      else if (min_bend==fabs(stub_register.at(5*j+4))
		       && min_off==int(stub_register.at(5*j)/2)
		       && min_chp==stub_register.at(5*j+1)
		       && min_add<stub_register.at(5*j+2)
		       ) // highest offset
	      {
		min_add=stub_register.at(5*j+2);
		idx=j;
	      }          
	    }
	  }
	  
	  //	  cout << idx << " /÷ " << i << endl;
            
	  // Here we exchange row i and idx, if needed;
            
	  if (idx!=int(i))
          {
	    tp1 = stub_register.at(5*i);
	    tp2 = stub_register.at(5*i+1);
	    tp3 = stub_register.at(5*i+2);
	    tp4 = stub_register.at(5*i+3);
	    tp5 = stub_register.at(5*i+4);
            
	    stub_register.at(5*i)   = stub_register.at(5*idx);
	    stub_register.at(5*i+1) = stub_register.at(5*idx+1);
	    stub_register.at(5*i+2) = stub_register.at(5*idx+2);
	    stub_register.at(5*i+3) = stub_register.at(5*idx+3);
	    stub_register.at(5*i+4) = stub_register.at(5*idx+4);
            
	    stub_register.at(5*idx)   = tp1;
	    stub_register.at(5*idx+1) = tp2;
	    stub_register.at(5*idx+2) = tp3;
	    stub_register.at(5*idx+3) = tp4;
	    stub_register.at(5*idx+4) = tp5;
	  }

	  //	  cout << i << " / " << stub_register.at(5*i) << " / " << stub_register.at(5*i+1) << endl;
          
	  std::bitset<3> bx  = stub_register.at(5*i);          // Offset in the sequence
	  std::bitset<3> chp = stub_register.at(5*i+1);        // Chip number
	  std::bitset<8> pos = stub_register.at(5*i+2);        // Stub position
	  std::bitset<4> col = stub_register.at(5*i+3);        // Z position (for PS module)
          
          
	  for (int j=0;j<3;++j) m_tri_data->push_back(bx[bit_bx-1-j]);
	  for (int j=0;j<3;++j) m_tri_data->push_back(chp[2-j]);
	  for (int j=0;j<8;++j) m_tri_data->push_back(pos[7-j]);
          
	  if (ps)
	  {
	    std::bitset<3> off = convert_bend_PS(lay,lad,stub_register.at(5*i+4));
	    for (int j=0;j<bend_bit_MPA;++j) m_tri_data->push_back(off[2-j]);
	  }
	  else
	  {
	    std::bitset<4> off = convert_bend_2S(lay,lad,stub_register.at(5*i+4));
	    for (int j=0;j<bend_bit_CBC;++j) m_tri_data->push_back(off[3-j]);
	  }
            
	  // Finally the column for MPA side
	  if (ps)
          {
	    for (int j=0;j<4;++j) m_tri_data->push_back(col[3-j]);
	  }
        }
    }
    
    
    // Potentially put a trailer here, only for the FE chips

    if (!conc)
    {
        if (!ps)
        {
            for (int j=m_tri_data->size();j<36;++j) m_tri_data->push_back(0);

            (chip_id>=0 && ovflow[chip_id]==1)         // BIT37
            ? m_tri_data->push_back(1)
            : m_tri_data->push_back(0);
      
            m_tri_data->push_back(0);  // OR254 (38)
            m_tri_data->push_back(0);  // Error (39)
            m_tri_data->push_back(1);  // Synchro bit (40)

            /// For the moment stubs are coded as A1B2A2B2A3B3, code them as A1A2A3B1B2B3
            
            for (int j=0;j<40;++j) CBC_word[j] = m_tri_data->at(j);
            
            for (int j=8;j<16;++j)  CBC_word[j] = m_tri_data->at(j+4); // A2 after A1
            for (int j=16;j<24;++j) CBC_word[j] = m_tri_data->at(j+8); // A3 after A2
            
            for (int j=24;j<28;++j) CBC_word[j] = m_tri_data->at(j-16); // B1 after A3
            for (int j=28;j<32;++j) CBC_word[j] = m_tri_data->at(j-8);  // B2 after B1
            
            for (int j=0;j<40;++j) m_tri_data->at(j) = CBC_word[j];
        }
        else
        {
            for (int j=m_tri_data->size();j<80;++j) m_tri_data->push_back(0);
        }
    }

    if (conc)
    {
        for (int j=0;j<8;++j)
        {
            if (ovflow[j]==1) m_tri_data->at(j+1) = 1; // FE chip ovflow bit in CIC header
        }
      
        if (ovflow_CIC==1)
        {
            m_tri_data->at(9) = 1; // Raise the CIC error bit if the number of stubs is in overflow
            std::bitset<4> nst = lim_stub_CIC;
            for (int j=0;j<4;++j) m_tri_data->at(22+j) = nst[3-j];
        }
        else
        {
            std::bitset<4> nst = nstubs_tot;
            for (int j=0;j<4;++j) m_tri_data->at(22+j) = nst[3-j];
        }

        for (int j=m_tri_data->size();j<256;++j)
        {
            m_tri_data->push_back(0); // padding
        }
    }

  m_tri_size=m_tri_data->size();
  m_tri_tree->Fill();

  for (unsigned int j=0;j<m_tri_data->size();++j)
  {
      if (!conc) FE_TRG_OUT << m_tri_data->at(j);
      if (conc) CIC_TRG_OUT << m_tri_data->at(j);
  }
 
  (conc)
  ?  CIC_TRG_OUT << "\n"
  :  FE_TRG_OUT << "\n";

}


void stim_builder::fill_CONC_TRG_header(int BXid,int MPA)
{
  // Format of the CONC TRIGGER word header
  //
  // CSSSSSSSSSCCCCCCCCCCCC : SS..SS (status) CC..CC (BX ID bet 0 and 3564)
  //

  m_tri_data->push_back(MPA); 

  for (int j=0;j<9;++j) m_tri_data->push_back(0); // Status bits

  std::bitset<12> BX_ID = BXid;

  for (int j=0;j<12;++j) m_tri_data->push_back(BX_ID[11-j]); // CC..CC

  for (int j=0;j<4;++j) m_tri_data->push_back(0); // Let room for the number of stubs
}


std::bitset<3> stim_builder::convert_bend_PS(int layer,int ring, float bend)
{
    // Format defined in these talks:
    //
    //  * Implementation of reduced bits for the bend information.
    // https://indico.cern.ch/event/350910/session/4/contribution/12/material/slides/0.pdf
    // https://indico.cern.ch/event/361219/session/7/contribution/17/material/slides/5.pdf
    // 3 bits for PS, 4 bits for 2S.
    //
    
    int code=0;
    
    if (layer<=10) // barrel
    {
        switch( layer )
        {
            case 5:
                if (fabs(bend)<=0.5) code=0;
                if (bend==-1)        code=5;
                if (bend==-1.5)      code=6;
                if (bend<=-2)        code=7;
                if (bend==1)         code=1;
                if (bend==1.5)       code=2;
                if (bend>=2)         code=3;
                break;
            case 6:
                if (fabs(bend)<=0.5) code=0;
                if (bend==-1)        code=5;
                if (bend==-1.5)      code=6;
                if (bend<=-2)        code=7;
                if (bend==1)         code=1;
                if (bend==1.5)       code=2;
                if (bend>=2)         code=3;
                break;
            case 7:
                if (fabs(bend)<=0.5)        code=0;
                if (bend==-1)               code=5;
                if (bend==-1.5 || bend==-2) code=6;
                if (bend<=-2.5)             code=7;
                if (bend==1)                code=1;
                if (bend==1.5 || bend==2)   code=2;
                if (bend<=2.5)              code=3;
                break;
            default:
                break;
        }
    }
    else
    {
        if (ring<=3)
        {
            if (fabs(bend)<=0.5)        code=0;
            if (bend==-1)               code=5;
            if (bend==-1.5 || bend==-2) code=6;
            if (bend<=-2.5)             code=7;
            if (bend==1)                code=1;
            if (bend==1.5 || bend==2)   code=2;
            if (bend<=2.5)              code=3;
        }
        else if (ring<=6)
        {
            if (fabs(bend)<=1.)  code=0;
            if (bend==-1.5)      code=5;
            if (bend==-2.)       code=6;
            if (bend<=-2.5)      code=7;
            if (bend==1.5)       code=1;
            if (bend==2.)        code=2;
            if (bend>=2.5)       code=3;
        }
        else if (ring==8)
        {
            if (fabs(bend)<=1.)         code=0;
            if (bend==-1.5 || bend==-2) code=5;
            if (bend==-2.5)             code=6;
            if (bend<=-3.)              code=7;
            if (bend==1.5 || bend==2)   code=1;
            if (bend==2.5)              code=2;
            if (bend<=3.)               code=3;
        }
        else
        {
            if (fabs(bend)<=1.)         code=0;
            if (bend==-1.5 || bend==-2) code=5;
            if (bend==-2.5 || bend==-3) code=6;
            if (bend<=-3.5)             code=7;
            if (bend==1.5 || bend==2)   code=1;
            if (bend==2.5 || bend==3)   code=2;
            if (bend>=3.5)              code=3;
        }
    }
    
    std::bitset<3> bendcode = code;
    
    return bendcode;
}

std::bitset<4> stim_builder::convert_bend_2S(int layer,int ring, float bend)
{
    // Format defined in these talks:
    //
    //  * Implementation of reduced bits for the bend information.
    // https://indico.cern.ch/event/350910/session/4/contribution/12/material/slides/0.pdf
    // https://indico.cern.ch/event/361219/session/7/contribution/17/material/slides/5.pdf
    // 3 bits for PS, 4 bits for 2S.
    //
    
    int code=0;
    
    if (layer<=10) // barrel
    {
        switch( layer )
        {
            case 8:
                if (fabs(bend)<=0.5) code=0;
                if (bend==-1)        code=9;
                if (bend==-1.5)      code=10;
                if (bend==-2)        code=11;
                if (bend==-2.5)      code=12;
                if (bend==-3)        code=13;
                if (bend==-3.5)      code=14;
                if (bend<=-4)        code=15;
                if (bend==1)         code=1;
                if (bend==1.5)       code=2;
                if (bend==2)         code=3;
                if (bend==2.5)       code=4;
                if (bend==3)         code=5;
                if (bend==3.5)       code=6;
                if (bend>=4)         code=7;
                break;
            case 9:
                if (fabs(bend)<=0.5) code=0;
                if (bend==-1)        code=9;
                if (bend==-1.5)      code=10;
                if (bend==-2)        code=10;
                if (bend==-2.5)      code=10;
                if (bend==-3)        code=11;
                if (bend==-3.5)      code=11;
                if (bend==-4)        code=12;
                if (bend==-4.5)      code=13;
                if (bend==-5)        code=14;
                if (bend<=-5.5)      code=15;
                if (bend==1)         code=1;
                if (bend==1.5)       code=2;
                if (bend==2)         code=2;
                if (bend==2.5)       code=2;
                if (bend==3)         code=3;
                if (bend==3.5)       code=3;
                if (bend==4)         code=4;
                if (bend==4.5)       code=5;
                if (bend==5)         code=6;
                if (bend>=5.5)       code=7;
                break;
            case 10:
                if (fabs(bend)<=0.5) code=0;
                if (bend==-1)        code=9;
                if (bend==-1.5)      code=10;
                if (bend==-2)        code=10;
                if (bend==-2.5)      code=10;
                if (bend==-3)        code=11;
                if (bend==-3.5)      code=11;
                if (bend==-4)        code=11;
                if (bend==-4.5)      code=12;
                if (bend==-5)        code=12;
                if (bend==-5.5)      code=13;
                if (bend==-6)        code=14;
                if (bend<=-6.5)      code=15;
                if (bend==1)         code=1;
                if (bend==1.5)       code=2;
                if (bend==2)         code=2;
                if (bend==2.5)       code=2;
                if (bend==3)         code=3;
                if (bend==3.5)       code=3;
                if (bend==4)         code=3;
                if (bend==4.5)       code=4;
                if (bend==5)         code=4;
                if (bend==5.5)       code=5;
                if (bend==6)         code=6;
                if (bend>=6.5)       code=7;
                break;
            default:
                break;
        }
    }
    else
    {
        if (fabs(bend)<=0.5) code=0;
        if (bend==-1)        code=9;
        if (bend==-1.5)      code=10;
        if (bend==-2)        code=10;
        if (bend==-2.5)      code=11;
        if (bend==-3)        code=11;
        if (bend==-3.5)      code=12;
        if (bend==-4)        code=12;
        if (bend==-4.5)      code=13;
        if (bend==-5)        code=14;
        if (bend<=-5.5)      code=15;
        if (bend==1)         code=1;
        if (bend==1.5)       code=2;
        if (bend==2)         code=2;
        if (bend==2.5)       code=3;
        if (bend==3)         code=3;
        if (bend==3.5)       code=4;
        if (bend==4)         code=4;
        if (bend==4.5)       code=5;
        if (bend==5)         code=6;
        if (bend>=5.5)       code=7;
    }
    
    std::bitset<4> bendcode = code;
    
    return bendcode;
}

