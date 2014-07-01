// Class for the rates determination
// For more info, look at the header file

#include "evtbuilder.h"

// Main constructor

// Constructor for the concentrator

// Meaning of the parameters
//
//
// filename : the input data files
// outfile  : name of the output root file
// npatt    : how many BXs should be produced (will be fitted correctly afterwards)
// rate     : input L1A rate (for the L1 raw block writing
// layer    : restrict info to a given layer, or not (give -1 in this case)
// sector   : the input file defining the modules
// RAW      : write the L1 block or not
// TRIG     : write the Trigger block or not
// npblock  : How many bits do we reserve to L1 in the concentrator block (between 8 and 320, should be a multiple of 8) 
// BMPA     : How many bits are used to code the stub bend on MPA concentrator level (0 to 5)
// BCBC     : How many bits are used to code the stub bend on CBC concentrator level (0 to 5)
// conc     : Are we building concentrator (true) or FE (false) sequences
// L1prop   : The proportion of PU4T events you want in the TRG flow (in %)


evtbuilder::evtbuilder(std::string filenameRAW, std::string filenameTRG, std::string outfile, int npatt, int rate, int layer, std::string sector, bool RAW, bool TRIG, int npblock, int BMPA, int BCBC,bool conc,int L1prop)
{
  m_rate       = rate;
  m_lay        = layer;
  m_npblock    = npblock;
  m_tri_data   = new  std::vector<int>;
  m_raw_data   = new  std::vector<int>;
  m_L1prop     = L1prop; 

  m_write_raw = RAW;
  m_write_trg = TRIG;


  bend_bit_MPA = BMPA;
  bend_bit_CBC = BCBC;


  evtbuilder::initVars();
  evtbuilder::convert(sector);
  evtbuilder::initTuple(filenameRAW,filenameTRG,outfile);
  evtbuilder::get_stores(npatt-npatt%8,conc);   
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

void evtbuilder::get_stores(int nevts, bool conc)
{
 
  int B_id; // The detector module IDs (defined in the header)

  int layer,ladder,module,strip,chip,seg,pt;

  int n_entries_TRG = L1TT->GetEntries();
  int n_entries_RAW = PIX->GetEntries();

  bool isPS=false;
  //  bool isCONC=true;

  int goodstub = -1;

  float pTgen=0.;
  float d0gen=0.;

  // Then loop over events

  std::vector<int> m_stub_list;
  std::vector<int> m_digi_list;

  m_data_trig.clear();
  m_data_raw.clear();

  //
  // We make a first loop over all entries available
  // in order to data stores, using the digis (L1raw block)
  // and the stubs (trigger block)
  // 
  // 2 data samples are defined for that, one with 50%/50%  minbias/physics events, 
  // for the trigger words, and one with pure physics events, for L1 data transmission 
  //
  // You can choose the proportion of complex events you want in the TRG block, in % 
  //

  int mixpar = 0;


  cout << "--> Entering loop 1, producing the big data stores for " 
       << std::min(n_entries_TRG,n_entries_RAW) << " events..." << endl;

  //  for (int j=0;j<std::min(n_entries_TRG,n_entries_RAW);++j)
  for (int j=0;j<200;++j)
  {    
    if (j%20==0)
      cout << "Processing event " <<  j << "/" << std::min(n_entries_TRG,n_entries_RAW) << endl;

    m_chip_trig.clear();
    m_chip_raw.clear();

    if (rand()%100>m_L1prop) // Mbias
    {
      mixpar=rand()%int(n_entries_TRG/2.); 
    }
    else // Phys
    {
      mixpar=rand()%int(n_entries_TRG/2.)+n_entries_TRG/2.; 
    }

    if (m_write_trg) L1TT->GetEntry(mixpar); 
    if (m_write_raw) PIX->GetEntry(j); 

    // Initialize the map with dummy values for all referenced modules
    //
    // Please keep in mind that this map is 8 times larger in the FE output case

    if (!conc)
    {
      for (unsigned int i=0;i<m_chips.size();++i) // All the chips if one asks the FE output
      {
	m_digi_list.clear();
	m_digi_list.push_back(-1);
	m_chip_raw.insert(std::make_pair(m_chips.at(i),m_digi_list));
	m_chip_trig.insert(std::make_pair(m_chips.at(i),m_digi_list));
      }
    }
    else
    {
      for (unsigned int i=0;i<m_concs.size();++i) // All the concentrators if one asks the CONC output
      {
	m_digi_list.clear();
	m_digi_list.push_back(-1);
	m_chip_raw.insert(std::make_pair(m_concs.at(i),m_digi_list));
	m_chip_trig.insert(std::make_pair(m_concs.at(i),m_digi_list));
      }
    }

    // First loop over the digis (raw block)
    if (m_write_raw)
    {
      for (int i=0;i<m_pix;++i)
      {
	isPS  = false;
	layer = m_pix_layer[i]; 
	
	if (layer!=m_lay && m_lay!=-1) continue;
	
	ladder= m_pix_ladder[i]-1;
	if (layer<8 || (layer>10 && ladder<9)) isPS = true; 
	
	module= static_cast<int>((m_pix_module[i]-1)/2); 
	seg   = m_pix_col[i]; 
	strip = m_pix_row[i]; 
	
	if (isPS && m_pix_module[i]%2==1)
	{
	  chip  = static_cast<int>(strip/127)+(seg/16)*8;      
	}
	else
	{
	  chip  = static_cast<int>(strip/127)+seg*8;
	}

	strip = strip%128+((m_pix_module[i]-1)%2)*128;

	B_id = layer*1000000 + ladder*10000 + module*100 + chip;

	if (conc) B_id = B_id - chip%8; // Here we get the concentrator ID
	
	// Look if this chip has already been touched in this event
	m_iter = m_chip_raw.find(B_id);
	
	if (m_iter == m_chip_raw.end()) // Unknown concentrator???
	{
	  continue;
	  cout << "Wow, this should not happen, the chip ref " << B_id << " is unknown!!!" << endl;
	}
	else // Otherwise complement
	{
	  m_digi_list.clear();
	  m_digi_list = m_iter->second;
	  m_digi_list.push_back(chip%8);
	  m_digi_list.push_back(strip);
	  
	  (isPS) 
	    ? m_digi_list.push_back(seg%16)
	    : m_digi_list.push_back(-1);

	  m_chip_raw.erase(m_iter->first);
	  m_chip_raw.insert(std::make_pair(B_id,m_digi_list));
	}
      } // End of digi collection
    }

    // Start the stub loop

    if (m_write_trg)
    {
      for (int i=0;i<m_stub;++i)
      {  
	// First of all we compute the ID of the stub's module
	isPS  = false;
	goodstub=-1;
	layer = m_stub_layer[i]; 
	
	if (layer!=m_lay && m_lay!=-1) continue;
	
	ladder= m_stub_ladder[i]; 
	
	if (layer<8 || (layer>10 && ladder<9)) isPS = true; 
	
	module= m_stub_module[i];
	seg   = m_stub_seg[i]; 
		
	if (isPS)
	{
	  chip  =  m_stub_chip[i]+(seg/16)*8;      
	}
	else
	{
	  chip  =  m_stub_chip[i]+seg*8;
	}

	if (isPS)
	{
	  strip = int(2*m_stub_strip[i])%256;  // Bet 0 and 255
	}
	else
	{
	  strip = int(2*m_stub_strip[i])%254+1; // Code between 1 and 254 to avoid 00000000 position
	}

	pTgen = sqrt(m_stub_pxGEN[i]*m_stub_pxGEN[i]+m_stub_pyGEN[i]*m_stub_pyGEN[i]);
	d0gen = sqrt(m_stub_X0[i]*m_stub_X0[i]+m_stub_Y0[i]*m_stub_Y0[i]);

	if (pTgen>2 && d0gen<0.5) goodstub=1;

	pt   = static_cast<int>(abs(2*m_stub_deltas[i]));
	B_id = layer*1000000 + ladder*10000 + module*100 + chip;
	
	if (conc) B_id = B_id - chip%8; // Here we get the concentrator ID
	
	// Look if this chip has already been touched in this event
	m_iter = m_chip_trig.find(B_id);
	
	if (m_iter == m_chip_trig.end()) // Unknown chip???
	{
	  continue;
	  cout << "Wow, this should not happen, the chip ref " << B_id << " is unknown!!!" << endl;
	}
	else // Otherwise complement
	{
	  m_stub_list.clear();
	  m_stub_list = m_iter->second;
	  m_stub_list.push_back(chip%8);
	  m_stub_list.push_back(strip);

	  (isPS) 
	    ? m_stub_list.push_back(seg%16)
	    : m_stub_list.push_back(-1);

	  m_stub_list.push_back(pt);
	  m_stub_list.push_back(goodstub);
	  m_chip_trig.erase(m_iter->first);
	  m_chip_trig.insert(std::make_pair(B_id,m_stub_list));
	}
      }  // End of digi collection
    }      

    // Finally we add both collection to the stores
    
    m_data_trig.push_back(m_chip_trig);
    m_data_raw.push_back(m_chip_raw);
  }

  // Now we have two collections one containing the raw data, and one the trigger data for all event at the concentrator level

  // Trigger data is sent at every BX
  // Raw data is sent every n BX, where n depends on the L1 rate

  // If the L1 rate is given in kHz, we have n=40000/rate
  //
  // For the moment we produce it only for a concentrator chip

  cout << "--> Entering loop 2, producing the random Trig/L1 sequence..." << endl;

  int delta    = static_cast<int>(40000/m_rate);

  int n_trig = m_data_trig.size();
  int n_raw  = m_data_raw.size();

  std::default_random_engine generator(time(NULL));
  std::normal_distribution<double> L1s(delta,delta/3);

  std::cout << "Average L1A period will be " << delta << " BXs with a sigma of " << delta/3 << " BXs" << std::endl;

  srand(time(NULL));

  int trg_evnum;
  //int raw_evnum;
  //int raw_rank;
  //int delta_ref=delta;

  std::vector<int> trig_seq;
  std::vector<int> raw_seq;

  trig_seq.clear();
  raw_seq.clear();

  // First we create random sequence for L1 satisfying the trigger rules

  // The first L1 is always sent at BX=0

  raw_seq.push_back(rand()%n_raw);

  int rank_r1 = 0; 
  int rank_r2 = 0; 
  int rank_r3 = 0; 

  int lim_r1 = 3; 
  int lim_r2 = 25; 
  int lim_r3 = 100; 

  bool L1_rule = true;

  cout << "... Generate the sequence of " << nevts << " by picking up events in the store..." << std::endl;

  while (rank_r1<nevts)
  {   
    L1_rule = false;

    // First we need to generate a space which is not enforcing the trigger rules
    //
    // Then check that the rules are not enforced
    // L1T rules are taken from:
    // 
    // http://www.hephy.at/project/cms/trigger/globalTrigger/trans/GT_status_dec03.pdf
    //

    while (!L1_rule)
    {
      delta = L1s(generator); // Compute a L1 space

      if (delta < lim_r1) continue;         // Rule 1 enforced
      if ((delta+rank_r1-rank_r2) < lim_r2 && rank_r2!=0) continue; // Rule 2 enforced
      if ((delta+rank_r1-rank_r3) < lim_r3 && rank_r3!=0) continue; // Rule 3 enforced

      L1_rule=true;
    }

    // Event satisfying the trigger rules update the ranks

    rank_r3 = rank_r2;
    rank_r2 = rank_r1;
    rank_r1 = rank_r2+delta;

    for (int i=0;i<delta-1;++i)
    {
      raw_seq.push_back(-1);
    }

    raw_seq.push_back(rand()%n_raw);
  }

  int n_l1 = 0;

  for (int i=0;i<nevts;++i)
  { 
    trg_evnum = rand()%n_trig;
    trig_seq.push_back(trg_evnum);

    if (raw_seq.at(i)!=-1) ++n_l1;
  }

  float av_rate = 1000./((float(nevts)/float(n_l1)*0.025));

  cout << "Got " << n_l1 << " L1A over " << nevts << " BXs" << endl; 
  cout << "This corresponds to av. L1A freq of " << av_rate << " kHz" << endl; 

  cout << "Sequence of L1 events is (BX,evt_ID) :" << endl; 

  for (unsigned int i=0;i<raw_seq.size();++i)
  { 
    if (raw_seq.at(i)==-1) continue;

    cout << "(" << i << "," << raw_seq.at(i) << "),"; 
  }

  cout << endl;


  // Now we have a sequence of nevts, with L1A distributed randomly 
  // following the trigger rules

  // Next stage consists in creating a ROOT/txt file containing this sequence 
  // of events, with the correspondind data for each Concentrator chip

  cout << "--> Entering loop 3, producing the final root file..." << endl;


  // Create the sequences, and write them on the fly 

  int L1_id = 0;

  std::vector<int> FIFO,FIFO_new;
  int FIFO_size;
  int last_rank;
  int last_BX;

  m_raw_FIFO.clear();

  if (conc)
  {
    cout << "Coding the bend over " << bend_bit_MPA 
	 << " bit(s) in the MPA and " << bend_bit_CBC 
	 << " bit(s) in the CBC in the concentrator output..." << endl;

    for (unsigned int i=0;i<m_concs.size();++i)
    {
      m_digi_list.clear();
      m_digi_list.push_back(-1);
      m_raw_FIFO.insert(std::make_pair(m_concs.at(i),m_digi_list));
    }
  }

  // First write the FE raw data (pretty simple 
  // because here we don't really bother with BX sequences)
  //
  // We just try to emulate the FIFO behaviour for the concentrator case

 
  if (m_write_raw)
  {
    for (int i=0;i<nevts;++i)
    { 
      if (i%100==0)
	cout << i << endl;

      if (nevts%3564==0) L1_id=0; // BCR, L1ID get back to 0

      evtbuilder::initVars();

      if (raw_seq.at(i)!=-1) // CBC/MPA received a L1 A, we write the raw block
      {
	++L1_id;
	
	m_chip_raw = m_data_raw.at(raw_seq.at(i)); // Pick up the event in the store
	
	for ( m_iter = m_chip_raw.begin(); m_iter != m_chip_raw.end();++m_iter )
	{

	  // First, look at the main info
	  isPS            = false;
	  m_raw_bx        = i;
	  m_raw_FIFO_FULL = 0;
	  m_raw_chip      = m_iter->first;
	  m_raw_lay       = m_raw_chip/1000000;
	  m_raw_lad       = (m_raw_chip-1000000*m_raw_lay)/10000;
	  m_raw_mod       = (m_raw_chip-1000000*m_raw_lay-10000*m_raw_lad)/100;
	  m_raw_np        = 0;
	  m_raw_ns        = 0;


	  if (m_raw_lay<8 || (m_raw_lay>10 && m_raw_lad<9)) isPS = true; 

	  // Now, we get the digis contained in the event, and fill the block
	  m_digi_list = m_iter->second;
	  
	  (conc)
	    ? evtbuilder::fill_CONC_RAW_block(m_digi_list,isPS,L1_id)
	    : evtbuilder::fill_RAW_block(m_digi_list,isPS,L1_id);
	  
	  m_raw_size  = m_raw_data->size();

	  // Concentrator FIFO size is 266*8*16 = 34048 bits
	  // You can store up to 17024 BX of data transmission

	  if (conc)
	  {
	    // Find the chip FIFO footprint
	    // defined as follows:

	    // < CHIP_ID, <-1,BXi(1),BXo(1), > >

	    m_iter2 = m_raw_FIFO.find(m_raw_chip);

	    FIFO_new.clear();
	    FIFO=m_iter2->second;
	    FIFO_size=(FIFO.size()-1)/2;	    
	    FIFO_new.push_back(-1);

	    if (FIFO_size==0) // Empty FIFO
	    {
	      FIFO_new.push_back(i);                            // BXi is the current BX
	      FIFO_new.push_back(i+m_raw_size/(m_npblock/8)+1); // BXi + numb of BXs to extract the event 
	      m_raw_FIFO_FULL = (FIFO_new.size()-1)/2;          // Number of events in the FIFO (1)
	    }
	    else  // FIFO not empty
	    {
	      last_BX=i;

	      for (unsigned int j=0;j<(FIFO.size()-1)/2;++j)	    
	      {
		if (FIFO.at(2*j+2)<i) continue; // Time to go out for this event (BXo(j)<i)

		FIFO_new.push_back(FIFO.at(2*j+1)); // Otherwise event stays there
		FIFO_new.push_back(FIFO.at(2*j+2));

		last_BX=FIFO.at(2*j+2); // The last BXo fixes the time where we can start to extract the next event
	      }

	      // Here we compute the size of the current FIFO

	      last_rank=0;
	      if (FIFO_new.size()>1) last_rank = (last_BX-FIFO_new.at(1))*(m_npblock/8);

	      if (last_rank+m_raw_size>34048) // The new event cannot fit in, we raise an error 
	      {
		m_raw_FIFO_FULL = -1;
	      }
	      else // OK, come in!
	      {
		FIFO_new.push_back(last_BX);
		FIFO_new.push_back(last_BX+m_raw_size/(m_npblock/8)+1);
		m_raw_FIFO_FULL = (FIFO_new.size()-1)/2;
	      }
	    }

	    m_raw_FIFO.erase(m_iter2->first);
	    m_raw_FIFO.insert(std::make_pair(m_raw_chip,FIFO_new)); // Update the FIFO
	  }

	  m_raw_tree->Fill();
	}
      }
    }
  }

  // Then write the CONC TRG data (8BX for both cases)
  
  if (m_write_trg)
  {
    if (conc) // For the concentrator we work on a 8BX basis
    {
      for (int i=0;i<nevts/8;++i)
      { 
	if (i%100==0)
	  cout << 8*i << endl;

	m_chip_trig   = m_data_trig.at(trig_seq.at(8*i)); // Pick up the event in the store
	
	for ( m_iter = m_chip_trig.begin(); m_iter != m_chip_trig.end();++m_iter )
	{
	  trig_sequence.clear(); 

	  // The main info

	  isPS            = false;
	  m_tri_bx        = 8*i;
	  m_tri_nstubs    = 0;
	  m_tri_nstubs_g  = 0;
	  m_tri_nstubs_s  = 0;
	  m_tri_nstubs_gs = 0;
	  m_tri_chip      = m_iter->first;
	  m_tri_lay       = m_tri_chip/1000000;
	  m_tri_lad       = (m_tri_chip-1000000*m_tri_lay)/10000;
	  m_tri_mod       = (m_tri_chip-1000000*m_tri_lay-10000*m_tri_lad)/100;
	  
	  if (m_tri_lay<8 || (m_tri_lay>10 && m_tri_lad<9)) isPS = true; 
	  
	  // Then we get all the stubs of the corresponding chip
	  // in trig_sequence

	  trig_sequence.push_back(m_iter->second);
	  
	  for (int j=1;j<8;++j)
	  {
	    trig_sequence.push_back((m_data_trig.at(trig_seq.at(8*i+j)).find(m_tri_chip))->second);
	  }
	
	  // And we write the word
	  // std::cout << m_tri_bx << " / " << m_tri_chip << " -> ";
	  evtbuilder::fill_TRG_block(trig_sequence,isPS,conc,8*i);
	}
      }
    }
    else // For the front end chips, we write on a 2BX basis
    {
      for (int i=0;i<nevts/2;++i)
      { 
	if (i%100==0)
	  cout << i << endl;

	m_chip_trig   = m_data_trig.at(trig_seq.at(2*i));
	
	for ( m_iter = m_chip_trig.begin(); m_iter != m_chip_trig.end();++m_iter )
	{
	  trig_sequence.clear(); 
	  
	  isPS        = false;
	  m_tri_bx    = 2*i;
	  m_tri_nstubs= 0;
	  m_tri_nstubs_g=0;
	  m_tri_nstubs_s=0;
	  m_tri_nstubs_gs=0;
	  m_tri_chip  = m_iter->first;
	  m_tri_lay   = m_tri_chip/1000000;
	  m_tri_lad   = (m_tri_chip-1000000*m_tri_lay)/10000;
	  m_tri_mod   = (m_tri_chip-1000000*m_tri_lay-10000*m_tri_lad)/100;

	  trig_sequence.push_back(m_iter->second);
	  
	  if (m_tri_lay<8 || (m_tri_lay>10 && m_tri_lad<9)) isPS = true; 
	  
	  if (isPS) 	  
	  {
	    m_tri_nstubs= 0;
	    m_tri_nstubs_g=0;
	    m_tri_nstubs_s=0;
	    m_tri_nstubs_gs=0;
	    trig_sequence.push_back((m_data_trig.at(trig_seq.at(2*i+1)).find(m_tri_chip))->second);
	    
	    std::cout << m_tri_bx << " / " << m_tri_chip << " -> ";
	    evtbuilder::fill_TRG_block(trig_sequence,isPS,conc,2*i);
	  }
	  else
	  {
	    m_tri_nstubs= 0;
	    m_tri_nstubs_g=0;
	    m_tri_nstubs_s=0;
	    m_tri_nstubs_gs=0;
	    std::cout << m_tri_bx << " / " << m_tri_chip << " -> ";
	    evtbuilder::fill_TRG_block(trig_sequence,isPS,conc,2*i);
	    trig_sequence.clear(); 
	    m_tri_bx    = 2*i+1;
	    m_tri_nstubs= 0;
	    m_tri_nstubs_g=0;
	    m_tri_nstubs_s=0;
	    m_tri_nstubs_gs=0;
	    trig_sequence.push_back((m_data_trig.at(trig_seq.at(2*i+1)).find(m_tri_chip))->second);
	    std::cout << m_tri_bx << " / " << m_tri_chip << " -> ";
	    evtbuilder::fill_TRG_block(trig_sequence,isPS,conc,2*i+1);
	  }
	}
      }
    }
  }
  
  m_outfile->Write();
  
  delete L1TT;
  delete PIX;
  delete m_outfile;
}




/////////////////////////////////////////////////////////////
//
// Basic methods, initializations,...
//
/////////////////////////////////////////////////////////////

void evtbuilder::initVars()
{
  m_tri_bx=0;
  m_raw_bx=0;
  m_tri_chip=0;
  m_raw_chip=0;
  m_tri_size=0;
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
}


void evtbuilder::initTuple(std::string inRAW,std::string inTRG,std::string out)
{
  m_outfile    = new TFile(out.c_str(),"recreate");

  m_tri_tree   = new TTree("Trigger_FE","L1Trigger words after FE");
  m_raw_tree   = new TTree("Raw_FE","Raw data words after FE");

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
  m_raw_tree->Branch("RAW_FIFULL",     &m_raw_FIFO_FULL,"RAW_FIFULL/I");

  m_raw_tree->Branch("RAW_NPCLUS",     &m_raw_np);
  m_raw_tree->Branch("RAW_NSCLUS",     &m_raw_ns);
 
  L1TT   = new TChain("TkStubs");
  PIX    = new TChain("Pixels"); 
 
  // Input RAW data file
  
  std::size_t found = inRAW.find(".root");
  
  // Case 1, it's a root file
  if (found!=std::string::npos)
  {
    PIX->Add(inRAW.c_str());
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
      if (found!=std::string::npos) PIX->Add(STRING.c_str());
    }
    
    in2.close();
  }

 
  // Input TRG data file
  
  found = inTRG.find(".root");
  
  // Case 1, it's a root file
  if (found!=std::string::npos)
  {
    L1TT->Add(inTRG.c_str());
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
      if (found!=std::string::npos) L1TT->Add(STRING.c_str());
    }
    
    in2.close();
  }


  
  pm_pix_layer=&m_pix_layer;
  pm_pix_ladder=&m_pix_ladder;
  pm_pix_module=&m_pix_module;
  pm_pix_row=&m_pix_row;
  pm_pix_col=&m_pix_col;

  PIX->SetBranchAddress("PIX_n",         &m_pix);
  PIX->SetBranchAddress("PIX_nPU",       &m_npu);
  PIX->SetBranchAddress("PIX_layer",     &pm_pix_layer);
  PIX->SetBranchAddress("PIX_ladder",    &pm_pix_ladder);
  PIX->SetBranchAddress("PIX_module",    &pm_pix_module);
  PIX->SetBranchAddress("PIX_row",       &pm_pix_row);
  PIX->SetBranchAddress("PIX_column",    &pm_pix_col);

  pm_stub_layer=&m_stub_layer;
  pm_stub_ladder=&m_stub_ladder;
  pm_stub_module=&m_stub_module;
  pm_stub_pt=&m_stub_pt;
  pm_stub_deltas=&m_stub_deltas;
  pm_stub_strip=&m_stub_strip;
  pm_stub_seg=&m_stub_seg;
  pm_stub_chip=&m_stub_chip;
  pm_stub_X0=&m_stub_X0;
  pm_stub_Y0=&m_stub_Y0;
  pm_stub_pxGEN=&m_stub_pxGEN;
  pm_stub_pyGEN=&m_stub_pyGEN;
 
  L1TT->SetBranchAddress("L1TkSTUB_n",         &m_stub);
  L1TT->SetBranchAddress("L1TkSTUB_layer",     &pm_stub_layer);
  L1TT->SetBranchAddress("L1TkSTUB_ladder",    &pm_stub_ladder);
  L1TT->SetBranchAddress("L1TkSTUB_module",    &pm_stub_module);
  L1TT->SetBranchAddress("L1TkSTUB_pt",        &pm_stub_pt);
  L1TT->SetBranchAddress("L1TkSTUB_deltas",    &pm_stub_deltas);
  L1TT->SetBranchAddress("L1TkSTUB_strip",     &pm_stub_strip);
  L1TT->SetBranchAddress("L1TkSTUB_seg",       &pm_stub_seg);
  L1TT->SetBranchAddress("L1TkSTUB_chip",      &pm_stub_chip);
  L1TT->SetBranchAddress("L1TkSTUB_pxGEN",     &pm_stub_pxGEN);
  L1TT->SetBranchAddress("L1TkSTUB_pyGEN",     &pm_stub_pyGEN);
  L1TT->SetBranchAddress("L1TkSTUB_X0",        &pm_stub_X0);
  L1TT->SetBranchAddress("L1TkSTUB_Y0",        &pm_stub_Y0);
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

bool evtbuilder::convert(std::string sectorfilename) 
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

      if (int(atoi(s.c_str())/10000)!=m_lay && m_lay!=-1) continue;

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
// the following presentation:
//
// https://indico.cern.ch/event/321805/session/5/contribution/13/material/slides/1.pdf

//
// 1: L1 raw block at the MPA/CBC level
//

void evtbuilder::fill_RAW_block(std::vector<int> digis,bool spars,int BXid)
{
  m_raw_data->clear();

  // First write the headers

  (spars)
    ? evtbuilder::fill_RAW_header_MPA(BXid)
    : evtbuilder::fill_RAW_header_CBC(BXid);

  // Then go for the data

  int ndata=(digis.size()-1)/3; // Default first item is -1
  int hsize=m_raw_data->size();

  if (!spars) // CBC (unsparsified)
  {
    for (int j=0;j<256;++j)   m_raw_data->push_back(0);
    for (int j=0;j<ndata;++j) m_raw_data->at(digis.at(3*j+2)+hsize) = 1; 
  }
  else // MPA
  {
    int data_MPA[128][17];
    int row,col;

    // The sparsification is done at the MPA level
    // Strips 0 to 127 belong to the pixel layer 
    // 128 to 253 belong to the strip layer

    for (int j=0;j<128;++j)
    {
      for (int i=0;i<17;++i) data_MPA[j][i] = 0;
    }

    for (int j=0;j<ndata;++j)
    {
      row = digis.at(3*j+2);
      col = digis.at(3*j+3);
    
      (row<128)
	? data_MPA[row][col]    = 1
	: data_MPA[row%128][16] = 1;
    }

    int np = 0;
    int ns = 0;

    int compt=0;

    int start_r = -1;

    std::vector<int> clus_p;
    std::vector<int> clus_s;

    clus_p.clear();
    clus_s.clear();

    for (int j=0;j<17;++j)
    {
      start_r = -1;
      compt   = 0;

      //      cout << j << " / " << clus_p.size()/3 << " / " << clus_s.size()/2 << endl;

      for (int i=0;i<128;++i)
      {
	if (data_MPA[i][j] == 1)
	{
	  if (start_r==-1) // new clus
	  {
	    start_r=i;
	    compt=0;
	  }
	  else 
	  {
	    ++compt;

	    if (compt==7) // Max width reached
	    {
	      if (j<16)
	      {
		clus_p.push_back(start_r);
		clus_p.push_back(j);
		clus_p.push_back(compt);
	      }
	      else
	      {
		clus_s.push_back(start_r);
		clus_s.push_back(compt);
	      }

	      start_r=-1;
	      compt=0;
	    }  
	  }
	}
	else
	{	
	  if (start_r!=-1)
	  {
	    if (j<16)
	    {
	      clus_p.push_back(start_r);
	      clus_p.push_back(j);
	      clus_p.push_back(compt);
	    }
	    else
	    {
	      clus_s.push_back(start_r);
	      clus_s.push_back(compt);
	    }

	    start_r=-1;
	    compt=0;
	  }
	}
      }

      if (start_r!=-1) // Just handle the case where rows ends up with 1
      {
	if (j<16)
	{
	  clus_p.push_back(start_r);
	  clus_p.push_back(j);
	  clus_p.push_back(compt);
	}
	else
	{
	  clus_s.push_back(start_r);
	  clus_s.push_back(compt);
	}
      }
    }

    //    std::cout << "This MPA chip contains " << clus_p.size()/3 << " pixels clusters and " 
    //	      << clus_s.size()/2 << " strip clusters" << std::endl; 

    // Now encode them 
    np = clus_p.size()/3 ;
    ns = clus_s.size()/2; 
    m_raw_np        = np;
    m_raw_ns        = ns;


    np = std::min(31,np);
    ns = std::min(31,ns);

    std::bitset<5> N_S = ns;
    std::bitset<5> N_P = np;

    for (int j=0;j<5;++j) m_raw_data->push_back(N_P[4-j]);
    for (int j=0;j<5;++j) m_raw_data->push_back(N_S[4-j]);
    m_raw_data->push_back(0);

    for (int j=0;j<np;++j)
    {
      std::bitset<7> row =  clus_p.at(3*j);
      std::bitset<4> col =  clus_p.at(3*j+1);
      std::bitset<3> wdt =  clus_p.at(3*j+2);

      for (int k=0;k<7;++k) m_raw_data->push_back(row[6-k]);
      for (int k=0;k<3;++k) m_raw_data->push_back(wdt[2-k]);
      for (int k=0;k<4;++k) m_raw_data->push_back(col[3-k]);
    }

    for (int j=0;j<ns;++j)
    {
      std::bitset<7> row =  clus_s.at(2*j);
      std::bitset<3> wdt =  clus_s.at(2*j+1);

      for (int k=0;k<7;++k) m_raw_data->push_back(row[6-k]);
      for (int k=0;k<3;++k) m_raw_data->push_back(wdt[2-k]);
    }
  }

  // Here we put a common trailer (if needed)


  //  std::cout << "Size of the L1 raw word / " << spars << " / " << m_raw_data->size() << std::endl;

  // for (int j=0;j<m_raw_data->size();++j)
  // {
  //  std::cout << m_raw_data->at(j);
  // }

  //  std::cout << std::endl;


  // Here we just count the maximum number of consecutive bits at 1

  int max_bits = 0;
  int nbits    = 0;

  for (unsigned int j=hsize;j<m_raw_data->size();++j)
  {
    if (m_raw_data->at(j)==1) ++nbits;
    if (m_raw_data->at(j)==0) nbits=0;
    if (nbits>=max_bits)      max_bits=nbits;
  }

  m_raw_mbits = max_bits;

}

//
// 2: L1 raw block at the concentrator level
//

void evtbuilder::fill_CONC_RAW_block(std::vector<int> digis,bool spars,int BXid)
{
  m_raw_data->clear();

  // First write the header

  evtbuilder::fill_CONC_RAW_header(BXid);

  // Then go for the data

  int np,ns;

  int ndata=(digis.size()-1)/3; // Default first item is -1
  int hsize=m_raw_data->size();

  if (!spars) // CBC case
  {
    int data_CBC[8][256];
    int row,chip;

    for (int j=0;j<256;++j)
    {
      for (int i=0;i<8;++i) data_CBC[i][j] = 0;
    }

    for (int j=0;j<ndata;++j)
    {
      chip = digis.at(3*j+1);
      row  = digis.at(3*j+2);
     
      //      std::cout << row << "," << chip << " / "; 
 
      data_CBC[chip][row] = 1;
    }
    //    if (ndata!=0) std::cout << std::endl;

    ns      = 0;
    int compt   = 0;
    int start_r = -1;

    std::vector<int> clus_s;

    clus_s.clear();

    for (int j=0;j<8;++j)
    {
      start_r = -1;
      compt   = 0;

      //      cout << j << " / " << clus_s.size()/3 << endl;

      for (int i=0;i<256;++i)
      {
	if (data_CBC[j][i] == 1)
	{
	  if (start_r==-1) // new clus
	  {
	    start_r=i;
	    compt=0;
	  }
	  else 
	  {
	    ++compt;

	    if (compt==7) // Max width reached
	    {
	      clus_s.push_back(start_r);
	      clus_s.push_back(compt);
	      clus_s.push_back(j);	      

	      start_r=-1;
	      compt=0;
	    }  
	  }
	}
	else
	{	
	  if (start_r!=-1)
	  {
	    clus_s.push_back(start_r);
	    clus_s.push_back(compt);
	    clus_s.push_back(j);	      

	    start_r=-1;
	    compt=0;
	  }
	}
      }

      if (start_r!=-1) // Just handle the case where rows ends up with 1
      {
	clus_s.push_back(start_r);
	clus_s.push_back(compt);
	clus_s.push_back(j);	      
      }
    }
    
    //    std::cout << "This CBC conc chip contains " << clus_s.size()/3 << " strip clusters" << std::endl; 

    // Now encode them 
    ns = clus_s.size()/3; 
    m_raw_ns        = ns;
    ns = std::min(31,ns);

    std::bitset<5> N_S = ns;

    for (int j=0;j<5;++j) m_raw_data->push_back(N_S[4-j]);
    m_raw_data->push_back(0);

    for (int j=0;j<ns;++j)
    {
      std::bitset<3> chp =  clus_s.at(3*j+2);
      std::bitset<8> row =  clus_s.at(3*j);
      std::bitset<3> wdt =  clus_s.at(3*j+1);

      for (int k=0;k<3;++k) m_raw_data->push_back(chp[2-k]);
      for (int k=0;k<8;++k) m_raw_data->push_back(row[7-k]);
      for (int k=0;k<3;++k) m_raw_data->push_back(wdt[2-k]);
    }
  }
  else // MPA case
  {
    int data_MPA[8][128][17];
    int row,col,chip;

    // The sparsification is done at the MPA level
    // Strips 0 to 127 belong to the pixel layer 
    // 128 to 253 belong to the strip layer
 
    for (int k=0;k<8;++k)
    {
      for (int j=0;j<128;++j)
      {
	for (int i=0;i<17;++i) data_MPA[k][j][i] = 0;
      }
    }

    for (int j=0;j<ndata;++j)
    {
      chip = digis.at(3*j+1);
      row  = digis.at(3*j+2);
      col  = digis.at(3*j+3);
     
      //      std::cout << chip << "," << row << "," << col << " / "; 
 
      (row<128)
	? data_MPA[chip][row][col]    = 1
	: data_MPA[chip][row%128][16] = 1;
    }
    //    if (ndata!=0) std::cout << std::endl;

    np = 0;
    ns = 0;

    int compt=0;

    int start_r = -1;

    std::vector<int> clus_p;
    std::vector<int> clus_s;

    clus_p.clear();
    clus_s.clear();

    for (int k=0;k<8;++k)
    {
      for (int j=0;j<17;++j)
      {
	start_r = -1;
	compt   = 0;

	//      cout << j << " / " << clus_p.size()/3 << " / " << clus_s.size()/2 << endl;

	for (int i=0;i<128;++i)
	{
	  if (data_MPA[k][i][j] == 1)
	  {
	    if (start_r==-1) // new clus
	    {
	      start_r=i;
	      compt=0;
	    }
	    else 
	    {
	      ++compt;

	      if (compt==7) // Max width reached
	      {
		if (j<16)
		{
		  clus_p.push_back(start_r);
		  clus_p.push_back(j);
		  clus_p.push_back(compt);
		  clus_p.push_back(k);
		}
		else
		{
		  clus_s.push_back(start_r);
		  clus_s.push_back(compt);
		  clus_s.push_back(k);
		}
		
		start_r=-1;
		compt=0;
	      }  
	    }
	  }
	  else
	  {	
	    if (start_r!=-1)
	    {
	      if (j<16)
	      {
		clus_p.push_back(start_r);
		clus_p.push_back(j);
		clus_p.push_back(compt);
		clus_p.push_back(k);
	      }
	      else
	      {
		clus_s.push_back(start_r);
		clus_s.push_back(compt);
		clus_s.push_back(k);
	      }
	      
	      start_r=-1;
	      compt=0;
	    }
	  }
	}

	if (start_r!=-1) // Just handle the case where rows ends up with 1
	{
	  if (j<16)
	  {
	    clus_p.push_back(start_r);
	    clus_p.push_back(j);
	    clus_p.push_back(compt);
	    clus_p.push_back(k);
	  }
	  else
	  {
	    clus_s.push_back(start_r);
	    clus_s.push_back(compt);
	    clus_s.push_back(k);
	  }
	}
      }
    }

    //    std::cout << "This MPA chip contains " << clus_p.size()/4 << " pixels clusters and " 
    //	      << clus_s.size()/3 << " strip clusters" << std::endl; 

    // Now encode them 
    np = clus_p.size()/4 ;
    ns = clus_s.size()/3; 
    m_raw_np        = np;
    m_raw_ns        = ns;
    np = std::min(31,np);
    ns = std::min(31,ns);

    std::bitset<5> N_S = ns;
    std::bitset<5> N_P = np;

    //    std::cout << "This MPA chip contains

    for (int j=0;j<5;++j) m_raw_data->push_back(N_P[4-j]);
    for (int j=0;j<5;++j) m_raw_data->push_back(N_S[4-j]);
    m_raw_data->push_back(0);

    for (int j=0;j<np;++j)
    {
      std::bitset<7> row =  clus_p.at(4*j);
      std::bitset<4> col =  clus_p.at(4*j+1);
      std::bitset<3> wdt =  clus_p.at(4*j+2);
      std::bitset<3> chp =  clus_p.at(4*j+3);

      for (int k=0;k<3;++k) m_raw_data->push_back(chp[2-k]);
      for (int k=0;k<7;++k) m_raw_data->push_back(row[6-k]);
      for (int k=0;k<3;++k) m_raw_data->push_back(wdt[2-k]);
      for (int k=0;k<4;++k) m_raw_data->push_back(col[3-k]);

    }

    for (int j=0;j<ns;++j)
    {
      std::bitset<7> row =  clus_s.at(3*j);
      std::bitset<3> wdt =  clus_s.at(3*j+1);
      std::bitset<3> chp =  clus_s.at(3*j+2);

      for (int k=0;k<3;++k) m_raw_data->push_back(chp[2-k]);
      for (int k=0;k<7;++k) m_raw_data->push_back(row[6-k]);
      for (int k=0;k<3;++k) m_raw_data->push_back(wdt[2-k]);
    }
  }



  // Here we put a common trailer (if needed)

  for (int k=0;k<4;++k) m_raw_data->push_back(0);


  int max_bits = 0;
  int nbits    = 0;

  for (unsigned int j=hsize;j<m_raw_data->size();++j)
  {
    if (m_raw_data->at(j)==1) ++nbits;
    if (m_raw_data->at(j)==0) nbits=0;
    if (nbits>=max_bits)      max_bits=nbits;
  }

  m_raw_mbits = max_bits;

  if (m_raw_mbits>17)
  {
    std::cout << "This chip contains " << np << " pixels clusters and " 
    	      << ns << " strip clusters" << std::endl; 

    std::cout << "Size of the L1 raw word / " << spars << " / " << m_raw_data->size() << std::endl;

    for (unsigned int j=0;j<m_raw_data->size();++j)
    {
      std::cout << m_raw_data->at(j);
    }

    std::cout << std::endl;
  }

}


//
// 3: Trigger block
//


void evtbuilder::fill_TRG_block(std::vector< std::vector<int> > stubs, bool spars, bool conc, int BXid)
{
  m_tri_data->clear();

  // Then go for the data

  int nstubs[8];
  int seqsize = static_cast<int>(stubs.size());

  int lim_CBC=3; // Max number of stubs passed by the CBC is 3/BX
  int lim_MPA=4; // Max number of stubs passed by the MPA is 4/2BXs

  if (bend_bit_MPA<4) lim_MPA=5;

  int gstub;

  m_tri_nstubs    = 0;
  m_tri_nstubs_s  = 0;
  m_tri_nstubs_g  = 0;
  m_tri_nstubs_gs = 0;

  std::vector<int> stublist;

  if (conc) evtbuilder::fill_CONC_TRG_header(BXid%3564);

  //std::cout << "Writing trigger block for " << seqsize << " BXs" << std::endl;
  //  
  //  (conc)
  //  ? std::cout << "Concentrator" << std::endl
  //  : std::cout << "FE" << std::endl;
  
  for (int i=0;i<seqsize;++i)
  {
    std::bitset<3> bx = i;
 
    stublist = stubs.at(i); // The 
    
    if (spars && i%2==0)
    {
      for (int j=0;j<8;++j) nstubs[j]=0;
    }

    if (!spars)
    {
      for (int j=0;j<8;++j) nstubs[j]=0;
    }

    //std::cout << i << "/" << nstubs << "/" << bend_bit_MPA << "/" << bend_bit_CBC << "/" ;
    
    // Loop over the stubs contained in the event

    for (unsigned int kk=0;kk<(stublist.size()-1)/5;++kk)
    {
      gstub=stublist.at(5*kk+5);
      
      ++m_tri_nstubs;
      if (gstub>0) ++m_tri_nstubs_g;

      if (nstubs[stublist.at(5*kk+1)]>=lim_CBC && !spars) continue; // No CHIP sorting for the moment
      if (nstubs[stublist.at(5*kk+1)]>=lim_MPA && spars) continue;  // No CHIP sorting for the moment

      ++nstubs[stublist.at(5*kk+1)];
      ++m_tri_nstubs_s;
      if (gstub>0) ++m_tri_nstubs_gs;

      std::bitset<3> chp = stublist.at(5*kk+1);        // Chip number
      std::bitset<8> pos = stublist.at(5*kk+2);        // Stub position
      std::bitset<4> col = stublist.at(5*kk+3);        // Z position (for PS module)
      std::bitset<5> off = abs(2*stublist.at(5*kk+4)); // Bend

      //std::cout <<  stublist.at(4*kk+1)<< "," 
      //	<<  stublist.at(4*kk+2)<< "," 
      //	<<  stublist.at(4*kk+3) << "," 
      //	<<  stublist.at(4*kk+4) << "/" ;
      
      if (spars && !conc) m_tri_data->push_back(i);
      if (conc) 
      {
	for (int j=0;j<3;++j) m_tri_data->push_back(bx[2-j]);
	for (int j=0;j<3;++j) m_tri_data->push_back(chp[2-j]);
      }

      for (int j=0;j<8;++j) m_tri_data->push_back(pos[7-j]);
      if (spars) for (int j=0;j<4;++j) m_tri_data->push_back(col[3-j]);
      
      if (conc)
      {
	if (spars)
	{
	  for (int j=0;j<bend_bit_MPA;++j) m_tri_data->push_back(off[4-j]);
	}
	else
	{
	  for (int j=0;j<bend_bit_CBC;++j) m_tri_data->push_back(off[4-j]);
	}
      }
      else
      {
	for (int j=0;j<5;++j) m_tri_data->push_back(off[4-j]);
      }
      
    }
    
  }
  //  std::cout << std::endl;

  // Potentially put a trailer here

  if (conc) 
  {
    for (int j=0;j<4;++j) m_tri_data->push_back(0);
  }

  m_tri_size=m_tri_data->size();
  m_tri_tree->Fill();


  //  std::cout << "Size of the trigger word / " << spars << " / " << conc << " / " << m_tri_size << std::endl;
  /*
  for (int j=0;j<m_tri_data->size();++j)
  {
    std::cout << m_tri_data->at(j);
  }
 
  std::cout << std::endl;
  */

}



void evtbuilder::fill_RAW_header_CBC(int L1id)
{
  // Format of the CBC L1 word header
  //
  // HHEECCCCCCCCCPPPPPPPPP : HHEE (header + error) CC..CC (L1 ID bet 0 and 512)
  //                          PP..PP CBC pipeline address

  for (int j=0;j<4;++j) m_raw_data->push_back(1); // HHEE

  if (L1id>512)
  {
    std::cout << "Too many L1ids, problem!!!" << std::endl;
  }

  std::bitset<9> L1_ID = L1id;

  for (int j=0;j<9;++j) m_raw_data->push_back(L1_ID[8-j]); // CC..CC
  for (int j=0;j<9;++j) m_raw_data->push_back(L1_ID[8-j]); // PP..PP
}


void evtbuilder::fill_RAW_header_MPA(int L1id)
{
  // Format of the CBC L1 word header
  //
  // 1111111111111111110EECCCCCCCCC0 : EE (error) CC..CC (L1 ID bet 0 and 512)
  //

  for (int j=0;j<18;++j) m_raw_data->push_back(1); 
  m_raw_data->push_back(0); 
  m_raw_data->push_back(1); 
  m_raw_data->push_back(1); 

  if (L1id>512)
  {
    std::cout << "Too many L1ids, problem!!!" << std::endl;
  }

  std::bitset<9> L1_ID = L1id;

  for (int j=0;j<9;++j) m_raw_data->push_back(L1_ID[8-j]); // CC..CC

  m_raw_data->push_back(0); 

}

void evtbuilder::fill_CONC_RAW_header(int L1id)
{
  // Format of the CONC L1 word header
  //
  // 1111111111111110EEEEEEEEECCCCCCCCC : E..E (error) CC..CC (L1 ID bet 0 and 512)
  //

  for (int j=0;j<15;++j) m_raw_data->push_back(1); 
  m_raw_data->push_back(0); 
  for (int j=0;j<9;++j) m_raw_data->push_back(0); 

  if (L1id>512)
  {
    std::cout << "Too many L1ids, problem!!!" << std::endl;
  }

  std::bitset<9> L1_ID = L1id;

  for (int j=0;j<9;++j) m_raw_data->push_back(L1_ID[8-j]); // CC..CC

  m_raw_data->push_back(0); 

}

void evtbuilder::fill_CONC_TRG_header(int BXid)
{
  // Format of the CONC TRIGGER word header
  //
  // CSSSSSSSSSCCCCCCCCCCCC : SS..SS (status) CC..CC (BX ID bet 0 and 3564)
  //

  for (int j=0;j<10;++j) m_tri_data->push_back(0); 

  std::bitset<12> BX_ID = BXid;

  for (int j=0;j<12;++j) m_tri_data->push_back(BX_ID[11-j]); // CC..CC
}
