// Base class for Phase II tracker FrontEnd tests

// Main constructor

// Constructor for the concentrator

// Meaning of the parameters
//
//
// filenameRAW : the input data files for QCD events (go to L1)
// filenameTRG : the input data files for Minbias events
// outfile     : name of the output root file
// npatt       : how many BXs should be produced (will be fitted correctly afterwards)
// rate        : input L1A rate (for the L1 raw block writing
// layer       : restrict info to a given layer, or not (give -1 in this case)
// ladder      : restrict info to a given ladder, or not (give -1 in this case)
// module      : restrict info to a given module, or not (give -1 in this case)
// sector      : the input file defining the modules (csv file)
// RAW         : write the L1 block or not
// TRIG        : write the Trigger block or not
// npblock     : How many bits do we reserve to L1 in the concentrator block (should be a multiple of 8) 
// conc        : Are we building concentrator (true) or CBC/MPA (false) data sequences
// L1prop      : The proportion of PU4T events you want in the TRG flow (in %)

#include "evtbuilder.h"

evtbuilder::evtbuilder(std::string filenameTRG, std::string filenameRAW, std::string outfile, 
		       int npatt, int rate, int layer, int ladder, int module, std::string sector, 
		       bool RAW, bool TRIG, int npblock, bool conc,float L1prop,bool dbg)
{
  m_rate          = rate;
  m_lay           = layer;
  m_lad           = ladder;
  m_mod           = module;
  m_npblock       = npblock;
  m_tri_data      = new  std::vector<int>;
  m_raw_data      = new  std::vector<int>;
  m_raw_chip_bx   = new  std::vector<int>;;
  m_raw_chip_fifo = new  std::vector<int>;;
  m_L1prop        = L1prop; 
  m_write_raw     = RAW;
  m_write_trg     = TRIG;
  m_dbg           = dbg;

  m_CICsize       = 8; // Size of the CIC block (in BX)
  bend_bit_MPA    = 3; // MPA bend encoding length
  bend_bit_CBC    = 4; // CBC bend encoding length

  m_write_out     = false;

  evtbuilder::initVars();                                 // Initialize everything
  evtbuilder::convert(sector);                            // Get the list of involved chips
  evtbuilder::initTuple(filenameRAW,filenameTRG,outfile); // Initialize ROOT stuff
  evtbuilder::get_stores(npatt-conc*npatt%8,conc);        // Prepare the data store where to pickup info
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

  int layer,ladder,module,strip,chip,seg,pt,nstrip,nseg;

  int n_entries_TRG = L1TT->GetEntries();
  int n_entries_RAW = PIX->GetEntries();

  bool  isPS        = false;
  int   goodstub    = -1;
  float pTgen       = 0.;
  float d0gen       = 0.;
  float x,z;

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

  if (m_write_raw) store_size = n_entries_RAW;
  if (m_write_trg) store_size = n_entries_TRG;

  cout << "--> Entering loop 1, producing the big data stores for " 
       << store_size << " events..." << endl;

  for (int j=0;j<store_size;++j)
  {    
    if (j%50==0)
      cout << "Processing event " <<  j << "/" << store_size << endl;

    m_chip_trig.clear();
    m_chip_raw.clear();

    if (m_write_trg) 
    {
      if (rand()%1000>static_cast<int>(10*m_L1prop)) // Mbias event in the trigger block
      {
	mixpar=rand()%(n_entries_TRG-m_PHYsize); 

	L1TT->GetEntry(m_PHYsize+mixpar); 
      }
      else // Phys event (stored in the s2nd half of the TRG sample)
      {
	mixpar=rand()%m_PHYsize; 
	L1TT->GetEntry(mixpar); 

	cout << "Event " << j << " is QCD " << endl;
      }
    }

    if (m_write_raw) PIX->GetEntry(rand()%n_entries_RAW); // For the L1 it's much more simple

    // Initialize the map with dummy values for all referenced modules
    //
    // Please keep in mind that this map is 8 times larger in the FE output case
    // Maps are quite heavy, so don't produce a lot of FE events

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
      for (unsigned int i=0;i<m_concs.size();++i) // All the concentrators if one asks the CONC output (8 times less)
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
      for (int i=0;i<m_pix;++i) // Loop over pixels
      {
	isPS  = false;
	layer = m_pix_layer[i]; 
	if (layer!=m_lay && m_lay!=-1) continue; // By default we loop over all layers (-1)
	
	ladder= m_pix_ladder[i];
	if (ladder!=m_lad && m_lad!=-1) continue; // By default we loop over all ladders (-1)

	int disk=-1;
	if (layer>10) disk = (layer-11)%7;

	if (layer<8) isPS=true;
	if (layer>10 && disk<2 && ladder<10) isPS = true; 
	if (layer>10 && disk>=2 && ladder<7) isPS = true; 
	
	module= static_cast<int>(m_pix_module[i]); 
	if (module!=m_mod && m_mod!=-1) continue; // By default we loop over all modules (-1)

	seg   = m_pix_col[i]; 
	strip = m_pix_row[i]; 
	nseg  = m_pix_ncol[i]; 
	nstrip= m_pix_nrow[i]; 
	z     = m_pix_z[i];

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

	B_id = layer*1000000 + ladder*10000 + module*100 + chip; // Finally get the module ID corresponding to the map

	if (conc) B_id = B_id - chip%8; // Here we get the concentrator ID

	// Look if this chip has already been touched in this event
	m_iter = m_chip_raw.find(B_id);
	
	if (m_iter == m_chip_raw.end()) // Unknown chip, this can happen because csv file is only considering chips involved in TRG
	{
	  //	  continue; // We don't consider them for the moment...
	  m_digi_list.clear();

	  (conc) 
	    ? m_concs.push_back(B_id)
	    : m_chips.push_back(B_id);

	  m_digi_list.push_back(-1);
	  m_chip_raw.insert(std::make_pair(B_id,m_digi_list));
	  m_chip_trig.insert(std::make_pair(B_id,m_digi_list));
	  m_iter = m_chip_raw.find(B_id);
	}

	m_digi_list.clear();
	m_digi_list = m_iter->second;
	m_digi_list.push_back(chip%8);
	m_digi_list.push_back(strip);
	  
	(isPS) 
	  ? m_digi_list.push_back(seg%16)
	  : m_digi_list.push_back(-1);

	m_chip_raw.erase(m_iter->first);
	m_chip_raw.insert(std::make_pair(B_id,m_digi_list));      
      } // End of digi collection
    } // End of raw event collection

    // Start the trigger store building

    if (m_write_trg)
    {
      for (int i=0;i<m_stub;++i) // Loop over stubs
      {  
	// First of all we compute the ID of the stub's module
	isPS     = false;
	goodstub =-1;
	layer    = m_stub_layer[i]; 
	
	if (layer!=m_lay && m_lay!=-1) continue;
	
	ladder= m_stub_ladder[i]; 
	if (ladder!=m_lad && m_lad!=-1) continue;

	int disk=-1;
	if (layer>10) disk = (layer-11)%7;
	
	if (layer<8) isPS=true;
	if (layer>10 && disk<2 && ladder<10) isPS = true; 
	if (layer>10 && disk>=2 && ladder<7) isPS = true; 
	
	module= m_stub_module[i];
	if (module!=m_mod && m_mod!=-1) continue; // By default we loop over all modules (-1)

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
	  strip = int(2*m_stub_strip[i])%240;  // Bet 0 and 255
	}
	else
	{
	  strip = int(2*m_stub_strip[i])%254+1; // Code between 1 and 254 to avoid 00000000 position
	}

	pTgen = sqrt(m_stub_pxGEN[i]*m_stub_pxGEN[i]+m_stub_pyGEN[i]*m_stub_pyGEN[i]);
	d0gen = sqrt(m_stub_X0[i]*m_stub_X0[i]+m_stub_Y0[i]*m_stub_Y0[i]);

	if (pTgen>2 && d0gen<0.5) goodstub=1; // Good stubs are the ones looked at by L1 tracking

	pt   = static_cast<int>(std::abs(2*m_stub_deltas[i]));
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
      }  // End of stub collection
    } // End of trigger event collection     

    // Finally we add both collection to the stores
    
    m_data_trig.push_back(m_chip_trig);
    m_data_raw.push_back(m_chip_raw);
  } // End of storage loop

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

  std::vector<int> cicseq;
  bool inseq;

  // We build the TRG sequence, taking care that a given event does not end
  // up twice in the same CIC train

  for (int i=0;i<nevts;++i)
  { 
    if (i%8==0) cicseq.clear();

    inseq=true;

    while (inseq)
    {
      trg_evnum = rand()%n_trig;
      inseq=false;

      for (unsigned int j=std::max(0,int(cicseq.size())-16);j<cicseq.size();++j)
      {
	if (trg_evnum==cicseq.at(j)) inseq=true;
      }
    }

    cicseq.push_back(trg_evnum);
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
  int last_rank;
  int last_BX;

  m_raw_FIFO.clear();
  m_chip_FIFOs.clear();

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
      m_digi_list.clear();
      m_digi_list.push_back(0);
      m_digi_list.push_back(0);
      m_chip_FIFOs.insert(std::make_pair(m_concs.at(i),m_digi_list));
    }
  }
  else
  {
    for (unsigned int i=0;i<m_chips.size();++i)
    {
      m_digi_list.clear();
      m_digi_list.push_back(-1);
      m_raw_FIFO.insert(std::make_pair(m_chips.at(i),m_digi_list));
      m_digi_list.clear();
      m_digi_list.push_back(0);
      m_digi_list.push_back(0);
      m_chip_FIFOs.insert(std::make_pair(m_chips.at(i),m_digi_list));
    }
  }


  //
  // First of all one writes the L1 word, either at the FE of the concentrator level
  // 
  //
  //

  if (m_write_raw)
  {
    for (int i=0;i<nevts;++i) 
    { 
      if (i%100==0)
	cout << i << endl;

      if (i%3564==0) L1_id=0; // BCR, L1ID get back to 0

      evtbuilder::initVars();

      if (raw_seq.at(i)!=-1) // Chip has received received a L1 A, we write the raw block
      {
	++L1_id;
	
	m_chip_raw = m_data_raw.at(raw_seq.at(i)); // Pick up the event in the store

	for ( m_iter = m_chip_raw.begin(); m_iter != m_chip_raw.end();++m_iter )
	{

	  // First, look at the main info
	  isPS            = false;
	  m_raw_bx        = i;
	  m_raw_FIFO_FULL = 0;
	  m_raw_FIFO_SIZE = 0;
	  m_raw_chip      = m_iter->first;
	  m_raw_lay       = m_raw_chip/1000000;
	  m_raw_lad       = (m_raw_chip-1000000*m_raw_lay)/10000;
	  m_raw_mod       = (m_raw_chip-1000000*m_raw_lay-10000*m_raw_lad)/100;
	  m_raw_np        = 0;
	  m_raw_ns        = 0;

	  if (m_raw_lay<8 || (m_raw_lay>10 && m_raw_lad<9)) isPS = true; 

	  // Now, we get the digis contained in the event, and fill the block
	  m_digi_list = m_iter->second;
	  
	  // We have the info, we now write the word, making a difference bet. CIC and FE chip words

	  (conc)
	    ? evtbuilder::fill_CONC_RAW_block(m_digi_list,isPS,L1_id)
	    : evtbuilder::fill_RAW_block(m_digi_list,isPS,L1_id);
	  

	  // The word is written, we now emulate the FIFO behavior for the CIC 

	  m_raw_size  = m_raw_data->size();

	  // Concentrator FIFO size is 266*8*16 = 34048 bits
	  
	  int FIFO_size = 0.;

	  int extracted_bit_per_BX = (m_npblock/m_CICsize);

	  if (!conc) extracted_bit_per_BX = 16/2;

	  if (conc || isPS) // It's only for the CIC and MPA chips
	  {
	    // Find the chip FIFO footprint
	    // defined as follows:

	    // < CHIP_ID, <-1,BXin(1),BXout(1),BXin(2),BXout(2),... > >

	    m_iter2 = m_raw_FIFO.find(m_raw_chip); // Find the chip

	    FIFO_new.clear();
	    FIFO=m_iter2->second;
	    FIFO_size=(FIFO.size()-1)/2; // How many L1 events are in the FIFO?	    
	    FIFO_new.push_back(-1);

	    if (FIFO_size==0) // Empty FIFO, initialize it!
	    {
	      FIFO_new.push_back(i);                                    // BXin is the current BX, the L1 event enters the CIC
	      FIFO_new.push_back(i+m_raw_size/(extracted_bit_per_BX)+1); // BXout = BXin + numb of BXs to extract the event (if FIFO is empty it's simple) 
	      m_raw_FIFO_FULL = (FIFO_new.size()-1)/2;                  // Number of events in the FIFO (1)
	      FIFO_size=m_raw_size;
	      m_raw_FIFO_SIZE=m_raw_size;

	      last_BX=i+m_raw_size/(extracted_bit_per_BX)+1;
	    }
	    else  // FIFO not empty, one need to increment
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
	      if (FIFO_new.size()>1) last_rank = (last_BX-i)*(extracted_bit_per_BX);

	      FIFO_size=m_raw_size+last_rank;
	      m_raw_FIFO_SIZE=m_raw_size+last_rank;

	      //if (last_rank+m_raw_size>34048)   // The new event cannot fit in, we raise an error... 
	      if (last_rank+m_raw_size>100000000) // ...or not, as here, as we want to study the FIFO divergence 
	      {
		m_raw_FIFO_FULL = -1;
	      }
	      else // OK, come in!
	      {
		FIFO_new.push_back(i);                                          // BXin is the current BX, the L1 event enters the CIC
		FIFO_new.push_back(last_BX+m_raw_size/(extracted_bit_per_BX)+1); // BXout tells us when the event really goes out
		m_raw_FIFO_FULL = (FIFO_new.size()-1)/2;

		last_BX=last_BX+m_raw_size/(extracted_bit_per_BX)+1;
	      }
	    }

	    m_raw_FIFO.erase(m_iter2->first);
	    m_raw_FIFO.insert(std::make_pair(m_raw_chip,FIFO_new)); // Update the FIFO

	    // Then update the chip FIFO list

	    m_iter2 = m_chip_FIFOs.find(m_raw_chip);
	    FIFO=m_iter2->second;
	    FIFO.push_back(i);
	    FIFO.push_back(last_BX-i);
	    m_chip_FIFOs.erase(m_iter2->first);
	    m_chip_FIFOs.insert(std::make_pair(m_raw_chip,FIFO)); // Update the FIFO list
	  } // End of the FIFO emulator

	  if (m_dbg) m_raw_tree->Fill();
	}
      }
    } // End of the loop on raw events
  } // End of the RAW block definition


  // Then write the trigger data 
  
  int limCIC=16;
  int disk=0;

  if (m_write_trg)
  {
    if (conc) // For the concentrator we work on a m_CICsize BXs basis (default is 8)
    {
      for (int i=0;i<nevts/m_CICsize;++i) // Loop over the number of sequences
      { 
          if (i%100==0) cout << m_CICsize*i << endl;

          for (int ii=0;ii<6;++ii)
          {
              for (int j=0;j<40;++j)
	  {
	    for (int k=0;k<40;++k)
	    {
	      for (int m=0;m<3;++m)
	      {
		for (int n=0;n<3;++n)
		{
		  m_ovflow_b_nstubs[ii][j][k][m][n]   = 0.;
		  m_tri_b_nstubs[ii][j][k][m][n]   = 0.;
		}
	      }
	    }
	  }
	}
	
	
	for (int ii=0;ii<5;++ii)
	{
	  for (int j=0;j<15;++j)
	  {
	    for (int k=0;k<40;++k)
	    {
	      for (int m=0;m<3;++m)
	      {
		for (int n=0;n<3;++n)
		{
		  m_ovflow_d_nstubs[ii][j][k][m][n]   = 0.;
		  m_tri_d_nstubs[ii][j][k][m][n]   = 0.;
		}
	      }
	    }
	  }
	}


	m_chip_trig = m_data_trig.at(trig_seq.at(m_CICsize*i)); // Pick up the event in the store
	
	for ( m_iter = m_chip_trig.begin(); m_iter != m_chip_trig.end();++m_iter )
	{
	  trig_sequence.clear(); 

	  // The main info

	  isPS            = false;
	  m_tri_bx        = m_CICsize*i;
	  m_tri_nstubs    = 0;
	  m_tri_nstubs_g  = 0;
	  m_tri_nstubs_s  = 0;
	  m_tri_nstubs_gs = 0;
	  m_tri_chip      = m_iter->first;
	  m_tri_lay       = m_tri_chip/1000000;
	  m_tri_lad       = (m_tri_chip-1000000*m_tri_lay)/10000;
	  m_tri_mod       = (m_tri_chip-1000000*m_tri_lay-10000*m_tri_lad)/100;
	  
	  disk=-1;
	  if (m_tri_lay>10) disk = (m_tri_lay-11)%7;

	  if (m_tri_lay<8) isPS=true;

	  if (m_tri_lay>10 && disk<2 && m_tri_lad<10) isPS = true; 
	  if (m_tri_lay>10 && disk>=2 && m_tri_lad<7) isPS = true; 
	  
	  if (isPS)
	  {
	    limCIC=17;
	    if (m_tri_lay==5) limCIC=35;
	    if (m_tri_lay>10)
	    {
	      if (disk<2 && m_tri_lad<3) limCIC=35;
	    }
	  }
	  // Then we get all the stubs of the corresponding chip
	  // in trig_sequence

	  trig_sequence.push_back(m_iter->second);
	  
	  for (int j=1;j<m_CICsize;++j)
	  {
	    trig_sequence.push_back((m_data_trig.at(trig_seq.at(m_CICsize*i+j)).find(m_tri_chip))->second);
	  }
	
	  // And we write the word
	  
	  if (m_write_out) std::cout << i << " / " << m_tri_bx << " / " << m_tri_chip << " -> ";

	  evtbuilder::fill_TRG_block(trig_sequence,isPS,conc,m_CICsize*i,limCIC);

	  //
	  // The matrix m_tri_nstubs_f contains, for the chip considered,
	  // info concerning multiplicities, for a given bend value
	  //
	  // m_tri_nstubs_f[A][B][C]
	  // 
	  // A: stub bend, in half strips (0 to 39)
	  // B: multplicity before FE (0) after FE cut (1) and after CIC cut (2)
	  // C: good stub (1) or not (0). A good stub is induced by a primary with pt>2GeV/c

	  //  m_ovflow_(b/d)_nstubs[***][A][B][C] : number of block with overflows 
	  //  m_tri_(b/d)_nstubs[***][A][B][C]    : number of stubs in the block


	  if (m_tri_lay<=10)
	  {
	    for (int j=0;j<40;++j) 
	    {
	      m_ovflow_b_nstubs[m_tri_lay-5][m_tri_mod][j][0][1]++;
	      m_ovflow_b_nstubs[m_tri_lay-5][m_tri_mod][j][0][0]++;

	      if (m_tri_nstubs_f[j][0][0]>m_tri_nstubs_f[j][1][0]) // FE loss
		m_ovflow_b_nstubs[m_tri_lay-5][m_tri_mod][j][1][0]++;
	      if (m_tri_nstubs_f[j][0][0]>m_tri_nstubs_f[j][2][0]) // CIC loss
		m_ovflow_b_nstubs[m_tri_lay-5][m_tri_mod][j][2][0]++;

	      if (m_tri_nstubs_f[j][0][1]>m_tri_nstubs_f[j][1][1]) // FE loss
		m_ovflow_b_nstubs[m_tri_lay-5][m_tri_mod][j][1][1]++;
	      if (m_tri_nstubs_f[j][0][1]>m_tri_nstubs_f[j][2][1]) // CIC loss
		m_ovflow_b_nstubs[m_tri_lay-5][m_tri_mod][j][2][1]++;

	      for (int k=0;k<3;++k) 
	      {
		for (int l=0;l<2;++l) 
		  m_tri_b_nstubs[m_tri_lay-5][m_tri_mod][j][k][l] += m_tri_nstubs_f[j][k][l];
	      }
	    }
	  }
	  else
	  {
	    for (int j=0;j<40;++j) 
	    {
	      m_ovflow_d_nstubs[disk][m_tri_lad][j][0][0]++;
	      m_ovflow_d_nstubs[disk][m_tri_lad][j][0][1]++;

	      if (m_tri_nstubs_f[j][0][0]>m_tri_nstubs_f[j][1][0]) // FE loss
		m_ovflow_d_nstubs[disk][m_tri_lad][j][1][0]++;
	      if (m_tri_nstubs_f[j][0][0]>m_tri_nstubs_f[j][2][0]) // CIC loss
		m_ovflow_d_nstubs[disk][m_tri_lad][j][2][0]++;
	      if (m_tri_nstubs_f[j][0][1]>m_tri_nstubs_f[j][1][1]) // FE loss
		m_ovflow_d_nstubs[disk][m_tri_lad][j][1][1]++;
	      if (m_tri_nstubs_f[j][0][1]>m_tri_nstubs_f[j][2][1]) // CIC loss
		m_ovflow_d_nstubs[disk][m_tri_lad][j][2][1]++;

	      for (int k=0;k<3;++k) 
	      {
		for (int l=0;l<2;++l) 
		  m_tri_d_nstubs[disk][m_tri_lad][j][k][l] += m_tri_nstubs_f[j][k][l];
	      }
	    }
	  }
	}
      
	// Finally estimate the losses

	for (int ii=0;ii<6;++ii)
	{
	  for (int j=0;j<40;++j)
	  {
	    for (int k=0;k<40;++k)
	    {
	      m_los_b_nstubs[ii][j][k][1][0] = 0;
	      m_los_b_nstubs[ii][j][k][2][0] = 0;
	      m_los_b_nstubs[ii][j][k][1][1] = 0;
	      m_los_b_nstubs[ii][j][k][2][1] = 0;

	      if (m_tri_b_nstubs[ii][j][k][0][0]==0) continue;

	      m_los_b_nstubs[ii][j][k][1][0] = (1-m_tri_b_nstubs[ii][j][k][1][0]/m_tri_b_nstubs[ii][j][k][0][0]);
	      m_los_b_nstubs[ii][j][k][2][0] = (1-m_tri_b_nstubs[ii][j][k][2][0]/m_tri_b_nstubs[ii][j][k][0][0]);
	      
	      if (m_tri_b_nstubs[ii][j][k][0][1]!=0)
	      {
		m_los_b_nstubs[ii][j][k][1][1] = (1-m_tri_b_nstubs[ii][j][k][1][1]/m_tri_b_nstubs[ii][j][k][0][1]);
		m_los_b_nstubs[ii][j][k][2][1] = (1-m_tri_b_nstubs[ii][j][k][2][1]/m_tri_b_nstubs[ii][j][k][0][1]);
	      }
	    }
	  }
	}

	for (int ii=0;ii<5;++ii)
	{
	  for (int j=0;j<15;++j)
	  {
	    for (int k=0;k<40;++k)
	    {	      
	      m_los_d_nstubs[ii][j][k][1][0] = 0;
	      m_los_d_nstubs[ii][j][k][2][0] = 0;
	      m_los_d_nstubs[ii][j][k][1][1] = 0;
	      m_los_d_nstubs[ii][j][k][2][1] = 0;

	      if (m_tri_d_nstubs[ii][j][k][0][0]==0) continue;
	      
	      m_los_d_nstubs[ii][j][k][1][0] = (1-m_tri_d_nstubs[ii][j][k][1][0]/m_tri_d_nstubs[ii][j][k][0][0]);
	      m_los_d_nstubs[ii][j][k][2][0] = (1-m_tri_d_nstubs[ii][j][k][2][0]/m_tri_d_nstubs[ii][j][k][0][0]);
	      
	      if (m_tri_d_nstubs[ii][j][k][0][1]!=0)
	      {
		m_los_d_nstubs[ii][j][k][1][1] = (1-m_tri_d_nstubs[ii][j][k][1][1]/m_tri_d_nstubs[ii][j][k][0][1]);
		m_los_d_nstubs[ii][j][k][2][1] = (1-m_tri_d_nstubs[ii][j][k][2][1]/m_tri_d_nstubs[ii][j][k][0][1]);
	      }
	    }
	  }
	}

	m_tri_summary->Fill();

      }      
    }
    else // For the front end chips, we write on a 2BX basis
    {
      for (int i=0;i<nevts/2;++i)
      { 
	if (i%100==0) cout << i << endl;

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
	  
	  disk=-1;
	  if (m_tri_lay>10) disk = (m_tri_lay-11)%7;

	  if (m_tri_lay<8) isPS=true;

	  if (m_tri_lay>10 && disk<2 && m_tri_lad<10) isPS = true; 
	  if (m_tri_lay>10 && disk>=2 && m_tri_lad<7) isPS = true; 
	  
	  if (isPS) // MPA case, data is sent over 2BXs	  
	  {
	    m_tri_nstubs= 0;
	    m_tri_nstubs_g=0;
	    m_tri_nstubs_s=0;
	    m_tri_nstubs_gs=0;
	    trig_sequence.push_back((m_data_trig.at(trig_seq.at(2*i+1)).find(m_tri_chip))->second);
	    
	    std::cout << m_tri_bx << " / " << m_tri_chip << " -> ";
	    evtbuilder::fill_TRG_block(trig_sequence,isPS,conc,2*i,0);
	  }
	  else // CBC case, one BX basis
	  {
	    m_tri_nstubs= 0;
	    m_tri_nstubs_g=0;
	    m_tri_nstubs_s=0;
	    m_tri_nstubs_gs=0;
	    std::cout << m_tri_bx << " / " << m_tri_chip << " -> ";
	    evtbuilder::fill_TRG_block(trig_sequence,isPS,conc,2*i,0);
	    trig_sequence.clear(); 
	    m_tri_bx    = 2*i+1;
	    m_tri_nstubs= 0;
	    m_tri_nstubs_g=0;
	    m_tri_nstubs_s=0;
	    m_tri_nstubs_gs=0;
	    trig_sequence.push_back((m_data_trig.at(trig_seq.at(2*i+1)).find(m_tri_chip))->second);
	    std::cout << m_tri_bx << " / " << m_tri_chip << " -> ";
	    evtbuilder::fill_TRG_block(trig_sequence,isPS,conc,2*i+1,0);
	  }
	}
      }
    }
  } // End of the trigger writing block


  // We evaluate the drift function : BX_out-BX_in = a*BX + b
  //                                        X      = a*Z  + b
  //
  // The idea here is quite simple, we look if the time to extract 
  // an event from the CIC is drifting or not. If it's drifting, it means that 
  // the FIFO is filling up slowly, which is bad, definitely...
  //
  // Along with the error on these params

  
  std::vector<int> FIFO_cnt;

  double parZX[2][2];
  double resZX[2];
  double invZX[2][2];
  double detZX = 0;
  
  float a,b;
  float da,db;
  
  if (m_write_raw)
  {    	
    for ( m_iter = m_chip_FIFOs.begin(); m_iter != m_chip_FIFOs.end();++m_iter )
    {
      a = 0;
      b = 0;
      da = 0;
      db = 0;

      for (int i=0;i<2;++i)
      {
	resZX[i] = 0.;
	for (int j=0;j<2;++j) parZX[i][j] = 0.;
	for (int j=0;j<2;++j) invZX[i][j] = 0.;
      }
      
      detZX = 0;

      m_raw_chip_bx->clear();
      m_raw_chip_fifo->clear();
      m_raw_chip_slope=0.;
      m_raw_chip_slope_err=0.;
      m_raw_chip_int=0.;
      m_raw_chip_int_err=0.;
   
      FIFO_cnt.clear();
      m_raw_chip  = m_iter->first;
      FIFO_cnt    = m_iter->second;

      // Chip FIFO footprint
      // defined as follows:
      // < m_raw_chip , FIFO >
      // FIFO.push_back(i); // BX where the L1 event enters the FIFO
      // FIFO.push_back(last_BX-i); // Time spent by L1 in the FIFO


      for (unsigned int j=1;j<(FIFO_cnt.size())/2;++j)	    
      {
	//	cout  << j << "/" << FIFO_cnt.at(2*j) << "/" << FIFO_cnt.at(2*j+1) << endl;

	m_raw_chip_bx->push_back(FIFO_cnt.at(2*j));
	m_raw_chip_fifo->push_back(FIFO_cnt.at(2*j+1));

	x = FIFO_cnt.at(2*j+1); // BXout-BXin
	z = FIFO_cnt.at(2*j);                    // BXin

	parZX[0][0] += z*z; // S_ZZ
	parZX[1][1] += 1;   // S
	parZX[1][0] += z;   // S_Z
	
	resZX[0] += x*z;    // S_XZ
	resZX[1] += x;      // S_X

	//	cout << m_raw_chip << " / " << z << " / " << x << endl;
      }


      detZX = parZX[0][0]*parZX[1][1]-parZX[1][0]*parZX[1][0];
 
      invZX[1][1] = parZX[1][1]/detZX; 
      invZX[1][0] = parZX[1][0]/detZX; 
      invZX[0][0] = parZX[0][0]/detZX; 

      b     = invZX[0][0]*resZX[1] - invZX[1][0]*resZX[0]; 
      a     = invZX[1][1]*resZX[0] - invZX[1][0]*resZX[1]; 
 
      //  cout << m_raw_chip << " / Fifo_cnt =  " << a << "*BX+" << b<< endl;
      
      db = sqrt(invZX[0][0]);
      da = sqrt(invZX[1][1]);
      
      //      cout <<" a =  " << a << "+/-" << da << endl;
      //      cout <<" b =  " << b << "+/-" << db << endl;

      m_raw_chip_slope     = a;
      m_raw_chip_slope_err = da;

      m_raw_chip_int     = b;
      m_raw_chip_int_err = db;

      m_raw_summary->Fill();
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
  m_raw_chip_bx->clear();
  m_raw_chip_fifo->clear();
  m_raw_chip_slope=0.;
  m_raw_chip_slope_err=0.;
  m_raw_chip_int=0.;
  m_raw_chip_int_err=0.;


  for (int i=0;i<6;++i)
  {
    for (int j=0;j<40;++j)
    {
      for (int k=0;k<40;++k)
      {
	for (int m=0;m<3;++m)
	{
	  for (int n=0;n<3;++n)
	  {
	    m_ovflow_b_nstubs[i][j][k][m][n]   = 0.;
	    m_tri_b_nstubs[i][j][k][m][n]   = 0.;
	  }
	}
      }
    }
  }

  
  for (int i=0;i<5;++i)
  {
    for (int j=0;j<15;++j)
    {
      for (int k=0;k<40;++k)
      {
	for (int m=0;m<3;++m)
	{
	  for (int n=0;n<3;++n)
	  {
	    m_ovflow_d_nstubs[i][j][k][m][n]   = 0.;
	    m_tri_d_nstubs[i][j][k][m][n]   = 0.;
	  }
	}
      }
    }
  }

  for (int i=0;i<40;++i)
  {
    m_sw[i] = 0.5*i;
  }

}


void evtbuilder::initTuple(std::string inRAW,std::string inTRG,std::string out)
{
  m_outfile    = new TFile(out.c_str(),"recreate");

  if (m_dbg)
  {
    m_tri_tree    = new TTree("Trigger_FE","L1Trigger words after FE");
    m_raw_tree    = new TTree("Raw_FE","Raw data words after FE");

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
    m_raw_tree->Branch("RAW_FISIZE",     &m_raw_FIFO_SIZE,"RAW_FISIZE/I");
    m_raw_tree->Branch("RAW_NPCLUS",     &m_raw_np);
    m_raw_tree->Branch("RAW_NSCLUS",     &m_raw_ns);
  }    


  m_tri_summary = new TTree("Trigger_SUM","L1Trigger summary info");
  m_raw_summary = new TTree("Raw_SUM","Raw data summary info");

  m_tri_summary->Branch("BAR_LOSS",        &m_los_b_nstubs,    "BAR_LOSS[6][40][40][3][2]/F");
  m_tri_summary->Branch("DIS_LOSS",        &m_los_d_nstubs,    "DIS_LOSS[5][15][40][3][2]/F");
  m_tri_summary->Branch("BAR_MULT",        &m_tri_b_nstubs,    "BAR_MULT[6][40][40][3][2]/F");
  m_tri_summary->Branch("DIS_MULT",        &m_tri_d_nstubs,    "DIS_MULT[5][15][40][3][2]/F");
  m_tri_summary->Branch("BAR_OVFLOW_I",    &m_ovflow_b_nstubs,     "BAR_OVFLOW[6][40][40][3][2]/F");
  m_tri_summary->Branch("DIS_OVFLOW_I",    &m_ovflow_d_nstubs,     "DIS_OVFLOW[5][15][40][3][2]/F");
  m_tri_summary->Branch("SW",              &m_sw,              "SW[40]/F");

  m_raw_summary->Branch("RAW_CHP",           &m_raw_chip,    "RAW_CHP/I");
  m_raw_summary->Branch("RAW_CHP_BX",        &m_raw_chip_bx);
  m_raw_summary->Branch("RAW_CHP_FIFO",      &m_raw_chip_fifo);
  m_raw_summary->Branch("RAW_CHP_SLOPE",     &m_raw_chip_slope);
  m_raw_summary->Branch("RAW_CHP_SLOPE_ERR", &m_raw_chip_slope_err);
  m_raw_summary->Branch("RAW_CHP_INT",       &m_raw_chip_int);
  m_raw_summary->Branch("RAW_CHP_INT_ERR",   &m_raw_chip_int_err);

  L1TT   = new TChain("TkStubs");
  PIX    = new TChain("Pixels"); 
 
  // Input RAW data file
  
  if (m_write_raw)
  {
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
  }
 
  // Input TRG data file
  
  if (m_write_trg)
  {
    std::size_t found = inRAW.find(".root");
  
    // Case 1, it's a root file
    if (found!=std::string::npos)
    {
      L1TT->Add(inRAW.c_str());
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
	if (found!=std::string::npos) L1TT->Add(STRING.c_str());
      }
      
      in2.close();
    }

    m_PHYsize = L1TT->GetEntries();

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


   
  }
  
  if (m_write_raw)
  {
    pm_pix_layer=&m_pix_layer;
    pm_pix_ladder=&m_pix_ladder;
    pm_pix_module=&m_pix_module;
    pm_pix_row=&m_pix_row;
    pm_pix_col=&m_pix_col;
    pm_pix_bot=&m_pix_bot;
    pm_pix_nrow=&m_pix_nrow;
    pm_pix_ncol=&m_pix_ncol;
    pm_pix_z=&m_pix_z;

    PIX->SetBranchAddress("PIX_n",         &m_pix);
    PIX->SetBranchAddress("PIX_nPU",       &m_npu);
    PIX->SetBranchAddress("PIX_layer",     &pm_pix_layer);
    PIX->SetBranchAddress("PIX_ladder",    &pm_pix_ladder);
    PIX->SetBranchAddress("PIX_module",    &pm_pix_module);
    PIX->SetBranchAddress("PIX_row",       &pm_pix_row);
    PIX->SetBranchAddress("PIX_column",    &pm_pix_col);
    PIX->SetBranchAddress("PIX_nrow",      &pm_pix_nrow);
    PIX->SetBranchAddress("PIX_ncolumn",   &pm_pix_ncol);
    PIX->SetBranchAddress("PIX_bottom",    &pm_pix_bot);
    PIX->SetBranchAddress("PIX_z",         &pm_pix_z);
  }    

  if (m_write_trg)
  {
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
    pm_stub_z=&m_stub_z;
 
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
    L1TT->SetBranchAddress("L1TkSTUB_z",         &pm_stub_z);
  }
}





/////////////////////////////////////////////////////////////////////////////////
//
// ==> convert(std::string sectorfilename) 
//
// Here we retrieve info from the TKLayout CSV file containing the tracker modules detids
//
// The role of this method is to create a list of modules involved in the FE analysis
//
/////////////////////////////////////////////////////////////////////////////////

bool evtbuilder::convert(std::string sectorfilename) 
{
  int modid,lay,lad,mod,disk,type;

  bool m_tilted=true;

  //  std::cout << "Starting the conversion" << std::endl;

  int m_sec_mult = 0;

  int limits[6][3];
  int n_tilted_rings[6];
  int n_flat_rings[6];

  for (int i=0; i < 6; ++i) n_tilted_rings[i]=0;
  for (int i=0; i < 6; ++i) n_flat_rings[i]=0;

  if (m_tilted)
  {
    n_tilted_rings[0]=12;
    n_tilted_rings[1]=12;
    n_tilted_rings[2]=12;
    n_flat_rings[0]=7;
    n_flat_rings[1]=11;
    n_flat_rings[2]=15;
  }

  for (int i=0; i < 6; ++i)
  {
    for (int j=0; j < 3; ++j)
    {
      limits[i][j]=0;

      if (n_tilted_rings[i]==0) continue;

      limits[i][j]=(j%2)*n_flat_rings[i]+(j>0)*n_tilted_rings[i];
    }
  }

  std::vector<int> module;

  m_modules.clear();

  for (unsigned int i=0;i<230000;++i)
  {
    module.clear();
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
  bool found;

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
      found=false;

      std::bitset<32> detid = modid; // Le detid

      int rmodid; // Le modid que l'on utilise

      if (detid[25]) // barrel;
      {
	lay  = 8*detid[23]+4*detid[22]+2*detid[21]+detid[20]+4;

	type = 2*detid[19]+detid[18];

	if (type==3) // Pas tilté
	{
	  lad  = 128*detid[17]+64*detid[16]+32*detid[15]+16*detid[14]+
	    8*detid[13]+4*detid[12]+2*detid[11]+detid[10]-1;
	  mod  = 128*detid[9]+64*detid[8]+32*detid[7]+16*detid[6]+
	    8*detid[5]+4*detid[4]+2*detid[3]+detid[2]-1+limits[lay-5][type-1];
	}
	else // tilté
	{
	  mod  = 128*detid[17]+64*detid[16]+32*detid[15]+16*detid[14]+
	    8*detid[13]+4*detid[12]+2*detid[11]+detid[10]-1+limits[lay-5][type-1];
	  lad  = 128*detid[9]+64*detid[8]+32*detid[7]+16*detid[6]+
	    8*detid[5]+4*detid[4]+2*detid[3]+detid[2]-1;
	}
      }
      else // endcap
      {
	disk  = 8*detid[21]+4*detid[20]+2*detid[19]+detid[18];
	lay   = 10+disk+abs(2-(2*detid[24]+detid[23]))*7;
	lad   = 32*detid[17]+16*detid[16]+8*detid[15]+4*detid[14]+2*detid[13]+detid[12]-1;
	mod   = 128*detid[9]+64*detid[8]+32*detid[7]+16*detid[6]+
	  8*detid[5]+4*detid[4]+2*detid[3]+detid[2]-1;
      }

      if (m_lay!=-1 && lay!=m_lay) continue;
      if (lad!=m_lad && m_lad!=-1) continue;
      if (mod!=m_mod && m_mod!=-1) continue;

      rmodid = 10000*lay+100*lad+mod;

      module = m_modules.at(rmodid);
      module.push_back(m_sec_mult-2);

      m_modules.at(rmodid) = module;
    }
  }

  for (int i=0;i<230000;++i)
  {
    if (m_modules.at(i).size()==0) continue;

    for (int j=0;j<16;++j) m_chips.push_back(100*i+j);
    m_concs.push_back(100*i);
    m_concs.push_back(100*i+8);
  }

  in.close();

  m_sec_mult -= 2;

  cout << m_chips.size() << " CBC/MPA chips"<< endl;
  cout << m_concs.size() << " CIC chips "<< endl;

  return true;
}

//
// List of method writing the data blocks, according to the format defined in
// the following document:
//
// https://espace.cern.ch/Tracker-Upgrade/Electronics/CIC/Shared%20Documents/Specifications/CIC_specs_v2p1.pdf
//


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

  if (!spars) // CBC (unsparsified), slide 9 of the presentation cited
  {
    for (int j=0;j<256;++j)   m_raw_data->push_back(0);
    
    for (int j=0;j<ndata;++j)
    {
      m_raw_data->at(digis.at(3*j+2)+hsize) = 1;
    }
  }
  else // MPA
  {
    int data_MPA[128][17];
    int row,col;

    // The sparsification is done at the MPA level
    // Strips 0 to 120 belong to the pixel layer 
    // 121 to 239 belong to the strip layer

    for (int j=0;j<128;++j)
    {
      for (int i=0;i<17;++i) data_MPA[j][i] = 0;
    }

    for (int j=0;j<ndata;++j)
    {
      row = digis.at(3*j+2);
      col = digis.at(3*j+3);
    
      (row<120)
	? data_MPA[row][col]    = 1
	: data_MPA[row%120][16] = 1;
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
	    start_r=i+1; // 0000000 forbiden
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


    // Now encode them 
    np = clus_p.size()/3 ;
    ns = clus_s.size()/2; 
    m_raw_np        = np;
    m_raw_ns        = ns;


    np = std::min(31,np); // One cannot pass more than
    ns = std::min(31,ns); // 31 clusters out of the MPA

    std::bitset<5> N_S = ns;
    std::bitset<5> N_P = np;

    for (int j=0;j<5;++j) m_raw_data->push_back(N_S[4-j]);
    for (int j=0;j<5;++j) m_raw_data->push_back(N_P[4-j]);
    m_raw_data->push_back(0);

    for (int j=0;j<ns;++j)
    {
      std::bitset<7> row =  clus_s.at(2*j);
      std::bitset<3> wdt =  clus_s.at(2*j+1);

      for (int k=0;k<7;++k) m_raw_data->push_back(row[6-k]);
      for (int k=0;k<3;++k) m_raw_data->push_back(wdt[2-k]);
      m_raw_data->push_back(0); // The MIP bit
    }

    for (int j=0;j<np;++j)
    {
      std::bitset<7> row =  clus_p.at(3*j);
      std::bitset<4> col =  clus_p.at(3*j+1);
      std::bitset<3> wdt =  clus_p.at(3*j+2);

      for (int k=0;k<7;++k) m_raw_data->push_back(row[6-k]);
      for (int k=0;k<3;++k) m_raw_data->push_back(wdt[2-k]);
      for (int k=0;k<4;++k) m_raw_data->push_back(col[3-k]);
    }

    m_raw_data->push_back(0); // Trailer
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
 
  // %%%%%%%%%%%%%%%%
  // !! Still to clarify (SV 26/08/14), number of bits in the CBC word (252 or 256)
  // %%%%%%%%%%%%%%%%

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

    ns          = 0;
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
	    start_r=i+1;
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
    ns = std::min(127,ns);

    std::bitset<7> N_S = ns;

    for (int j=0;j<7;++j) m_raw_data->push_back(N_S[6-j]);

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
    // Strips 0 to 120 belong to the pixel layer 
    // 120 to 239 belong to the strip layer
 
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
 
      (row<120)
	? data_MPA[chip][row][col]    = 1
	: data_MPA[chip][row%120][16] = 1;
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
	      start_r=i+1;
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
    np = std::min(127,np);
    ns = std::min(127,ns);

    std::bitset<7> N_S = ns;
    std::bitset<7> N_P = np;

    //    std::cout << "This MPA chip contains

    for (int j=0;j<7;++j) m_raw_data->push_back(N_S[6-j]);
    for (int j=0;j<7;++j) m_raw_data->push_back(N_P[6-j]);

    for (int j=0;j<ns;++j)
    {
      std::bitset<7> row =  clus_s.at(3*j);
      std::bitset<3> wdt =  clus_s.at(3*j+1);
      std::bitset<3> chp =  clus_s.at(3*j+2);

      for (int k=0;k<3;++k) m_raw_data->push_back(chp[2-k]);
      for (int k=0;k<7;++k) m_raw_data->push_back(row[6-k]);
      for (int k=0;k<3;++k) m_raw_data->push_back(wdt[2-k]);
      m_raw_data->push_back(0);
    }

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
  }

  // Here we put a common trailer

  m_raw_data->push_back(0);


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


void evtbuilder::fill_TRG_block(std::vector< std::vector<int> > stubs, bool spars, bool conc, int BXid, int CIClim)
{
  m_tri_data->clear();

  // Then go for the data

  int nstubs[8][40]; // Number of stubs stored in each chip 
  int seqsize = static_cast<int>(stubs.size());

  int lim_CBC=3; // Max number of stubs passed by the CBC is 3/BX
  int lim_MPA=4; // Max number of stubs passed by the MPA is 4/2BXs

  if (bend_bit_MPA<4) lim_MPA=5;

  int gstub;

  std::vector<int> stublist;
  std::vector<float> stub_register;

  int sw=0;
  int position=0;
  int idx;
  std::vector<float> tmp_register;

  m_tri_nstubs    = 0;
  m_tri_nstubs_s  = 0;
  m_tri_nstubs_g  = 0;
  m_tri_nstubs_gs = 0;

  stub_register.clear();
  stublist.clear();
 
  for (int j=0;j<40;++j)  
  { 
    for (int i=0;i<3;++i)  
    { 
      m_tri_nstubs_f[j][i][0] = 0;
      m_tri_nstubs_f[j][i][1] = 0;
    }
  }


  if (conc) 
  {
    (spars)
      ? evtbuilder::fill_CONC_TRG_header(BXid%3564,1)
      : evtbuilder::fill_CONC_TRG_header(BXid%3564,0);
  }


  if (!conc) // Put a specific header for FE words  
  {
    m_tri_data->push_back(1); // Always start with synchro bit 1

    if (spars) // For the PS one puts the number of stubs in the first BX
    {          // will be completed at the end 
      for (int j=0;j<3;++j) m_tri_data->push_back(0);
    }
  }

  //  std::cout << "Writing trigger block for " << seqsize << " BXs" << std::endl;
    
  //  (conc)
  //    ? std::cout << "Concentrator" << std::endl
  //    : std::cout << "FE" << std::endl;
  
  int bit_bx = 3; // Number of bits necessary to precise the offset number in the CIC word
  if (m_CICsize>8) bit_bx=4;

  std::vector<int> sequence;
  sequence.clear();

  for (int i=0;i<seqsize;++i) // Loop over the block of events 
  {
    sequence.push_back(0);
  }

  for (int i=0;i<seqsize;++i) // Loop over the block of events 
  {
    //std::bitset<4> bx = i;
 
    stublist = stubs.at(i); // Get the stubs for event i
    
    // Here we initialize the stub per chip counters

    if (spars && i%2==0) // For the MPA, it's every 2BXs
    {
      for (int j=0;j<8;++j) 
      {
	for (int k=0;k<40;++k) nstubs[j][k]=0;
      }
    }

    if (spars && i%2==1 && !conc) // Here we have the number of stubs in the first BX for the MPA
    {
      std::bitset<3> cnt= m_tri_nstubs_s;
      for (int j=0;j<3;++j) m_tri_data->at(j+1)=cnt[2-j];
    }

    if (!spars) // For the CBC, it's every BX
    {
      for (int j=0;j<8;++j) 
      {
	for (int k=0;k<40;++k) nstubs[j][k]=0;
      }
    }

    // Then loop over the stubs contained in the event

    for (unsigned int kk=0;kk<(stublist.size()-1)/5;++kk)
    {
      gstub=stublist.at(5*kk+5); // Is the stub interesting or not
      sw   = std::min(stublist.at(5*kk+4),39);
      idx  = 100*kk+i;
      ++m_tri_nstubs;
      if (gstub>0) ++m_tri_nstubs_g;

      for (int j=sw;j<40;++j)
      {
	++m_tri_nstubs_f[j][0][0];
	if (gstub>0) ++m_tri_nstubs_f[j][0][1];
      }


      // Here we apply the chip limits

      if (nstubs[stublist.at(5*kk+1)][int(sw)]>=lim_CBC && !spars)      
      {
	m_tri_data->at(stublist.at(5*kk+1)+1)=1;
	continue; // No CHIP sorting for the moment
      }

      if (nstubs[stublist.at(5*kk+1)][int(sw)]>=lim_MPA && spars)
      {
	m_tri_data->at(stublist.at(5*kk+1)+1)=1;
	continue;  // No CHIP sorting for the moment
      }

      for (int j=sw;j<40;++j)
	++nstubs[stublist.at(5*kk+1)][j];

      ++m_tri_nstubs_s;
      if (gstub>0) ++m_tri_nstubs_gs;

      for (int j=sw;j<40;++j)
      {
	++m_tri_nstubs_f[j][1][0];
	if (gstub>0) ++m_tri_nstubs_f[j][1][1];
      }

      // Stubs are ordered by bend

      position = -1;
      tmp_register.clear();

      for (unsigned int kl=0;kl<stub_register.size()/3;++kl)
      {
	if (stub_register.at(3*kl+1)>=sw) 
	{
	  //	  std::cout << "?? " << kl << "/" << sw << std::endl; 

	  position = 3*kl;
	  break;
	}
      }

      if (position==-1)
      {
	stub_register.push_back(gstub);
	stub_register.push_back(sw);
	stub_register.push_back(idx);
      }
      else
      {
	for (unsigned int kl=0;kl<position;++kl) tmp_register.push_back(stub_register.at(kl));
	tmp_register.push_back(gstub);
	tmp_register.push_back(sw);
	tmp_register.push_back(idx);
	
	for (unsigned int kl=position;kl<stub_register.size();++kl) 
	  tmp_register.push_back(stub_register.at(kl));
	
	stub_register.clear();
	
	for (unsigned int kl=0;kl<tmp_register.size();++kl) 
	  stub_register.push_back(tmp_register.at(kl));
      }

      ++sequence.at(i);
      if (!conc)
      {
	//std::bitset<3> chp = stublist.at(5*kk+1);        // Chip number
	std::bitset<8> pos = stublist.at(5*kk+2);        // Stub position
	std::bitset<4> col = stublist.at(5*kk+3);        // Z position (for PS module)
	std::bitset<5> off = stublist.at(5*kk+4);        // Bend

	// Write the position
	for (int j=0;j<8;++j) 
	{
	  if (m_tri_data->size()==40 && spars) m_tri_data->push_back(0); // Second synchro bit
	  m_tri_data->push_back(pos[7-j]);
	}

	// The column for MPA side 
	if (spars)
	{ 
	  for (int j=0;j<4;++j)
	  {
	    if (m_tri_data->size()==40) m_tri_data->push_back(0); // Second synchro bit
	    m_tri_data->push_back(col[3-j]);
	  }
	}

	// Finally the bend
	for (int j=0;j<5;++j)
	{
	  if (m_tri_data->size()==40 && spars) m_tri_data->push_back(0); // Second synchro bit
	  m_tri_data->push_back(off[4-j]);
	}            
      }
    }
  }

  // For CIC we wrote after bend reordering

  if (conc)
  {
    for (unsigned int kl=0;kl<stub_register.size()/3;++kl)
    {
      sw =  stub_register.at(3*kl+1);
      sw = std::min(sw,39);

      for (int j=sw;j<40;++j)
      {
	if (m_tri_nstubs_f[j][2][0]>=CIClim) continue;
	++m_tri_nstubs_f[j][2][0];
	if (stub_register.at(3*kl)>0) ++m_tri_nstubs_f[j][2][1];
      }
      

      idx = stub_register.at(3*kl+2);

      int bxid = idx%100;
      idx = (idx-bxid)/100;

      std::bitset<3> bx  = bxid;
      std::bitset<3> chp = stubs.at(bxid).at(5*idx+1);        // Chip number
      std::bitset<8> pos = stubs.at(bxid).at(5*idx+2);        // Stub position
      std::bitset<4> col = stubs.at(bxid).at(5*idx+3);        // Z position (for PS module)
      std::bitset<5> off = stub_register.at(3*kl+1);          // Bend

      //      std::cout <<  kl << "," <<  std::abs(2*stub_register.at(3*kl+1)) << std::endl;
 
      // For the CIC, start with the offset and chip number

      for (int j=0;j<3;++j) m_tri_data->push_back(bx[bit_bx-1-j]);
      for (int j=0;j<3;++j) m_tri_data->push_back(chp[2-j]);
      

      // Then the position
      for (int j=0;j<8;++j) m_tri_data->push_back(pos[7-j]);

      // The column for MPA side 
      if (spars) 
      {
	for (int j=0;j<4;++j) m_tri_data->push_back(col[3-j]);
      }

      // Finally the bend
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
    }
  }
  //  std::cout << std::endl;

  // Potentially put padding bits here, only for the FE chips

  if (!conc) 
  {
    if (!spars)
    {
      for (int j=m_tri_data->size();j<40;++j) m_tri_data->push_back(0);
    }
    else
    {
      for (int j=m_tri_data->size();j<80;++j) m_tri_data->push_back(0);
    }
  }

  int nzeros=0;

  for (int i=0;i<seqsize;++i) // Loop over the block of events 
  {
    if (sequence.at(i)==0) ++nzeros;
  }

  if (conc) 
  {
    if (m_tri_nstubs_s>CIClim) 
    {
      m_tri_data->at(9) = 1; // Raise the CIC error bit if the number of stubs is in overflow

      for (int j=0;j<6;++j) m_tri_data->at(22+j) = 1;
    }
    else
    {
      std::bitset<6> nst = m_tri_nstubs_s;
      for (int j=0;j<6;++j) m_tri_data->at(22+j) = nst[5-j];
    }
  }


  m_tri_size=m_tri_data->size();
  if (m_dbg) m_tri_tree->Fill();


  //  std::cout << "Size of the trigger word / " << spars << " / " << conc << " / " << m_tri_size << std::endl;
  if (m_write_out && conc) 
  {
    std::cout << std::endl;
    std::cout << "HEADER  : ";

    for (unsigned int j=0;j<28;++j)
    {
      std::cout << m_tri_data->at(j);
    }
    std::cout << std::endl;


    std::cout << "PAYLOAD : "<< std::endl;

    for (unsigned int i=0;i<m_tri_nstubs_s;++i)
    {
      std::cout << "STUB " << i << "-> ";

      for (unsigned int j=28+i*(18+3*spars);j<28+(i+1)*(18+3*spars);++j)
      {
	std::cout << m_tri_data->at(j);
      }
      std::cout << std::endl;
    }

  }

}



void evtbuilder::fill_RAW_header_CBC(int L1id)
{
  // Format of the CBC L1 word header
  //
  // HHEEPPPPPPPPPLLLLLLLLL : HHEE (header + error) LL..LL (L1 ID bet 0 and 512)
  //                          PP..PP CBC pipeline address 

  for (int j=0;j<4;++j) m_raw_data->push_back(1); // HHEE

  if (L1id>512)
  {
    std::cout << "Too many L1ids, problem!!!" << std::endl;
  }

  std::bitset<9> L1_ID = L1id;

  for (int j=0;j<9;++j) m_raw_data->push_back(L1_ID[8-j]); // PP..PP
  for (int j=0;j<9;++j) m_raw_data->push_back(L1_ID[8-j]); // CC..CC
}


void evtbuilder::fill_RAW_header_MPA(int L1id)
{
  // Format of the MPA L1 word header
  //
  // 1111111111111111110EECCCCCCCCC0 : EE (error) CC..CC (L1 ID bet 0 and 512)
  //

  for (int j=0;j<18;++j) m_raw_data->push_back(1); 
  m_raw_data->push_back(0); 
  m_raw_data->push_back(0); 
  m_raw_data->push_back(0); 

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

void evtbuilder::fill_CONC_TRG_header(int BXid,int MPA)
{
  // Format of the CONC TRIGGER word header
  //
  // CSSSSSSSSSCCCCCCCCCCCC : SS..SS (status) CC..CC (BX ID bet 0 and 3564)
  //

  m_tri_data->push_back(MPA); 

  for (int j=0;j<9;++j) m_tri_data->push_back(0); // Status bits

  std::bitset<12> BX_ID = BXid;

  for (int j=0;j<12;++j) m_tri_data->push_back(BX_ID[11-j]); // CC..CC

  for (int j=0;j<6;++j) m_tri_data->push_back(0); // Let room for the number of stubs
}
