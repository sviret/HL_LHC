// Base class for Phase II tracker FrontEnd tests

// Main constructor

// Constructor for the concentrator

// Meaning of the parameters
//
//
// filenameRAW : the input data files for L1 block generation
// filenameTRG : the input data files for trigger block generation
// outfile     : name of the output root file
// npatt       : how many BXs should be produced (will be fitted correctly afterwards)
// rate        : input L1A rate (for the L1 raw block writing
// layer       : restrict info to a given layer, or not (give -1 in this case)
// sector      : the input file defining the modules
// RAW         : write the L1 block or not
// TRIG        : write the Trigger block or not
// npblock     : How many bits do we reserve to L1 in the concentrator block (between 8 and 320, should be a multiple of 8) 
// BMPA        : How many bits are used to code the stub bend on MPA concentrator level (0 to 5)
// BCBC        : How many bits are used to code the stub bend on CBC concentrator level (0 to 5)
// conc        : Are we building concentrator (true) or FE (false) data sequences
// L1prop      : The proportion of PU4T events you want in the TRG flow (in %)
// CICsize     : The size of the concentrator block, in BX

#include "losses.h"

losses::losses(std::string filenameRAW, std::string filenameTRG, std::string outfile, 
		       int npatt, int rate, int layer, std::string sector, 
		       bool RAW, bool TRIG, int npblock, int BMPA, int BCBC,
		       bool conc,int L1prop, int CICsize)
{
  m_rate       = rate;
  m_lay        = layer;
  m_npblock    = npblock;
  m_tri_data   = new  std::vector<int>;
  m_raw_data   = new  std::vector<int>;
  m_raw_chip_bx= new  std::vector<int>;;
  m_raw_chip_fifo= new  std::vector<int>;;
  m_L1prop     = L1prop; 
  m_CICsize    = CICsize; 
  m_write_raw  = RAW;
  m_write_trg  = TRIG;
  bend_bit_MPA = BMPA;
  bend_bit_CBC = BCBC;

  m_write_out  = false;
  int m_tower      = -1;
  
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

	    m_los_b_nstubs[i][j][k][m][n]   = 0.;
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

	    m_los_d_nstubs[i][j][k][m][n]   = 0.;
	  }
	}
      }
    }
  }

  losses::initVars();                                 // Initialize everything
  losses::convert(sector,m_tower);                    // Get the trigger tower module info
  losses::initTuple(filenameRAW,filenameTRG,outfile); // Initialize ROOT stuff
  losses::get_stores(npatt-conc*npatt%8,conc);        // Prepare the data store where to pickup info
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

void losses::get_stores(int nevts, bool conc)
{
 
  int B_id; // The detector module IDs (defined in the header)

  int layer,ladder,module,strip,chip,seg,pt;

  int n_entries_TRG = L1TT->GetEntries();
  int n_entries_RAW = PIX->GetEntries();

  bool  isPS        = false;
  int   goodstub    = -1;
  float pTgen       = 0.;
  float d0gen       = 0.;

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

//  for (int j=0;j<100;++j)
  for (int j=0;j<store_size;++j)
  {    
    if (j%50==0)
      cout << "Processing event " <<  j << "/" << store_size << endl;

    m_chip_trig.clear();
    

    if (m_write_trg) 
    {
      if (rand()%100>m_L1prop) // Mbias event in the trigger block
      {
	mixpar=rand()%(n_entries_TRG-m_PHYsize); 

	L1TT->GetEntry(m_PHYsize+mixpar); 
      }
      else // Phys event (stored in the s2nd half of the TRG sample)
      {
	mixpar=rand()%m_PHYsize; 
	L1TT->GetEntry(mixpar); 
      }
    }

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
	
	int disk=-1;
	if (layer>10) disk = (layer-11)%7;
	
	if (layer<8) isPS=true;
	if (layer>10 && disk<2 && ladder<10) isPS = true; 
	if (layer>10 && disk>=2 && ladder<7) isPS = true; 
	
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
	  strip = int(2*m_stub_strip[i])%240;  // Bet 0 and 255
	}
	else
	{
	  strip = int(2*m_stub_strip[i])%254+1; // Code between 1 and 254 to avoid 00000000 position
	}

	pTgen = sqrt(m_stub_pxGEN[i]*m_stub_pxGEN[i]+m_stub_pyGEN[i]*m_stub_pyGEN[i]);
	d0gen = sqrt(m_stub_X0[i]*m_stub_X0[i]+m_stub_Y0[i]*m_stub_Y0[i]);

	if (pTgen>2 && d0gen<0.5) goodstub=1; // Good stubs are the ones looked at by L1 tracking

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
      }  // End of stub collection
    } // End of trigger event collection     

    // Finally we add both collection to the stores
    
    m_data_trig.push_back(m_chip_trig);
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
  int n_raw  = m_data_trig.size();

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

  for (int i=0;i<nevts;++i)
  { 
    if (i%8==0) cicseq.clear();

    inseq=true;

    while (inseq)
    {
      trg_evnum = rand()%n_trig;
      inseq=false;

      //      cout << trg_evnum << "//" << inseq << endl;

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

  std::vector<int> FIFO,FIFO_new;

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


  // Then write the trigger data 

  int limCIC=16;
  int disk=0;

  if (m_write_trg)
  {
    if (conc) // For the concentrator we work on a m_CICsize BXs basis (default is 8)
    {
      for (int i=0;i<nevts/m_CICsize;++i) // Loop over the number of sequences
      { 
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

	if (i%100==0) cout << m_CICsize*i << endl;

	losses::initVars();     

	m_chip_trig = m_data_trig.at(trig_seq.at(m_CICsize*i)); // Pick up the event in the store
	
	for ( m_iter = m_chip_trig.begin(); m_iter != m_chip_trig.end();++m_iter )
	{
	  trig_sequence.clear(); 

	  // The main info

	  isPS            = false;
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
	 
	  //	  cout << m_tri_chip << " / " << limCIC << endl;
	  losses::fill_TRG_block(trig_sequence,isPS,limCIC);

	  if (m_tri_lay<=10)
	  {
	    for (int j=0;j<40;++j) 
	    {
	      m_ovflow_b_nstubs[m_tri_lay-5][m_tri_mod][j][0][0]++;
	      if (m_tri_nstubs[j][0][0]>m_tri_nstubs[j][1][0]) // FE loss
		m_ovflow_b_nstubs[m_tri_lay-5][m_tri_mod][j][1][0]++;
	      if (m_tri_nstubs[j][0][0]>m_tri_nstubs[j][2][0]) // CIC loss
		m_ovflow_b_nstubs[m_tri_lay-5][m_tri_mod][j][2][0]++;

	      m_ovflow_b_nstubs[m_tri_lay-5][m_tri_mod][j][0][1]++;
	      if (m_tri_nstubs[j][0][1]>m_tri_nstubs[j][1][1]) // FE loss
		m_ovflow_b_nstubs[m_tri_lay-5][m_tri_mod][j][1][1]++;
	      if (m_tri_nstubs[j][0][1]>m_tri_nstubs[j][2][1]) // CIC loss
		m_ovflow_b_nstubs[m_tri_lay-5][m_tri_mod][j][2][1]++;

	      for (int k=0;k<3;++k) 
	      {
		for (int l=0;l<2;++l) 
		  m_tri_b_nstubs[m_tri_lay-5][m_tri_mod][j][k][l] += m_tri_nstubs[j][k][l];
	      }
	    }
	  }
	  else
	  {
	    for (int j=0;j<40;++j) 
	    {
	      m_ovflow_d_nstubs[disk][m_tri_lad][j][0][0]++;
	      if (m_tri_nstubs[j][0][0]>m_tri_nstubs[j][1][0]) // FE loss
		m_ovflow_d_nstubs[disk][m_tri_lad][j][1][0]++;
	      if (m_tri_nstubs[j][0][0]>m_tri_nstubs[j][2][0]) // CIC loss
		m_ovflow_d_nstubs[disk][m_tri_lad][j][2][0]++;

	      m_ovflow_d_nstubs[disk][m_tri_lad][j][0][1]++;
	      if (m_tri_nstubs[j][0][1]>m_tri_nstubs[j][1][1]) // FE loss
		m_ovflow_d_nstubs[disk][m_tri_lad][j][1][1]++;
	      if (m_tri_nstubs[j][0][1]>m_tri_nstubs[j][2][1]) // CIC loss
		m_ovflow_d_nstubs[disk][m_tri_lad][j][2][1]++;

	      for (int k=0;k<3;++k) 
	      {
		for (int l=0;l<2;++l) 
		  m_tri_d_nstubs[disk][m_tri_lad][j][k][l] += m_tri_nstubs[j][k][l];


		//		cout << i << " / " << m_tri_d_nstubs[disk][m_tri_lad][j][k][0] << endl;
	      }
	    }
	  }
	}
      
	for (int i=0;i<6;++i)
	{
	  for (int j=0;j<40;++j)
	  {
	    for (int k=0;k<40;++k)
	    {
	      m_los_b_nstubs[i][j][k][1][0] = 0;
	      m_los_b_nstubs[i][j][k][2][0] = 0;
	      m_los_b_nstubs[i][j][k][1][1] = 0;
	      m_los_b_nstubs[i][j][k][2][1] = 0;
	      m_lp_b_nstubs[i][j][k][1][0] = 0;
	      m_lp_b_nstubs[i][j][k][2][0] = 0;
	      m_lp_b_nstubs[i][j][k][1][1] = 0;
	      m_lp_b_nstubs[i][j][k][2][1] = 0;

	      if (m_tri_b_nstubs[i][j][k][0][0]==0) continue;

	      m_los_b_nstubs[i][j][k][1][0] = (1-m_tri_b_nstubs[i][j][k][1][0]/m_tri_b_nstubs[i][j][k][0][0]);
	      m_los_b_nstubs[i][j][k][2][0] = (1-m_tri_b_nstubs[i][j][k][2][0]/m_tri_b_nstubs[i][j][k][0][0]);
	      
	      if (m_tri_b_nstubs[i][j][k][0][1]!=0)
	      {
		m_los_b_nstubs[i][j][k][1][1] = (1-m_tri_b_nstubs[i][j][k][1][1]/m_tri_b_nstubs[i][j][k][0][1]);
		m_los_b_nstubs[i][j][k][2][1] = (1-m_tri_b_nstubs[i][j][k][2][1]/m_tri_b_nstubs[i][j][k][0][1]);
	      }

	      m_lp_b_nstubs[i][j][k][1][0] = (m_los_b_nstubs[i][j][k][1][0]>0);
	      m_lp_b_nstubs[i][j][k][2][0] = (m_los_b_nstubs[i][j][k][2][0]>0);
	      m_lp_b_nstubs[i][j][k][1][1] = (m_los_b_nstubs[i][j][k][1][1]>0);
	      m_lp_b_nstubs[i][j][k][2][1] = (m_los_b_nstubs[i][j][k][2][1]>0);
	    }
	  }
	}

	//	std::cout << (1-m_tri_b_nstubs[0][19][39][2][1]/m_tri_b_nstubs[0][19][39][0][1]) << std::endl;


	for (int i=0;i<5;++i)
	{
	  for (int j=0;j<15;++j)
	  {
	    for (int k=0;k<40;++k)
	    {
	      /*
	      if (m_tri_d_nstubs[i][j][k][0][0]==0) continue;
	      
	      m_los_d_nstubs[i][j][k][1][0] += 1/fact*(1-m_tri_d_nstubs[i][j][k][1][0]/m_tri_d_nstubs[i][j][k][0][0]);
	      m_los_d_nstubs[i][j][k][2][0] += 1/fact*(1-m_tri_d_nstubs[i][j][k][2][0]/m_tri_d_nstubs[i][j][k][0][0]);
	      
	      if (m_tri_d_nstubs[i][j][k][0][1]!=0)
	      {
		m_los_d_nstubs[i][j][k][1][1] += 1/fact*(1-m_tri_d_nstubs[i][j][k][1][1]/m_tri_d_nstubs[i][j][k][0][1]);
		m_los_d_nstubs[i][j][k][2][1] += 1/fact*(1-m_tri_d_nstubs[i][j][k][2][1]/m_tri_d_nstubs[i][j][k][0][1]);
	      }
	      */
	      
	      m_los_d_nstubs[i][j][k][1][0] = 0;
	      m_los_d_nstubs[i][j][k][2][0] = 0;
	      m_los_d_nstubs[i][j][k][1][1] = 0;
	      m_los_d_nstubs[i][j][k][2][1] = 0;
	      m_lp_d_nstubs[i][j][k][1][0] = 0;
	      m_lp_d_nstubs[i][j][k][2][0] = 0;
	      m_lp_d_nstubs[i][j][k][1][1] = 0;
	      m_lp_d_nstubs[i][j][k][2][1] = 0;

	      if (m_tri_d_nstubs[i][j][k][0][0]==0) continue;
	      
	      m_los_d_nstubs[i][j][k][1][0] = (1-m_tri_d_nstubs[i][j][k][1][0]/m_tri_d_nstubs[i][j][k][0][0]);
	      m_los_d_nstubs[i][j][k][2][0] = (1-m_tri_d_nstubs[i][j][k][2][0]/m_tri_d_nstubs[i][j][k][0][0]);
	      
	      if (m_tri_d_nstubs[i][j][k][0][1]!=0)
	      {
		m_los_d_nstubs[i][j][k][1][1] = (1-m_tri_d_nstubs[i][j][k][1][1]/m_tri_d_nstubs[i][j][k][0][1]);
		m_los_d_nstubs[i][j][k][2][1] = (1-m_tri_d_nstubs[i][j][k][2][1]/m_tri_d_nstubs[i][j][k][0][1]);
	      }


	      m_lp_d_nstubs[i][j][k][1][0] = (m_los_d_nstubs[i][j][k][1][0]>0);
	      m_lp_d_nstubs[i][j][k][2][0] = (m_los_d_nstubs[i][j][k][2][0]>0);
	      m_lp_d_nstubs[i][j][k][1][1] = (m_los_d_nstubs[i][j][k][1][1]>0);
	      m_lp_d_nstubs[i][j][k][2][1] = (m_los_d_nstubs[i][j][k][2][1]>0);
	    }
	  }
	}

	m_tri_tree->Fill();

      }
    }
  } // End of the trigger writing block





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

void losses::initVars()
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


void losses::initTuple(std::string inRAW,std::string inTRG,std::string out)
{
  m_outfile    = new TFile(out.c_str(),"recreate");

  m_tri_tree    = new TTree("Trigger_Loss","L1Trigger words after FE");
  m_raw_tree    = new TTree("Raw_FE","Raw data words after FE");
  m_raw_summary = new TTree("Raw_SUM","Raw data summary info");


  m_tri_tree->Branch("BAR_LOSS",        &m_los_b_nstubs,    "BAR_LOSS[6][40][40][3][2]/F");
  m_tri_tree->Branch("DIS_LOSS",        &m_los_d_nstubs,    "DIS_LOSS[5][15][40][3][2]/F");
  m_tri_tree->Branch("BAR_MULT",        &m_tri_b_nstubs,    "BAR_MULT[6][40][40][3][2]/F");
  m_tri_tree->Branch("DIS_MULT",        &m_tri_d_nstubs,    "DIS_MULT[5][15][40][3][2]/F");
  m_tri_tree->Branch("BAR_OVFLOW_I",    &m_ovflow_b_nstubs,     "BAR_OVFLOW[6][40][40][3][2]/F");
  m_tri_tree->Branch("DIS_OVFLOW_I",    &m_ovflow_d_nstubs,     "DIS_OVFLOW[5][15][40][3][2]/F");
  m_tri_tree->Branch("SW",              &m_sw,              "SW[40]/F");


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
// Here we retrieve info from the TKLayout CSV file containing the sector definition
//
// This file contains, for each sector, the ids of the modules and chips contained in the sector sec_num
//
// The role of this method is to create the opposite, ie a vector containing, for every module the list of sectors belonging to it
//
/////////////////////////////////////////////////////////////////////////////////

bool losses::convert(std::string sectorfilename, int sec_num) 
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

  while (!in.eof()) 
  {
    ++m_sec_mult;
    getline(in,STRING);

    if (m_sec_mult<2) continue;
    if (m_sec_mult-2!=sec_num && sec_num!=-1) continue;

    std::istringstream ss(STRING);
    npar = 0;

    while (ss)
    {
      std::string s;
      if (!getline( ss, s, ',' )) break;

      ++npar;
      if (npar<=2) continue;

      modid = atoi(s.c_str());

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

	//	if (disk>=3 && m_tilted) lad += 2;

	mod  = 128*detid[9]+64*detid[8]+32*detid[7]+16*detid[6]+
	  8*detid[5]+4*detid[4]+2*detid[3]+detid[2]-1;
      }

      if (m_lay!=-1 && lay!=m_lay) continue;

      rmodid = 10000*lay+100*lad+mod;

      //      if (m_sec_mult-2==0)
      //      std::cout << modid << " / " << rmodid << std::endl; 

      module = m_modules.at(rmodid);
      module.push_back(m_sec_mult-2);

      m_modules.at(rmodid) = module;
    }
  }

  //  std::cout << "Found " << m_modules.size() << " modules" << endl;


  for (int i=0;i<230000;++i)
  {
    if (m_modules.at(i).size()==0) continue;

    //    std::cout << i << endl;

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
// 3: Trigger block
//


void losses::fill_TRG_block(std::vector< std::vector<int> > stubs, bool ps, int lim_CIC)
{
  int seqsize = static_cast<int>(stubs.size());

  int lim_CBC=3; // Max number of stubs passed by the CBC is 3/BX
  int lim_MPA=5; // Max number of stubs passed by the MPA is 5/2BXs



  int gstub;

  std::vector<int> stublist;
  std::vector<float> stub_register;

  int sw=0;
  int position=0;
  std::vector<float> tmp_register;

  stub_register.clear();
  stublist.clear();
 
  for (int j=0;j<40;++j)  
  { 
    for (int i=0;i<3;++i)  
    { 
      m_tri_nstubs[j][i][0] = 0;
      m_tri_nstubs[j][i][1] = 0;
    }
  }

  int nstubs[8][40];

  for (int i=0;i<seqsize;++i) // Loop over the block of events 
  { 
    stublist = stubs.at(i); // Get the stubs for event i
 

    if (ps && i%2==0) // For the MPA, it's every 2BXs
    {
      for (int j=0;j<8;++j) 
      {
	for (int k=0;k<40;++k) nstubs[j][k]=0;
      }
    }

    if (!ps) // For the CBC, it's every BX
    {
      for (int j=0;j<8;++j) 
      {
	for (int k=0;k<40;++k) nstubs[j][k]=0;
      }
    }

    // Then loop over the stubs contained in the i-th event

    for (unsigned int kk=0;kk<(stublist.size()-1)/5;++kk)
    {
      gstub=stublist.at(5*kk+5); // Is the stub interesting or not
      sw   = std::min(stublist.at(5*kk+4),39);

      for (int j=sw;j<40;++j)
      {
	++m_tri_nstubs[j][0][0];
	if (gstub>0) ++m_tri_nstubs[j][0][1];
      }

      // Here we apply the chip limits

      if (nstubs[stublist.at(5*kk+1)][int(sw)]>=lim_CBC && !ps) continue; // No CHIP sorting for the moment
      if (nstubs[stublist.at(5*kk+1)][int(sw)]>=lim_MPA && ps) continue;  // No CHIP sorting for the moment


      //      if (nstubs[stublist.at(5*kk+1)][39]>=lim_CBC && !ps) continue; // No CHIP sorting for the moment
      //      if (nstubs[stublist.at(5*kk+1)][39]>=lim_MPA && ps) continue;  // No CHIP sorting for the moment


      for (int j=sw;j<40;++j)
	++nstubs[stublist.at(5*kk+1)][j];

      for (int j=sw;j<40;++j)
      {
	++m_tri_nstubs[j][1][0];
	if (gstub>0) ++m_tri_nstubs[j][1][1];
      }

      position = -1;
      tmp_register.clear();

      for (unsigned int kl=0;kl<stub_register.size()/2;++kl)
      {
	if (stub_register.at(2*kl+1)>=sw) 
	{
	  position = 2*kl;
	  break;
	}
      }

      if (position==-1)
      {
	stub_register.push_back(gstub);
	stub_register.push_back(sw);
      }
      else
      {
	for (int kl=0;kl<position;++kl) tmp_register.push_back(stub_register.at(kl));
	tmp_register.push_back(gstub);
	tmp_register.push_back(sw);
	for (int kl=position;kl<static_cast<int>(stub_register.size());++kl) tmp_register.push_back(stub_register.at(kl));
	stub_register.clear();
	for (int kl=0;kl<static_cast<int>(tmp_register.size());++kl) stub_register.push_back(tmp_register.at(kl));
      }
    }
  }

  for (unsigned int kl=0;kl<stub_register.size()/2;++kl)
  {
    // if (kl>=lim_CIC) continue;
    sw =  stub_register.at(2*kl+1);
    sw = std::min(sw,39);

    for (int j=sw;j<40;++j)
    {
      if (m_tri_nstubs[j][2][0]>=lim_CIC) continue;
      ++m_tri_nstubs[j][2][0];
      if (stub_register.at(2*kl)>0) ++m_tri_nstubs[j][2][1];
    }
  }
  
  
  /*
  for (int j=0;j<40;++j)  
  { 
    if (m_tri_nstubs[j][0][0]!=m_tri_nstubs[j][1][0] && m_tri_nstubs[j][1][0]<=lim_MPA)
    {
      std::cout << j << " : " << m_tri_nstubs[j][0][0] << "//" 
		<< m_tri_nstubs[j][1][0] << "//" 
		<< m_tri_nstubs[j][2][0] << "%%%" 
		<< m_tri_nstubs[j][0][1] << "//" 
		<< m_tri_nstubs[j][1][1] << "//" 
		<< m_tri_nstubs[j][2][1] << "//" 
		<< stub_register.size()/2<< std::endl;
    }
  }
  */
  //  m_tri_tree->Fill();
}


void losses::fill_CONC_TRG_header(int BXid,int MPA)
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
