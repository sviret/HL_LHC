// Base class for Phase II tracker FrontEnd tests

// Main constructor

// Constructor for the concentrator

// Meaning of the parameters
//
//
// filenameRAW : the input data files for L1 block generation
// outfile     : name of the output root file
// sectorfile  : name of the csv file containing the module info (to build chip lists)
// npatt       : how many BXs should be produced (will be fitted correctly afterwards)
// rate        : input L1A rate (for the L1 raw block writing), given in kHz
// layer       : restrict info to a given layer, or not (give -1 in this case)
// ladder      : restrict info to a given ladder, or not (give -1 in this case)
// module      : restrict info to a given module, or not (give -1 in this case)
// npblock     : How many bits do we reserve to L1 in the concentrator block (between 8 and 320, should be a multiple of 8) 
//             : in FE chips this value is always 8


#include "l1_builder.h"

l1_builder::l1_builder(std::string filenameRAW, std::string outfile, std::string sectorfile,
		       int npatt, int rate, int layer, int ladder, int module, int npblock, int error)
{
    m_lay        = layer;
    m_lad        = ladder;
    m_mod        = module;
    m_rate       = rate;
    m_npblock    = npblock;
    
    m_MPA_L1_delay = 7;  // The minimal delay between L1A reception by the MPA and transmission to the CIC
    m_CBC_L1_delay = 20; // The minimal delay between L1A reception by the CBC and transmission to the CIC
    m_CIC_L1_delay = 0;  // The minimal delay between L1A reception by the CIC and transmission to the GBT

    m_MPA_FIFO_depth = 4;
    m_CIC_FIFO_depth = 1000;
    
    m_MPA_FIFO_size  = 100000;
    m_CIC_FIFO_size  = 100000;
    
    m_raw_data     = new  std::vector<int>;
    m_raw_chip_fifo= new  std::vector<int>;
    m_raw_chip_bx  = new  std::vector<int>;
    m_raw_conc_fifo= new  std::vector<int>;
    m_raw_conc_bx  = new  std::vector<int>;
    
    // ERROR types : ABCDE
    // A: random errors (flip bits) out of the L1 data (SEU)  
    // B: header errors (flip bits)   
    // C: cluster size errors (flip bits)   
    // D: not used
    // E: not used

    m_error_FE     = error;
    m_error_prop   = 0.001; // Proportion of errors

    m_error_rdm    = m_error_FE/10000;
    m_error_hdr    = (m_error_FE-10000*m_error_rdm)/1000;
    m_error_siz    = (m_error_FE-10000*m_error_rdm-1000*m_error_hdr)/100;

    unsp = false;

    if (error==2) unsp=true; // Turn on unsparsification for the CBC


    l1_builder::initVars();                       // Initialize everything
    l1_builder::convert(sectorfile);              // Get the trigger tower module info
    l1_builder::initTuple(filenameRAW,outfile);   // Initialize ROOT stuff
    l1_builder::get_stores(npatt);                // Prepare the data store where to pickup info

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

void l1_builder::get_stores(int nevts)
{
 
    int B_id, B_id_conc; // The detector module IDs (defined in the header)

    int layer,ladder,module,strip,chip,seg;

    int n_entries_RAW = PIX->GetEntries();

    bool  isPS        = false;

    // Then loop over events

    std::vector<int> m_digi_list;

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

    store_size = n_entries_RAW;

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
    
    //for (int j=0;j<10;++j)
    for (int j=0;j<store_size;++j)
    {
        if (j%1==20)
            cout << "Processing event " <<  j << "/" << store_size << endl;

        m_chip_raw.clear();
        m_conc_raw.clear();
        
        mixpar = rand()%n_entries_RAW;

        PIX->GetEntry(mixpar); // For the L1 it's much more simple

        // Initialize the map with dummy values for all referenced modules
        //
        // Please keep in mind that this map is 8 times larger in the FE output case
        // Maps are quite heavy, so don't produce a lot of FE events
        
        for (unsigned int i=0;i<m_chips.size();++i) // All the chips if one asks the FE output
        {
            m_digi_list.clear();
            m_digi_list.push_back(mixpar);
            m_chip_raw.insert(std::make_pair(m_chips.at(i),m_digi_list));
        }
    
        for (unsigned int i=0;i<m_concs.size();++i) // All the concentrators if one asks the CONC output (8 times less)
        {
            m_digi_list.clear();
            m_digi_list.push_back(mixpar);
            m_conc_raw.insert(std::make_pair(m_concs.at(i),m_digi_list));
        }

        
        for (int i=0;i<m_pix;++i) // Loop over pixels
        {
            isPS  = false;
            layer = m_pix_layer[i];
	
            if (layer!=m_lay && m_lay!=-1) continue; // By default we loop over all layers (-1)
	
            ladder= m_pix_ladder[i]-1;
            
            if (ladder!=m_lad && m_lad!=-1) continue; // By default we loop over all ladders (-1)
            if (layer<8 || (layer>10 && ladder<9)) isPS = true;
            
            module= static_cast<int>((m_pix_module[i]-1)/2);
            
            if (module!=m_mod && m_mod!=-1) continue; // By default we loop over all modules (-1)
            
            seg   = m_pix_col[i];
            strip = m_pix_row[i];
	
            
            ///////////////////////
            //////// SV 27/08/15: SHOULD BE CROSSCHECKED !!!!
            ///////////////////////
            
            if (isPS && m_pix_module[i]%2==1) // Ger the chip number for the PS-P
            {
                chip  = static_cast<int>(strip/120)+(seg/16)*8;
                strip = strip%120+((m_pix_module[i]-1)%2)*120;
            }
            else if (isPS && m_pix_module[i]%2==0) // Ger the chip number for the PS-S
            {
                chip  = static_cast<int>(strip/120)+seg*8;
                strip = strip%120+((m_pix_module[i]-1)%2)*120;
            }
            else // For the 2S
            {
                chip  = static_cast<int>(strip/127)+seg*8;
                strip = strip%127+((m_pix_module[i]-1)%2)*127;
            }

            ///////////////////////
            ///////////////////////
            
            B_id      = layer*1000000 + ladder*10000 + module*100 + chip; // Finally get the module ID corresponding to the map
            B_id_conc = B_id - chip%8; // Here we get the concentrator ID

            // Look if this chip has already been touched in this event
            m_iter  = m_chip_raw.find(B_id);
	    m_iter2 = m_conc_raw.find(B_id_conc);
            
            if (m_iter == m_chip_raw.end()) // Unknown chip, this can happen because csv file is only considering chips involved in TRG
            {
                m_digi_list.clear();
                m_digi_list.push_back(-1);
                m_chip_raw.insert(std::make_pair(B_id,m_digi_list));
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
            
            if (m_iter2 == m_conc_raw.end()) // Unknown chip, this can happen because csv file is only considering chips involved in TRG
            {
                m_digi_list.clear();
                m_digi_list.push_back(-1);
                m_conc_raw.insert(std::make_pair(B_id_conc,m_digi_list));
                m_iter = m_conc_raw.find(B_id_conc);
            }
            
            m_digi_list.clear();
            m_digi_list = m_iter2->second;
            m_digi_list.push_back(chip%8);
            m_digi_list.push_back(strip);
            
            (isPS)
            ? m_digi_list.push_back(seg%16)
            : m_digi_list.push_back(-1);
            
            m_conc_raw.erase(m_iter2->first);
            m_conc_raw.insert(std::make_pair(B_id_conc,m_digi_list));

        } // End of digi collection
        
        m_data_raw.push_back(m_chip_raw);
        m_data_raw.push_back(m_conc_raw);
        
    } // End of raw event collection

    // Now we have the collection containing the raw data
    // Raw data is sent every n BX, where n depends on the L1 rate
    // If the L1 rate is given in kHz, we have n=40000/rate
    //

    cout << "--> Entering loop 2, producing the random L1 sequence..." << endl;

    cout << "Size of the store: " << m_data_raw.size()/2 << endl;

    int delta    = static_cast<int>(40000/m_rate);

    int n_raw  = m_data_raw.size()/2;

    std::default_random_engine generator(time(NULL));
    std::normal_distribution<double> L1s(delta,delta/3);

    std::cout << "Average L1A period will be " << delta
            << " BXs with a sigma of " << delta/3 << " BXs" << std::endl;

    srand(time(NULL));

    std::vector<int> raw_seq;
    raw_seq.clear();

    // First we create random sequence for L1 satisfying the trigger rules

    // Rules are following the ones described in part 2.3 of following ref (2009):
    //
    // http://arxiv.org/pdf/0911.5422.pdf
    //
    // The last rule (4 in 240) is discarded as it prevents 1MHz L1A rate.
    
    
    // The first L1 is always sent at BX=0

    raw_seq.push_back(rand()%n_raw);

    int rank_r1 = 0;
    int rank_r2 = 0;
    int rank_r3 = 0;


    int lim_r1 = 3;
    int lim_r2 = 25;
    int lim_r3 = 100;

    //    int lim_r1 = 0;
    //    int lim_r2 = 0;
    //    int lim_r3 = 0;

    bool L1_rule = true;

    cout << "... Generate the sequence of " << nevts << " by picking up events in the store..." << std::endl;

    while (rank_r1<nevts)
    {
        L1_rule = false;

        // First we need to generate a space which is not enforcing the trigger rules
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
    // of events, with the correspondind data for each FE/Concentrator chip
    //
    // L1A signal is received at the same time in all the chips.
    //
    // FIFO behavior is the key here
    // FIFOs have a size and depth. Depth is the number of events which can be stored, and size is the size allocated for each event (sparsified)
    //
    // When an L1A is received, there are two possibilities:
    //
    // 1/ FIFO is not full: the sparsified event is stored and waits for its extraction
    // 2/ FIFO is full    : event is lost and an empty word is sent
    //
    // For each event in the FIFO one stores the time of arrival, and the time when the extraction is finished (therefore freeing a FIFO slots)
    
    
    
    
    cout << "--> Entering loop 3, producing the final root file..." << endl;

    
    // Write the L1 data
    // Everything is written down on a 1BX basis
    
    char buffer1[200];
    char buffer2[200];
    char buffer3[200];
    char buffer4[200];
    

    if (m_lay!=-1 && m_lad!=-1 && m_mod!=-1)
    {
        sprintf(buffer1, "l1_FE_IN_L%dR%dM%d_%d_%d.txt",m_lay,m_lad,m_mod,m_rate,nevts);
        sprintf(buffer2, "l1_FE_OUT_L%dR%dM%d_%d_%d.txt",m_lay,m_lad,m_mod,m_rate,nevts);
        sprintf(buffer4, "l1_FE_OUT_L%dR%dM%d_%d_%d_err_%d%d%d.txt",m_lay,m_lad,m_mod,m_rate,nevts,m_error_rdm,m_error_hdr,m_error_siz);
        sprintf(buffer3, "l1_CIC_OUT_L%dR%dM%d_%d_%d.txt",m_lay,m_lad,m_mod,m_rate,nevts);
    }
    else if (m_lay!=-1 && m_lad!=-1)
    {
        sprintf(buffer1, "l1_FE_IN_L%dR%d_%d_%d.txt",m_lay,m_lad,m_rate,nevts);
        sprintf(buffer2, "l1_FE_OUT_L%dR%d_%d_%d.txt",m_lay,m_lad,m_rate,nevts);
        sprintf(buffer4, "l1_FE_OUT_L%dR%d_%d_%d_err_%d%d%d.txt",m_lay,m_lad,m_rate,nevts,m_error_rdm,m_error_hdr,m_error_siz);
        sprintf(buffer3, "l1_CIC_OUT_L%dR%d_%d_%d.txt",m_lay,m_lad,m_rate,nevts);
    }
    else if (m_lay!=-1 && m_mod!=-1)
    {
        sprintf(buffer1, "l1_FE_IN_L%dM%d_%d_%d.txt",m_lay,m_mod,m_rate,nevts);
        sprintf(buffer2, "l1_FE_OUT_L%dM%d_%d_%d.txt",m_lay,m_mod,m_rate,nevts);
        sprintf(buffer4, "l1_FE_OUT_L%dM%d_%d_%d_err_%d%d%d.txt",m_lay,m_mod,m_rate,nevts,m_error_rdm,m_error_hdr,m_error_siz);
        sprintf(buffer3, "l1_CIC_OUT_L%dM%d_%d_%d.txt",m_lay,m_mod,m_rate,nevts);
    }
    else if (m_lad!=-1 && m_mod!=-1)
    {
        sprintf(buffer1, "l1_FE_IN_R%dM%d_%d_%d.txt",m_lad,m_mod,m_rate,nevts);
        sprintf(buffer2, "l1_FE_OUT_R%dM%d_%d_%d.txt",m_lad,m_mod,m_rate,nevts);
        sprintf(buffer4, "l1_FE_OUT_R%dM%d_%d_%d_err_%d%d%d.txt",m_lad,m_mod,m_rate,nevts,m_error_rdm,m_error_hdr,m_error_siz);
        sprintf(buffer3, "l1_CIC_OUT_R%dM%d_%d_%d.txt",m_lad,m_mod,m_rate,nevts);
    }
    else if (m_lay!=-1)
    {
        sprintf(buffer1, "l1_FE_IN_L%dR%d_%d_%d.txt",m_lay,m_lad,m_rate,nevts);
        sprintf(buffer2, "l1_FE_OUT_L%dR%d_%d_%d.txt",m_lay,m_lad,m_rate,nevts);
        sprintf(buffer4, "l1_FE_OUT_L%dR%d_%d_%d_err_%d%d%d.txt",m_lay,m_lad,m_rate,nevts,m_error_rdm,m_error_hdr,m_error_siz);
        sprintf(buffer3, "l1_CIC_OUT_L%dR%d_%d_%d.txt",m_lay,m_lad,m_rate,nevts);
    }
    else if (m_mod!=-1)
    {
        sprintf(buffer1, "l1_FE_IN_L%dM%d_%d_%d.txt",m_lay,m_mod,m_rate,nevts);
        sprintf(buffer2, "l1_FE_OUT_L%dM%d_%d_%d.txt",m_lay,m_mod,m_rate,nevts);
        sprintf(buffer4, "l1_FE_OUT_L%dM%d_%d_%d_err_%d%d%d.txt",m_lay,m_mod,m_rate,nevts,m_error_rdm,m_error_hdr,m_error_siz);
        sprintf(buffer3, "l1_CIC_OUT_L%dM%d_%d_%d.txt",m_lay,m_mod,m_rate,nevts);
    }
    else if (m_lad!=-1)
    {
        sprintf(buffer1, "l1_FE_IN_R%dM%d_%d_%d.txt",m_lad,m_mod,m_rate,nevts);
        sprintf(buffer2, "l1_FE_OUT_R%dM%d_%d_%d.txt",m_lad,m_mod,m_rate,nevts);
        sprintf(buffer4, "l1_FE_OUT_R%dM%d_%d_%d_err_%d%d%d.txt",m_lad,m_mod,m_rate,nevts,m_error_rdm,m_error_hdr,m_error_siz);
        sprintf(buffer3, "l1_CIC_OUT_R%dM%d_%d_%d.txt",m_lad,m_mod,m_rate,nevts);
    }
    else
    {
        sprintf(buffer1, "l1_FE_IN_ALL_%d_%d.txt",m_rate,nevts);
        sprintf(buffer2, "l1_FE_OUT_ALL_%d_%d.txt",m_rate,nevts);
        sprintf(buffer4, "l1_FE_OUT_ALL_%d_%d_err_%d%d%d.txt",m_rate,nevts,m_error_rdm,m_error_hdr,m_error_siz);
        sprintf(buffer3, "l1_CIC_OUT_ALL_%d_%d.txt",m_rate,nevts);
    }
    
    
    FE_L1_IN.open(buffer1);
    FE_L1_OUT.open(buffer2);
    FE_L1_OUT_E.open(buffer4);
    CIC_L1_OUT.open(buffer3);
    
    FE_L1_IN << "Digital input to the FE chip.\n";
    FE_L1_IN << "\n";
    FE_L1_IN << "Provides a list of all the L1 data contained in the considered modules\n";
    FE_L1_IN << "\n";
    
    
    FE_L1_OUT << "Digital L1 output of the FE chip.\n";
    FE_L1_OUT << "\n";
    FE_L1_OUT << "Format defined in:\n";
    FE_L1_OUT << "https://espace.cern.ch/Tracker-Upgrade/Electronics/CIC/Shared_Documents/Data_formats/CIC_IO_Formats_v2.pdf\n";
    FE_L1_OUT << "\n";
    FE_L1_OUT << "For each FE chip required, the list of L1 fragments is given for a period of " << nevts << " bunch crossings\n";
    FE_L1_OUT << "Errors in the stream (outside L1 data)      : 0\n";
    FE_L1_OUT << "Errors in the cluster size fields (MPA only): 0\n";
    FE_L1_OUT << "Errors in the header fields                 : 0\n";
    
    FE_L1_OUT_E<< "Digital L1 output of the FE chip (with random errors).\n";
    FE_L1_OUT_E << "\n";
    FE_L1_OUT_E << "Format defined in:\n";
    FE_L1_OUT_E << "https://espace.cern.ch/Tracker-Upgrade/Electronics/CIC/Shared_Documents/Data_formats/CIC_IO_Formats_v2.pdf\n";
    FE_L1_OUT_E << "\n";
    FE_L1_OUT_E << "For each FE chip required, the list of L1 fragments is given for a period of " << nevts << " bunch crossings\n";
    FE_L1_OUT_E << "Errors in the stream (outside L1 data)      : "<< m_error_rdm <<"\n";
    FE_L1_OUT_E << "Errors in the cluster size fields (MPA only): "<< m_error_siz <<"\n";
    FE_L1_OUT_E << "Errors in the header fields                 : "<< m_error_hdr <<"\n";

    CIC_L1_OUT << "Digital L1 output of the CIC chip.\n";
    CIC_L1_OUT << "\n";
    CIC_L1_OUT << "Format defined in:\n";
    CIC_L1_OUT << "https://espace.cern.ch/Tracker-Upgrade/Electronics/CIC/Shared_Documents/Data_formats/CIC_IO_Formats_v2.pdf\n";
    CIC_L1_OUT << "\n";
    CIC_L1_OUT << "For each CIC chip required, the list of L1 fragments is given for a period of " << nevts << " bunch crossings\n";

    
    
    // Create the sequences, and write them on the fly

    int L1_id = 0;

    std::vector<int> FIFO,FIFO_new,word;
    int last_rank;
    int last_BX;

    int extracted_bit_per_BX = 8;

    int nevts_tmp = nevts+10000;

    m_words.clear();
    m_raw_FIFO.clear();
    m_chip_FIFOs.clear();
    
    m_c_words.clear();
    m_c_raw_FIFO.clear();
    m_c_chip_FIFOs.clear();

    for (unsigned int i=0;i<m_concs.size();++i)
    {
      m_digi_list.clear();

      for (int j=0;j<m_npblock*nevts_tmp;++j)
      {
          m_digi_list.push_back(0);
      }
      m_c_words.insert(std::make_pair(m_concs.at(i),m_digi_list));
        
      m_digi_list.clear();
      m_digi_list.push_back(-1);
      m_c_raw_FIFO.insert(std::make_pair(m_concs.at(i),m_digi_list));
      m_digi_list.clear();
      m_digi_list.push_back(-1);
      m_digi_list.push_back(-1);
      m_c_chip_FIFOs.insert(std::make_pair(m_concs.at(i),m_digi_list));
    }

    for (unsigned int i=0;i<m_chips.size();++i)
    {
      m_digi_list.clear();
    
      for (int j=0;j<extracted_bit_per_BX*nevts_tmp;++j)
      {
          m_digi_list.push_back(m_error_rdm*badbit());
      }
      m_words.insert(std::make_pair(m_chips.at(i),m_digi_list));
        
      m_digi_list.clear();
      m_digi_list.push_back(-1);
      m_raw_FIFO.insert(std::make_pair(m_chips.at(i),m_digi_list));
      m_digi_list.clear();
      m_digi_list.push_back(-1);
      m_digi_list.push_back(-1);
      m_digi_list.push_back(-1);
      m_chip_FIFOs.insert(std::make_pair(m_chips.at(i),m_digi_list));
    }



    //
    // First of all one writes the L1 word, either at the FE of the concentrator level
    //
    //
    //

    for (int i=0;i<nevts;++i) 
    { 
        if (i%1000==0)
            cout << i << endl;

        if (i%3564==0) L1_id=0; // BCR, L1ID get back to 0

        l1_builder::initVars();
        
        if (raw_seq.at(i)!=-1) // Chip has received received a L1 A, we write the raw block
        {
            ++L1_id;
	
            
            FE_L1_IN << "\n";
            FE_L1_IN << "##################################################\n";
            FE_L1_IN << "#\n";
            FE_L1_IN << "# L1A signal number " << L1_id << "\n";
            FE_L1_IN << "# Received at BX ID " << i << "\n";
            FE_L1_IN << "# FE chip contents:\n";
            
            // First of all for the FE chips
            
            m_chip_raw = m_data_raw.at(2*raw_seq.at(i)); // Pick up the event in the store
            
            for ( m_iter = m_chip_raw.begin(); m_iter != m_chip_raw.end();++m_iter )
            {

                // First, look at the main info
                isPS            = false;
                m_raw_cic       = 0;
                m_raw_bx        = i;
                m_raw_FIFO_FULL = 0;
                m_raw_FIFO_SIZE = 0;
                m_raw_chip      = m_iter->first;
                m_raw_lay       = m_raw_chip/1000000;
                m_raw_lad       = (m_raw_chip-1000000*m_raw_lay)/10000;
                m_raw_mod       = (m_raw_chip-1000000*m_raw_lay-10000*m_raw_lad)/100;
                m_raw_np        = 0;
                m_raw_ns        = 0;
                
                FE_L1_IN << "\n";
                FE_L1_IN << "__________________________________________________\n";
                FE_L1_IN << "Layer/Disk  " << m_raw_lay << "\n";
                FE_L1_IN << "Ladder/Ring " << m_raw_lad << "\n";
                FE_L1_IN << "Module      " << m_raw_mod << "\n";
                FE_L1_IN << "Chip        " << m_raw_chip%100 << "\n";
                
                if (m_raw_lay<8 || (m_raw_lay>10 && m_raw_lad<9)) isPS = true;

                // Now, we get the digis contained in the event, and fill the block
                m_digi_list = m_iter->second;
	  
                // We have the info, we now write the word, making a difference bet. CIC and FE chip words

                l1_builder::fill_RAW_block(m_digi_list,isPS,L1_id);
                
                // Start the printing loop
                
                if (m_digi_list.size()-1 != 0)
                {
                    for (unsigned int k=0;k<(m_digi_list.size()-1)/3;++k) // Loop over all digis
                    {
                        (isPS)
                            ? FE_L1_IN << "digi " << std::setw(5) << m_digi_list.at(3*k+2) << " " << std::setw(3) << m_digi_list.at(3*k+3)  << "\n"
                            : FE_L1_IN << "digi " << std::setw(5) << m_digi_list.at(3*k+2) << "\n";
                    }
                }
                    
                // The sparsified word is written in the m_raw_data vector. We now store it in the FIFO, adding it's
                // BX_in and BX_out coordinates

                m_raw_size  = m_raw_data->size();
	  
                int FIFO_size = 0.;
                
                // Number of L1 bits extracted per bunch crossing is 8 for all the chips (one diff line at 320 MHz)
                // could be 16 for the CIC in the 10G-GBT scenario
                

                int delay                = 0; // Delay is the minimal period, in BX, between L1A reception and data emission by the chip
                
                (isPS)
                ? delay = m_MPA_L1_delay+m_raw_np+m_raw_ns // From D.Ceresa
                : delay = m_CBC_L1_delay;
                
                // Find the chip FIFO footprint
                // defined as follows:

                // < CHIP_ID, <-1,1,BXin(1),BXout(1),2,BXin(2),BXout(2),... > >

                m_iter2 = m_raw_FIFO.find(m_raw_chip); // Find the chip

                FIFO_new.clear();
                FIFO      = m_iter2->second; // Recover the FIFO info for the chip
                FIFO_size = (FIFO.size()-1)/3; // How many L1 events are in the FIFO (depth)?
                FIFO_new.push_back(-1);

		// We will fill FIFO_new with the new FIFO content

                if (FIFO_size==0) // Empty FIFO, initialize it!
                {
                    // BXin is BX for for the event starts to be extractible from the chip
                    // BXout = BXin + numb of BXs to extract the event (if FIFO is empty it's simple)
                    
                    FIFO_new.push_back(L1_id);
                    FIFO_new.push_back(i+delay);  // BXin
                    FIFO_new.push_back(i+delay+m_raw_size/(extracted_bit_per_BX)+1); // BXout
                    
                    m_raw_FIFO_FULL = (FIFO_new.size()-1)/3; // New FIFO size
                    m_raw_FIFO_SIZE = i+delay;  
                    last_BX         = i+delay+m_raw_size/(extracted_bit_per_BX)+1;
                    
                    m_iter3 = m_words.find(m_raw_chip);
                    word    = m_iter3->second; // Full data stream for the chip
                    
                    for (int k=extracted_bit_per_BX*(i+delay);k<extracted_bit_per_BX*(i+delay)+m_raw_size;++k)
                    {
                        if (k>=extracted_bit_per_BX*nevts_tmp) continue;

			if (word.at(k)==1) std::cout << "E/STRANGE " << m_raw_chip << " / " << k/extracted_bit_per_BX << std::endl;
                        word.at(k) = m_raw_data->at(k-extracted_bit_per_BX*(i+delay));
                    }
                    
                    m_words.erase(m_iter3->first);
                    m_words.insert(std::make_pair(m_raw_chip,word));
                        
                }
                else  // FIFO not empty, one need to increment
                {
                    last_BX=i+delay;
                    
		    //		    std::cout << m_raw_chip << " / " << (FIFO.size()-1)/3 << std::endl;

                    for (unsigned int j=0;j<(FIFO.size()-1)/3;++j)
                    {
                        if (FIFO.at(3*j+3)<i) continue; // Time to go out for this event (BXo(j)<i)

                        FIFO_new.push_back(FIFO.at(3*j+1)); // Otherwise event stays there
                        FIFO_new.push_back(FIFO.at(3*j+2)); 
                        FIFO_new.push_back(FIFO.at(3*j+3));
                        
                        // The last BXo fixes the time where we can start to extract the next event
                        last_BX=std::max(FIFO.at(3*j+3),i+delay);

			//			std::cout << m_raw_chip << " / " << i << " / " << j << " / " << FIFO.at(3*j+1)  << " / " << FIFO.at(3*j+2)  << " / " << FIFO.at(3*j+3)  << " / " << last_BX  << std::endl;
                    }

                    // FIFO has been updated, check if it is FULL or not
                    
                    if (static_cast<int>((FIFO_new.size()-1)/3)>=1000000*m_MPA_FIFO_depth)
                    {
		      cout << "SIZE ERROR -> " << static_cast<int>((FIFO_new.size()-1)/3) << " > " << m_MPA_FIFO_depth << endl;
                    }
                    
                    // Put the latest event
                    
		    FIFO_new.push_back(L1_id);
                    FIFO_new.push_back(i+delay);
                    FIFO_new.push_back(last_BX+m_raw_size/(extracted_bit_per_BX)+1);
                    
		    //		    std::cout << m_raw_chip << " / " << " / " << L1_id  << " / " << i+delay  << " / " << last_BX+m_raw_size/(extracted_bit_per_BX)+1  << std::endl;

                    // Here we compute the size of the current FIFO

                    m_raw_FIFO_SIZE = last_BX;
                    m_raw_FIFO_FULL = (FIFO_new.size()-1)/3;

                    m_iter3 = m_words.find(m_raw_chip);
                    word=m_iter3->second;
                    
                   // cout << last_BX << endl;
                    
                    for (int k=extracted_bit_per_BX*last_BX;k<extracted_bit_per_BX*last_BX+m_raw_size;++k)
                    {
                        if (k>=extracted_bit_per_BX*nevts_tmp) continue;
			if (word.at(k)==1) std::cout << "F/STRANGE " << m_raw_chip << " / " << k/extracted_bit_per_BX << std::endl;
                        word.at(k) = m_raw_data->at(k-extracted_bit_per_BX*last_BX);
                    }
                    
                    m_words.erase(m_iter3->first);
                    m_words.insert(std::make_pair(m_raw_chip,word));
                    
                    last_BX=last_BX+m_raw_size/(extracted_bit_per_BX)+1;
                }

		//		std::cout << m_raw_chip << " / " << (FIFO_new.size()-1)/3 << std::endl;

                m_raw_FIFO.erase(m_iter2->first);
                m_raw_FIFO.insert(std::make_pair(m_raw_chip,FIFO_new)); // Update the FIFO for this chip

                // Then update the chip FIFO list

                m_iter2 = m_chip_FIFOs.find(m_raw_chip);
                FIFO=m_iter2->second;
                
		FIFO.push_back(L1_id);
                FIFO.push_back(i+delay);
                FIFO.push_back(last_BX-i-delay);
                m_chip_FIFOs.erase(m_iter2->first);
                m_chip_FIFOs.insert(std::make_pair(m_raw_chip,FIFO)); // Update the FIFO list
                
                m_raw_tree->Fill();
                
            } // End of the FE chips loop

            // Then for the corresponding CIC chips
        
            m_chip_raw = m_data_raw.at(2*raw_seq.at(i)+1); // Pick up the event in the store
            
            for ( m_iter = m_chip_raw.begin(); m_iter != m_chip_raw.end();++m_iter )
            {
                // First, look at the main info
                isPS            = false;
                m_raw_cic       = 1;
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
                
                l1_builder::fill_CONC_RAW_block(m_digi_list,isPS,L1_id,unsp);
            
                // The sparsified word is written in the m_raw_data vector. We now store it in the FIFO, adding it's
                // BX_in and BX_out coordinates
                // For the CIC, the BX_in corresponds to the last BX_out of the corresponding FE chips
                //
                // We therefore first have to retrieve them
                
                int bx_in = i;
                
                //std::cout << "L1A was received at BX " << i << " in CIC chip " << m_raw_chip << std::endl;
                
                for (int ii=0;ii<8;++ii)
                {
                    m_iter2 = m_chip_FIFOs.find(m_raw_chip+ii);
                    FIFO=m_iter2->second;

                    for (unsigned int ij=0;ij<FIFO.size()/3;++ij)
                    {
                        if (FIFO.at(3*ij)==L1_id)
                        {
			  //std::cout << "FE chip " << m_raw_chip+ii << " ends up extraction at BX " << FIFO.at(3*ij+2)+FIFO.at(3*ij+1) << std::endl;
                            if (FIFO.at(3*ij+2)+FIFO.at(3*ij+1)>=bx_in) bx_in = FIFO.at(3*ij+2)+FIFO.at(3*ij+1);
                            break;
                        }
                    }
                }
                
                
		//                std::cout << "L1 event will start to be extractible at BX " << bx_in << std::endl;
                
                m_raw_size  = m_raw_data->size();
            
                int FIFO_size = 0.;
            
                // Number of L1 bits extracted per bunch crossing is 8 for all the chips (one diff line at 320 MHz)
                // could be 16 for the CIC in the 10G-GBT scenario
            
                int extracted_bit_per_BX = m_npblock;
            
                // Find the chip FIFO footprint
                // defined as follows:
            
                // < CHIP_ID, <-1,BXin(1),BXout(1),BXin(2),BXout(2),... > >
            
                m_iter2 = m_c_raw_FIFO.find(m_raw_chip); // Find the chip
            
                FIFO_new.clear();
                FIFO=m_iter2->second;
                FIFO_size=(FIFO.size()-1)/2; // How many L1 events are in the FIFO (depth)?
                FIFO_new.push_back(-1);
            
                if (FIFO_size==0) // Empty FIFO, initialize it!
                {
                    FIFO_new.push_back(bx_in);                                     // BXin is the current BX, the L1 event enters the CIC/MPA/CBC chip
                    FIFO_new.push_back(bx_in+m_raw_size/(extracted_bit_per_BX)+1); // BXout = BXin + numb of BXs to extract the event (if FIFO is empty it's simple)
                    m_raw_FIFO_FULL = (FIFO_new.size()-1)/2;                       // Number of events in the FIFO (1)
                    FIFO_size=m_raw_size;
                    m_raw_FIFO_SIZE=bx_in;
                    last_BX=bx_in+m_raw_size/(extracted_bit_per_BX)+1;
                    
                    m_iter3 = m_c_words.find(m_raw_chip);
                    word=m_iter3->second;
                    
                    for (int k=m_npblock*bx_in;k<m_npblock*bx_in+m_raw_size;++k)
                    {
                        if (k>=m_npblock*nevts_tmp) continue;
                        word.at(k) = m_raw_data->at(k-m_npblock*bx_in);
                    }
                    
                    m_c_words.erase(m_iter3->first);
                    m_c_words.insert(std::make_pair(m_raw_chip,word));
                    
                }
                else  // FIFO not empty, one need to increment
                {
                    last_BX=bx_in;
                
                    for (unsigned int j=0;j<(FIFO.size()-1)/2;++j)
                    {
                        if (FIFO.at(2*j+2)<i) continue; // Time to go out for this event (BXo(j)<i)
                    
                        FIFO_new.push_back(FIFO.at(2*j+1)); // Otherwise event stays there
                        FIFO_new.push_back(FIFO.at(2*j+2));
                    
                        last_BX=std::max(FIFO.at(2*j+2),bx_in); // The last BXo fixes the time where we can start to extract the next event
                    }
                
                    // Here we compute the size of the current FIFO
                
                    last_rank=0;
                    if (FIFO_new.size()>1) last_rank = (last_BX-i)*(extracted_bit_per_BX);
                
                    FIFO_size=m_raw_size+last_rank;
                    m_raw_FIFO_SIZE=last_BX;
                
                    FIFO_new.push_back(bx_in);                                          // BXin is the current BX, the L1 event enters the CIC
                    FIFO_new.push_back(last_BX+m_raw_size/(extracted_bit_per_BX)+1); // BXout tells us when the event really goes out
                    m_raw_FIFO_FULL = (FIFO_new.size()-1)/2;
                
                    
                    m_iter3 = m_c_words.find(m_raw_chip);
                    word=m_iter3->second;
                    
                    for (int k=m_npblock*last_BX;k<m_npblock*last_BX+m_raw_size;++k)
                    {
                        if (k>=m_npblock*nevts_tmp) continue;
                        word.at(k) = m_raw_data->at(k-m_npblock*last_BX);
                    }
                    
                    m_c_words.erase(m_iter3->first);
                    m_c_words.insert(std::make_pair(m_raw_chip,word));

                    
                    last_BX=last_BX+m_raw_size/(extracted_bit_per_BX)+1;
                }
            
                m_c_raw_FIFO.erase(m_iter2->first);
                m_c_raw_FIFO.insert(std::make_pair(m_raw_chip,FIFO_new)); // Update the FIFO for this chip
            
                // Then update the chip FIFO list
            
                m_iter2 = m_c_chip_FIFOs.find(m_raw_chip);
                FIFO=m_iter2->second;
                FIFO.push_back(i);
                FIFO.push_back(bx_in);
                FIFO.push_back(last_BX-bx_in);
                m_c_chip_FIFOs.erase(m_iter2->first);
                m_c_chip_FIFOs.insert(std::make_pair(m_raw_chip,FIFO)); // Update the FIFO list
                
                m_raw_tree->Fill();
                
            } // End of the CIC chip loop
        }
    } // End of the loop on raw events

    
    // Print the L1 data fragments
    
    int n_err;

    for (unsigned int i=0;i<m_chips.size();++i)
    {
        L1_id = 0;
        m_iter = m_words.find(m_chips.at(i));
        FE_L1_OUT << "     CHIP ID |    BX ID | L1 ID | L1 fragment " << "\n";
        FE_L1_OUT_E << "     CHIP ID |    BX ID | L1 ID | L1 fragment " << "\n";
        
        for (unsigned int j=0;j<m_iter->second.size();++j)
        {
            if (j%8==0)
            {
              n_err=0;
	      if ((j/8)%3564==0) L1_id=0; // BCR, L1ID get back to 0
	      FE_L1_OUT << std::setw(12 - (j!=0)) << m_iter->first << " | " << std::setw(8) << j/8 << " | ";
	      FE_L1_OUT_E << std::setw(12 - (j!=0)) << m_iter->first << " | " << std::setw(8) << j/8 << " | ";
	      
	      if (j/8<nevts)
	      {
		if (raw_seq.at(j/8)!=-1) // Chip has received received a L1 A, we write the raw block
		{
		  ++L1_id;
		  FE_L1_OUT << std::setw(5) << L1_id << " | ";
		  FE_L1_OUT_E << std::setw(5) << L1_id << " | ";
		}
		else
		{
		  FE_L1_OUT << "      | ";
		  FE_L1_OUT_E << "      | ";
		}
	      }
	      else
	      {
		FE_L1_OUT << "      | ";
		FE_L1_OUT_E << "      | ";
	      }
	    }
            
            FE_L1_OUT << m_iter->second.at(j)%2;

	    if (m_iter->second.at(j)>1)
	    {
	      ++n_err;
	      FE_L1_OUT_E << 1-m_iter->second.at(j)%2;
	    }
	    else
	    {
	      FE_L1_OUT_E << m_iter->second.at(j)%2;
	    }

            if (j%8==7) FE_L1_OUT << " | |\n ";
            if (j%8==7 && n_err==0) FE_L1_OUT_E << " | |\n ";
            if (j%8==7 && n_err>0) FE_L1_OUT_E << " |E|\n ";
        }
        
        FE_L1_OUT << "\n";
	FE_L1_OUT_E << "\n";
    }
    
    for (unsigned int i=0;i<m_concs.size();++i)
    {
        L1_id = 0;
        m_iter = m_c_words.find(m_concs.at(i));
        CIC_L1_OUT << "      CIC ID |    BX ID | L1 ID | L1 fragment " << "\n";
        
        for (unsigned int j=0;j<m_iter->second.size();++j)
        {
            if (j%m_npblock==0)
            {
                if ((j/m_npblock)%3564==0) L1_id=0; // BCR, L1ID get back to 0
                CIC_L1_OUT << std::setw(12 - (j!=0)) << m_iter->first << " | " << std::setw(8) << j/m_npblock << " | ";

		if (j/m_npblock<nevts)
		{
		  if (raw_seq.at(j/m_npblock)!=-1) // Chip has received received a L1 A, we write the raw block
		  {
                    ++L1_id;
                    CIC_L1_OUT << std::setw(5) << L1_id << " | ";
		  }
		  else
		  {
                    CIC_L1_OUT << "      | ";
		  }
		}
		else
		{
		  CIC_L1_OUT << "      | ";
		}

            }
            CIC_L1_OUT << m_iter->second.at(j);
            if (int(j)%m_npblock==m_npblock-1) CIC_L1_OUT << " \n ";
        }
        
        CIC_L1_OUT << "\n";
    }
    
    FE_L1_IN.close();
    FE_L1_OUT.close();
    FE_L1_OUT_E.close();
    CIC_L1_OUT.close();
    
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
  
  float x,z;
	
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
  

 
    m_outfile->Write();
  
    delete PIX;
    delete m_outfile;
}




/////////////////////////////////////////////////////////////
//
// Basic methods, initializations,...
//
/////////////////////////////////////////////////////////////

void l1_builder::initVars()
{
  m_raw_cic=0;
  m_raw_bx=0;
  m_raw_chip=0;
  m_raw_size=0;
  m_raw_mbits=0;
  m_raw_np=0;
  m_raw_ns=0;

  m_raw_lay=0;
  m_raw_lad=0;
  m_raw_mod=0;
  m_raw_FIFO_FULL=0;

  m_raw_data->clear();

  m_raw_chip=0;
  m_raw_chip_bx->clear();
  m_raw_chip_fifo->clear();
  m_raw_conc_bx->clear();
  m_raw_conc_fifo->clear();
  m_raw_chip_slope=0.;
  m_raw_chip_slope_err=0.;
  m_raw_chip_int=0.;
  m_raw_chip_int_err=0.;
}


void l1_builder::initTuple(std::string inRAW,std::string out)
{
    m_outfile    = new TFile(out.c_str(),"recreate");

    m_raw_tree    = new TTree("Raw_FE","Raw data words after FE");
    m_raw_summary = new TTree("Raw_SUM","Raw data summary info");

    m_raw_tree->Branch("RAW_CIC",        &m_raw_cic,     "RAW_CIC/I");
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

bool l1_builder::convert(std::string sectorfilename)
{
  std::vector<int> module;

  std::cout << "Convert in" << std::endl;

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
        

        if (mlay!=m_lay && m_lay!=-1) continue;
        if (mlad!=m_lad && m_lad!=-1) continue;
        if (mmod!=m_mod && m_mod!=-1) continue;

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
// https://espace.cern.ch/Tracker-Upgrade/Electronics/CIC/Shared%20Documents/Data%20formats/CIC_IO_Formats_v3.pdf
//


//
// 1: L1 raw block at the MPA/CBC level
//

void l1_builder::fill_RAW_block(std::vector<int> digis,bool spars,int BXid)
{
    m_raw_data->clear();

    // First write the headers
  
    (spars)
        ? l1_builder::fill_RAW_header_MPA(BXid)
        : l1_builder::fill_RAW_header_CBC(BXid);

    // Then go for the data

    int ndata=(digis.size()-1)/3; // Default first item is -1
    int hsize=m_raw_data->size();

    if (!spars) // CBC (unsparsified)
    {
        for (int j=0;j<254;++j)   m_raw_data->push_back(0);
    
        for (int j=0;j<ndata;++j)
        {
//            cout << "digiz " << std::setw(5) << digis.at(3*j+1) << " " << std::setw(5) << digis.at(3*j+2) << " " << std::setw(3) << digis.at(3*j+3)  << endl;
            m_raw_data->at(digis.at(3*j+2)+hsize) = 1;
        }
    }
    else // MPA
    {
        int data_MPA[128][17];
        int row,col;

        // The sparsification is done at the MPA level
        // Strips 0 to 120 belong to the pixel layer
        // 120 to 239 belong to the strip layer

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

            for (int i=0;i<128;++i)
            {
                if (data_MPA[i][j] == 1)
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

        if (np>31 || ns>28)
        {
            cout << "OVERFLOW ERROR IN MPA" << endl;
            m_raw_data->at(19) = 1;
        }

        np = std::min(31,np); // One cannot pass more than
        ns = std::min(28,ns); // 31 pixel clusters and 28 strips clusters out of the MPA

        std::bitset<5> N_S = ns;
        std::bitset<5> N_P = np;

	for (int j=0;j<5;++j) m_raw_data->push_back(m_error_siz*badbit()+N_S[4-j]);
        for (int j=0;j<5;++j) m_raw_data->push_back(m_error_siz*badbit()+N_P[4-j]);
	m_raw_data->push_back(0);

        for (int j=0;j<ns;++j)
        {
            std::bitset<7> row =  clus_s.at(2*j);
            std::bitset<3> wdt =  clus_s.at(2*j+1);

            for (int k=0;k<7;++k) m_raw_data->push_back(row[6-k]);
            for (int k=0;k<3;++k) m_raw_data->push_back(wdt[2-k]);

	    //////////////
	    // SV 17/12/2015
	    // Adding the MIP flag (for the moment at 0 by default)
	    // Neew to update this when the new digitizer will be in CMSSW
	    m_raw_data->push_back(0);
	    //
	    //////////////
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


    //std::cout << m_raw_chip << " / size of the L1 raw word / " << spars << " / " << m_raw_data->size() << std::endl;
    
    //for (int j=0;j<m_raw_data->size();++j) std::cout << m_raw_data->at(j);
    
    //std::cout << std::endl;
    
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

void l1_builder::fill_CONC_RAW_block(std::vector<int> digis,bool MPA,int BXid,bool unsparsified)
{
//    cout << "Into fill_CONC_RAW_block" << endl;
    
    m_raw_data->clear();

    // First write the header

    l1_builder::fill_CONC_RAW_header(BXid);

    if (unsparsified && !MPA) m_raw_data->clear();

    // Then go for the data

    int np,ns;

    int ndata=(digis.size()-1)/3; // Default first item is -1
    int hsize=m_raw_data->size();
 
    if (!MPA) // CBC case
    {
        int data_CBC[8][254];
        int row,chip;

        for (int j=0;j<254;++j)
        {
            for (int i=0;i<8;++i) data_CBC[i][j] = 0;
        }

        for (int j=0;j<ndata;++j)
        {
            chip = digis.at(3*j+1);
            row  = digis.at(3*j+2);
     
            //std::cout << row << "," << chip << " / ";
 
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

            //cout << j << " / " << clus_s.size()/3 << endl;

            for (int i=0;i<254;++i)
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
    
        //std::cout << "This CBC conc chip contains " << clus_s.size()/3 << " strip clusters" << std::endl;

        // Now encode them
        ns = clus_s.size()/3;
        m_raw_ns        = ns;

	if (unsparsified)
	{
	  for (int j=0;j<15;++j) m_raw_data->push_back(1);
	  m_raw_data->push_back(0);

	  for (int j=0;j<8;++j)
	  {
	    l1_builder::fill_RAW_header_CBC(BXid);
            for (int i=0;i<254;++i) m_raw_data->push_back(data_CBC[j][i]);
	  }
	}
	else
	{
	  if (ns>127)
	  {
            cout << "OVERFLOW ERROR IN CIC 2S" << endl;
            m_raw_data->at(24) = 1;
	  }
	  
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
    }
    else // MPA case
    {
      //  cout << "Start the MPA treatment with " << ndata << " digis " << endl;
        
        int data_MPA[8][120][17];
        int row,col,chip;

        // The sparsification is done at the MPA level
        // Strips 0 to 119 belong to the pixel layer
        // 120 to 239 belong to the strip layer
 
        for (int k=0;k<8;++k)
        {
            for (int j=0;j<120;++j)
            {
                for (int i=0;i<17;++i) data_MPA[k][j][i] = 0;
            }
        }

        for (int j=0;j<ndata;++j)
        {
            chip = digis.at(3*j+1);
            row  = digis.at(3*j+2);
            col  = digis.at(3*j+3);
     
            //std::cout << chip << "," << row << "," << col << " / ";

            (row<120)
            ? data_MPA[chip][row][col]    = 1
            : data_MPA[chip][row%120][16] = 1;
        }
        
//        if (ndata!=0) std::cout << std::endl;

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

          //cout << j << " / " << clus_p.size()/4 << " / " << clus_s.size()/3 << endl;

          for (int i=0;i<120;++i)
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
//                << clus_s.size()/3 << " strip clusters" << std::endl;

    // Now encode them 
    np = clus_p.size()/4 ;
    ns = clus_s.size()/3;
        
    if (np>127 || ns>127)
    {
        cout << "OVERFLOW ERROR IN CIC/PS" << endl;
        m_raw_data->at(24) = 1;
    }
        
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

      //////////////
      // SV 17/12/2015
      // Adding the MIP flag (for the moment at 0 by default)
      // Neew to update this when the new digitizer will be in CMSSW
      m_raw_data->push_back(0);
      //
      //////////////
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
//    std::cout << "This chip contains " << np << " pixels clusters and "
//    	      << ns << " strip clusters" << std::endl;

//    std::cout << "Size of the L1 raw word / " << spars << " / " << m_raw_data->size() << std::endl;

    for (unsigned int j=0;j<m_raw_data->size();++j)
    {
      std::cout << m_raw_data->at(j);
    }

    std::cout << std::endl;
  }

  //  cout << "Out of fill_CONC_RAW_block" << endl;
}



void l1_builder::fill_RAW_header_CBC(int L1id)
{
  // Format of the CBC L1 word header
  //
  // HHEEPPPPPPPPPLLLLLLLLL : HHEE (header + error) LL..LL (L1 ID bet 0 and 512)
  //                          PP..PP CBC pipeline address 

    for (int j=0;j<2;++j) m_raw_data->push_back(m_error_hdr*badbit()+1); // HH
    for (int j=0;j<2;++j) m_raw_data->push_back(m_error_hdr*badbit()+0); // EE

    if (L1id>512)
    {
        std::cout << "Too many L1ids, problem!!!" << std::endl;
    }

    std::bitset<9> L1_ID = L1id;

    for (int j=0;j<9;++j) m_raw_data->push_back(m_error_hdr*badbit()+L1_ID[8-j]); // PP..PP
    for (int j=0;j<9;++j) m_raw_data->push_back(m_error_hdr*badbit()+L1_ID[8-j]); // CC..CC
}


void l1_builder::fill_RAW_header_MPA(int L1id)
{
  // Format of the MPA L1 word header
  //
  // 1111111111111111110EECCCCCCCCC0 : EE (error) CC..CC (L1 ID bet 0 and 512)
  //
  
  // if we add random errors we add 2 to the bits 

  for (int j=0;j<18;++j) m_raw_data->push_back(m_error_hdr*badbit()+1);
  m_raw_data->push_back(m_error_hdr*badbit());
  m_raw_data->push_back(m_error_hdr*badbit());
  m_raw_data->push_back(m_error_hdr*badbit());
  
  if (L1id>512)
  {
    std::cout << "Too many L1ids, problem!!!" << std::endl;
 }

  std::bitset<9> L1_ID = L1id;
    
  for (int j=0;j<9;++j) m_raw_data->push_back(m_error_hdr*badbit()+L1_ID[8-j]); // CC..CC

  m_raw_data->push_back(m_error_hdr*badbit());
}

void l1_builder::fill_CONC_RAW_header(int L1id)
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
}

int l1_builder::badbit()
{
  return 2*(rand()%(static_cast<int>(10/m_error_prop))<10);
}
