
#include "../interface/PixelExtractor.h"


PixelExtractor::PixelExtractor(edm::EDGetTokenT< edm::DetSetVector< Phase2TrackerDigi> > pixToken, edm::EDGetTokenT< edm::DetSetVector< PixelDigiSimLink> > pixslToken, edm::EDGetTokenT< std::vector<PileupSummaryInfo> > puToken, bool doTree, bool doMatch)
{
  m_OK = false;
  m_pixToken   = pixToken;
  m_pixslToken = pixslToken;
  m_puToken    = puToken;
  

  m_matching = doMatch;

  // Tree definition

  m_pixclus_x        = new std::vector<float>;    
  m_pixclus_y        = new std::vector<float>;  
  m_pixclus_z        = new std::vector<float>;    
  m_pixclus_e        = new std::vector<float>;    
  m_pixclus_row      = new std::vector<int>;      
  m_pixclus_column   = new std::vector<int>;      
  m_pixclus_simhit   = new std::vector<int>;      
  m_pixclus_simhitID = new std::vector< std::vector<int> >;  
  m_pixclus_evtID    = new std::vector< std::vector<int> >;  
  m_pixclus_bcID     = new std::vector< std::vector<int> >;  

  m_pixclus_layer    = new std::vector<int>;    
  m_pixclus_module   = new std::vector<int>; 
  m_pixclus_ladder   = new std::vector<int>;    
  m_pixclus_nrow     = new std::vector<int>;    
  m_pixclus_ncolumn  = new std::vector<int>;  
  m_pixclus_bot      = new std::vector<int>;  
  m_pixclus_type     = new std::vector<int>;  
  m_pixclus_pitchx   = new std::vector<float>;  
  m_pixclus_pitchy   = new std::vector<float>; 

  PixelExtractor::reset();

  if (doTree)
  {
    m_OK = true;

    m_tree      = new TTree("Pixels","RECO Pixel info") ;

    // Branches definition

    m_tree->Branch("PIX_n",         &m_pclus,    "PIX_n/I");
    m_tree->Branch("PIX_nPU",       &m_nPU,      "PIX_nPU/I");

    m_tree->Branch("PIX_x",         &m_pixclus_x);
    m_tree->Branch("PIX_y",         &m_pixclus_y);
    m_tree->Branch("PIX_z",         &m_pixclus_z);
    m_tree->Branch("PIX_charge",    &m_pixclus_e);
    m_tree->Branch("PIX_row",       &m_pixclus_row);
    m_tree->Branch("PIX_column",    &m_pixclus_column);
    m_tree->Branch("PIX_simhit",    &m_pixclus_simhit);
    m_tree->Branch("PIX_simhitID",  &m_pixclus_simhitID);
    m_tree->Branch("PIX_evtID",     &m_pixclus_evtID);
    m_tree->Branch("PIX_bcID",      &m_pixclus_bcID);
    m_tree->Branch("PIX_layer",     &m_pixclus_layer);
    m_tree->Branch("PIX_module",    &m_pixclus_module);
    m_tree->Branch("PIX_ladder",    &m_pixclus_ladder);
    m_tree->Branch("PIX_nrow",      &m_pixclus_nrow);
    m_tree->Branch("PIX_ncolumn",   &m_pixclus_ncolumn);
    m_tree->Branch("PIX_bottom",    &m_pixclus_bot);
    m_tree->Branch("PIX_type",      &m_pixclus_type);
    m_tree->Branch("PIX_pitchx",    &m_pixclus_pitchx);
    m_tree->Branch("PIX_pitchy",    &m_pixclus_pitchy);
  }
}

PixelExtractor::PixelExtractor(TFile *a_file)
{
  std::cout << "PixelExtractor object is retrieved" << std::endl;

  // Tree definition
  m_OK = false;

  m_pixclus_x        = new std::vector<float>;    
  m_pixclus_y        = new std::vector<float>;  
  m_pixclus_z        = new std::vector<float>;    
  m_pixclus_e        = new std::vector<float>;    
  m_pixclus_row      = new std::vector<int>;      
  m_pixclus_column   = new std::vector<int>;      
  m_pixclus_simhit   = new std::vector<int>;      
  m_pixclus_simhitID = new std::vector< std::vector<int> >;  
  m_pixclus_evtID    = new std::vector< std::vector<int> >; 
  m_pixclus_bcID     = new std::vector< std::vector<int> >;  
  m_pixclus_layer    = new std::vector<int>;    
  m_pixclus_module   = new std::vector<int>; 
  m_pixclus_ladder   = new std::vector<int>;    
  m_pixclus_nrow     = new std::vector<int>;    
  m_pixclus_ncolumn  = new std::vector<int>;  
  m_pixclus_bot      = new std::vector<int>;
  m_pixclus_type     = new std::vector<int>;    
  m_pixclus_pitchx   = new std::vector<float>;  
  m_pixclus_pitchy   = new std::vector<float>; 

  PixelExtractor::reset();

  m_tree = dynamic_cast<TTree*>(a_file->Get("Pixels"));

  if (!m_tree)
  {
    std::cout << "This tree (Pixels) doesn't exist!!!" << std::endl;
    return;
  }

  m_OK = true;
  m_matching = false;

  m_n_events = m_tree->GetEntries();

  std::cout << "This file contains " << m_n_events << " events..." << std::endl;

  m_tree->SetBranchAddress("PIX_n",         &m_pclus);
  m_tree->SetBranchAddress("PIX_nPU",       &m_nPU);

  m_tree->SetBranchAddress("PIX_x",         &m_pixclus_x);
  m_tree->SetBranchAddress("PIX_y",         &m_pixclus_y);
  m_tree->SetBranchAddress("PIX_z",         &m_pixclus_z);
  m_tree->SetBranchAddress("PIX_charge",    &m_pixclus_e);
  m_tree->SetBranchAddress("PIX_row",       &m_pixclus_row);
  m_tree->SetBranchAddress("PIX_column",    &m_pixclus_column);
  m_tree->SetBranchAddress("PIX_simhit",    &m_pixclus_simhit);
  m_tree->SetBranchAddress("PIX_simhitID",  &m_pixclus_simhitID);
  m_tree->SetBranchAddress("PIX_evtID",     &m_pixclus_evtID);
  m_tree->SetBranchAddress("PIX_bcID",      &m_pixclus_bcID);
  m_tree->SetBranchAddress("PIX_layer",     &m_pixclus_layer);
  m_tree->SetBranchAddress("PIX_module",    &m_pixclus_module);
  m_tree->SetBranchAddress("PIX_ladder",    &m_pixclus_ladder);
  m_tree->SetBranchAddress("PIX_nrow",      &m_pixclus_nrow);
  m_tree->SetBranchAddress("PIX_ncolumn",   &m_pixclus_ncolumn);
  m_tree->SetBranchAddress("PIX_bottom",    &m_pixclus_bot);
  m_tree->SetBranchAddress("PIX_type",      &m_pixclus_type);
}



PixelExtractor::~PixelExtractor()
{}

void PixelExtractor::init(const edm::EventSetup *setup, bool isFlat)
{
  setup->get<TrackerDigiGeometryRecord>().get(theTrackerGeometry);
  setup->get<TrackerTopologyRcd>().get(theTrackerTopology);

  m_tilted = (!isFlat);


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
}

//
// Method filling the main particle tree
//

void PixelExtractor::writeInfo(const edm::Event *event) 
{
  PixelExtractor::reset();

  event->getByToken(m_pixToken, pDigiColl);


  const TrackerTopology* const tTopo = theTrackerTopology.product();
  const TrackerGeometry* const theTrackerGeom = theTrackerGeometry.product();

  if (m_matching)
  {
    event->getByToken(m_puToken,pileupinfos);

    std::vector<PileupSummaryInfo>::const_iterator PVI;
    
    for(PVI = pileupinfos->begin(); PVI != pileupinfos->end(); ++PVI) 
    {
      if (PVI->getBunchCrossing()==0) m_nPU = PVI->getPU_NumInteractions();      
    }

    event->getByToken(m_pixslToken, pDigiLinkColl);
  }
  

  bool barrel;
  bool endcap;

  int disk;  

  int cols;    
  int rows;    
  float pitchX;
  float pitchY;

  LocalPoint  clustlp;
  GlobalPoint pos;    

  for (auto gd=theTrackerGeom->dets().begin(); gd != theTrackerGeom->dets().end(); gd++) 
  {
    DetId detid = (*gd)->geographicalId();
    if(detid.subdetId()!=StripSubdetector::TOB && detid.subdetId()!=StripSubdetector::TID ) continue; // only run on OT

    if (pDigiColl->find( detid ) == pDigiColl->end() ) continue;

    barrel = (detid.subdetId()==StripSubdetector::TOB);
    endcap = (detid.subdetId()==StripSubdetector::TID);

    if (!barrel && !endcap) continue;

    /// Get the DetSets of the Digis
    const edm::DetSet< Phase2TrackerDigi > digis = (*pDigiColl)[ detid ];
    const GeomDetUnit* det0 = theTrackerGeom->idToDetUnit( detid );
    const PixelGeomDetUnit* theGeomDet = dynamic_cast< const PixelGeomDetUnit* >( det0 );
    const PixelTopology* topol = dynamic_cast< const PixelTopology* >( &(theGeomDet->specificTopology()) );
    
    cols     = topol->ncolumns();
    rows     = topol->nrows();
    pitchX   = topol->pitch().first;
    pitchY   = topol->pitch().second;
    
    if (m_matching)
    {
      if (pDigiLinkColl->find( detid ) == pDigiLinkColl->end() ) continue;
      pDigiLinks = (*pDigiLinkColl)[ detid ];	        
    }
   

    for ( auto iter = digis.begin();iter != digis.end();++iter ) 
    {

      clustlp = topol->localPosition( MeasurementPoint(float((*iter).row())+0.5,float((*iter).column())+0.5));
      pos     =  theGeomDet->surface().toGlobal(clustlp);

      the_ids.clear();
      the_eids.clear();
      the_bids.clear();

      m_pixclus_x->push_back(pos.x());
      m_pixclus_y->push_back(pos.y());
      m_pixclus_z->push_back(pos.z());
      m_pixclus_e->push_back(1);
      m_pixclus_row->push_back((*iter).row()); 
      m_pixclus_column->push_back((*iter).column());
      m_pixclus_bot->push_back(static_cast<int>(tTopo->isLower(detid)));
      m_pixclus_type->push_back(static_cast<int>(tTopo->tobSide(detid))); // Tilt-/Tilt+/Flat <-> 1/2/3

      if (m_matching)
      {
	if (pDigiLinks.data.size() != 0)
	{
	  // Loop over DigisSimLink in this det unit

	  for(auto it = pDigiLinks.data.begin();  it != pDigiLinks.data.end(); it++) 
	  {         
	    if (static_cast<int>(it->channel())!=static_cast<int>((*iter).channel()))
	      continue;

	    //	    std::cout << it->SimTrackId() << " ##### " << it->eventId().event() << std::endl;

	    the_ids.push_back(it->SimTrackId()); 
	    the_eids.push_back(it->eventId().event()); 
	    the_bids.push_back(it->eventId().bunchCrossing()); 
	  }
	}
      }
      

      //      std::cout << the_ids.size() << " #///# " << the_eids.size() << std::endl;

      m_pixclus_simhit->push_back(the_ids.size());
      m_pixclus_simhitID->push_back(the_ids);
      m_pixclus_evtID->push_back(the_eids);
      m_pixclus_bcID->push_back(the_bids);

      if (barrel)
      {
	m_pixclus_layer->push_back(static_cast<int>(tTopo->layer(detid))+4); 

	if (tTopo->tobSide(detid)==3)
	{ 
	  m_pixclus_ladder->push_back(static_cast<int>(tTopo->tobRod(detid))-1); 
	  m_pixclus_module->push_back(static_cast<int>(tTopo->module(detid))-1
				      +limits[static_cast<int>(tTopo->layer(detid))-1][tTopo->tobSide(detid)-1]); 
	}
	else
	{ 
	  m_pixclus_ladder->push_back(static_cast<int>(tTopo->module(detid))-1); 
	  m_pixclus_module->push_back(static_cast<int>(tTopo->tobRod(detid))-1
				      +limits[static_cast<int>(tTopo->layer(detid))-1][tTopo->tobSide(detid)-1]); 
	}
      }

      if (endcap)
      {
	disk = 10+static_cast<int>(tTopo->tidWheel(detid))+abs(2-static_cast<int>(tTopo->side(detid)))*7;
	m_pixclus_layer->push_back(disk); 
	m_pixclus_ladder->push_back(tTopo->tidRing(detid)-1); 
	m_pixclus_module->push_back(tTopo->module(detid)-1); 
      }

      m_pixclus_nrow->push_back(rows);
      m_pixclus_ncolumn->push_back(cols);
      m_pixclus_pitchx->push_back(pitchX);
      m_pixclus_pitchy->push_back(pitchY);
      
      m_pclus++;
    }
  }

  PixelExtractor::fillTree();
}


//
// Method getting the info from an input file
//

void PixelExtractor::getInfo(int ievt) 
{
  m_tree->GetEntry(ievt); 
}

// Method initializing everything (to do for each event)

void PixelExtractor::reset()
{
  m_pclus = 0;
  m_nPU = 0;

  m_pixclus_x->clear(); 
  m_pixclus_y->clear(); 
  m_pixclus_z->clear(); 
  m_pixclus_e->clear(); 
  m_pixclus_row->clear();   
  m_pixclus_column->clear();
  m_pixclus_simhit->clear();
  m_pixclus_simhitID->clear();
  m_pixclus_evtID->clear();
  m_pixclus_bcID->clear();
  m_pixclus_layer->clear();  
  m_pixclus_module->clear(); 
  m_pixclus_ladder->clear(); 
  m_pixclus_nrow->clear();   
  m_pixclus_ncolumn->clear();
  m_pixclus_bot->clear();
  m_pixclus_type->clear();
  m_pixclus_pitchx->clear(); 
  m_pixclus_pitchy->clear(); 
}


void PixelExtractor::fillTree()
{
  m_tree->Fill(); 
}
 
void PixelExtractor::fillSize(int size)
{
  m_pclus=size;
}

int  PixelExtractor::getSize()
{
  return m_pclus;
}

