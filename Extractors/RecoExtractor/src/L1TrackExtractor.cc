#include "../interface/L1TrackExtractor.h"


L1TrackExtractor::L1TrackExtractor(edm::EDGetTokenT< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > STUB_tag, edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > PATT_tag,  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TC_tag,  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TRCK_tag, bool doTree)
{
  m_OK = false;
  n_tot_evt=0;

  m_STUB_tag   = STUB_tag;
  m_PATT_tag   = PATT_tag;
  m_TC_tag     = TC_tag;
  m_TRCK_tag   = TRCK_tag;

  // Tree definition
 
  m_patt_links   = new  std::vector< std::vector<int> >;
  m_patt_secid   = new  std::vector<int>;
  m_patt_pattid  = new  std::vector<int>;
  m_patt_miss    = new  std::vector<int>;

  m_tc_pt        = new  std::vector<float>;
  m_tc_eta       = new  std::vector<float>;
  m_tc_phi       = new  std::vector<float>;
  m_tc_z         = new  std::vector<float>;
  m_tc_links     = new  std::vector< std::vector<int> >;
  m_tc_secid     = new  std::vector<int>;
  m_tc_pattid    = new  std::vector<int>;

  m_trk_pt       = new  std::vector<float>;
  m_trk_eta      = new  std::vector<float>;
  m_trk_phi      = new  std::vector<float>;
  m_trk_z        = new  std::vector<float>;
  m_trk_chi      = new  std::vector<float>;
  m_trk_links    = new  std::vector< std::vector<int> >;
  m_trk_secid    = new  std::vector<int>;
  m_trk_pattid   = new  std::vector<int>;

  L1TrackExtractor::reset();

  if (doTree)
  {
    m_OK = true;

    m_tree      = new TTree("L1tracks","Official L1-AM tracks info") ;

    // Branches definition

    m_tree->Branch("L1evt", &n_tot_evt); // Simple evt number or event ID

    m_tree->Branch("L1PATT_n",           &m_patt);
    m_tree->Branch("L1PATT_links",       &m_patt_links);
    m_tree->Branch("L1PATT_secid",       &m_patt_secid);
    m_tree->Branch("L1PATT_pattid",      &m_patt_pattid);
    m_tree->Branch("L1PATT_nmiss",       &m_patt_miss);

    m_tree->Branch("L1TC_n",             &m_tc);
    m_tree->Branch("L1TC_links",         &m_tc_links);
    m_tree->Branch("L1TC_secid",         &m_tc_secid);
    m_tree->Branch("L1TC_pattid",        &m_tc_pattid);
    m_tree->Branch("L1TC_pt",            &m_tc_pt);
    m_tree->Branch("L1TC_phi",           &m_tc_phi);
    m_tree->Branch("L1TC_z",             &m_tc_z);
    m_tree->Branch("L1TC_eta",           &m_tc_eta);

    m_tree->Branch("L1TRK_n",            &m_trk);
    m_tree->Branch("L1TRK_links",        &m_trk_links);
    m_tree->Branch("L1TRK_secid",        &m_trk_secid);
    m_tree->Branch("L1TRK_pattid",       &m_trk_pattid);
    m_tree->Branch("L1TRK_pt",           &m_trk_pt);
    m_tree->Branch("L1TRK_phi",          &m_trk_phi);
    m_tree->Branch("L1TRK_z",            &m_trk_z);
    m_tree->Branch("L1TRK_eta",          &m_trk_eta);
    m_tree->Branch("L1TRK_chi2",         &m_trk_chi);
  }
}

L1TrackExtractor::L1TrackExtractor(TFile *a_file)
{
  std::cout << "L1TrackExtractor object is retrieved" << std::endl;
 
  m_patt_links   = new  std::vector< std::vector<int> >;
  m_patt_secid   = new  std::vector<int>;
  m_patt_pattid  = new  std::vector<int>;
  m_patt_miss    = new  std::vector<int>;
   
  m_tc_pt        = new  std::vector<float>;
  m_tc_eta       = new  std::vector<float>;
  m_tc_phi       = new  std::vector<float>;
  m_tc_z         = new  std::vector<float>;
  m_tc_links     = new  std::vector< std::vector<int> >;
  m_tc_secid     = new  std::vector<int>;
  m_tc_pattid    = new  std::vector<int>;

  m_trk_pt       = new  std::vector<float>;
  m_trk_eta      = new  std::vector<float>;
  m_trk_phi      = new  std::vector<float>;
  m_trk_z        = new  std::vector<float>;
  m_trk_chi      = new  std::vector<float>;
  m_trk_links    = new  std::vector< std::vector<int> >;
  m_trk_secid    = new  std::vector<int>;
  m_trk_pattid   = new  std::vector<int>;

  // Tree definition
  m_OK = false;


  L1TrackExtractor::reset();

  m_tree = dynamic_cast<TTree*>(a_file->Get("L1tracks"));

  if (!m_tree)
  {
    std::cout << "This tree (L1tracks) doesn't exist!!!" << std::endl;
    return;
  }

  m_OK = true;

  m_n_events = m_tree->GetEntries();

  m_tree->SetBranchAddress("L1evt", &n_tot_evt); // Simple evt number or event ID
  
  m_tree->SetBranchAddress("L1PATT_n",           &m_patt);
  m_tree->SetBranchAddress("L1PATT_links",       &m_patt_links);
  m_tree->SetBranchAddress("L1PATT_secid",       &m_patt_secid);
  m_tree->SetBranchAddress("L1PATT_pattid",      &m_patt_pattid);
  m_tree->SetBranchAddress("L1PATT_miss",        &m_patt_miss);
 
  m_tree->SetBranchAddress("L1TC_n",             &m_tc);
  m_tree->SetBranchAddress("L1TC_links",         &m_tc_links);
  m_tree->SetBranchAddress("L1TC_secid",         &m_tc_secid);
  m_tree->SetBranchAddress("L1TC_pattid",        &m_tc_pattid);
  m_tree->SetBranchAddress("L1TC_pt",            &m_tc_pt);
  m_tree->SetBranchAddress("L1TC_phi",           &m_tc_phi);
  m_tree->SetBranchAddress("L1TC_z",             &m_tc_z);
  m_tree->SetBranchAddress("L1TC_eta",           &m_tc_eta);
 
  m_tree->SetBranchAddress("L1TRK_n",            &m_trk);
  m_tree->SetBranchAddress("L1TRK_links",        &m_trk_links);
  m_tree->SetBranchAddress("L1TRK_secid",        &m_trk_secid);
  m_tree->SetBranchAddress("L1TRK_pattid",       &m_trk_pattid);
  m_tree->SetBranchAddress("L1TRK_pt",           &m_trk_pt);
  m_tree->SetBranchAddress("L1TRK_phi",          &m_trk_phi);
  m_tree->SetBranchAddress("L1TRK_z",            &m_trk_z);
  m_tree->SetBranchAddress("L1TRK_eta",          &m_trk_eta);
  m_tree->SetBranchAddress("L1TRK_chi2",         &m_trk_chi);
  
  std::cout << "This file contains " << m_n_events << " events..." << std::endl;
}



L1TrackExtractor::~L1TrackExtractor()
{}


void L1TrackExtractor::init(const edm::EventSetup *setup, bool isFlat)
{
  setup->get<TrackerTopologyRcd>().get(tTopoHandle);
  setup->get<TrackerDigiGeometryRecord>().get(tGeomHandle);

  m_tilted = (!isFlat);

  int n_tilted_rings[6];
  int n_flat_rings[6];

  for (int i=0; i < 6; ++i) n_tilted_rings[i]=0;
  for (int i=0; i < 6; ++i) n_flat_rings[i]=0;

  if (m_tilted)
  {
    n_tilted_rings[0]=11;
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

  /// Magnetic Field
  //  edm::ESHandle< MagneticField > magneticFieldHandle;
  // setup->get< IdealMagneticFieldRecord >().get(magneticFieldHandle);
  // const MagneticField* theMagneticField = magneticFieldHandle.product();
  // mMagneticFieldStrength = theMagneticField->inTesla(GlobalPoint(0,0,0)).z();

  mMagneticFieldStrength = 3.8;

}

//
// Method filling the main particle tree
//

void L1TrackExtractor::writeInfo(const edm::Event *event, StubExtractor *stub) 
{

  const TrackerTopology* const tTopo = tTopoHandle.product();
  const TrackerGeometry* const theTrackerGeom = tGeomHandle.product();

  L1TrackExtractor::reset();
  ++n_tot_evt;    

  /// Get the Stubs already stored away
  edm::Handle< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > TTStubHandle;
  edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TTPatternHandle;
  edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TTTCHandle;
  edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TTTrackHandle;

  event->getByToken( m_STUB_tag, TTStubHandle );
  event->getByToken( m_PATT_tag, TTPatternHandle );
  event->getByToken( m_TC_tag,   TTTCHandle );
  event->getByToken( m_TRCK_tag, TTTrackHandle );
  
  int layer  = 0;
  int ladder = 0;
  int module = 0;
  int type;
  std::vector<int> stub_list;

  /// STEP 1
  /// Loop over patterns

  std::vector< TTTrack< Ref_Phase2TrackerDigi_ > >::const_iterator iterTTTrack;
  
  /// Go on only if there are Patterns from PixelDigis
  if ( TTPatternHandle->size() > 0 )
  {
    /// Loop over Patterns
    unsigned int patCnt = 0;

    for ( iterTTTrack = TTPatternHandle->begin();
	  iterTTTrack != TTPatternHandle->end();
	  ++iterTTTrack )
    {
      edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > tempTrackPtr( TTPatternHandle, patCnt++ );

      /// Get everything is relevant

      std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_  > >, TTStub< Ref_Phase2TrackerDigi_  > > > trackStubs = tempTrackPtr->getStubRefs();

      // Loop over stubs contained in the pattern to recover the info
      // and match them with the stubs in the STUB_extractor container

      m_patt_secid->push_back(tempTrackPtr->getSector()%100);
      m_patt_pattid->push_back(tempTrackPtr->getSector()/100);
      m_patt_miss->push_back(tempTrackPtr->getWedge());

      stub_list.clear();

      for(unsigned int i=0;i<trackStubs.size();i++)
      {

	edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > tempStubRef = trackStubs.at(i);

	DetId stackDetid = tempStubRef->getDetId();
	DetId detid = stackDetid;

	// DIRTY!! but need this loop to get back the geographicalId from the detid 
	// Is there a method in TrackerTopology.h for that????
	for (auto gd=tGeomHandle->dets().begin(); gd != tGeomHandle->dets().end(); gd++) 
	{
	    DetId detidg = (*gd)->geographicalId();
	    if(detidg.subdetId()!=StripSubdetector::TOB && detidg.subdetId()!=StripSubdetector::TID ) continue; // only run on OT
	    if(tTopoHandle->stack(detidg)!=stackDetid) continue; 

	    detid=detidg;
	    break;
	}

	const GeomDetUnit* det0 = tGeomHandle->idToDetUnit( detid );
	const PixelGeomDetUnit* theGeomDet = dynamic_cast< const PixelGeomDetUnit* >( det0 );
	const PixelTopology* topol = dynamic_cast< const PixelTopology* >( &(theGeomDet->specificTopology()) );

	/// Calculate average coordinates col/row for inner/outer Cluster
	/// These are already corrected for being at the center of each pixel
	MeasurementPoint coords = tempStubRef->getClusterRef(0)->findAverageLocalCoordinatesCentered();
	LocalPoint clustlp      = topol->localPosition(coords);
	GlobalPoint posStub     =  theGeomDet->surface().toGlobal(clustlp);


	if ( detid.subdetId()==StripSubdetector::TOB )
	{
	  type   = static_cast<int>(tTopo->tobSide(detid)); // Tilt-/Tilt+/Flat <-> 1/2/3
	  layer  = static_cast<int>(tTopo->layer(detid))+4;
	  ladder = static_cast<int>(tTopo->tobRod(detid))-1;
	  module = static_cast<int>(tTopo->module(detid))-1+limits[layer-5][type-1];

	  if (type<3)
	  {
	    ladder = static_cast<int>(tTopo->module(detid))-1;
	    module = static_cast<int>(tTopo->tobRod(detid))-1+limits[layer-5][type-1];
	  }
	}
	else if ( detid.subdetId()==StripSubdetector::TID )
	{	
	  layer  = 10+static_cast<int>(tTopo->tidWheel(detid))+abs(2-static_cast<int>(tTopo->side(detid)))*7;
	  ladder = static_cast<int>(tTopo->tidRing(detid))-1;
	  module = static_cast<int>(tTopo->module(detid))-1;
	  type   = 0;
	}


	double SW = tempStubRef->getTriggerDisplacement();

	stub_list.push_back(stub->get_id(layer,ladder,module,posStub.x(),posStub.y(),posStub.z(),SW));


      } /// End of loop over track stubs

      m_patt_links->push_back(stub_list);


    } // End of loop over patterns

    m_patt =  patCnt;
  }

  /// STEP 2
  /// Loop over TCs

  if ( TTTCHandle->size() > 0 )
  {
    /// Loop over TCs
    unsigned int tkCnt = 0;
    
    //      std::cout << TTTCHandle->size() << std::endl;
    
    for ( iterTTTrack = TTTCHandle->begin();
	  iterTTTrack != TTTCHandle->end();
	  ++iterTTTrack )
    {
      edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > tempTrackPtr( TTTCHandle, tkCnt++ );
      
      /// Get everything is relevant
	
      std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_  > >, TTStub< Ref_Phase2TrackerDigi_  > > > trackStubs = tempTrackPtr->getStubRefs();

      // Loop over stubs contained in the pattern to recover the info
      // and match them with the stubs in the STUB_extractor container
	
      m_tc_secid->push_back(tempTrackPtr->getSector()%100);
      m_tc_pattid->push_back(tempTrackPtr->getSector()/100);
      m_tc_pt->push_back(tempTrackPtr->getMomentum(5).perp() );
      m_tc_eta->push_back(tempTrackPtr->getMomentum(5).eta());
      m_tc_phi->push_back(tempTrackPtr->getMomentum(5).phi());
      m_tc_z->push_back(tempTrackPtr->getPOCA(5).z());
      
      stub_list.clear();
	
      for(unsigned int i=0;i<trackStubs.size();i++)
      {
	edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > tempStubRef = trackStubs.at(i);


	DetId stackDetid = tempStubRef->getDetId();
	DetId detid = stackDetid;

	// DIRTY!! but need this loop to get back the geographicalId from the detid 
	// Is there a method in TrackerTopology.h for that????
	for (auto gd=tGeomHandle->dets().begin(); gd != tGeomHandle->dets().end(); gd++) 
	{
	    DetId detidg = (*gd)->geographicalId();
	    if(detidg.subdetId()!=StripSubdetector::TOB && detidg.subdetId()!=StripSubdetector::TID ) continue; // only run on OT
	    if(tTopoHandle->stack(detidg)!=stackDetid) continue; 

	    detid=detidg;
	    break;
	}

	const GeomDetUnit* det0 = tGeomHandle->idToDetUnit( detid );
	const PixelGeomDetUnit* theGeomDet = dynamic_cast< const PixelGeomDetUnit* >( det0 );
	const PixelTopology* topol = dynamic_cast< const PixelTopology* >( &(theGeomDet->specificTopology()) );

	/// Calculate average coordinates col/row for inner/outer Cluster
	/// These are already corrected for being at the center of each pixel
	MeasurementPoint coords = tempStubRef->getClusterRef(0)->findAverageLocalCoordinatesCentered();
	LocalPoint clustlp   = topol->localPosition(coords);
	GlobalPoint posStub  =  theGeomDet->surface().toGlobal(clustlp);

	if ( detid.subdetId()==StripSubdetector::TOB )
	{
	  type   = static_cast<int>(tTopo->tobSide(detid)); // Tilt-/Tilt+/Flat <-> 1/2/3
	  layer  = static_cast<int>(tTopo->layer(detid))+4;
	  ladder = static_cast<int>(tTopo->tobRod(detid))-1;
	  module = static_cast<int>(tTopo->module(detid))-1+limits[layer-5][type-1];

	  if (type<3)
	  {
	    ladder = static_cast<int>(tTopo->module(detid))-1;
	    module = static_cast<int>(tTopo->tobRod(detid))-1+limits[layer-5][type-1];
	  }
	}
	else if ( detid.subdetId()==StripSubdetector::TID )
	{	
	  layer  = 10+static_cast<int>(tTopo->tidWheel(detid))+abs(2-static_cast<int>(tTopo->side(detid)))*7;
	  ladder = static_cast<int>(tTopo->tidRing(detid))-1;
	  module = static_cast<int>(tTopo->module(detid))-1;
	  type   = 0;
	}

	double SW = tempStubRef->getTriggerDisplacement();

	stub_list.push_back(stub->get_id(layer,ladder,module,posStub.x(),posStub.y(),posStub.z(),SW));

      } /// End of loop over track stubs	
	
      m_tc_links->push_back(stub_list);

    } // End of loop over patterns
      
    m_tc =  tkCnt;
  }      

  /// STEP 3
  /// Loop over fitted tracks

  if ( TTTrackHandle->size() > 0 )
  {
    /// Loop over Patterns
    unsigned int tkCnt = 0;

    for ( iterTTTrack = TTTrackHandle->begin();
	  iterTTTrack != TTTrackHandle->end();
	  ++iterTTTrack )
    {
      edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > tempTrackPtr( TTTrackHandle, tkCnt++ );
      
      /// Get everything is relevant
		
      std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_  > >, TTStub< Ref_Phase2TrackerDigi_  > > > trackStubs = tempTrackPtr->getStubRefs();
      
      // Loop over stubs contained in the pattern to recover the info
      // and match them with the stubs in the STUB_extractor container
      
      m_trk_secid->push_back(tempTrackPtr->getSector()%100);
      m_trk_pattid->push_back(tempTrackPtr->getSector()/100);
      m_trk_pt->push_back(tempTrackPtr->getMomentum(5).perp() );
      m_trk_eta->push_back(tempTrackPtr->getMomentum(5).eta());
      m_trk_phi->push_back(tempTrackPtr->getMomentum(5).phi());
      m_trk_z->push_back(tempTrackPtr->getPOCA(5).z());
      m_trk_chi->push_back(tempTrackPtr->getChi2(5));
      
      stub_list.clear();
      
      for(unsigned int i=0;i<trackStubs.size();i++)
      {
	edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > tempStubRef = trackStubs.at(i);


	DetId stackDetid = tempStubRef->getDetId();
	DetId detid = stackDetid;

	// DIRTY!! but need this loop to get back the geographicalId from the detid 
	// Is there a method in TrackerTopology.h for that????
	for (auto gd=tGeomHandle->dets().begin(); gd != tGeomHandle->dets().end(); gd++) 
	{
	    DetId detidg = (*gd)->geographicalId();
	    if(detidg.subdetId()!=StripSubdetector::TOB && detidg.subdetId()!=StripSubdetector::TID ) continue; // only run on OT
	    if(tTopoHandle->stack(detidg)!=stackDetid) continue; 

	    detid=detidg;
	    break;
	}

	const GeomDetUnit* det0 = tGeomHandle->idToDetUnit( detid );
	const PixelGeomDetUnit* theGeomDet = dynamic_cast< const PixelGeomDetUnit* >( det0 );
	const PixelTopology* topol = dynamic_cast< const PixelTopology* >( &(theGeomDet->specificTopology()) );

	/// Calculate average coordinates col/row for inner/outer Cluster
	/// These are already corrected for being at the center of each pixel
	MeasurementPoint coords = tempStubRef->getClusterRef(0)->findAverageLocalCoordinatesCentered();
	LocalPoint clustlp   = topol->localPosition(coords);
	GlobalPoint posStub  =  theGeomDet->surface().toGlobal(clustlp);


	/// Find pixel pitch and topology related information

	if ( detid.subdetId()==StripSubdetector::TOB )
	{
	  type   = static_cast<int>(tTopo->tobSide(detid)); // Tilt-/Tilt+/Flat <-> 1/2/3
	  layer  = static_cast<int>(tTopo->layer(detid))+4;
	  ladder = static_cast<int>(tTopo->tobRod(detid))-1;
	  module = static_cast<int>(tTopo->module(detid))-1+limits[layer-5][type-1];

	  if (type<3)
	  {
	    ladder = static_cast<int>(tTopo->module(detid))-1;
	    module = static_cast<int>(tTopo->tobRod(detid))-1+limits[layer-5][type-1];
	  }
	}
	else if ( detid.subdetId()==StripSubdetector::TID )
	{	
	  layer  = 10+static_cast<int>(tTopo->tidWheel(detid))+abs(2-static_cast<int>(tTopo->side(detid)))*7;
	  ladder = static_cast<int>(tTopo->tidRing(detid))-1;
	  module = static_cast<int>(tTopo->module(detid))-1;
	  type   = 0;
	}

	double SW = tempStubRef->getTriggerDisplacement();
	
	stub_list.push_back(stub->get_id(layer,ladder,module,posStub.x(),posStub.y(),posStub.z(),SW));
      } /// End of loop over track stubs	
	
      m_trk_links->push_back(stub_list);

    } // End of loop over patterns
      
    m_trk =  tkCnt;
  }      

  L1TrackExtractor::fillTree();
}


//
// Method getting the info from an input file
//

void L1TrackExtractor::getInfo(int ievt) 
{
  m_tree->GetEntry(ievt); 
}

// Method initializing everything (to do for each event)

void L1TrackExtractor::reset()
{
  m_patt = 0;
  m_trk  = 0;
  m_tc  = 0; 

  m_patt_links->clear(); 
  m_patt_secid->clear(); 
  m_patt_pattid->clear(); 
  m_patt_miss->clear(); 
  
  m_tc_pt->clear();    
  m_tc_eta->clear();   
  m_tc_phi->clear();   
  m_tc_z->clear();     
  m_tc_links->clear(); 
  m_tc_secid->clear();
  m_tc_pattid->clear();
 
  m_trk_pt->clear();    
  m_trk_eta->clear();   
  m_trk_phi->clear();   
  m_trk_z->clear();     
  m_trk_links->clear();
  m_trk_secid->clear(); 
  m_trk_pattid->clear(); 
  m_trk_chi->clear();
}


void L1TrackExtractor::fillTree()
{
  m_tree->Fill(); 
}
 
void L1TrackExtractor::fillSize(int size)
{
  m_patt=size;
}

int  L1TrackExtractor::getSize()
{
  return m_patt;
}

