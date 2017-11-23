#include "../interface/StubExtractor.h"


StubExtractor::StubExtractor(edm::EDGetTokenT< edmNew::DetSetVector< TTCluster< Ref_Phase2TrackerDigi_ > > > ctoken,edm::EDGetTokenT< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > stoken, edm::EDGetTokenT< TTClusterAssociationMap< Ref_Phase2TrackerDigi_ > > cttoken, edm::EDGetTokenT< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > sttoken, edm::EDGetTokenT< std::vector< TrackingParticle > > tptoken, edm::EDGetTokenT< std::vector< TrackingVertex > > tvtoken, edm::EDGetTokenT< edm::SimTrackContainer > simttoken, edm::EDGetTokenT< edm::SimVertexContainer > simvtoken, bool doTree)
{

  m_ctoken = ctoken;
  m_cttoken= cttoken;
  m_stoken = stoken;
  m_sttoken= sttoken;
  m_tptoken= tptoken;
  m_tvtoken= tvtoken;
  m_simttoken=simttoken;
  m_simvtoken=simvtoken;


  m_OK = false;
  n_tot_evt=0;

  // Tree definition
 
  m_clus_x       = new  std::vector<float>;
  m_clus_y       = new  std::vector<float>;
  m_clus_z       = new  std::vector<float>; 
  m_clus_e       = new  std::vector<float>;
  m_clus_layer   = new  std::vector<int>;
  m_clus_module  = new  std::vector<int>;  
  m_clus_ladder  = new  std::vector<int>;  
  m_clus_seg     = new  std::vector<int>; 
  m_clus_type    = new  std::vector<int>; 
  m_clus_strip   = new  std::vector<float>; 
  m_clus_used    = new  std::vector<int>;  
  m_clus_sat     = new  std::vector<int>; 
  m_clus_nstrips = new  std::vector<int>; 
  m_clus_matched = new  std::vector<int>;  
  m_clus_PS      = new  std::vector<int>;  
  m_clus_nrows   = new  std::vector<int>;
  m_clus_bot     = new  std::vector<int>;
  m_clus_pid     = new  std::vector<int>;  
  m_clus_pdgID   = new  std::vector<int>;  
  m_clus_ptGEN   = new  std::vector<float>; 
  m_clus_tp      = new  std::vector<int>;  
  m_clus_pix     = new  std::vector<std::vector<int> >; 
  m_stub_pt      = new  std::vector<float>;
  m_stub_ptMC    = new  std::vector<float>; 
  m_stub_pxGEN   = new  std::vector<float>; 
  m_stub_pyGEN   = new  std::vector<float>; 
  m_stub_etaGEN  = new  std::vector<float>; 
  m_stub_X0      = new  std::vector<float>; 
  m_stub_Y0      = new  std::vector<float>; 
  m_stub_Z0      = new  std::vector<float>; 
  m_stub_PHI0    = new  std::vector<float>; 
  m_stub_layer   = new  std::vector<int>; 
  m_stub_module  = new  std::vector<int>;  
  m_stub_ladder  = new  std::vector<int>; 
  m_stub_detid   = new  std::vector<int>;  
  m_stub_seg     = new  std::vector<int>;  
  m_stub_type    = new  std::vector<int>; 
  m_stub_chip    = new  std::vector<int>;  
  m_stub_strip   = new  std::vector<float>; 
  m_stub_x       = new  std::vector<float>;  
  m_stub_y       = new  std::vector<float>;  
  m_stub_z       = new  std::vector<float>;  
  m_stub_clust1  = new  std::vector<int>;  
  m_stub_clust2  = new  std::vector<int>;  
  m_stub_cw1     = new  std::vector<int>;  
  m_stub_cw2     = new  std::vector<int>;  
  m_stub_deltas  = new  std::vector<float>;  
  m_stub_deltasf = new  std::vector<float>;  
  m_stub_deltash = new  std::vector<float>;  
  m_stub_cor     = new  std::vector<float>;  
  m_stub_corf    = new  std::vector<float>;  
  m_stub_tp      = new  std::vector<int>;  
  m_stub_pdg     = new  std::vector<int>;  
  m_stub_pid     = new  std::vector<int>;  
  m_stub_rank    = new  std::vector<int>;  

  StubExtractor::reset();

  if (doTree)
  {
    m_OK = true;

    m_tree      = new TTree("TkStubs","Official stub info") ;

    // Branches definition

    m_tree->Branch("L1Tkevt", &n_tot_evt); // Simple evt number or event ID

    // If we don't request only matched stubs, we keep all the info
    // otherwise we skim the data file (useful for BANK generation)

    m_tree->Branch("L1TkCLUS_n",         &m_clus);
    m_tree->Branch("L1TkCLUS_x",         &m_clus_x);
    m_tree->Branch("L1TkCLUS_y",         &m_clus_y);
    m_tree->Branch("L1TkCLUS_z",         &m_clus_z);
    m_tree->Branch("L1TkCLUS_charge",    &m_clus_e);
    m_tree->Branch("L1TkCLUS_layer",     &m_clus_layer);
    m_tree->Branch("L1TkCLUS_module",    &m_clus_module);
    m_tree->Branch("L1TkCLUS_ladder",    &m_clus_ladder);
    m_tree->Branch("L1TkCLUS_seg",       &m_clus_seg);
    m_tree->Branch("L1TkCLUS_type",      &m_clus_type);
    m_tree->Branch("L1TkCLUS_strip",     &m_clus_strip);
    m_tree->Branch("L1TkCLUS_nstrip",    &m_clus_nstrips);
    m_tree->Branch("L1TkCLUS_nsat",      &m_clus_sat);
    m_tree->Branch("L1TkCLUS_match",     &m_clus_matched);
    m_tree->Branch("L1TkCLUS_PS",        &m_clus_PS);
    m_tree->Branch("L1TkCLUS_nrows",     &m_clus_nrows);
    m_tree->Branch("L1TkCLUS_bottom",    &m_clus_bot);
    m_tree->Branch("L1TkCLUS_pdgID",     &m_clus_pdgID);
    m_tree->Branch("L1TkCLUS_ptGEN",     &m_clus_ptGEN);
    m_tree->Branch("L1TkCLUS_process",   &m_clus_pid);
    m_tree->Branch("L1TkCLUS_PIX",       &m_clus_pix);
    m_tree->Branch("L1TkCLUS_tp",        &m_clus_tp);

    m_tree->Branch("L1TkSTUB_clust1",    &m_stub_clust1);
    m_tree->Branch("L1TkSTUB_clust2",    &m_stub_clust2);
    m_tree->Branch("L1TkSTUB_cor",       &m_stub_cor);
    m_tree->Branch("L1TkSTUB_corf",      &m_stub_corf);
    m_tree->Branch("L1TkSTUB_PHI0",      &m_stub_PHI0);
    m_tree->Branch("L1TkSTUB_tp",        &m_stub_tp);
    m_tree->Branch("L1TkSTUB_pdgID",     &m_stub_pdg);
    m_tree->Branch("L1TkSTUB_process",   &m_stub_pid);
    m_tree->Branch("L1TkSTUB_n",         &m_stub);
    m_tree->Branch("L1TkSTUB_pt",        &m_stub_pt);
    m_tree->Branch("L1TkSTUB_pxGEN",     &m_stub_pxGEN);
    m_tree->Branch("L1TkSTUB_pyGEN",     &m_stub_pyGEN);
    m_tree->Branch("L1TkSTUB_etaGEN",    &m_stub_etaGEN);
    m_tree->Branch("L1TkSTUB_layer",     &m_stub_layer);
    m_tree->Branch("L1TkSTUB_module",    &m_stub_module);
    m_tree->Branch("L1TkSTUB_ladder",    &m_stub_ladder);
    m_tree->Branch("L1TkSTUB_seg",       &m_stub_seg);
    m_tree->Branch("L1TkSTUB_type",      &m_stub_type);
    m_tree->Branch("L1TkSTUB_chip",      &m_stub_chip);
    m_tree->Branch("L1TkSTUB_strip",     &m_stub_strip);
    m_tree->Branch("L1TkSTUB_detid",     &m_stub_detid);
    m_tree->Branch("L1TkSTUB_x",         &m_stub_x);
    m_tree->Branch("L1TkSTUB_y",         &m_stub_y);
    m_tree->Branch("L1TkSTUB_z",         &m_stub_z);
    m_tree->Branch("L1TkSTUB_deltas",    &m_stub_deltas);
    m_tree->Branch("L1TkSTUB_deltasf",   &m_stub_deltasf);
    m_tree->Branch("L1TkSTUB_deltash",   &m_stub_deltash);
    m_tree->Branch("L1TkSTUB_X0",        &m_stub_X0);
    m_tree->Branch("L1TkSTUB_Y0",        &m_stub_Y0);
    m_tree->Branch("L1TkSTUB_Z0",        &m_stub_Z0);
    m_tree->Branch("L1TkSTUB_rank",      &m_stub_rank);
  }
}

StubExtractor::StubExtractor(TFile *a_file)
{
  std::cout << "StubExtractor object is retrieved" << std::endl;
 
  m_clus_x       = new  std::vector<float>;
  m_clus_y       = new  std::vector<float>;
  m_clus_z       = new  std::vector<float>; 
  m_clus_e       = new  std::vector<float>;
  m_clus_layer   = new  std::vector<int>;
  m_clus_module  = new  std::vector<int>;  
  m_clus_ladder  = new  std::vector<int>;  
  m_clus_seg     = new  std::vector<int>; 
  m_clus_type    = new  std::vector<int>; 
  m_clus_strip   = new  std::vector<float>; 
  m_clus_used    = new  std::vector<int>;  
  m_clus_sat     = new  std::vector<int>; 
  m_clus_bot     = new  std::vector<int>; 
  m_clus_nstrips = new  std::vector<int>; 
  m_clus_matched = new  std::vector<int>;  
  m_clus_PS      = new  std::vector<int>;  
  m_clus_nrows   = new  std::vector<int>;
  m_clus_pid     = new  std::vector<int>;  
  m_clus_pdgID   = new  std::vector<int>;  
  m_clus_ptGEN   = new  std::vector<float>; 
  m_clus_tp      = new  std::vector<int>;  
  m_clus_pix     = new  std::vector<std::vector<int> >; 

  m_stub_pt      = new  std::vector<float>;
  m_stub_ptMC    = new  std::vector<float>; 
  m_stub_pxGEN   = new  std::vector<float>; 
  m_stub_pyGEN   = new  std::vector<float>; 
  m_stub_etaGEN  = new  std::vector<float>; 
  m_stub_X0      = new  std::vector<float>; 
  m_stub_Y0      = new  std::vector<float>; 
  m_stub_Z0      = new  std::vector<float>; 
  m_stub_PHI0    = new  std::vector<float>; 
  m_stub_layer   = new  std::vector<int>; 
  m_stub_module  = new  std::vector<int>;  
  m_stub_ladder  = new  std::vector<int>; 
  m_stub_seg     = new  std::vector<int>;  
  m_stub_type    = new  std::vector<int>;  
  m_stub_strip   = new  std::vector<float>; 
  m_stub_chip    = new  std::vector<int>; 
  m_stub_detid   = new  std::vector<int>; 
  m_stub_x       = new  std::vector<float>;  
  m_stub_y       = new  std::vector<float>;  
  m_stub_z       = new  std::vector<float>;  
  m_stub_clust1  = new  std::vector<int>;  
  m_stub_clust2  = new  std::vector<int>;  
  m_stub_cw1     = new  std::vector<int>;  
  m_stub_cw2     = new  std::vector<int>;  
  m_stub_deltas  = new  std::vector<float>;  
  m_stub_deltasf = new  std::vector<float>;  
  m_stub_deltash = new  std::vector<float>;  
  m_stub_cor     = new  std::vector<float>;  
  m_stub_corf    = new  std::vector<float>;  
  m_stub_tp      = new  std::vector<int>;  
  m_stub_pdg     = new  std::vector<int>;  
  m_stub_pid     = new  std::vector<int>; 
  m_stub_rank    = new  std::vector<int>; 

  // Tree definition
  m_OK = false;


  StubExtractor::reset();

  m_tree = dynamic_cast<TTree*>(a_file->Get("TkStubs"));

  if (!m_tree)
  {
    std::cout << "This tree (TkStubs) doesn't exist!!!" << std::endl;
    return;
  }

  m_OK = true;
  m_matching = false;

  m_n_events = m_tree->GetEntries();

  m_tree->SetBranchAddress("L1Tkevt", &n_tot_evt); // Simple evt number or event ID
  
  // If we don't request only matched stubs, we keep all the info
  // otherwise we skim the data file (useful for BANK generation)
  
  m_tree->SetBranchAddress("L1TkCLUS_n",         &m_clus);
  m_tree->SetBranchAddress("L1TkCLUS_x",         &m_clus_x);
  m_tree->SetBranchAddress("L1TkCLUS_y",         &m_clus_y);
  m_tree->SetBranchAddress("L1TkCLUS_z",         &m_clus_z);
  m_tree->SetBranchAddress("L1TkCLUS_charge",    &m_clus_e);
  m_tree->SetBranchAddress("L1TkCLUS_layer",     &m_clus_layer);
  m_tree->SetBranchAddress("L1TkCLUS_module",    &m_clus_module);
  m_tree->SetBranchAddress("L1TkCLUS_ladder",    &m_clus_ladder);
  m_tree->SetBranchAddress("L1TkCLUS_seg",       &m_clus_seg);
  m_tree->SetBranchAddress("L1TkCLUS_type",      &m_clus_type);
  m_tree->SetBranchAddress("L1TkCLUS_strip",     &m_clus_strip);
  m_tree->SetBranchAddress("L1TkCLUS_nstrip",    &m_clus_nstrips);
  m_tree->SetBranchAddress("L1TkCLUS_nsat",      &m_clus_sat);
  m_tree->SetBranchAddress("L1TkCLUS_match",     &m_clus_matched);
  m_tree->SetBranchAddress("L1TkCLUS_PS",        &m_clus_PS);
  m_tree->SetBranchAddress("L1TkCLUS_nrows",     &m_clus_nrows);
  m_tree->SetBranchAddress("L1TkCLUS_bottom",    &m_clus_bot);
  m_tree->SetBranchAddress("L1TkCLUS_pdgID",     &m_clus_pdgID);
  m_tree->SetBranchAddress("L1TkCLUS_ptGEN",     &m_clus_ptGEN);
  m_tree->SetBranchAddress("L1TkCLUS_process",   &m_clus_pid);
  m_tree->SetBranchAddress("L1TkCLUS_tp",        &m_clus_tp);

  m_tree->SetBranchAddress("L1TkSTUB_clust1",    &m_stub_clust1);
  m_tree->SetBranchAddress("L1TkSTUB_clust2",    &m_stub_clust2);
  m_tree->SetBranchAddress("L1TkSTUB_cor",       &m_stub_cor);
  m_tree->SetBranchAddress("L1TkSTUB_corf",      &m_stub_corf);
  m_tree->SetBranchAddress("L1TkSTUB_PHI0",      &m_stub_PHI0);
  m_tree->SetBranchAddress("L1TkSTUB_tp",        &m_stub_tp);
  m_tree->SetBranchAddress("L1TkSTUB_pdgID",     &m_stub_pdg);
  m_tree->SetBranchAddress("L1TkSTUB_process",   &m_stub_pid);
  m_tree->SetBranchAddress("L1TkSTUB_n",         &m_stub);
  m_tree->SetBranchAddress("L1TkSTUB_pt",        &m_stub_pt);
  m_tree->SetBranchAddress("L1TkSTUB_pxGEN",     &m_stub_pxGEN);
  m_tree->SetBranchAddress("L1TkSTUB_pyGEN",     &m_stub_pyGEN);
  m_tree->SetBranchAddress("L1TkSTUB_etaGEN",    &m_stub_etaGEN);
  m_tree->SetBranchAddress("L1TkSTUB_layer",     &m_stub_layer);
  m_tree->SetBranchAddress("L1TkSTUB_module",    &m_stub_module);
  m_tree->SetBranchAddress("L1TkSTUB_ladder",    &m_stub_ladder);
  m_tree->SetBranchAddress("L1TkSTUB_seg",       &m_stub_seg);
  m_tree->SetBranchAddress("L1TkSTUB_type",      &m_stub_type);
  m_tree->SetBranchAddress("L1TkSTUB_strip",     &m_stub_strip);
  m_tree->SetBranchAddress("L1TkSTUB_chip",      &m_stub_chip);
  m_tree->SetBranchAddress("L1TkSTUB_detid",     &m_stub_detid);
  m_tree->SetBranchAddress("L1TkSTUB_x",         &m_stub_x);
  m_tree->SetBranchAddress("L1TkSTUB_y",         &m_stub_y);
  m_tree->SetBranchAddress("L1TkSTUB_z",         &m_stub_z);
  m_tree->SetBranchAddress("L1TkSTUB_deltas",    &m_stub_deltas);
  m_tree->SetBranchAddress("L1TkSTUB_deltasf",   &m_stub_deltasf);
  m_tree->SetBranchAddress("L1TkSTUB_deltash",   &m_stub_deltash);
  m_tree->SetBranchAddress("L1TkSTUB_X0",        &m_stub_X0);
  m_tree->SetBranchAddress("L1TkSTUB_Y0",        &m_stub_Y0);
  m_tree->SetBranchAddress("L1TkSTUB_Z0",        &m_stub_Z0);
  m_tree->SetBranchAddress("L1TkSTUB_rank",      &m_stub_rank);

  std::cout << "This file contains " << m_n_events << " events..." << std::endl;

}



StubExtractor::~StubExtractor()
{}


void StubExtractor::init(const edm::EventSetup *setup, bool isFlat)
{
  setup->get<TrackerTopologyRcd>().get(tTopoHandle);
  setup->get<TrackerDigiGeometryRecord>().get(tGeomHandle);

  m_tilted=(!isFlat);



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

void StubExtractor::writeInfo(const edm::Event *event, MCExtractor *mc, bool MCinfo) 
{
  StubExtractor::reset();
  ++n_tot_evt;    


  const TrackerTopology* const tTopo = tTopoHandle.product();
  const TrackerGeometry* const theTrackerGeom = tGeomHandle.product();

  if (MCinfo) mc->clearTP(0.001,10000000.0);
  /// Sim Tracks and Vtx
  event->getByToken(m_simttoken,SimTrackHandle);  
  event->getByToken(m_simvtoken,SimVtxHandle);  

  /// Track Trigger
  event->getByToken( m_ctoken, PixelDigiL1TkClusterHandle );
  event->getByToken( m_stoken, PixelDigiL1TkStubHandle );

  /// Track Trigger MC Truth
  event->getByToken( m_cttoken, MCTruthTTClusterHandle );
  event->getByToken( m_sttoken, MCTruthTTStubHandle );


  /// TrackingParticles

  event->getByToken(m_tptoken,TrackingParticleHandle);  
  event->getByToken(m_tvtoken,TrackingVertexHandle);  


  int layer  = 0;
  int ladder = 0;
  int module = 0;
  int type   = 0;
  int segs   = 0;
  int rows   = 0;

  std::vector<int> clus_coords;

  std::vector<int> clus_rows;
  std::vector<int> clus_cols;

  LocalPoint  clustlp;
  GlobalPoint posClu;  
  GlobalPoint posStub;

  //  std::cout << "A" << std::endl;

  /// Go on only if there are L1TkCluster from PixelDigis
  if ( PixelDigiL1TkClusterHandle->size() > 0 )
  {
    /// Loop over L1TkClusters
    for (auto gd=theTrackerGeom->dets().begin(); gd != theTrackerGeom->dets().end(); gd++) 
    {
      DetId detid = (*gd)->geographicalId();
      if(detid.subdetId()!=StripSubdetector::TOB && detid.subdetId()!=StripSubdetector::TID ) continue; // only run on OT

      if (PixelDigiL1TkClusterHandle->find( detid ) == PixelDigiL1TkClusterHandle->end() ) continue;

      /// Get the DetSets of the Clusters
      edmNew::DetSet< TTCluster< Ref_Phase2TrackerDigi_ > > clusters = (*PixelDigiL1TkClusterHandle)[ detid ];
      const GeomDetUnit* det0 = theTrackerGeom->idToDetUnit( detid );
      const PixelGeomDetUnit* theGeomDet = dynamic_cast< const PixelGeomDetUnit* >( det0 );
      const PixelTopology* topol = dynamic_cast< const PixelTopology* >( &(theGeomDet->specificTopology()) );

      for ( auto clusterIter = clusters.begin();clusterIter != clusters.end();++clusterIter ) 
      {

	edm::Ref< edmNew::DetSetVector< TTCluster< Ref_Phase2TrackerDigi_  > >, TTCluster< Ref_Phase2TrackerDigi_  > > tempCluRef = 
	  edmNew::makeRefTo( PixelDigiL1TkClusterHandle, clusterIter );
       	bool genuineClu = MCTruthTTClusterHandle->isGenuine( tempCluRef );

	MeasurementPoint coords = tempCluRef->findAverageLocalCoordinatesCentered();

	clustlp = topol->localPosition(coords);
	posClu  =  theGeomDet->surface().toGlobal(clustlp);

	clus_rows = tempCluRef->getRows();
	clus_cols = tempCluRef->getCols();
	clus_coords.clear();

	for (unsigned int i=0;i<clus_rows.size();++i) 
	{
	  clus_coords.push_back(clus_rows.at(i));
	  clus_coords.push_back(clus_cols.at(i));
	}

	++m_clus;

	m_clus_pix->push_back(clus_coords);

	m_clus_x->push_back(posClu.x());
	m_clus_y->push_back(posClu.y());
	m_clus_z->push_back(posClu.z());
	
	m_clus_seg->push_back(coords.y());
	m_clus_strip->push_back(coords.x());
	m_clus_nstrips->push_back(tempCluRef->findWidth());
	
	segs= topol->ncolumns();
	rows= topol->nrows();;

	m_clus_bot->push_back(static_cast<int>(tTopo->isLower(detid)));

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
	
	m_clus_layer->push_back(layer);
	m_clus_ladder->push_back(ladder);
	m_clus_module->push_back(module);
	m_clus_type->push_back(type);
	m_clus_PS->push_back(segs);

	if (genuineClu)
	{
	  edm::Ptr< TrackingParticle > tpPtr = MCTruthTTClusterHandle->findTrackingParticlePtr( tempCluRef );

	  if ( tpPtr.isNull() )  
	  {
	    m_clus_matched->push_back(0);
	    m_clus_ptGEN->push_back(0); 
	    m_clus_pdgID->push_back(0);
	    m_clus_tp->push_back(-1);
	  }
	  else
	  {
	    m_clus_matched->push_back(1);
	    m_clus_ptGEN->push_back(tpPtr->p4().pt()); 
	    m_clus_pdgID->push_back(tpPtr->pdgId());

	    if (MCinfo) m_clus_tp->push_back(mc->getMatchingTP(tpPtr->vertex().x(),tpPtr->vertex().y(),tpPtr->vertex().z(),
							       tpPtr->p4().px(),tpPtr->p4().py(),tpPtr->p4().pz()));
	  }
	}
	else
	{
	  m_clus_matched->push_back(0);
	  m_clus_ptGEN->push_back(0); 
	  m_clus_pdgID->push_back(0);
	  m_clus_tp->push_back(-1);
	}

	m_clus_sat->push_back(0);
	m_clus_nrows->push_back(rows);
	
      } /// End of Loop over L1TkClusters in the detId
    }
  } /// End of if ( PixelDigiL1TkClusterHandle->size() > 0 )

  // Clusters are built, now look at the stubs
    
  int clust1     = -1;
  int clust2     = -1;
 
  /// Go on only if there are L1TkStub from PixelDigis
  if ( PixelDigiL1TkStubHandle->size() > 0 )
  {
    /// Loop over L1TkStubs

    /// Loop over L1TkClusters
    for (auto gd=theTrackerGeom->dets().begin(); gd != theTrackerGeom->dets().end(); gd++) 
    {
      DetId detid = (*gd)->geographicalId();
      if(detid.subdetId()!=StripSubdetector::TOB && detid.subdetId()!=StripSubdetector::TID ) continue; // only run on OT
      if(!tTopo->isLower(detid) ) continue; // loop on the stacks: choose the lower arbitrarily
      DetId stackDetid = tTopo->stack(detid); // Stub module detid

      if (PixelDigiL1TkStubHandle->find( stackDetid ) == PixelDigiL1TkStubHandle->end() ) continue;

      /// Get the DetSets of the Clusters
      edmNew::DetSet< TTStub< Ref_Phase2TrackerDigi_ > > stubs = (*PixelDigiL1TkStubHandle)[ stackDetid ];
      const GeomDetUnit* det0 = theTrackerGeom->idToDetUnit( detid );
      const PixelGeomDetUnit* theGeomDet = dynamic_cast< const PixelGeomDetUnit* >( det0 );
      const PixelTopology* topol = dynamic_cast< const PixelTopology* >( &(theGeomDet->specificTopology()) );

      for ( auto stubIter = stubs.begin();stubIter != stubs.end();++stubIter ) 
      {
	edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_  > >, TTStub< Ref_Phase2TrackerDigi_  > > tempStubPtr = edmNew::makeRefTo( PixelDigiL1TkStubHandle, stubIter );

      
	// The stub position is the inner cluster position
	//

	MeasurementPoint coords = tempStubPtr->getClusterRef(0)->findAverageLocalCoordinatesCentered();

	clustlp = topol->localPosition(coords);
	posStub =  theGeomDet->surface().toGlobal(clustlp);

	++m_stub;

	double displStub    = tempStubPtr->getTriggerDisplacement();
	double offsetStub   = tempStubPtr->getTriggerOffset();
	double offsetRStub  = tempStubPtr->getRealTriggerOffset();

	bool genuineStub    = MCTruthTTStubHandle->isGenuine( tempStubPtr );
	
	segs= topol->ncolumns();
	rows= topol->nrows();;
	
	clust1 = StubExtractor::getClust1Idx(posStub.x(),posStub.y(),posStub.z());
	clust2 = StubExtractor::getClust2Idx(clust1,m_clus_strip->at(clust1)+displStub);
	
	m_stub_x->push_back(posStub.x());
	m_stub_y->push_back(posStub.y());
	m_stub_z->push_back(posStub.z());
	m_stub_clust1->push_back(clust1);
	m_stub_clust2->push_back(clust2);
	m_stub_seg->push_back(m_clus_seg->at(clust1));
	m_stub_chip->push_back(m_clus_strip->at(clust1)/(rows/8));
	m_stub_strip->push_back(m_clus_strip->at(clust1));
	m_stub_detid->push_back(static_cast<int>(tTopo->stack(detid)));
	m_stub_deltas->push_back(displStub-offsetStub);
	m_stub_deltash->push_back(tempStubPtr->getHardwareBend());
	m_stub_cor->push_back(offsetStub);
	m_stub_deltasf->push_back(displStub-offsetRStub);
	m_stub_corf->push_back(offsetRStub);

	m_stub_pt->push_back(0);

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

	m_stub_layer->push_back(layer);
	m_stub_ladder->push_back(ladder);
	m_stub_module->push_back(module);
	m_stub_type->push_back(type);
	m_stub_rank->push_back(0);

	if ( genuineStub )
	{
	  edm::Ptr< TrackingParticle > tpPtr = MCTruthTTStubHandle->findTrackingParticlePtr( tempStubPtr );
	  
	  m_stub_pxGEN->push_back(tpPtr->p4().px());
	  m_stub_pyGEN->push_back(tpPtr->p4().py());
	  m_stub_etaGEN->push_back(tpPtr->momentum().eta());
	  m_stub_pdg->push_back(tpPtr->pdgId());
	  m_stub_pid->push_back(0);
	  m_stub_X0->push_back(tpPtr->vertex().x());
	  m_stub_Y0->push_back(tpPtr->vertex().y());
	  m_stub_Z0->push_back(tpPtr->vertex().z());
	  m_stub_PHI0->push_back(tpPtr->momentum().phi());
	  if (MCinfo) m_stub_tp->push_back(mc->getMatchingTP(tpPtr->vertex().x(),tpPtr->vertex().y(),tpPtr->vertex().z(),
							     tpPtr->p4().px(),tpPtr->p4().py(),tpPtr->p4().pz()));
	}
	else
	{
	  m_stub_pxGEN->push_back(0);
	  m_stub_pyGEN->push_back(0);
	  m_stub_etaGEN->push_back(0);
	  m_stub_X0->push_back(0);
	  m_stub_Y0->push_back(0);
	  m_stub_Z0->push_back(0);
	  m_stub_PHI0->push_back(0);
	  m_stub_pdg->push_back(0);
	  m_stub_pid->push_back(0);
	  m_stub_tp->push_back(-1);
	}
      } /// End of loop over L1TkStubs
    } 
  } /// End of if ( PixelDigiL1TkStubHandle->size() > 0 ) 

  StubExtractor::fillTree();
}


//
// Method getting the info from an input file
//

void StubExtractor::getInfo(int ievt) 
{
  m_tree->GetEntry(ievt); 
}

// Method initializing everything (to do for each event)

void StubExtractor::reset()
{
  m_clus = 0;
  m_stub = 0;

  m_clus_x->clear(); 
  m_clus_y->clear(); 
  m_clus_z->clear(); 
  m_clus_e->clear(); 
  m_clus_layer->clear(); 
  m_clus_module->clear();
  m_clus_ladder->clear();
  m_clus_seg->clear();   
  m_clus_type->clear(); 
  m_clus_strip->clear(); 
  m_clus_sat->clear();   
  m_clus_nstrips->clear();
  m_clus_used->clear();   
  m_clus_matched->clear();
  m_clus_PS->clear();
  m_clus_nrows->clear();
  m_clus_bot->clear();
  m_clus_pid->clear();
  m_clus_pdgID->clear();
  m_clus_ptGEN->clear();
  m_clus_pix->clear(); 
  m_clus_tp->clear(); 

  m_stub_X0->clear();     
  m_stub_Y0->clear();     
  m_stub_Z0->clear();     
  m_stub_PHI0->clear();     
  m_stub_tp->clear();     
  m_stub_pt->clear();     
  m_stub_ptMC->clear();   
  m_stub_pxGEN->clear();  
  m_stub_pyGEN->clear();  
  m_stub_etaGEN->clear();  
  m_stub_layer->clear();  
  m_stub_module->clear(); 
  m_stub_ladder->clear(); 
  m_stub_seg->clear();  
  m_stub_type->clear();
  m_stub_chip->clear();   
  m_stub_strip->clear(); 
  m_stub_detid->clear(); 
  m_stub_x->clear(); 
  m_stub_y->clear(); 
  m_stub_z->clear(); 
  m_stub_clust1->clear(); 
  m_stub_clust2->clear(); 
  m_stub_cw1->clear(); 
  m_stub_cw2->clear(); 
  m_stub_deltas->clear(); 
  m_stub_cor->clear(); 
  m_stub_deltasf->clear(); 
  m_stub_deltash->clear(); 
  m_stub_corf->clear();
  m_stub_pdg->clear();
  m_stub_pid->clear();
  m_stub_rank->clear(); 

}


void StubExtractor::fillTree()
{
  m_tree->Fill(); 
}
 
void StubExtractor::fillSize(int size)
{
  m_clus=size;
}

int  StubExtractor::getSize()
{
  return m_clus;
}

int  StubExtractor::getClust1Idx(float x, float y, float z)
{
  for (int i=0;i<m_clus;++i) // Loop over clusters
  { 
    if (fabs(m_clus_x->at(i)-x) > 0.001) continue;
    if (fabs(m_clus_y->at(i)-y) > 0.001) continue;
    if (fabs(m_clus_z->at(i)-z) > 0.001) continue;

    return i;
  }

  return -1;
}

int  StubExtractor::getClust2Idx(int idx1, float dist)
{
  float dmax = 20;
  int idx2 = -1;

  for (int i=0;i<m_clus;++i) // Loop over clusters
  { 
    if (i==idx1) continue;
    if (m_clus_layer->at(i)!=m_clus_layer->at(idx1)) continue;
    if (m_clus_ladder->at(i)!=m_clus_ladder->at(idx1)) continue;
    if (m_clus_module->at(i)!=m_clus_module->at(idx1)) continue;

    if (fabs(m_clus_strip->at(i)-dist)<dmax)
    {  
      dmax = fabs(m_clus_strip->at(i)-dist);
      idx2 = i;
    }
  }

  return idx2;
}

int  StubExtractor::get_id(int lay,int lad,int mod,float x,float y,float z,float sw)
{
  int idx = -1;

  for (int i=0;i<m_stub;++i) // Loop over stubs
  { 
    if (idx!=-1) return idx;

    if (m_stub_layer->at(i)  !=lay) continue;
    if (m_stub_ladder->at(i) !=lad) continue;
    if (m_stub_module->at(i) !=mod) continue;
    if (m_stub_deltas->at(i)+m_stub_cor->at(i) !=sw) continue;
    if (fabs(m_stub_x->at(i)-x)>0.001) continue;
    if (fabs(m_stub_y->at(i)-y)>0.001) continue;
    if (fabs(m_stub_z->at(i)-z)>0.001) continue;

    idx=i;
  }

  if (idx!=-1) return idx;

  std::cout << "StubExtractor::get_id: you are not supposed to get here!" << std::endl;

  return idx;
}
