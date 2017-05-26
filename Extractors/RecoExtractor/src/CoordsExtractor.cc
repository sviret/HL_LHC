
#include "../interface/CoordsExtractor.h"


CoordsExtractor::CoordsExtractor(bool doTree, bool skim)
{
  if (doTree)
  {
    m_tree      = new TTree("StripCoords","") ;

    m_tree->Branch("layer",     &m_c_layer);
    m_tree->Branch("module",    &m_c_module);
    m_tree->Branch("ladder",    &m_c_ladder);
    m_tree->Branch("modid",     &m_c_modid);
    m_tree->Branch("xm",        &m_c_xm);
    m_tree->Branch("ym",        &m_c_ym);
    m_tree->Branch("zm",        &m_c_zm);
    m_tree->Branch("strip_dx",  &m_c_stdx);
    m_tree->Branch("strip_dy",  &m_c_stdy);
    m_tree->Branch("strip_dz",  &m_c_stdz);
    m_tree->Branch("seg_dx",    &m_c_sedx);
    m_tree->Branch("seg_dy",    &m_c_sedy);
    m_tree->Branch("seg_dz",    &m_c_sedz);

    if (!skim)
    {
      m_tree->Branch("strip",     &m_c_row);
      m_tree->Branch("segment",   &m_c_column);
      m_tree->Branch("x",         &m_c_x);
      m_tree->Branch("y",         &m_c_y);
      m_tree->Branch("z",         &m_c_z);
    }
  }

  m_skim = skim;
}

CoordsExtractor::~CoordsExtractor()
{}


void CoordsExtractor::init(const edm::EventSetup *setup, bool isFlat)
{
  setup->get<TrackerTopologyRcd>().get(tTopoHandle);
  setup->get<TrackerDigiGeometryRecord>().get(tGeomHandle);

  const TrackerTopology* const tTopo = tTopoHandle.product();
  const TrackerGeometry* const theTrackerGeom = tGeomHandle.product();

  bool barrel;
  bool endcap;

  int layer  = 0;
  int ladder = 0;
  int module = 0;
  int type   = 0;
  int segs   = 0;
  int rows   = 0;

  int limits[6][3];

  int n_tilted_rings[6];
  int n_flat_rings[6];

  bool m_tilt=(!isFlat);

  for (int i=0; i < 6; ++i) n_tilted_rings[i]=0;
  for (int i=0; i < 6; ++i) n_flat_rings[i]=0;

  if (m_tilt)
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

  LocalPoint  clustlp,clustlp2;
  GlobalPoint pos,pos2;    
  float phi,r,z;

  float x0,y0,z0,x1,y1,z1;

  for (auto id=theTrackerGeom->dets().begin(); id != theTrackerGeom->dets().end(); id++) 
  { 
    DetId detid = (*id)->geographicalId();
    if(detid.subdetId()==1 || detid.subdetId()==2 ) continue; // only run on OT
    if(!tTopo->isLower(detid)) continue; // loop on the stacks: choose the lower arbitrarily
    DetId lowerDetid = detid;

    barrel = (detid.subdetId()==StripSubdetector::TOB);
    endcap = (detid.subdetId()==StripSubdetector::TID);

    if (!barrel && !endcap) continue;

    /// Get chip size information
    const GeomDetUnit* det0 = theTrackerGeom->idToDetUnit( lowerDetid );
    const PixelGeomDetUnit* theGeomDet = dynamic_cast< const PixelGeomDetUnit* >( det0 );
    const PixelTopology* topol = dynamic_cast< const PixelTopology* >( &(theGeomDet->specificTopology()) );

    segs     = topol->ncolumns();
    rows     = topol->nrows();
    
    r   = det0->position().perp();
    z   = det0->position().z();
    phi = atan2(det0->position().y(),det0->position().x()) ;
    
    m_c_xm = r*cos(phi);
    m_c_ym = r*sin(phi);
    m_c_zm = z;

    std::string binary = std::bitset<32>(static_cast<int>(tTopo->stack(detid))).to_string(); //to binary


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

    m_c_layer   = layer;
    m_c_ladder  = ladder;
    m_c_module  = module;
    m_c_modid   = 10000*m_c_layer+100*m_c_ladder+m_c_module;

    // Get the strip direction

    clustlp  = topol->localPosition( MeasurementPoint(0,0));
    clustlp2 = topol->localPosition( MeasurementPoint(1,0));

    pos     =  theGeomDet->surface().toGlobal(clustlp);
    pos2    =  theGeomDet->surface().toGlobal(clustlp2);

    x0      = pos.x();
    y0      = pos.y();
    z0      = pos.z();
    x1      = pos2.x();
    y1      = pos2.y();
    z1      = pos2.z();

    m_c_stdx = (x1-x0);
    m_c_stdy = (y1-y0);
    m_c_stdz = (z1-z0);

    // Get the segment direction
    clustlp  = topol->localPosition( MeasurementPoint(0,0));
    clustlp2 = topol->localPosition( MeasurementPoint(0,1));

    pos     =  theGeomDet->surface().toGlobal(clustlp);
    pos2    =  theGeomDet->surface().toGlobal(clustlp2);

    x0      = pos.x();
    y0      = pos.y();
    z0      = pos.z();
    x1      = pos2.x();
    y1      = pos2.y();
    z1      = pos2.z();

    m_c_sedx = (x1-x0);
    m_c_sedy = (y1-y0);
    m_c_sedz = (z1-z0);

    std::cout << layer << " / " << ladder << " / " << module;
    std::cout << " / " << (static_cast<float>(rows)-1)/2 << " / " << (static_cast<float>(segs)-1)/2;
    std::cout << " / " << std::setw(7) << std::fixed << std::setprecision(4) << m_c_xm; 
    std::cout << " / " << std::setw(7) << std::fixed << std::setprecision(4) << m_c_ym; 
    std::cout << " / " << std::setw(7) << std::fixed << std::setprecision(3) << m_c_zm; 
    std::cout << " / " << std::setw(7) << std::fixed << std::setprecision(4) << m_c_stdx; 
    std::cout << " / " << std::setw(7) << std::fixed << std::setprecision(4) << m_c_stdy; 
    std::cout << " / " << std::setw(7) << std::fixed << std::setprecision(3) << m_c_stdz; 
    std::cout << " / " << std::setw(7) << std::fixed << std::setprecision(4) << m_c_sedx; 
    std::cout << " / " << std::setw(7) << std::fixed << std::setprecision(4) << m_c_sedy; 
    std::cout << " / " << std::setw(7) << std::fixed << std::setprecision(3) << m_c_sedz << std::endl;
    //    std::cout << " / " << std::setw(7) << std::fixed << m_c_modid;
    //    std::cout << " / " << binary;
    //    std::cout << " / " << static_cast<int>(tTopo->stack(detid)) << std::endl;


    if (!m_skim)
    {
      for (int i=0; i < segs; ++i)
      {
	for (int j=0; j < rows; ++j)
	{
	  clustlp = topol->localPosition( MeasurementPoint(j+0.5,i+0.5));
	  pos     =  theGeomDet->surface().toGlobal(clustlp);

	  m_c_row    = j;
	  m_c_column = i;
	  m_c_x      = pos.x();
	  m_c_y      = pos.y();
	  m_c_z      = pos.z();
	  
	  m_tree->Fill(); 
	}
      }
    }   
    else
    {
      m_tree->Fill(); 
    }
  }
}


