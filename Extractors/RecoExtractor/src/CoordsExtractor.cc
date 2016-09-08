
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

  int cols;    
  int rows;    

  int limits[6][3];

  int n_tilted_rings[6];
  int n_flat_rings[6];

  bool m_tilt=(!isFlat);

  if (!m_tilt)
  {
    for (int i=0; i < 6; ++i) n_tilted_rings[i]=0;
    for (int i=0; i < 6; ++i) n_flat_rings[i]=0;
  }
  else
  {
    n_tilted_rings[0]=11;
    n_tilted_rings[1]=12;
    n_tilted_rings[2]=13;
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

  LocalPoint  clustlp;
  GlobalPoint pos;    
  float phi,r,z;

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

    cols     = topol->ncolumns();
    rows     = topol->nrows();
    
    r   = det0->position().perp();
    z   = det0->position().z();
    phi = atan2(det0->position().y(),det0->position().x()) ;
    
    std::string binary = std::bitset<32>(static_cast<int>(tTopo->stack(detid))).to_string(); //to binary

    if (barrel) // barrel module
    {
      m_c_layer   = static_cast<int>(tTopo->layer(detid))+4;
      m_c_ladder  = static_cast<int>(tTopo->tobRod(detid))-1;
      m_c_module  = static_cast<int>(tTopo->module(detid))-1;

      if (tTopo->tobSide(detid)==3) 
      {	
	std::cout << "Barrel FLAT   / ";
	m_c_modid = 10000*m_c_layer+100*m_c_ladder+limits[m_c_layer-5][tTopo->tobSide(detid)-1]+m_c_module;
      }
      else
      {
	std::cout << "Barrel TILTED / ";
	m_c_modid = 10000*m_c_layer+100*m_c_module+limits[m_c_layer-5][tTopo->tobSide(detid)-1]+m_c_ladder;
      }
    }
    else //Disk
    {
      m_c_layer  = 10+static_cast<int>(tTopo->tidWheel(detid))+abs(2-static_cast<int>(tTopo->side(detid)))*7;
      m_c_ladder = static_cast<int>(tTopo->tidRing(detid))-1;
      m_c_module = static_cast<int>(tTopo->module(detid))-1;

      std::cout << "Disk          / ";

      m_c_modid = 10000*m_c_layer+100*m_c_ladder+m_c_module;
    }

    std::cout << std::setw(7) << std::fixed << std::setprecision(2) << r; 
    std::cout << " / " << std::setw(7) << std::fixed << std::setprecision(2) << z; 
    std::cout << " / " << std::setw(7) << std::fixed << std::setprecision(2) << phi; 
    std::cout << " / " << std::setw(7) << std::fixed << m_c_modid;
    std::cout << " / " << binary << std::endl;

    if (!m_skim)
    {
      for (int i=0; i < cols; ++i)
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


