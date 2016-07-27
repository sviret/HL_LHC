
#include "../interface/CoordsExtractor.h"


CoordsExtractor::CoordsExtractor(bool doTree)
{
  if (doTree)
  {
    m_tree      = new TTree("StripCoords","") ;

    m_tree->Branch("layer",     &m_c_layer);
    m_tree->Branch("module",    &m_c_module);
    m_tree->Branch("ladder",    &m_c_ladder);
    m_tree->Branch("strip",     &m_c_row);
    m_tree->Branch("segment",   &m_c_column);
    m_tree->Branch("x",         &m_c_x);
    m_tree->Branch("y",         &m_c_y);
    m_tree->Branch("z",         &m_c_z);
  }
}

CoordsExtractor::~CoordsExtractor()
{}


void CoordsExtractor::init(const edm::EventSetup *setup)
{
  setup->get<TrackerTopologyRcd>().get(tTopoHandle);
  setup->get<TrackerDigiGeometryRecord>().get(tGeomHandle);

  const TrackerTopology* const tTopo = tTopoHandle.product();
  const TrackerGeometry* const theTrackerGeom = tGeomHandle.product();

  bool barrel;
  bool endcap;

  int disk;  

  int cols;    
  int rows;    

  int limits[3][3];

  int n_tilted_rings[3];
  int n_flat_rings[3];

  bool m_tilt=true;

  if (!m_tilt)
  {
    for (int i=0; i < 3; ++i) n_tilted_rings[i]=0;
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

  for (int i=0; i < 3; ++i)
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

  int tilt, layer, module, ladder;

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

    int modid;
    
    std::string binary = std::bitset<32>(static_cast<int>(tTopo->stack(detid))).to_string(); //to binary

    // DetID of the module containing the stub
    std::bitset<32> MODID = static_cast<int>(tTopo->stack(detid));

    if (barrel) // barrel module
    {
      tilt  = tTopo->tobSide(detid);
      layer = tTopo->tobLayer(detid);

      if (tilt==3) std::cout << "Barrel FLAT   / ";
      if (tilt<=2) std::cout << "Barrel TILTED / ";

      if (layer>3) // 2S barrel module
      {
	modid = 10000*(layer+4)+100*(static_cast<int>(tTopo->tobRod(detid))-1)+static_cast<int>(tTopo->module(detid));
      }
      else // PS barrel module
      {
	modid = 10000*(layer+4);

	if (tilt==3) modid += 100*(static_cast<int>(tTopo->tobRod(detid))-1)
		       + limits[layer-1][tilt-1]+static_cast<int>(tTopo->module(detid)); // Flat 

	if (tilt<=2) modid += 100*(static_cast<int>(tTopo->module(detid))-1)
		       + limits[layer-1][tilt-1]+static_cast<int>(tTopo->tobRod(detid)); // Tilted
      }
    }
    else //Disk
    {
      //      std::cout << tTopo->print(detid) << std::endl;

      std::cout << "Disk          / ";
      layer= 4*MODID[20]+2*MODID[19]+MODID[18];

      if (MODID[24]) layer+=10;
      if (MODID[23]) layer+=17;

      modid = 10000*(layer)+100*(static_cast<int>(tTopo->tidRing(detid))-1)+static_cast<int>(tTopo->module(detid));
    }


    std::cout << std::setw(7) << std::fixed << std::setprecision(2) << r; 
    std::cout << " / " << std::setw(7) << std::fixed << std::setprecision(2) << z; 
    std::cout << " / " << std::setw(7) << std::fixed << std::setprecision(2) << phi; 
    std::cout << " / " << std::setw(7) << std::fixed << modid-1;
    std::cout << " / " << binary << std::endl;

    if (barrel)
    {
      m_c_layer   = static_cast<int>(tTopo->layer(detid));
      m_c_module  = static_cast<int>(tTopo->module(detid));
      m_c_ladder  = static_cast<int>(tTopo->tobRod(detid));
 
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
    
    if (endcap)
    {
      disk = (static_cast<int>(tTopo->side(detid))*2-3)*static_cast<int>(tTopo->tidWheel(detid));

      if (disk>5 || disk<-5) continue; 
      m_c_layer = disk; 
      m_c_module = static_cast<int>(tTopo->module(detid)); 
      m_c_ladder = static_cast<int>(tTopo->tidRing(detid)); 

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
  }
}


