
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
  setup->get<TrackerDigiGeometryRecord>().get(theTrackerGeometry);

  bool barrel;
  bool endcap;

  int disk;  

  int cols;    
  int rows;    

  LocalPoint  clustlp;
  GlobalPoint pos;    

  typedef std::vector<DetId>                 DetIdContainer;

  DetIdContainer allIds = theTrackerGeometry->detIds();

  for( DetIdContainer::const_iterator id = allIds.begin(), detUnitIdEnd = allIds.end(); id != detUnitIdEnd; ++id ) 
  {
    barrel = (id->subdetId() == static_cast<int>(PixelSubdetector::PixelBarrel));
    endcap = (id->subdetId() == static_cast<int>(PixelSubdetector::PixelEndcap));

    if (!barrel && !endcap) continue;
      
    const PixelGeomDetUnit* theGeomDet = dynamic_cast<const PixelGeomDetUnit*>(theTrackerGeometry->idToDet(*id));      
    const PixelTopology* topol         = &(theGeomDet->specificTopology());

    cols     = topol->ncolumns();
    rows     = topol->nrows();
            
    if (barrel)
    {
      PXBDetId bdetid(*id);
      
      m_c_layer   = static_cast<int>(bdetid.layer());

      if (m_c_layer<=4) continue;

      m_c_module  = static_cast<int>(bdetid.module());

      if (m_c_module%2==0) continue;

      m_c_module = (m_c_module-1)/2;

      m_c_ladder  = static_cast<int>(bdetid.ladder());
 
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
      PXFDetId fdetid(*id);
      
      // Disk 1 to 10 are Pixels, 11 to 15 Tracker

      disk = (static_cast<int>(fdetid.side())*2-3)*static_cast<int>(fdetid.disk());

      if (disk<11 && disk>-11) continue; 

      if (disk>=11)  m_c_layer = disk; 
      if (disk<=-11) m_c_layer = 7-disk; 

      m_c_module = static_cast<int>(fdetid.module()); 
  
      if (m_c_module%2==0) continue;

      m_c_module = (m_c_module-1)/2;
      m_c_ladder = static_cast<int>(fdetid.ring()); 

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


