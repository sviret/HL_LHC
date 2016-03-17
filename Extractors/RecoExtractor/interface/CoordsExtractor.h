#ifndef COORDSEXTRACTOR_H
#define COORDSEXTRACTOR_H

/**
 * PixelExtractor
 * \brief: Base class for extracting tracker coordinates info
 */


//Include RECO inf
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"

//Include std C++
#include <iostream>
#include <vector>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

class CoordsExtractor
{

 public:

  CoordsExtractor(bool doTree);
  ~CoordsExtractor();

  void init(const edm::EventSetup *setup);

 private:
  
  TTree* m_tree;

  edm::ESHandle<TrackerGeometry> theTrackerGeometry;

  int    m_c_layer;
  int    m_c_module;
  int    m_c_ladder;
  int    m_c_row;
  int    m_c_column;
  float  m_c_x;
  float  m_c_y;
  float  m_c_z;
};

#endif 
