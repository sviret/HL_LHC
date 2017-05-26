#ifndef L1TRACKEXTRACTOR_H
#define L1TRACKEXTRACTOR_H

/**
 * L1TrackExtractor
 * \brief: Base class for extracting pixel info
 */



#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"


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
#include <boost/shared_ptr.hpp>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "StubExtractor.h"

class L1TrackExtractor
{

 public:

  L1TrackExtractor(edm::EDGetTokenT< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > STUB_tag, edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > PATT_tag,  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TC_tag,  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TRCK_tag, bool doTree);
  L1TrackExtractor(TFile *a_file);
  ~L1TrackExtractor();


  void init(const edm::EventSetup *setup, bool isFlat);
  void writeInfo(const edm::Event *event, StubExtractor *stub); 
  void getInfo(int ievt); 

  void reset();
  void fillTree(); 
  void fillSize(int size);
  int  getSize();
  int  n_events() {return m_n_events;}

  bool isOK() {return m_OK;}

  //Getters

 private:
  

  edm::EDGetTokenT< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > m_STUB_tag;
  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > m_PATT_tag;
  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > m_TC_tag;
  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > m_TRCK_tag;

  TTree* m_tree;
  bool m_OK;
  double mMagneticFieldStrength ;
  int m_n_events;

  bool m_tilted;

  int limits[6][3];

  /// Geometry handles etc

  edm::ESHandle<TrackerTopology> tTopoHandle;
  edm::ESHandle<TrackerGeometry> tGeomHandle;

  int n_tot_evt;

  int m_patt;  // Number of patterns in the event

  // Size of the following vectors is m_patt
  std::vector< std::vector<int> > *m_patt_links; // Links to the stubs making the patterns in L1TkStubs tree
  std::vector<int>                *m_patt_secid; // Sector number
  std::vector<int>                *m_patt_pattid;// Pattern rank in the bank
  std::vector<int>                *m_patt_miss;  // Number of missed SStrips


  // Size of the following vectors is m_tc

  int m_tc;  // Number of L1TCs in the event
  
  std::vector<float>               *m_tc_pt;    // pt calculated for track i (in GeV/c) 
  std::vector<float>               *m_tc_eta;   // eta calculated for track i 
  std::vector<float>               *m_tc_phi;   // phi calculated for track i (in rad)
  std::vector<float>               *m_tc_z;     // z0 calculated for track i (in mm)
  std::vector< std::vector<int> >  *m_tc_links; // Links to the stubs making the tracks in L1TkStubs tree
  std::vector<int>                 *m_tc_secid; // Sector number
  std::vector<int>                 *m_tc_pattid;// Pattern leading to the TC rank in the bank


  // Size of the following vectors is m_trk

  int m_trk;  // Number of L1Tracks in the event
  
  std::vector<float>               *m_trk_pt;    // pt calculated for track i (in GeV/c) 
  std::vector<float>               *m_trk_eta;   // eta calculated for track i 
  std::vector<float>               *m_trk_phi;   // phi calculated for track i (in rad)
  std::vector<float>               *m_trk_z;     // z0 calculated for track i (in mm)
  std::vector< std::vector<int> >  *m_trk_links; // Links to the stubs making the tracks in L1TkStubs tree
  std::vector<int>                 *m_trk_secid; // Sector number
  std::vector<int>                 *m_trk_pattid;// Pattern leading to the Track rank in the bank
  std::vector<float>               *m_trk_chi;   // Chi2/dof of the track 

};

#endif 
