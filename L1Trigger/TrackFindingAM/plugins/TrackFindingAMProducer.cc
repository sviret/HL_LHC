/*! \class   TrackFindingAMProducer
 *
 *  \author Nicola Pozzobon
 *  \date   2013, Jul 18
 *
 *  \update by S.Viret 
 *  \date   2014, Feb 17
 *
 */

#ifndef TRACK_BUILDER_AM_H
#define TRACK_BUILDER_AM_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerDetUnit.h"

#include "L1Trigger/TrackFindingAM/interface/CMSPatternLayer.h"
#include "L1Trigger/TrackFindingAM/interface/PatternFinder.h"
#include "L1Trigger/TrackFindingAM/interface/SectorTree.h"
#include "L1Trigger/TrackFindingAM/interface/Hit.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/shared_ptr.hpp>
#include <memory>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>

#ifndef __APPLE__
BOOST_CLASS_EXPORT_IMPLEMENT(CMSPatternLayer)
#endif

class TrackFindingAMProducer : public edm::EDProducer
{
  public:
    /// Constructor
    explicit TrackFindingAMProducer( const edm::ParameterSet& iConfig );

    /// Destructor;
    ~TrackFindingAMProducer();

  private:

  /// Data members
  double                       mMagneticField;
  std::string                  nBKName;
  int                          nThresh;
  SectorTree                   m_st;
  PatternFinder                *m_pf;
  const StackedTrackerGeometry *theStackedTracker;
  edm::InputTag                TTStubsInputTag;
  edm::InputTag                TTClustersInputTag;
  std::string                  TTPatternOutputTag;
  std::string                  TTStubOutputTag;
  std::string                  TTClusOutputTag;
  std::vector<int>             stored_IDs;

  /// Mandatory methods
  virtual void beginRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void endRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void produce( edm::Event& iEvent, const edm::EventSetup& iSetup );

  bool inPattern(int j);

}; /// Close class

/*! \brief   Implementation of methods
 */

/// Constructors
TrackFindingAMProducer::TrackFindingAMProducer( const edm::ParameterSet& iConfig )
{
  TTStubsInputTag    = iConfig.getParameter< edm::InputTag >( "TTInputStubs" );
  TTClustersInputTag = iConfig.getParameter< edm::InputTag >( "TTInputClusters" );
  TTPatternOutputTag = iConfig.getParameter< std::string >( "TTPatternName" );
  TTStubOutputTag    = iConfig.getParameter< std::string >( "TTFiltStubsName" );
  TTClusOutputTag    = iConfig.getParameter< std::string >( "TTFiltClustersName" );
  nBKName            = iConfig.getParameter< std::string >("inputBankFile");
  nThresh            = iConfig.getParameter< int >("threshold");

  std::cout << "Loading pattern bank file : " << std::endl;
  std::cout << nBKName << std::endl;

  std::ifstream ifs(nBKName.c_str());
  boost::archive::text_iarchive ia(ifs);

  ia >> m_st;
  m_pf = new PatternFinder( m_st.getSuperStripSize(), nThresh, &m_st, "", "" );

  produces< std::vector< TTTrack< Ref_PixelDigi_ > > >( TTPatternOutputTag );
  produces<  edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > > >( TTStubOutputTag );
  produces<  edmNew::DetSetVector< TTCluster< Ref_PixelDigi_ > > >( TTClusOutputTag );
}

/// Destructor
TrackFindingAMProducer::~TrackFindingAMProducer() {}

/// Begin run
void TrackFindingAMProducer::beginRun( const edm::Run& run, const edm::EventSetup& iSetup )
{
  /// Get the geometry references
  edm::ESHandle< StackedTrackerGeometry > StackedTrackerGeomHandle;
  iSetup.get< StackedTrackerGeometryRecord >().get( StackedTrackerGeomHandle );
  theStackedTracker = StackedTrackerGeomHandle.product();

  /// Get magnetic field
  edm::ESHandle<MagneticField> magneticFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticFieldHandle);
  const MagneticField* theMagneticField = magneticFieldHandle.product();
  double mMagneticFieldStrength = theMagneticField->inTesla(GlobalPoint(0,0,0)).z();
  mMagneticField = (floor(mMagneticFieldStrength*10.0 + 0.5))/10.0;
}

/// End run
void TrackFindingAMProducer::endRun( const edm::Run& run, const edm::EventSetup& iSetup ) {}

/// Implement the producer
void TrackFindingAMProducer::produce( edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  /// Prepare output
  /// The temporary collection is used to store tracks
  /// before removal of duplicates
  std::auto_ptr< std::vector< TTTrack< Ref_PixelDigi_ > > > TTTracksForOutput( new std::vector< TTTrack< Ref_PixelDigi_ > > );

  // The container for filtered stubs / clusters
  std::auto_ptr< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > > > TTStubsForOutputAccepted( new edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > > );
  std::auto_ptr< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > > > TTStubsForOutput( new edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > > );
  std::auto_ptr< edmNew::DetSetVector< TTCluster< Ref_PixelDigi_ > > > TTClustsForOutput( new edmNew::DetSetVector< TTCluster< Ref_PixelDigi_ > > );


  /// Get the Stubs/Cluster already stored
  edm::Handle< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > > > TTStubHandle;
  iEvent.getByLabel( TTStubsInputTag, TTStubHandle );

  edm::Handle< edmNew::DetSetVector< TTCluster< Ref_PixelDigi_ > > > TTClusterHandle;
  iEvent.getByLabel( TTClustersInputTag, TTClusterHandle );

  /// STEP 0
  /// Prepare output

  TTTracksForOutput->clear();

  int layer  = 0;
  int ladder = 0;
  int module = 0;
  int n_active =  m_st.getAllSectors().at(0)->getNbLayers();

  /// STEP 1
  /// Loop over input stubs

  std::vector< Hit* > m_hits;
  for(unsigned int i=0;i<m_hits.size();i++) delete m_hits[i];
  m_hits.clear();

  unsigned int j = 0;
  std::map< unsigned int, edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > > stubMap;

  edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >::const_iterator inputIter;
  edmNew::DetSetVector< TTCluster< Ref_PixelDigi_ > >::const_iterator inputclusIter;
  edmNew::DetSet< TTStub< Ref_PixelDigi_ > >::const_iterator stubIter;
  edmNew::DetSet< TTCluster< Ref_PixelDigi_ > >::iterator clusterIter;

  for ( inputIter = TTStubHandle->begin(); inputIter != TTStubHandle->end(); ++inputIter )
  {
    for ( stubIter = inputIter->begin(); stubIter != inputIter->end(); ++stubIter )
    {
      /// Increment the counter
      j++;

      /// Make the Ref to be put in the Track
      edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > tempStubRef = makeRefTo( TTStubHandle, stubIter );

      stubMap.insert( std::make_pair( j, tempStubRef ) );

      /// Calculate average coordinates col/row for inner/outer Cluster
      /// These are already corrected for being at the center of each pixel
      MeasurementPoint mp0 = tempStubRef->getClusterRef(0)->findAverageLocalCoordinates();
      GlobalPoint posStub  = theStackedTracker->findGlobalPosition( &(*tempStubRef) );

      StackedTrackerDetId detIdStub( tempStubRef->getDetId() );

      const GeomDetUnit* det0 = theStackedTracker->idToDetUnit( detIdStub, 0 );
      const GeomDetUnit* det1 = theStackedTracker->idToDetUnit( detIdStub, 1 );

      /// Find pixel pitch and topology related information
      const PixelGeomDetUnit* pix0 = dynamic_cast< const PixelGeomDetUnit* >( det0 );
      const PixelGeomDetUnit* pix1 = dynamic_cast< const PixelGeomDetUnit* >( det1 );
      const PixelTopology* top0    = dynamic_cast< const PixelTopology* >( &(pix0->specificTopology()) );
      const PixelTopology* top1    = dynamic_cast< const PixelTopology* >( &(pix1->specificTopology()) );

      /// Find the z-segment
      int cols0   = top0->ncolumns();
      int cols1   = top1->ncolumns();
      int ratio   = cols0/cols1; /// This assumes the ratio is integer!
      int segment = floor( mp0.y() / ratio );

      // Here we rearrange the number in order to be compatible with the AM emulator
      if ( detIdStub.isBarrel() )
      {
        layer  = detIdStub.iLayer()+4;
        ladder = detIdStub.iPhi()-1;
        module = detIdStub.iZ()-1;
      }
      else if ( detIdStub.isEndcap() )
      {
        layer  = 10+detIdStub.iZ()+abs((int)(detIdStub.iSide())-2)*7;
        ladder = detIdStub.iRing()-1;
        module = detIdStub.iPhi()-1;
      }

      module = CMSPatternLayer::getModuleCode(layer,module);

      // the stub is on the third Z position on the other side of the tracker -> out of range
      if ( module < 0 )  continue;

      ladder = CMSPatternLayer::getLadderCode(layer, ladder);

      float x    = posStub.x();
      float y    = posStub.y();
      float z    = posStub.z();

      Hit* h = new Hit(layer,ladder, module, segment, mp0.x(), j, -1, 0, 0, 0, 0, x, y, z, 0, 0, 0);
      m_hits.push_back(h);
    } /// End of loop over input stubs
  } /// End of loop over DetSetVector

  /// STEP 2
  /// PAssing the superstrips into the AM chip

  std::vector< Sector* > patternsSectors = m_pf->find(m_hits); // AM PR is done here....

  /// STEP 3
  /// Collect the info and store the track seed stuff

  std::vector< Hit* > hits;
  bool found;

  stored_IDs.clear();
  std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > > tempVec;

  for ( unsigned int i = 0; i < patternsSectors.size(); i++ )
  {
    std::vector< GradedPattern* > pl = patternsSectors[i]->getPatternTree()->getLDPatterns();

    if ( pl.size() == 0 ) continue; // No patterns

    int secID = patternsSectors[i]->getOfficialID();

    //    std::cout<<"Found "<<pl.size()<<" patterns in sector " << secID<<std::endl;
    //    std::cout<<"containing "<<n_active<<" layers " << secID<<std::endl;
  
    //delete the GradedPattern objects

    for ( unsigned j = 0; j < pl.size(); j++ )
    {
      hits.clear();
      hits = pl[j]->getHits();

      /// Create the Seed in the form of a Track and store it in the output
      tempVec.clear();

      for(unsigned k = 0; k < hits.size(); k++ )
      {
        tempVec.push_back( stubMap[ hits[k]->getID() ] );

	found = false;

	for(unsigned l = 0; l < stored_IDs.size(); ++l )
	{
	  if (found) continue;
	  if (stored_IDs.at(l)==hits[k]->getID()) found=true;
	}

	if (!found) stored_IDs.push_back(hits[k]->getID());
      }

      TTTrack< Ref_PixelDigi_ > tempTrack( tempVec );
      tempTrack.setSector( secID );
      tempTrack.setWedge( n_active );
      TTTracksForOutput->push_back( tempTrack );

      delete pl[j];
    }

    //delete the Sectors
    delete patternsSectors[i];
  }


  /// STEP 4
  /// Fill the DetSetVector of filtered stubs/clusters
  
  unsigned int j2 = 0;

  for ( inputIter = TTStubHandle->begin(); inputIter != TTStubHandle->end(); ++inputIter )
  {
    DetId thisStackedDetId = inputIter->id();

    /// Get its DetUnit

    const StackedTrackerDetUnit* thisUnit = theStackedTracker->idToStack( thisStackedDetId );

    DetId id0 = thisUnit->stackMember(0);
    DetId id1 = thisUnit->stackMember(1);
 
    /// Check that everything is ok in the maps
    if ( theStackedTracker->findPairedDetector( id0 ) != id1 ||
         theStackedTracker->findPairedDetector( id1 ) != id0 )
    {
      std::cerr << "A L E R T! error in detector association within Pt module (detector-to-detector)" << std::endl;
      continue;
    }

    if ( theStackedTracker->findStackFromDetector( id0 ) != thisStackedDetId ||
         theStackedTracker->findStackFromDetector( id1 ) != thisStackedDetId )
    {
      std::cerr << "A L E R T! error in detector association within Pt module (detector-to-module)" << std::endl;
      continue;
    }

    /// Go on only if both detectors have Clusters
    if ( TTClusterHandle->find( id0 ) == TTClusterHandle->end() ||
         TTClusterHandle->find( id1 ) == TTClusterHandle->end() )
      continue;

    /// Get the DetSets of the Clusters
    // edmNew::DetSet< TTCluster< T > > innerClusters = (*TTClusterHandle)[ id0 ];
    // edmNew::DetSet< TTCluster< T > > outerClusters = (*TTClusterHandle)[ id1 ];

    /// Create the vectors of objects to be passed to the FastFillers
    std::vector< TTCluster< Ref_PixelDigi_ > > *tempInner = new std::vector< TTCluster< Ref_PixelDigi_ > >();
    std::vector< TTCluster< Ref_PixelDigi_ > > *tempOuter = new std::vector< TTCluster< Ref_PixelDigi_ > >();
    std::vector< TTStub< Ref_PixelDigi_ > > *tempOutput = new std::vector< TTStub< Ref_PixelDigi_ > >();
    tempInner->clear();
    tempOuter->clear();
    tempOutput->clear();

    for ( stubIter = inputIter->begin(); stubIter != inputIter->end(); ++stubIter )
    {
      ++j2;

      if (!TrackFindingAMProducer::inPattern(j2)) continue;

      TTStub< Ref_PixelDigi_ > tempTTStub( stubIter->getDetId() );

      tempTTStub.addClusterRef(stubIter->getClusterRef(0));
      tempTTStub.addClusterRef(stubIter->getClusterRef(1));
      tempTTStub.setTriggerDisplacement( stubIter->getTriggerDisplacement() );
      tempTTStub.setTriggerOffset( stubIter->getTriggerOffset() );

      tempOutput->push_back( tempTTStub );

      //Get the clusters 
      edm::Ref< edmNew::DetSetVector< TTCluster< Ref_PixelDigi_ > >, TTCluster< Ref_PixelDigi_ > > 
      	clus_i = stubIter->getClusterRef(0);
      edm::Ref< edmNew::DetSetVector< TTCluster< Ref_PixelDigi_ > >, TTCluster< Ref_PixelDigi_ > > 
      	clus_o = stubIter->getClusterRef(1);

      TTCluster< Ref_PixelDigi_ > tempTTInnerClus(clus_i->getHits(),clus_i->getDetId(),clus_i->getStackMember(),true);
      TTCluster< Ref_PixelDigi_ > tempTTOuterClus(clus_o->getHits(),clus_o->getDetId(),clus_o->getStackMember(),true);

      tempInner->push_back( tempTTInnerClus );
      tempOuter->push_back( tempTTOuterClus );
    }


    /// Create the FastFillers
    if ( tempInner->size() > 0 )
    {
      typename edmNew::DetSetVector< TTCluster< Ref_PixelDigi_ > >::FastFiller innerOutputFiller( *TTClustsForOutput, id0 );
      for ( unsigned int m = 0; m < tempInner->size(); m++ )
      {
        innerOutputFiller.push_back( tempInner->at(m) );
      }
      if ( innerOutputFiller.empty() )
        innerOutputFiller.abort();
    }

    if ( tempOuter->size() > 0 )
    {
      typename edmNew::DetSetVector< TTCluster< Ref_PixelDigi_ > >::FastFiller outerOutputFiller( *TTClustsForOutput, id1 );
      for ( unsigned int m = 0; m < tempOuter->size(); m++ )
      {
        outerOutputFiller.push_back( tempOuter->at(m) );
      }
      if ( outerOutputFiller.empty() )
        outerOutputFiller.abort();
    }

    if ( tempOutput->size() > 0 )
    {
      typename edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >::FastFiller tempOutputFiller( *TTStubsForOutput, thisStackedDetId );
      for ( unsigned int m = 0; m < tempOutput->size(); m++ )
      {
        tempOutputFiller.push_back( tempOutput->at(m) );
      }
      if ( tempOutputFiller.empty() )
        tempOutputFiller.abort();
    }
  }

  /// Get also the OrphanHandle of the accepted clusters
  edm::OrphanHandle< edmNew::DetSetVector< TTCluster< Ref_PixelDigi_ > > > TTClusterAcceptedHandle = iEvent.put( TTClustsForOutput, TTClusOutputTag );

  // Do a second loop in order to link Cluster and Stub collections

  /// Now, correctly reset the output

  for ( inputIter = TTStubsForOutput->begin();
        inputIter != TTStubsForOutput->end();
        ++inputIter )
  {
    /// Get the DetId and prepare the FastFiller
    DetId thisStackedDetId = inputIter->id();
    edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >::FastFiller acceptedOutputFiller( *TTStubsForOutputAccepted, thisStackedDetId );
    
    /// Get its DetUnit
    const StackedTrackerDetUnit* thisUnit = theStackedTracker->idToStack( thisStackedDetId );
    DetId id0 = thisUnit->stackMember(0);
    DetId id1 = thisUnit->stackMember(1);

    /// Check that everything is ok in the maps
    /// Redundant up to (*)
    if ( theStackedTracker->findPairedDetector( id0 ) != id1 ||
         theStackedTracker->findPairedDetector( id1 ) != id0 )
    {
      std::cerr << "A L E R T! error in detector association within Pt module (detector-to-detector)" << std::endl;
      continue;
    }

    if ( theStackedTracker->findStackFromDetector( id0 ) != thisStackedDetId ||
         theStackedTracker->findStackFromDetector( id1 ) != thisStackedDetId )
    {
      std::cerr << "A L E R T! error in detector association within Pt module (detector-to-module)" << std::endl;
      continue;
    }

    /// Go on only if both detectors have clusters
    if ( TTClusterAcceptedHandle->find( id0 ) == TTClusterAcceptedHandle->end() ||
         TTClusterAcceptedHandle->find( id1 ) == TTClusterAcceptedHandle->end() )
      continue;
    

    /// 
    
    /// Get the DetSets of the clusters
    edmNew::DetSet< TTCluster< Ref_PixelDigi_ > > innerClusters = (*TTClusterAcceptedHandle)[ id0 ];
    edmNew::DetSet< TTCluster< Ref_PixelDigi_ > > outerClusters = (*TTClusterAcceptedHandle)[ id1 ];

    /// Get the DetSet of the stubs
    edmNew::DetSet< TTStub< Ref_PixelDigi_ > > theseStubs = (*TTStubsForOutput)[ thisStackedDetId ];

    /// Prepare the new DetSet to replace the current one
    /// Loop over the stubs
    edmNew::DetSet< TTCluster< Ref_PixelDigi_ > >::iterator clusterIter;
    edmNew::DetSet< TTStub< Ref_PixelDigi_ > >::iterator stubIter;  

    for ( stubIter = theseStubs.begin();
          stubIter != theseStubs.end();
          ++stubIter )
    {
      /// Create a temporary stub
      TTStub< Ref_PixelDigi_ > tempTTStub( stubIter->getDetId() );

      /// Compare the clusters stored in the stub with the ones of this module
      edm::Ref< edmNew::DetSetVector< TTCluster< Ref_PixelDigi_ > >, TTCluster< Ref_PixelDigi_ > > innerClusterToBeReplaced = stubIter->getClusterRef(0);
      edm::Ref< edmNew::DetSetVector< TTCluster< Ref_PixelDigi_ > >, TTCluster< Ref_PixelDigi_ > > outerClusterToBeReplaced = stubIter->getClusterRef(1);

      bool innerOK = false;
      bool outerOK = false;
      
      for ( clusterIter = innerClusters.begin();
            clusterIter != innerClusters.end() && !innerOK;
            ++clusterIter )
      {
	
        if ( clusterIter->getHits() == innerClusterToBeReplaced->getHits() )
        {
          tempTTStub.addClusterRef( edmNew::makeRefTo( TTClusterAcceptedHandle, clusterIter ) );
          innerOK = true;
        }
       
      }

      for ( clusterIter = outerClusters.begin();
            clusterIter != outerClusters.end() && !outerOK;
            ++clusterIter )
      {
	
        if ( clusterIter->getHits() == outerClusterToBeReplaced->getHits() )
        {
          tempTTStub.addClusterRef( edmNew::makeRefTo( TTClusterAcceptedHandle, clusterIter ) );
          outerOK = true;
        }
      }
      
      /// If no compatible clusters were found, skip to the next one
      if ( !innerOK || !outerOK )
        continue;

      tempTTStub.setTriggerDisplacement( stubIter->getTriggerDisplacement() );
      tempTTStub.setTriggerOffset( stubIter->getTriggerOffset() );

      acceptedOutputFiller.push_back( tempTTStub );

    } /// End of loop over stubs of this module

    if ( acceptedOutputFiller.empty() )
      acceptedOutputFiller.abort();
  
  } /// End of loop over stub DetSetVector
  
  /// Put in the event content
  iEvent.put( TTTracksForOutput, TTPatternOutputTag);
  iEvent.put( TTStubsForOutputAccepted, TTStubOutputTag);

}


bool TrackFindingAMProducer::inPattern(int j)
{
  for(unsigned l = 0; l < stored_IDs.size(); ++l )
  {    
    if (stored_IDs.at(l)==j) return true;
  }

  return false;
}

// DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(TrackFindingAMProducer);

#endif

