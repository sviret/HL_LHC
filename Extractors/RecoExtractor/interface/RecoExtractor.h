#ifndef RecoExtractor_h
#define RecoExtractor_h

/** \class RecoExtractor
 *  \author by sviret
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "DataFormats/Common/interface/DetSetVectorNew.h"

#include "../interface/CoordsExtractor.h"
#include "../interface/PixelExtractor.h"
#include "../interface/StubExtractor.h"
#include "../interface/L1TrackExtractor.h"
#include "../interface/MCExtractor.h"
//#include "../interface/L1TrackTrigger_analysis.h"
#include "../interface/StubTranslator.h"
#include "../interface/AnalysisSettings.h"

#include "TFile.h"
#include "TRFIOFile.h"

class RecoExtractor : public edm::EDAnalyzer{
 public:
  /// Constructor
  RecoExtractor(const edm::ParameterSet& pset);

  /// Destructor
  virtual ~RecoExtractor(){ }
   int nevent;
  int nevent_tot;

  /// Method called before the event loop
  void beginJob();
  void endJob();

  void beginRun(edm::Run const&, edm::EventSetup const&);
  void endRun(edm::Run const&, edm::EventSetup const&);

  /// Method called once per event
  void analyze(const edm::Event&, const edm::EventSetup& );

  void fillInfo(const edm::Event *event);
  void getInfo(int ievent);
  void initialize();
  void retrieve();
  void doAna();
  


 private:
 
  bool do_fill_;


  bool do_COORDS_;
  bool do_PIX_;
  bool do_MC_;
  bool do_STUB_;
  bool do_L1TRK_;
  bool do_BANK_;
  bool do_MATCH_;
  bool do_L1tt_;

  bool use_flat_;
  bool fullinfo_;

  int  nevts_;
  int  skip_;

  edm::EDGetTokenT< edmNew::DetSetVector< TTCluster< Ref_Phase2TrackerDigi_ > > > clustersToken_;
  edm::EDGetTokenT< TTClusterAssociationMap< Ref_Phase2TrackerDigi_ > > clustersTToken_;
  edm::EDGetTokenT< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > stubsToken_;
  edm::EDGetTokenT< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > stubsTToken_;
  edm::EDGetTokenT< reco::GenParticleCollection > gpToken_;
  edm::EDGetTokenT< TrackingParticleCollection > tpToken_;
  edm::EDGetTokenT< edm::DetSetVector< Phase2TrackerDigi> > pixToken_;
  edm::EDGetTokenT< edm::DetSetVector< PixelDigiSimLink> > pixslToken_;
  edm::EDGetTokenT< std::vector<PileupSummaryInfo> > puToken_;
  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > pattToken_;
  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > tcToken_;
  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > trkToken_;

  edm::EDGetTokenT< std::vector< TrackingParticle > > tpToken2_;
  edm::EDGetTokenT< std::vector< TrackingVertex > > tvToken_;
  edm::EDGetTokenT< edm::SimTrackContainer > simtToken_;
  edm::EDGetTokenT< edm::SimVertexContainer > simvToken_;


  //
  // Definition of root-tuple :
  //

  std::string outFilename_;
  std::string inFilename_;

  std::vector<std::string> m_settings_;

  TFile* m_dummyfile;
  TFile* m_infile;
  TFile* m_outfile;

  CoordsExtractor*          m_COORDS;
  PixelExtractor*           m_PIX;
  MCExtractor*              m_MC;
  MCExtractor*              m_dummy_MC;
  StubExtractor*            m_STUB;
  L1TrackExtractor*         m_L1TRK;
  StubTranslator*           m_BK;
  AnalysisSettings*         m_ana_settings;
  // L1TrackTrigger_analysis*  m_L1TT_analysis;

};


#endif
