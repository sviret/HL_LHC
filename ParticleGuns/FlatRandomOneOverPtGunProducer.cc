#include <ostream>

#include "IOMC/ParticleGuns/interface/FlatRandomOneOverPtGunProducer.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CLHEP/Units/GlobalSystemOfUnits.h"
using namespace edm;

FlatRandomOneOverPtGunProducer::FlatRandomOneOverPtGunProducer(const edm::ParameterSet& pset) : 
  BaseFlatGunProducer(pset) {


  edm::ParameterSet defpset ;
  edm::ParameterSet pgun_params = 
    pset.getParameter<ParameterSet>("PGunParameters") ;
  
  fMinOneOverPt = pgun_params.getParameter<double>("MinOneOverPt");
  fMaxOneOverPt = pgun_params.getParameter<double>("MaxOneOverPt");
  fXFlatSpread  = pgun_params.getParameter<double>("XFlatSpread");
  fYFlatSpread  = pgun_params.getParameter<double>("YFlatSpread");
  fZFlatSpread  = pgun_params.getParameter<double>("ZFlatSpread");
  towerID       = pgun_params.getParameter<int>("towerID");

  produces<HepMCProduct>();
  produces<GenEventInfoProduct>();

  edm::LogInfo("ParticleGun") << "FlatRandomOneOverPtGunProducer: initialized with minimum and maximum 1/pt " << fMinOneOverPt << ":" << fMaxOneOverPt;
}

FlatRandomOneOverPtGunProducer::~FlatRandomOneOverPtGunProducer() {
  // no need to cleanup GenEvent memory - done in HepMCProduct
}

void FlatRandomOneOverPtGunProducer::produce(Event &e, const EventSetup& es) {

  LogDebug("ParticleGun") << " FlatRandomOneOverPtGunProducer : Begin New Event Generation"; 

  // event loop (well, another step in it...)
          
  // no need to clean up GenEvent memory - done in HepMCProduct
  // 
   
  // here re-create fEvt (memory)
  //
  fEvt = new HepMC::GenEvent() ;
   
  double t_eta_min[6] = {-2.4,-1.7,-1.1,-0.4,0.4,1.2};
  double t_eta_max[6] = {-1.2,-0.4,0.4,1.1,1.7,2.4};
  double t_phi_min[8] = {-0.5,0.3,1.1,1.9,2.7,-2.9,-2.1,-1.3};
  double t_phi_max[8] = {1.3,2.1,2.9,3.7,4.5,-1.1,-0.3,0.5};
  
  int eta_sec;
  int phi_sec;
  
  // Special treatment to produce events in a given trigger tower
  // according to the 8x6 config defined here:
  //
  // http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCData
  //
  
  if (towerID>=0 && towerID<=47)
  {
    phi_sec=towerID%8;
    eta_sec=(towerID-phi_sec)/8;
    
    fMinEta = t_eta_min[eta_sec];
    fMaxEta = t_eta_max[eta_sec];
    fMinPhi = t_phi_min[phi_sec];
    fMaxPhi = t_phi_max[phi_sec];
  }

  // now actualy, cook up the event from PDGTable and gun parameters
  //
  // 1st, primary vertex
  //
  double xpos     = fRandomGenerator->fire(-fXFlatSpread,fXFlatSpread);
  double ypos     = fRandomGenerator->fire(-fYFlatSpread,fYFlatSpread);
  double zpos     = fRandomGenerator->fire(-fZFlatSpread,fZFlatSpread);

  double charge   = fRandomGenerator->fire(0,1);

  HepMC::GenVertex* Vtx = new HepMC::GenVertex(HepMC::FourVector(xpos,ypos,zpos));

  // loop over particles
  //
  int barcode = 1 ;
  for (unsigned int ip=0; ip<fPartIDs.size(); ++ip) {

    double pt     = fRandomGenerator->fire(fMinOneOverPt,fMaxOneOverPt);
    double eta    = fRandomGenerator->fire(fMinEta, fMaxEta) ;
    double phi    = fRandomGenerator->fire(fMinPhi, fMaxPhi) ;
    if (pt != 0) pt = 1./pt;
    int PartID = fPartIDs[ip] ;

    if (charge<0.5) PartID = -PartID;

    const HepPDT::ParticleData* 
      PData = fPDGTable->particle(HepPDT::ParticleID(abs(PartID))) ;
    double mass   = PData->mass().value() ;
    double theta  = 2.*atan(exp(-eta)) ;
    double mom    = pt/sin(theta) ;
    double px     = pt*cos(phi) ;
    double py     = pt*sin(phi) ;
    double pz     = mom*cos(theta) ;
    double energy2= mom*mom + mass*mass ;
    double energy = sqrt(energy2) ; 
    HepMC::FourVector p(px,py,pz,energy) ;
    HepMC::GenParticle* Part = new HepMC::GenParticle(p,PartID,1);
    Part->suggest_barcode( barcode ) ;
    barcode++ ;
    Vtx->add_particle_out(Part);
    LogDebug("ParticleGun") << "FlatRandomOneOverPtGunProducer: Event generated with pt:eta:phi " << pt << " " << eta << " " << phi << " (" << theta/CLHEP::deg << ":" << phi/CLHEP::deg << ")";

    if ( fAddAntiParticle ) {
      HepMC::FourVector ap(-px,-py,-pz,energy) ;
      int APartID = -PartID ;
      if ( PartID == 22 || PartID == 23 ) {
	APartID = PartID ;
      }	  
      HepMC::GenParticle* APart = new HepMC::GenParticle(ap,APartID,1);
      APart->suggest_barcode( barcode ) ;
      barcode++ ;
      Vtx->add_particle_out(APart) ;
    }

  }

  fEvt->add_vertex(Vtx) ;
  fEvt->set_event_number(e.id().event()) ;
  fEvt->set_signal_process_id(20) ; 
        
  if ( fVerbosity > 0 ) {
    fEvt->print() ;  
  }

  std::auto_ptr<HepMCProduct> BProduct(new HepMCProduct()) ;
  BProduct->addHepMCData( fEvt );
  e.put(BProduct);

  std::auto_ptr<GenEventInfoProduct> genEventInfo(new GenEventInfoProduct(fEvt));
  e.put(genEventInfo);
    
  LogDebug("ParticleGun") << " FlatRandomOneOverPtGunProducer : Event Generation Done ";
}
