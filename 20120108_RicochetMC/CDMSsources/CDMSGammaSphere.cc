// $Id: CDMSGammaSphere.cc,v 1.7 2011/07/07 05:04:47 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//  File:        CDMSGammaSphere.cc                                   //
//  Description: Inward spherical source for detector response        //
//                                                                    //
//  Author:      Michael Kelsey                                       //
//  Date:        23 May 2011                                          //
//                                                                    //
//  20110524  M. Kelsey -- Add configuration reporting, messenger     //
//  20110705  M. Kelsey -- Use CDMSMultiGenerator for production. Use //
//		new (unreleased) RadioactiveDecay hook to collimate   //
//		daughters, if CPP flag is set.                        //
//  20110707  M. Kelsey -- Drop AddXXX() functions; no longer needed  //
////////////////////////////////////////////////////////////////////////

#include "CDMSsources/CDMSGammaSphere.hh"
#include "CDMSsources/CDMSGammaSphereMessenger.hh"
#include "G4Gamma.hh"
#include "G4GenericIon.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4RadioactiveDecay.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include <iostream>


// Constructor and destructor

CDMSGammaSphere::CDMSGammaSphere()
  : CDMSVSourceConstruction("GammaSphere"), generator("GammaSphere"),
    inwardGammas(true), radius(1.*meter), halfAngle(0.*deg),
    particlesPerEvent(1), messenger(new CDMSGammaSphereMessenger(this)) {}

CDMSGammaSphere::~CDMSGammaSphere() {
  delete messenger;
}


// Generate user requested gammas from among configured sources

void CDMSGammaSphere::GeneratePrimaries(G4Event* event) {
  generator.SetVerboseLevel(verboseLevel);	// Pass on verbosity

  if (verboseLevel) G4cout << "CDMSGammaSphere::GeneratePrimaries" << G4endl;

  if (generator.GetNSources() == 0) return;	// No decays?  No work to do

  for (int i=0; i<particlesPerEvent; i++) GeneratePrimary(event);
}


// Generate single gamma or radioisotope vertex

void CDMSGammaSphere::GeneratePrimary(G4Event* event) {
  if (verboseLevel>1) G4cout << "CDMSGammaSphere::GeneratePrimary" << G4endl;

  setGunPosition();	// Must choose vertex first; used for radial vector

  generator.SetDirection(getGunDirection(), halfAngle);
  
  G4ParticleDefinition* particle;	// Buffers for generated event
  G4ThreeVector dir;
  G4double ekin;
  generator.shoot(ekin, dir, particle);

#ifdef CDMS_RDM_COLLIMATION
  // Configure collimation for radioactive decay products
  if (dynamic_cast<G4GenericIon*>(particle)) {
    static G4RadioactiveDecay* rdmProcess = findRadioactiveDecayProcess();
    if (rdmProcess) {
      rdmProcess->SetDecayCollimation(generator.GetDirection(),
				      generator.GetHalfAngle());
      if (verboseLevel>2)
	G4cout << " decay daughters collimated to "
	       << generator.GetDirection() << G4endl;
    }
  }
#endif	/* CDMS_RDM_COLLIMATION */

  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleEnergy(ekin);
  particleGun->SetParticleMomentumDirection(dir);

  // Add gamma production or radioactive decay to event
  particleGun->GeneratePrimaryVertex(event);
}


// Generate random position on source sphere

void CDMSGammaSphere::setGunPosition() {
  if (verboseLevel>1) G4cout << "CDMSGammaSphere::setGunPosition" << G4endl;

  G4double phi = 2.*pi*G4UniformRand();
  G4double theta = std::acos(2.*G4UniformRand() - 1.);

  if (verboseLevel>2)
    G4cout << " choose R " << radius << " phi " << phi << " theta "
	   << theta << G4endl;

  G4ThreeVector pos;
  pos.setRThetaPhi(radius, theta, phi);
  pos += GetPosition();	        // Apply offset of sphere from origin

  particleGun->SetParticlePosition(pos);
}


// Retrieve source position and use to build radial direction vector

G4ThreeVector CDMSGammaSphere::getGunDirection() const {
  G4ThreeVector dir = particleGun->GetParticlePosition();	// Outward!
  dir -= GetPosition();			// Remove offset of sphere center
  dir *= (inwardGammas?-1.:1.);

  return dir.unit();
}


// Search process list for RadioactiveDecay (there should be only one)

G4RadioactiveDecay* CDMSGammaSphere::findRadioactiveDecayProcess() const {
  if (verboseLevel > 1)
    G4cout << "CDMSGammaSphere::findRadioactiveDecayProcess" << G4endl;

  G4RadioactiveDecay* rdm = 0;		// Buffer for return value

  static G4ParticleDefinition* allIons = G4GenericIon::Definition();

  G4ProcessManager* procMan = allIons->GetProcessManager();
  if (procMan) {
    G4ProcessVector* ionProcs = procMan->GetProcessList();
    if (ionProcs) {
      for (G4int ip=0; ip<ionProcs->size(); ip++) {
	rdm = dynamic_cast<G4RadioactiveDecay*>((*ionProcs)[ip]);
	if (rdm) break;			// Found the right process
      }
    }
  }

  if (rdm && verboseLevel > 1)
    G4cout << " found RadioactiveDecay associated with GenericIon" << G4endl;

  return rdm;		// Will be null if conditionals above failed
}


// Report configuration using CDMSVDetectorGeometry interface

void CDMSGammaSphere::PrintParameters(std::ostream& os) const {
  os << GetName() << " parameters"
     << "\n radius " << radius << " mm"
     << "\n center " << GetPosition() << " mm"
     << "\n " << particlesPerEvent << " gammas/event,"
     << " radially " << (inwardGammas?"in":"out") << "ward" << G4endl;

  generator.PrintSources(os);
}
