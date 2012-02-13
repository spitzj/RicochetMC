// $Id: CDMSVParticleGenerator.cc,v 1.2 2011/07/06 04:08:21 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSVParticleGenerator.hh                            //
//  Description: base class for generating background particles with  //
//               specific characteristics (radiogenic or cosmogenic)  //
//                                                                    //
//  Client code: Configure particle type in ctor or with SetParticle; //
//               On each event, get back 4-vector and particle type.  //
//                                                                    //
//  Subclasses:  Implement functions to throw energy, direction.      //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        16 June 2011                                         //
//                                                                    //
//  20110705  M. Kelsey -- Add call to new shootParticle().           //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSsources/CDMSVParticleGenerator.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"


// Null vector for convenience

const G4ThreeVector CDMSVParticleGenerator::origin(0.,0.,0.);


// Constructor and destructor

CDMSVParticleGenerator::CDMSVParticleGenerator(const G4String& name,
					       G4ParticleDefinition* type)
  : sourceName(name), particle(type ? type : G4Gamma::Definition()),
    verboseLevel(0) {}

CDMSVParticleGenerator::~CDMSVParticleGenerator() {}


// Use subclass functions to generate particle dynamics

const G4LorentzVector& CDMSVParticleGenerator::shoot() const {
  // Make sure we have a particle definition, to construct four vector
  particle = shootParticle();			// Subclass generates this
  if (!particle) particle = G4Gamma::Definition();

  G4double ekin = shootKineticEnergy();		// Subclass generates this
  G4double m = particle->GetPDGMass();
  G4double p = std::sqrt(ekin*(ekin + 2.*m));	// Minimize calculations

  G4cout << "ekin = " << ekin << G4endl;

  lastShoot.setVectM(p*shootDirection(), m);	// Save value in buffer
  return lastShoot;
}


// Fill in user buffers after generating particle

void CDMSVParticleGenerator::shoot(G4LorentzVector& mom,
				   G4ParticleDefinition*& type) const {
  mom = shoot();
  type = GetParticle();
}

void CDMSVParticleGenerator::shoot(G4double& ekin, G4ThreeVector& dir,
				   G4ParticleDefinition*& type) const {
  shoot();				// Fills "lastShoot" buffer

  ekin = lastShoot.e()-lastShoot.m(); 
  dir = lastShoot.vect().unit();
  type = GetParticle();
}

