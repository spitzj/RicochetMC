// $Id: CDMSNeutronWall.cc,v 1.4 2010/12/11 05:29:34 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSNeutronWall.cc                                    //
//  Description: Example G4PrimaryGeneratorAction for CDMS Simulation //
//                                                                    //
//  Author:      Dennis Wright                                        //
//  Date:        16 July 2010                                         //
//                                                                    //
//  20101026  M. Kelsey -- Use CDMSVSourceConstruction as base, add   //
//            data members for all gun parameters (settable?)         //
//  20101210  M. Kelsey -- Follow geometry class change to dimensions //
////////////////////////////////////////////////////////////////////////

#include "CDMSsources/CDMSNeutronWall.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

// FIXME:  Want to have parameters configurable from macro
CDMSNeutronWall::CDMSNeutronWall() 
  : CDMSVSourceConstruction(), particlesPerEvent(1), particleName("neutron") {
  direction.set(0.,1.,0.);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName);

  particleGun->SetNumberOfParticles(particlesPerEvent);
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(direction);
}

CDMSNeutronWall::~CDMSNeutronWall() {}


void CDMSNeutronWall::GeneratePrimaries(G4Event* anEvent)
{
  position.set(2*(G4UniformRand()-0.5)*m, -3.0*m, 2*(G4UniformRand()-0.5)*m);
  particleGun->SetParticlePosition(position);
  
  // Deal with the weird spectrum
  G4double randGen = G4UniformRand();
  G4double randEnergy;
  if(randGen < 2./3.)
    randEnergy = G4UniformRand()*MeV;
  else
    randEnergy = (2 - sqrt(1 - G4UniformRand()))*MeV;

  particleGun->SetParticleEnergy(randEnergy);
  particleGun->GeneratePrimaryVertex(anEvent);
}

