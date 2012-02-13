// $Id: CDMSDemoSource.cc,v 1.4 2010/12/11 05:29:34 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSDemoSource.cc                                    //
//  Description: Example G4PrimaryGeneratorAction for CDMS Simulation //
//                                                                    //
//  Author:      Dennis Wright                                        //
//  Date:        16 July 2010                                         //
//                                                                    //
//  20101026  M. Kelsey -- Use CDMSVSourceConstruction as base, add   //
//            data members for all gun parameters (settable?)         //
//  20101210  M. Kelsey -- Follow geometry class change to dimensions //
////////////////////////////////////////////////////////////////////////

#include "CDMSsources/CDMSDemoSource.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"


// FIXME:  Want to have parameters configurable from macro
CDMSDemoSource::CDMSDemoSource() 
  : CDMSVSourceConstruction(), particlesPerEvent(1), particleName("e-") {
  position.set(0.,0.,-2.0*cm);
  direction.set(0.,0.,1.);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName);

  particleGun->SetNumberOfParticles(particlesPerEvent);
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(direction);
  particleGun->SetParticlePosition(position);
  particleGun->SetParticleEnergy(1*MeV);
}

CDMSDemoSource::~CDMSDemoSource() {}


void CDMSDemoSource::GeneratePrimaries(G4Event* anEvent)
{
  particleGun->GeneratePrimaryVertex(anEvent);
}

