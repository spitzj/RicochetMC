////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCVSourceConstruction.hh                            //
//  Description: base class for geometrically based RMC sources       //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Adapted from: Michael Kelsey (SLAC)                  //
//  Date:        29 November 2010                                     //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "RMCg4base/RMCVSourceConstruction.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4Geantino.hh"

// Constructor and destructor

RMCVSourceConstruction::RMCVSourceConstruction(const G4String& nameString)
  : G4VUserPrimaryGeneratorAction(), RMCVDetectorGeometry(nameString),
    particleGun(new G4ParticleGun) {}


RMCVSourceConstruction::~RMCVSourceConstruction() {
  delete particleGun;
}


// Set beam type, for convenience of subclasses

void RMCVSourceConstruction::SetParticleType(G4ParticleDefinition* beam) {
  if (!beam) beam = G4Geantino::Definition();
  particleGun->SetParticleDefinition(beam);
}
