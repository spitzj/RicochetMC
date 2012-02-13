// $Id: CDMSVSourceConstruction.cc,v 1.4 2011/06/29 22:24:54 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSVSourceConstruction.hh                           //
//  Description: base class for geometrically based CDMS sources      //
//               (radiogenic or cosmogenic)                           //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        29 November 2010                                     //
//                                                                    //
//  20110426  M. Kelsey -- Add interface to set beam type             //
//  20110629  M. Kelsey -- Fix constness of ctor arg                  //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSg4base/CDMSVSourceConstruction.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4Geantino.hh"

// Constructor and destructor

CDMSVSourceConstruction::CDMSVSourceConstruction(const G4String& nameString)
  : G4VUserPrimaryGeneratorAction(), CDMSVDetectorGeometry(nameString),
    particleGun(new G4ParticleGun) {}


CDMSVSourceConstruction::~CDMSVSourceConstruction() {
  delete particleGun;
}


// Set beam type, for convenience of subclasses

void CDMSVSourceConstruction::SetParticleType(G4ParticleDefinition* beam) {
  if (!beam) beam = G4Geantino::Definition();
  particleGun->SetParticleDefinition(beam);
}
