// $Id: Am241Source.cc,v 1.8 2011/07/22 21:08:36 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        Am241Source.cc                                       //
//  Description: specific source class for CDMS                       //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        2 September 2010                                     //
//                                                                    //
//  20101206  M. Kelsey -- Don't create or delete particleGun here.   //
//  20101210  M. Kelsey -- Follow geometry class change to dimensions //
//  20110427  M. Kelsey -- Move gamma lines to separate class.        //
//  20110722  M. Kelsey -- Discard matTable, using singleton instead. //
////////////////////////////////////////////////////////////////////////

#include "CDMSsources/Am241Source.hh"
#include "CDMSsources/Am241Messenger.hh"
#include "CDMSsources/Am241Lines.hh"
#include "CDMSg4base/CDMSMaterialTable.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4VisAttributes.hh"
#include "Randomize.hh"


Am241Source::Am241Source()
  : CDMSVSourceConstruction(), R_active(3.0*cm), L_active(0.2*cm),
    gammaLines(new Am241Lines), sourceMessenger(new Am241Messenger(this)) {
  direction.set(0.,0.,1.0);		// Shoot in +ve Z direction

  // Only gammas are used
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("gamma");
  particleGun->SetParticleDefinition(particle);
}


Am241Source::~Am241Source() {
  delete sourceMessenger;
}


G4LogicalVolume* Am241Source::BuildGeometry() {
  if (verboseLevel) G4cout << " Am241Source::BuildGeometry" << G4endl;
  if (verboseLevel>1) PrintParameters(G4cout);

  G4Tubs* sourceDisk = new G4Tubs("SourceDisk", 0.0, R_active, L_active/2,
                                  0, 360*deg);
  G4Material* srcMat = CDMSMaterialTable::GetMaterial("G4_Galactic");
  G4LogicalVolume* logicalSource = 
     new G4LogicalVolume(sourceDisk, srcMat, "LogicalSource", 0,0,0);

  // Visualization attributes
  G4VisAttributes* VisAttSourceDisk = new G4VisAttributes(G4Colour(1,1,0));
  VisAttSourceDisk->SetForceWireframe(true);
  logicalSource->SetVisAttributes(VisAttSourceDisk);

  return logicalSource;
}


void Am241Source::GeneratePrimaries(G4Event* anEvent)
{
  if (verboseLevel > 1)
    G4cout << " Am241Source::GeneratePrimaries (gammas)" << G4endl;

  // Generate initial particle position (random within active area of source)
  G4double radius = R_active*std::sqrt(G4UniformRand());
  G4double theta = 2.*pi*G4UniformRand();
  G4double x0 = position.x() + radius*std::cos(theta);
  G4double y0 = position.y() + radius*std::sin(theta);
  G4double z0 = position.z();
  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

  if (verboseLevel > 1)
    G4cout << " position " << x0 << " " << y0 << " " << z0 << G4endl;

  // Set direction (if zero, isotropic random distribuion)
  gammaLines->SetDirection(direction);

  G4LorentzVector theGamma = gammaLines->shoot();
  if (verboseLevel > 1) {
    G4cout << " direction " << theGamma.vect().unit() << G4endl
	   << " energy " << theGamma.e()/keV << " keV" << G4endl;
  }

  particleGun->SetParticleMomentum(theGamma.vect());

  particleGun->GeneratePrimaryVertex(anEvent);
}


void Am241Source::PrintParameters(std::ostream& os) const {
  os << " Am241Source parameters"
     << "\n Radius " << R_active << " mm, Thickness " << L_active << " mm"
     << "\n Position " << position << " mm"
     << "\n Direction " << direction << " mm";
  os << std::endl;
}
  
