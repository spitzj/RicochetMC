////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CosmogenicSource.cc                                  //
//  Description: Primary track generator using CosmicMuon generator   //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Adapted from: Dennis Wright (SLAC)                   //
//  Date:        13 February 2012                                     //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "RMCsources/CosmogenicSource.hh"
#include "RMCsources/CosmogenicSourceMessenger.hh"
#include "RMCsources/CosmicMuon.hh"
#include "RMCg4base/RMCMaterialTable.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include <math.h>


CosmogenicSource::CosmogenicSource(G4String name)
 : RMCVSourceConstruction(name), R_active(20.*m),
   eventTime(0.), beamEnergy(0.), paramChange(false),
   smearingXwidth(0.0*m), smearingYwidth(0.0*m),
   muonSpectrum(new CosmicMuon),
   theMessenger(new CosmogenicSourceMessenger(this) )
{ 
  Esum = 0.0;
  NEsum = 0;
}


CosmogenicSource::~CosmogenicSource()
{
  // Report cumulative statistics
  if (verboseLevel) 
    G4cout << " CosmogenicSource: Mean energy (GeV) = " << Esum/G4double(NEsum);

  G4cout << "CosmogenicSource::~CosmogenicSource()" << G4endl;

  delete muonSpectrum;
  delete theMessenger;
}


G4LogicalVolume* CosmogenicSource::BuildGeometry()
{
  if (verboseLevel) G4cout << " SurfaceCosmogenicSource::BuildGeometry" << G4endl;

  G4Sphere* sourceSphere = new G4Sphere("CosmicSourceSphere", 0.0, R_active,
					0*deg, 360*deg, 0*deg, 90*deg);
  G4Material* srcMat = RMCMaterialTable::GetMaterial("G4_Galactic");
  G4LogicalVolume* logicalSource =
    new G4LogicalVolume(sourceSphere, srcMat, "LogicalCosmicSource", 0,0,0);

  // Visualization attributes
  G4VisAttributes* VisAttSourceDisk = new G4VisAttributes(G4Color(1,1,0));
  VisAttSourceDisk->SetForceWireframe(true);
  logicalSource->SetVisAttributes(VisAttSourceDisk);

  return logicalSource;
}


void CosmogenicSource::GeneratePrimaries(G4Event* anEvent)
{ 
  G4int cycles = 0;
  
  GenerateBeam();
  cycles++;
  GenerateTime();	// Don't increment time until valid particle made
  
  if (verboseLevel > 2) {
    G4cout << "Generated valid muon in " << cycles << " cycles." << G4endl;
  }	


  Esum += beamEnergy/GeV;
  NEsum += 1;

  if (verboseLevel > 1) {
    G4cout << "Generated Muon @ " << beamPos/m << " m"
	   << "\n  Phi " << beamDir.phi() << " Theta " << beamDir.theta()
	   << " Energy " << beamEnergy/GeV << " GeV"
	   << "\nNEsum = " << NEsum << " , Running mean = "
	   << Esum/G4double(NEsum) << G4endl;
  }

  FillParticleGun();
  particleGun->GeneratePrimaryVertex(anEvent);
}


// NOTE: distributions take theta as angle from downward direction, but
// particle gun momentum direction requires angle from upwards direction

void CosmogenicSource::FillParticleGun() {
  beamDir.set(-beamDir.x(), -beamDir.y(), -beamDir.z());
  particleGun->SetParticleMomentumDirection(beamDir);

  particleGun->SetParticleEnergy(beamEnergy);

  particleGun->SetParticleDefinition(partDef);
  particleGun->SetParticlePosition(beamPos);

  particleGun->SetParticleTime(eventTime);
}

// Generates event time (cumulative since start of run)
void CosmogenicSource::GenerateTime()
{
  G4double area = 0.;
  eventTime += muonSpectrum->RandomTime(area);
}


// Generate a random position in the top plane of the lab volume
void CosmogenicSource::GeneratePosition()
{
  beamPos.set(R_active*sin(beamDir.getTheta())*cos(beamDir.getPhi()),
	      R_active*sin(beamDir.getTheta())*sin(beamDir.getPhi()),
	      R_active*cos(beamDir.getTheta()));
  
  //beamPos += position;		// Convert from local to global coordinates

  // Apply smearing over 5m by 5m square
  G4double xSmear = (G4UniformRand() - 0.5) * smearingXwidth;
  G4double ySmear = (G4UniformRand() - 0.5) * smearingYwidth;
  beamPos.set(beamPos.getX() + xSmear,
	      beamPos.getY() + ySmear, beamPos.getZ());

  beamPos += (position + position);		// Convert from local to global coordinates
}


// Generate random direction, energy, and particle
void CosmogenicSource::GenerateBeam()
{
  muonSpectrum->shoot(beamEnergy, beamDir, partDef);
  GeneratePosition();
}


