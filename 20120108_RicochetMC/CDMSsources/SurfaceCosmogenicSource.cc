// $Id: CosmogenicSource.cc,v 1.4 2011/07/22 21:08:36 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CosmogenicSource.cc                                  //
//  Description: Primary track generator based on the Cassiday        //
//               (Groom) spectrum parameterization                    //
//               (Phys. Rev. D7, 2022, 1973)                          //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        27 May 2011                                          //
//                                                                    //
//  20110722  M. Kelsey -- Minor rewrite of generating functions to   //
//		make better use of CDMSVParticleGenerator interface.  //
//              Discard matTable, using singleton instead.            //
////////////////////////////////////////////////////////////////////////
//

// $Id: CosmogenicSource.cc,v 1.4 2011/07/22 21:08:36 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CosmogenicSource.cc                                  //
//  Description: Primary track generator based on the Cassiday        //
//               (Groom) spectrum parameterization                    //
//               (Phys. Rev. D7, 2022, 1973)                          //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        27 May 2011                                          //
//                                                                    //
//  20110722  M. Kelsey -- Minor rewrite of generating functions to   //
//		make better use of CDMSVParticleGenerator interface.  //
//              Discard matTable, using singleton instead.            //
////////////////////////////////////////////////////////////////////////
//

#include "CDMSsources/SurfaceCosmogenicSource.hh"
#include "CDMSsources/CosmogenicSourceMessenger.hh"
#include "CDMSsources/SurfaceCosmicMuon.hh"
#include "CDMSg4base/CDMSMaterialTable.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include <math.h>


SurfaceCosmogenicSource::SurfaceCosmogenicSource(G4String name)
 : CDMSVSourceConstruction(name), R_active(1*m), L_active(1*mm),
   eventTime(0), beamEnergy(0.), 
   hallX(4.*m), hallY(4.*m), hallZ(4.*m), rockDensity(1*g/cm3),
   overBurden(10.*m), muEffectRad(10.*m), paramChange(false),
   muonSpectrum(new SurfaceCosmicMuon)
{ 
  GenerationArea = pi*R_active*R_active;
  EffectiveXRange = hallX/2. + muEffectRad;
  EffectiveYRange = hallY/2. + muEffectRad;
  EffectiveZTravel = hallZ + overBurden/2. + muEffectRad;
  Esum = 0.0;
  NEsum = 0;
}


SurfaceCosmogenicSource::~SurfaceCosmogenicSource()
{
  // Report cumulative statistics
  if (verboseLevel) 
    G4cout << " CosmogenicSource: Mean energy (GeV) = " << Esum/G4double(NEsum);

  delete muonSpectrum;
  delete theMessenger;
}


G4LogicalVolume* SurfaceCosmogenicSource::BuildGeometry()
{
  if (verboseLevel) G4cout << " CosmogenicSource::BuildGeometry" << G4endl;

  G4Tubs* sourceDisk = new G4Tubs("CosmicSourceDisk", 0.0, R_active, L_active/2,
                                  0, 360*deg);
  G4Material* srcMat = CDMSMaterialTable::GetMaterial("G4_Galactic");
  G4LogicalVolume* logicalSource =
    new G4LogicalVolume(sourceDisk, srcMat, "LogicalCosmicSource", 0,0,0);

  // Visualization attributes
  G4VisAttributes* VisAttSourceDisk = new G4VisAttributes(G4Color(1,1,0));
  VisAttSourceDisk->SetForceSolid(true);
  logicalSource->SetVisAttributes(VisAttSourceDisk);

  return logicalSource;
}


void SurfaceCosmogenicSource::ResetGenerator()
{
  muonSpectrum->Initialize();
  if (verboseLevel > 1)
    G4cout << " Generator parameters changed: reset " << G4endl;
}

void SurfaceCosmogenicSource::TestGenerator()
{
  muonSpectrum->Test();
}

void SurfaceCosmogenicSource::SetLabDepth(G4double depth)
{
  G4double mweDepth = (depth - overBurden - hallZ)*rockDensity/(g/cm3);
  muonSpectrum->SetMWEDepth(mweDepth);
  paramChange = true;

  if (verboseLevel > 1)
    G4cout << " SetLabDepth: rock density = " << rockDensity/(g/cm3)
           << " g/cm3, depth = " << depth/m
           << " m , MWEDepth = " << mweDepth/m
           << " m " << G4endl;  
}


void SurfaceCosmogenicSource::GeneratePrimaries(G4Event* anEvent)
{
  if (paramChange) {
    ResetGenerator();
    paramChange = false;
  }

  G4int cycles = 0;
  
  do {
    GeneratePosition();
    GenerateBeam();
    cycles++;
  } while (InvalidDirection());

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
	   << Esum/G4double(NEsum) << " partdef = " 
	   << partDef->GetPDGMass() << " eventTime = " << eventTime << G4endl;
  }

  FillParticleGun();
  G4cout << "test" << G4endl;
  particleGun->GeneratePrimaryVertex(anEvent);
  G4cout << "generated primary vertex" << G4endl;
}


// NOTE: distributions take theta as angle from downward direction, but
// particle gun momentum direction requires angle from upwards direction

void SurfaceCosmogenicSource::FillParticleGun() {
  beamDir.setZ(-beamDir.z());
  particleGun->SetParticleMomentumDirection(beamDir);

  particleGun->SetParticleEnergy(beamEnergy);

  particleGun->SetParticleDefinition(partDef);

  particleGun->SetParticlePosition(beamPos);

  particleGun->SetParticleTime(eventTime);
}


// Generates event time (cumulative since start of run)
void SurfaceCosmogenicSource::GenerateTime()
{
  eventTime += muonSpectrum->RandomTime(GenerationArea);
}


// Generate a random position in the top plane of the lab volume
void SurfaceCosmogenicSource::GeneratePosition()
{
  G4double r = R_active*std::sqrt(G4UniformRand());
  G4double phi = twopi*G4UniformRand();

  beamPos.setRhoPhiZ(r,phi,0.);
  beamPos += position;		// Convert from local to global coordinates
}


// Generate random direction, energy, and particle
void SurfaceCosmogenicSource::GenerateBeam()
{
  muonSpectrum->shoot(beamEnergy, beamDir, partDef);
}


// Coordinates x,y,z,phi,theta are invalid if the current location
// and direction will not pass near the Soudan2 Hall
G4bool SurfaceCosmogenicSource::InvalidDirection()
{
  // Extract coordinates from position and direction vectors
  G4double x_pos = beamPos.x(), y_pos = beamPos.y();
  G4double theta_dir = beamDir.theta(), phi_dir = beamDir.phi();

  G4double phi_1;	// First angle in ccw direction from phi = 0
  G4double phi_2;	// Second angle in ccw direction from phi = 0
  G4double min_xy_travel;
  G4double max_xy_travel;
  
  // Maximum distance muon will travel in xy plane before passing outside
  // of vertically acceptable region
  max_xy_travel = EffectiveZTravel * tan(theta_dir);
	
  // Check -x location
  if (x_pos < (0 - EffectiveXRange)) {
    // Non-detector lies between phi_1 and phi_2
		
    phi_1 = atan2(0 + EffectiveYRange - y_pos,
		  0 - EffectiveXRange - x_pos);
    if (phi_1 < 0) phi_1 = 0;
    if (phi_dir > phi_1) return true;
		
    phi_2 = atan2(0 - EffectiveYRange - y_pos,
                  0 - EffectiveXRange - x_pos);
    if (phi_2 > 0) phi_2 = 0;
    if (phi_dir < phi_2) return true;
		
    // Check if muon is too vertical to pass near Soudan2 Hall
    min_xy_travel = fabs((0 - EffectiveXRange - x_pos) / cos(phi_dir));
    if (max_xy_travel < min_xy_travel) return true;
  }
	
  // Check +x location
  if (x_pos > (0 + EffectiveXRange)) {
    // Detector lies between phi_1 and phi_2
		
    phi_1 = atan2(0 + EffectiveYRange - y_pos,
                  0 + EffectiveXRange - x_pos);
    if (phi_1 < 0) phi_1 = pi;
		
    phi_2 = atan2(0 - EffectiveYRange - y_pos,
                  0 + EffectiveXRange - x_pos);
    if (phi_2 > 0) phi_2 = pi;
		
    if ((phi_dir < phi_1) & (phi_dir > phi_2)) return true;
		
    // Check if muon is too vertical to pass near Soudan2 Hall
    min_xy_travel = fabs((0 + EffectiveXRange - x_pos) / cos(phi_dir));
    if (max_xy_travel < min_xy_travel) return true;
  }
	
  // Check -y location
  if (y_pos < (0 - EffectiveYRange)) {
    // Detector lies between phi_1 and phi_2
		
    phi_1 = atan2(0 - EffectiveYRange - y_pos,
                  0 + EffectiveXRange - x_pos);
    if (phi_1 > halfpi) phi_1 = halfpi;
    if (phi_dir < phi_1) return true;
		
    phi_2 = atan2(0 - EffectiveYRange - y_pos,
                  0 - EffectiveXRange - x_pos);
    if (phi_2 < halfpi) phi_2 = halfpi;
    if (phi_dir > phi_2) return true;
		
    // Check if muon is too vertical to pass near Soudan2 Hall
    min_xy_travel = fabs((0 - EffectiveYRange - y_pos) / sin(phi_dir));
    if (max_xy_travel < min_xy_travel) return true;
  }
	
  // Check +y location
  if (y_pos > (0 + EffectiveYRange)) {
    // Detector lies between phi_1 and phi_2
		
    phi_1 = atan2(0 + EffectiveYRange - y_pos,
                  0 - EffectiveXRange - x_pos);
    if (phi_1 > 0 - halfpi) phi_1 = 0 - halfpi;
    if (phi_dir < phi_1) return true;
		
    phi_2 = atan2(0 + EffectiveYRange - y_pos,
                  0 + EffectiveXRange - x_pos);
    if (phi_2 < 0 - halfpi) phi_2 = 0 - halfpi;
    if (phi_dir > phi_2) return true;
		
    // Check if muon is too vertical to pass near Soudan2 Hall
    min_xy_travel = fabs((0 + EffectiveYRange - y_pos) / sin(phi_dir));
    if (max_xy_travel < min_xy_travel) return true;
  }
	
  return false;
}


/*#include "CDMSsources/SurfaceCosmogenicSource.hh"
#include "CDMSsources/CosmogenicSourceMessenger.hh"
#include "CDMSsources/SurfaceCosmicMuon.hh"
#include "CDMSg4base/CDMSMaterialTable.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include <math.h>


SurfaceCosmogenicSource::SurfaceCosmogenicSource(G4String name)
 : CDMSVSourceConstruction(name), R_active(10*m),
   eventTime(0), beamEnergy(0.),  paramChange(false),
   smearingXwidth(5.0*m), smearingYwidth(5.0*m),
   muonSpectrum(new SurfaceCosmicMuon)
{ 
  Esum = 0.0;
  NEsum = 0;
}


SurfaceCosmogenicSource::~SurfaceCosmogenicSource()
{
  delete muonSpectrum;
  delete theMessenger;
}


G4LogicalVolume* SurfaceCosmogenicSource::BuildGeometry()
{
  if (verboseLevel) G4cout << " SurfaceCosmogenicSource::BuildGeometry" << G4endl;

  G4Sphere* sourceSphere = new G4Sphere("CosmicSourceSphere", 0.0, R_active,
					0*deg, 360*deg, 0*deg, 90*deg);
  G4Material* srcMat = CDMSMaterialTable::GetMaterial("G4_Galactic");
  G4LogicalVolume* logicalSource =
    new G4LogicalVolume(sourceSphere, srcMat, "LogicalCosmicSource", 0,0,0);

  // Visualization attributes
  G4VisAttributes* VisAttSourceDisk = new G4VisAttributes(G4Color(1,1,0));
  VisAttSourceDisk->SetForceSolid(true);
  logicalSource->SetVisAttributes(VisAttSourceDisk);

  return logicalSource;
}


void SurfaceCosmogenicSource::GeneratePrimaries(G4Event* anEvent)
{ 
  G4int cycles = 0;
  
G4cout << "primaries..." << G4endl;
  GenerateBeam();
  cycles++;
  GenerateTime();	// Don't increment time until valid particle made
  
  if (verboseLevel > 2) {
    G4cout << "Generated valid muon in " << cycles << " cycles." << G4endl;
  }	


  Esum += beamEnergy/GeV;
  NEsum += 1;

  //if (verboseLevel > 1) {
    G4cout << "Generated Muon @ " << beamPos/m << " m"
	   << "\n  Phi " << beamDir.phi() << " Theta " << beamDir.theta()
	   << " Energy " << beamEnergy/GeV << " GeV"
	   << "\nNEsum = " << NEsum << " , Running mean = "
	   << Esum/G4double(NEsum) << G4endl;
    //}

  FillParticleGun();
  particleGun->GeneratePrimaryVertex(anEvent);
}


// NOTE: distributions take theta as angle from downward direction, but
// particle gun momentum direction requires angle from upwards direction

void SurfaceCosmogenicSource::FillParticleGun() {
  beamDir.setZ(-beamDir.z());
  particleGun->SetParticleMomentumDirection(beamDir);

  particleGun->SetParticleEnergy(beamEnergy);

  particleGun->SetParticleDefinition(partDef);

  particleGun->SetParticlePosition(beamPos);

  particleGun->SetParticleTime(eventTime);
}


// Generates event time (cumulative since start of run)
void SurfaceCosmogenicSource::GenerateTime()
{
  G4double area = 0.;
  eventTime += muonSpectrum->RandomTime(area);
}


// Generate a random position in the top plane of the lab volume
void SurfaceCosmogenicSource::GeneratePosition()
{
  beamPos.set(R_active*sin(beamDir.getTheta())*cos(beamDir.getPhi()),
	      R_active*sin(beamDir.getTheta())*sin(beamDir.getPhi()),
	      R_active*cos(beamDir.getTheta()));
  
  //beamPos += position;		// Convert from local to global coordinates

  // Apply smearing over 5m by 5m square
  G4double xSmear = (G4UniformRand() - 0.5) * smearingXwidth;
  G4double ySmear = (G4UniformRand() - 0.5) * smearingYwidth;
  G4cout << "xSmear = " << smearingXwidth << ", ySmear = " << smearingYwidth << G4endl;
  beamPos.set(beamPos.getX() + xSmear,
	      beamPos.getY() + ySmear, beamPos.getZ());
}


// Generate random direction, energy, and particle
void SurfaceCosmogenicSource::GenerateBeam()
{
  muonSpectrum->shoot(beamEnergy, beamDir, partDef);
  GeneratePosition();
}*/


