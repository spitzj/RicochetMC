// $Id: CDMSGeomConstructor.cc,v 1.18 2011/07/22 21:09:17 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//  File:        CDMSVDetectorGeometry.hh                             //
//  Description: Wrapper class to build any user-defined geometry     //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        26 October 2010                                      //
//                                                                    //
//  20101210  M. Kelsey -- GetLocation() changed to GetPosition(),    //
//		GetWorldRadius() and GetWorldLength() dimensions.     //
//  20110105  M. Kelsey -- Drop "shield" -- move to MultiTower det.   //
//  20110215  M. Kelsey -- Don't delete detector components.          //
//  20110429  M. Kelsey -- Move material table to data member.        //
//  20110520  M. Kelsey -- Implement SetVerboseLevel passing to child //
//  20110525  M. Kelsey -- Allow non-geometric sources to set radius  //
//  20110722  M. Kelsey -- Discard matTable, using singleton instead. //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSgeometry/interface/CDMSGeomConstructor.hh"

#include "CDMSg4base/CDMSVDetectorGeometry.hh"
#include "CDMSg4base/CDMSVLabConstruction.hh"
#include "CDMSg4base/CDMSMaterialTable.hh"
#include "CDMSgeometry/interface/CDMSGeomMessenger.hh"
#include "CDMSgeometry/interface/CDMSGeometryManager.hh"
#include "CDMSsources/CosmogenicSource.hh"
#include "CDMSsources/SurfaceCosmogenicSource.hh"
#include "G4Box.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RegionStore.hh"
#include "G4RunManager.hh"
#include "G4SolidStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include <cmath>
#include <sstream>

const G4ThreeVector CDMSGeomConstructor::theOrigin(0.,0.,0.);


// Constructor and destructor

CDMSGeomConstructor::CDMSGeomConstructor()
  : verboseLevel(0), messenger(new CDMSGeomMessenger(this)),
    lab(0), detector(0), detectorOffset(G4ThreeVector(3.5*m, 0., 1.5*m)) {
  lab = CDMSGeometryManager::Instance()->GetNoLab();
}

CDMSGeomConstructor::~CDMSGeomConstructor() {
  // NOTE: Detector components are not owned here (GeometryManager)
  // NOTE: Sources are not owned here (RunManager)
  delete messenger;
}


// Pass verbosity through to detector components and sources

void CDMSGeomConstructor::SetVerboseLevel(G4int verbose) {
  verboseLevel = verbose;
  if (lab)      lab->SetVerboseLevel(verboseLevel);
  if (detector) detector->SetVerboseLevel(verboseLevel);

  for (size_t i=0; i<sources.size(); i++) {
    if (sources[i]) sources[i]->SetVerboseLevel(verboseLevel);
  }
}


// Assign or replace detector components as directed

void CDMSGeomConstructor::SetLab(CDMSVLabConstruction* theLab) {
  if (theLab) theLab->SetVerboseLevel(verboseLevel);
  lab = theLab;
}

void CDMSGeomConstructor::SetDetector(CDMSVDetectorGeometry* theDet) {
  if (theDet) theDet->SetVerboseLevel(verboseLevel);
  detector = theDet;
}

// Users may have multiple sources simultaneously

void CDMSGeomConstructor::AddSource(CDMSVDetectorGeometry* aSource) {
  if (aSource) {
    aSource->SetVerboseLevel(verboseLevel);
    sources.push_back(aSource);			// Don't add null pointers!
  }
}


// Register geometry changes with RunManager

void CDMSGeomConstructor::UpdateGeometry() {
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}


// Discard any existing geometry before construction

void CDMSGeomConstructor::ClearGeometry() {
  if (verboseLevel) G4cout << "CDMSGeomConstructor::ClearGeometry" << G4endl;

  G4GeometryManager::GetInstance()->OpenGeometry();
  // G4RegionStore::GetInstance()->Clean();
  G4PhysicalVolumeStore::Clean();
  G4LogicalVolumeStore::Clean();
  G4SolidStore::Clean();

  theWorld = 0;
  theOutside = 0;
}


// Build complete geometry model from current set of components

G4VPhysicalVolume* CDMSGeomConstructor::Construct() {
  if (verboseLevel) G4cout << "CDMSGeomConstructor::Construct" << G4endl;

  ClearGeometry();		// Discard any existing model first!

  ConstructWorld();
  ConstructLab();
  ConstructDetector();
  ConstructSources();

  return theWorld;	// This was set by ConstructWorld()
}

void CDMSGeomConstructor::ConstructWorld() {
  if (verboseLevel > 1)
    G4cout << "CDMSGeomConstructor::ConstructWorld" << G4endl;

  G4Box* worldBox = new G4Box("WorldBox", GetWorldRadius(), GetWorldRadius(),
			      GetWorldLength()/2.);

  // The Universe is empty but not empty...
  G4Material* worldMaterial = CDMSMaterialTable::GetMaterial("G4_Galactic");

  // This will be used when building lab, detectors, etc.
  theOutside = new G4LogicalVolume(worldBox, worldMaterial,"World");
  theOutside->SetVisAttributes(G4VisAttributes::Invisible);

  // This will be returned by base ::Construct() method
  theWorld = new G4PVPlacement(0, theOrigin, theOutside, "World", 0, false, 0);
}

void CDMSGeomConstructor::ConstructLab() {
  if (!lab) return;		// No laboratory defined; just use world

  if (verboseLevel > 1)
    G4cout << "CDMSGeomConstructor::ConstructLab" << G4endl;

  G4LogicalVolume* logicalLab = lab->BuildGeometry();

  // NOTE:  We don't need to capture pointer; handled automatically by GEANT4
  new G4PVPlacement(0, lab->GetPosition()+detectorOffset, logicalLab, "Lab", theOutside, false, 0);

  theOutside = logicalLab;	// This replaces world for subsequent builds
}


void CDMSGeomConstructor::ConstructDetector() {
  if (!detector) return;

  if (verboseLevel > 1)
    G4cout << "CDMSGeomConstructor::ConstructDetector" << G4endl;

  if (verboseLevel > 2)
    G4cout << " Placing detector at " << detector->GetPosition()
	   << " mm" << G4endl;

  detector->SetPosition(detectorOffset);

  // NOTE:  We don't need to capture pointer; handled automatically by GEANT4
  new G4PVPlacement(0, detector->GetPosition(), detector->BuildGeometry(),
		    "Detector", theOutside, false, 0);
}

void CDMSGeomConstructor::ConstructSources() {
  if (sources.empty()) return;

  if (verboseLevel > 1)
    G4cout << "CDMSGeomConstructor::ConstructSources" << G4endl;

  for (size_t i=0; i<sources.size(); i++) ConstructSource(i);
}

void CDMSGeomConstructor::ConstructSource(size_t i) {
  if (verboseLevel > 1)
    G4cout << "CDMSGeomConstructor::ConstructSource(" << i << ")" << G4endl;

  if (i >= sources.size()) {
    G4cerr << "CDMSGeomConstructor::ConstructSource(" << i << ") out of range"
	   << G4endl;
    return;
  }

  CDMSVDetectorGeometry* src = sources[i];	// Only the geometry interface!
  if (!src) {
    G4cerr << " CDMSGeomConstructor has no source [" << i << "]" << G4endl;
    return;
  }

  G4LogicalVolume* srcVol = src->BuildGeometry();
  if (srcVol) {					// Source has geometric extent
    std::ostringstream sourceName;
    sourceName << "Source" << i+1;

    if (src->GetName() == "cosmoSource") {
      G4cout << src->GetName() << " chosen " << G4endl;
      CosmogenicSource* cosmoSrc = (CosmogenicSource*)src;
      //G4double zPos = lab->GetLength()/2. - cosmoSrc->GetLength()/2.;
      //cosmoSrc->SetRadius(lab->GetRadius() );
      //cosmoSrc->SetRockDensity(lab->GetOverBurdenDensity() );
      //cosmoSrc->SetOverBurden(lab->GetOverBurden() );
      //cosmoSrc->SetCavernLength(lab->GetCavernLength() );
      //cosmoSrc->SetCavernWidth(lab->GetCavernWidth() );
      //cosmoSrc->SetCavernHeight(lab->GetCavernHeight() );
      src = cosmoSrc;
      src->SetPosition(detectorOffset);
    }

    if (verboseLevel > 2)
      G4cout << " Placing " << sourceName << " at " << src->GetPosition()
	     << " mm" << G4endl;

    // NOTE:  Don't need to capture pointer; handled automatically by GEANT4
    new G4PVPlacement(0, src->GetPosition(), srcVol, sourceName.str(),
		      theOutside, false, 0);
  }
}


// Evaluate nested objects to determine needed radius of universe
// NOTE:  Sources are assumed to be small compared to lab and detectors

G4double CDMSGeomConstructor::GetWorldRadius() const {
  G4double dx=0., dy=0., rad=1*m;	// Dummy value for empty universe

  CDMSVDetectorGeometry* mother = lab ? lab : detector ? detector : 0;
  if (mother) {
    dx = std::abs(mother->GetPosition().x());
    dy = std::abs(mother->GetPosition().y());
    rad = mother->GetRadius();
  }

    
    G4cout << dx << "  " << dy << "  " << rad << G4endl;
  return (rad + std::max(dx,dy));
}

G4double CDMSGeomConstructor::GetWorldLength() const {
  G4double dz=0., len=1*m;	// Dummy value for empty universe

  CDMSVDetectorGeometry* mother = lab ? lab : detector ? detector : 0;
  if (mother) {
    dz = std::abs(mother->GetPosition().z());
    len = mother->GetLength();
  }

  return (len + dz);
}
