////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCZipConstruction.hh                               //
//                                                                    //
//  Description: Construction of single ZIP with enclosure            //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        3 November 2010                                      //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "RMCgeometry/detectors/RMCZipConstruction.hh"
#include "RMCgeometry/detectors/RMCZipSD.hh"
#include "RMCgeometry/interface/RMCZipMessenger.hh"
#include "RMCg4base/RMCMaterialTable.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Polyhedra.hh"
#include "G4SDManager.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include <cmath>
#include <iostream>


RMCZipConstruction::RMCZipConstruction() 
  : RMCVDetectorGeometry("Zip"), makeHousingForZip(false),
    ZipMaterial("G4_Si"), HousingMaterial("G4_Cu"), DetectorName("RMCZip"),
    ZipRad(50*mm), ZipThick(25.4*mm), ZipAxis1Len(100*mm), ZipAxis2Len(100*mm),
    ZipClearanceR(0.67*mm), ZipClearanceZ(1*mm), 
    HousingSides(12), HousingThickness(4*mm), housingHeight(0.),
    housingRinner(0.), housingRouter(0.),
    messenger(new RMCZipMessenger(this)) {}

RMCZipConstruction::~RMCZipConstruction() {
  delete messenger; messenger=0;
}


G4double RMCZipConstruction::GetRadius() const {
  return (makeHousingForZip ? housingRouter : ZipRad);
}

G4double RMCZipConstruction::GetLength() const {
  return (makeHousingForZip ? housingHeight : ZipThick);
}

G4Material* RMCZipConstruction::GetZipG4Material() const {
  return RMCMaterialTable::GetMaterial(ZipMaterial);
}

G4Material* RMCZipConstruction::GetHousingG4Material() const {
  return RMCMaterialTable::GetMaterial(HousingMaterial);
}


void RMCZipConstruction::FillExtraParameters() {
  housingHeight = ZipThick + ZipClearanceZ;
  housingRinner = ZipRad + ZipClearanceR;
  housingRouter = housingRinner + HousingThickness;
}


G4LogicalVolume* RMCZipConstruction::BuildGeometry() {
  if (verboseLevel)
    G4cout << "RMCZipConstruction::BuildGeometry()" << G4endl;

  FillExtraParameters();	// Ensure constituents are up to date

  if (verboseLevel > 1) PrintParameters(G4cout);

  G4LogicalVolume* logicalZip = 0;

  if (makeHousingForZip) {
    G4Material* vacuum = RMCMaterialTable::GetMaterial("G4_Galactic");
    logicalZip =
      new G4LogicalVolume(BuildHousingShape(true), vacuum, GetName(),
			  0,0,0);
    
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), BuildZip(), "PhysicalZip",
		      logicalZip, false, copyNumber);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), BuildHousing(),
		      "PhysicalHousing", logicalZip, false, copyNumber);
  } else {
    logicalZip = BuildZip();
  }

  if (verboseLevel) 
    G4cout << "ZIP mass " << logicalZip->GetMass()/kg << " kg" << G4endl;

  return logicalZip;
}

G4LogicalVolume* RMCZipConstruction::BuildZip() {
  if (verboseLevel>1) G4cout << "RMCZipConstruction::BuildZip()" << G4endl;

  G4Material* mat = GetZipG4Material();
  if (!mat) return 0;		// Can't build from non-existent material

  // Start with a cylindrical crystal
  G4VSolid* zipShape = new G4Tubs(GetName(), 0, ZipRad, ZipThick/2., 0, 360*deg);

  // Make side cuts if necessary
  if (ZipAxis1Len > 0. && ZipAxis2Len > 0.) {
    if (verboseLevel>2) {
      G4cout << " Making " << ZipMaterial << " cylinder with clipped sides"
	     << G4endl;
    }

    G4Box* zipCutBox = new G4Box("ZipCutBox", ZipAxis1Len/2., ZipAxis2Len/2.,
				 (ZipThick+tolerance)/2.);
    zipShape = new G4IntersectionSolid(GetName(), zipShape, zipCutBox);
  }

  G4LogicalVolume* logicalZip =
    new G4LogicalVolume(zipShape, mat, GetName(), 0, 0, 0);

  G4VisAttributes* zipAtt;
  if (ZipMaterial == "G4_Si") {
    zipAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  } else {
    zipAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  }
  zipAtt->SetForceSolid(true);
  logicalZip->SetVisAttributes(zipAtt);

  // Attach data recording system
  logicalZip->SetSensitiveDetector(BuildSensitiveDetector());

  return logicalZip;
}


G4LogicalVolume* RMCZipConstruction::BuildHousing() {
  if (verboseLevel>1) G4cout << "RMCZipConstruction::BuildHousing()" << G4endl;

  G4Material* mat = GetHousingG4Material();
  if (!mat) return 0;		// Can't build from non-existent material

  G4LogicalVolume* logicalHousing =
    new G4LogicalVolume(BuildHousingShape(), mat, "LogicalHousing", 0,0,0);
  G4VisAttributes* housingAtt = new G4VisAttributes(G4Colour(1.0,0.6,0.4));
  housingAtt->SetForceSolid(true);
  logicalHousing->SetVisAttributes(housingAtt);

  return logicalHousing;
}


G4VSolid* RMCZipConstruction::BuildHousingShape(G4bool solid) {
  if (verboseLevel>1)
    G4cout << "RMCZipConstruction::BuildHousingShape()" << G4endl;

  G4double zEnds[2] = { housingHeight/2., -housingHeight/2. };
  G4double rInner[2] = { solid?0:housingRinner, solid?0:housingRinner };
  G4double rOuter[2] = { housingRouter, housingRouter };

  return new G4Polyhedra("HousingShape", 0, 360*deg, HousingSides, 2, 
			 zEnds, rInner, rOuter);
}


G4VSensitiveDetector* RMCZipConstruction::BuildSensitiveDetector() {
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4VSensitiveDetector* azipSD =
    SDman->FindSensitiveDetector(DetectorName, false);	// Search quietly

  if (!azipSD) {
    if (verboseLevel > 1)
      G4cout << " Creating sensitive detector " << DetectorName << G4endl;

    azipSD = new RMCZipSD(DetectorName);
    SDman->AddNewDetector(azipSD);
  }

  return azipSD;
}


void RMCZipConstruction::PrintParameters(std::ostream& os) const {
  os << " RMCZipConstruction parameters"
     << "\n makeHousingForZip " << (makeHousingForZip?"true":"false")
     << "\n ZipMaterial " << ZipMaterial
     << "\n HousingMaterial " << HousingMaterial
     << "\n ZipRad " << ZipRad << " mm"
     << "\n ZipThick " << ZipThick << " mm"
     << "\n ZipAxis1Len " << ZipAxis1Len << " mm"
     << "\n ZipAxis2Len " << ZipAxis2Len << " mm"
     << "\n ZipClearanceR " << ZipClearanceR << " mm"
     << "\n ZipClearanceZ " << ZipClearanceZ << " mm"
     << "\n HousingSides " << HousingSides
     << "\n HousingThickness " << HousingThickness << " mm"    
     << std::endl;
  
  if (makeHousingForZip) {
    os << " housingHeight " << housingHeight << " mm"
       << "\n housingRinner " << housingRinner << " mm"
       << "\n housingRouter " << housingRouter << " mm"
       << std::endl;
  }
}
