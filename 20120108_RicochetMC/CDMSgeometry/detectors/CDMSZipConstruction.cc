// $Id: CDMSZipConstruction.cc,v 1.14 2011/07/22 21:09:17 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSZipConstruction.hh                               //
//                                                                    //
//  Description: Construction of single ZIP with enclosure            //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        3 November 2010                                      //
//                                                                    //
//  20101130  M. Kelsey -- Add verbosity reporting                    //
//  20101203  M. Kelsey -- Add flag to suppress housing construction  //
//  20101208  M. Kelsey -- Use FillExtraParameters() for housing pars //
//  20101209  M. Kelsey -- Add support for sensitive detector         //
//  20101210  M. Kelsey -- Follow base class change to dimensions     //
//  20110204  M. Kelsey -- Add accessors to get material properties   //
//  20110215  M. Kelsey -- Report mass at end of construction.        //
//  20110421  M. Kelsey -- Housing radius doesn't need cos adjust;    //
//		reduce radial clearance to 0.67 mm.                   //
//  20110502  M. Kelsey -- Mother volume with housing must be solid   //
//  20110629  M. Kelsey -- Initialize and use name arg in base class  //
//  20110722  M. Kelsey -- Discard matTable, using singleton instead. //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSgeometry/detectors/CDMSZipConstruction.hh"
#include "CDMSgeometry/detectors/CDMSZipSD.hh"
#include "CDMSgeometry/interface/CDMSZipMessenger.hh"
#include "CDMSg4base/CDMSMaterialTable.hh"
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


CDMSZipConstruction::CDMSZipConstruction() 
  : CDMSVDetectorGeometry("Zip"), makeHousingForZip(false),
    ZipMaterial("G4_Si"), HousingMaterial("G4_Cu"), DetectorName("CDMSZip"),
    ZipRad(50*mm), ZipThick(25.4*mm), ZipAxis1Len(100*mm), ZipAxis2Len(100*mm),
    ZipClearanceR(0.67*mm), ZipClearanceZ(1*mm), 
    HousingSides(12), HousingThickness(4*mm), housingHeight(0.),
    housingRinner(0.), housingRouter(0.),
    messenger(new CDMSZipMessenger(this)) {}

CDMSZipConstruction::~CDMSZipConstruction() {
  delete messenger; messenger=0;
}


G4double CDMSZipConstruction::GetRadius() const {
  return (makeHousingForZip ? housingRouter : ZipRad);
}

G4double CDMSZipConstruction::GetLength() const {
  return (makeHousingForZip ? housingHeight : ZipThick);
}

G4Material* CDMSZipConstruction::GetZipG4Material() const {
  return CDMSMaterialTable::GetMaterial(ZipMaterial);
}

G4Material* CDMSZipConstruction::GetHousingG4Material() const {
  return CDMSMaterialTable::GetMaterial(HousingMaterial);
}


void CDMSZipConstruction::FillExtraParameters() {
  housingHeight = ZipThick + ZipClearanceZ;
  housingRinner = ZipRad + ZipClearanceR;
  housingRouter = housingRinner + HousingThickness;
}


G4LogicalVolume* CDMSZipConstruction::BuildGeometry() {
  if (verboseLevel)
    G4cout << "CDMSZipConstruction::BuildGeometry()" << G4endl;

  FillExtraParameters();	// Ensure constituents are up to date

  if (verboseLevel > 1) PrintParameters(G4cout);

  G4LogicalVolume* logicalZip = 0;

  if (makeHousingForZip) {
    G4Material* vacuum = CDMSMaterialTable::GetMaterial("G4_Galactic");
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

G4LogicalVolume* CDMSZipConstruction::BuildZip() {
  if (verboseLevel>1) G4cout << "CDMSZipConstruction::BuildZip()" << G4endl;

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


G4LogicalVolume* CDMSZipConstruction::BuildHousing() {
  if (verboseLevel>1) G4cout << "CDMSZipConstruction::BuildHousing()" << G4endl;

  G4Material* mat = GetHousingG4Material();
  if (!mat) return 0;		// Can't build from non-existent material

  G4LogicalVolume* logicalHousing =
    new G4LogicalVolume(BuildHousingShape(), mat, "LogicalHousing", 0,0,0);
  G4VisAttributes* housingAtt = new G4VisAttributes(G4Colour(1.0,0.6,0.4));
  housingAtt->SetForceSolid(true);
  logicalHousing->SetVisAttributes(housingAtt);

  return logicalHousing;
}


G4VSolid* CDMSZipConstruction::BuildHousingShape(G4bool solid) {
  if (verboseLevel>1)
    G4cout << "CDMSZipConstruction::BuildHousingShape()" << G4endl;

  G4double zEnds[2] = { housingHeight/2., -housingHeight/2. };
  G4double rInner[2] = { solid?0:housingRinner, solid?0:housingRinner };
  G4double rOuter[2] = { housingRouter, housingRouter };

  return new G4Polyhedra("HousingShape", 0, 360*deg, HousingSides, 2, 
			 zEnds, rInner, rOuter);
}


G4VSensitiveDetector* CDMSZipConstruction::BuildSensitiveDetector() {
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4VSensitiveDetector* azipSD =
    SDman->FindSensitiveDetector(DetectorName, false);	// Search quietly

  if (!azipSD) {
    if (verboseLevel > 1)
      G4cout << " Creating sensitive detector " << DetectorName << G4endl;

    azipSD = new CDMSZipSD(DetectorName);
    SDman->AddNewDetector(azipSD);
  }

  return azipSD;
}


void CDMSZipConstruction::PrintParameters(std::ostream& os) const {
  os << " CDMSZipConstruction parameters"
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
