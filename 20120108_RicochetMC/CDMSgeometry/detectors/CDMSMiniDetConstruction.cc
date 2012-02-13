// $Id: CDMSMiniDetConstruction.cc,v 1.9 2011/07/22 21:09:17 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSMiniDetConstruction.cc                           //
//  Description: single zip detector construction                     //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        23 June 2010                                         //
//                                                                    //
//  20101026  M. Kelsey -- Add GetMaximumSize() function, follow name //
//		changes in base class                                 //
//  20101210  M. Kelsey -- Follow base class change to dimensions,    //
//		make SetX() functions inline, set location with Z     //
//  20110722  M. Kelsey -- Discard matTable, using singleton instead. //
////////////////////////////////////////////////////////////////////////
 
#include "CDMSgeometry/detectors/CDMSMiniDetConstruction.hh"

#include "CDMSgeometry/detectors/CDMSMiniGeomMessenger.hh"
#include "CDMSgeometry/detectors/CDMSZipSD.hh"
#include "CDMSg4base/CDMSMaterialTable.hh"
#include "G4Region.hh"
#include "G4Tubs.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"

 
CDMSMiniDetConstruction::CDMSMiniDetConstruction()
  : CDMSVDetectorGeometry("Zip"), Zip_Rout(3.81*cm), Zip_z(1.0*cm),
    DetBoxShim(0.1*micrometer), detectorRegion(0),
    geomMessenger(new CDMSMiniGeomMessenger(this)) {
  position.setZ(Zip_z/2.);
}

 
CDMSMiniDetConstruction::~CDMSMiniDetConstruction() {
  delete geomMessenger;
}


G4LogicalVolume* CDMSMiniDetConstruction::BuildGeometry() {
  if (verboseLevel)
    G4cout << " CDMSMiniDetConstruction::BuildGeometry()" << G4endl;

  // Clear out regions is previously defined
  if (!detectorRegion) delete detectorRegion;

  // Detector Box
  G4double DetBox_Rout = Zip_Rout + DetBoxShim;
  G4double DetBox_Z = Zip_z + DetBoxShim;

  if (verboseLevel > 1) {
    G4cout << " virtual cylinder " << DetBox_Rout << " x " << DetBox_Z
	   << " mm" << G4endl;
  }

  G4Tubs* solidDetectorBox = new G4Tubs("DetBox_S", 0.0, DetBox_Rout, DetBox_Z, 0, 360*deg);
  G4Material* vacuum = CDMSMaterialTable::GetMaterial("G4_Galactic");
  G4LogicalVolume* logicalDetectorBox = new G4LogicalVolume(solidDetectorBox, vacuum, "DetBox_L", 0,0,0);
  logicalDetectorBox->SetVisAttributes(G4VisAttributes::Invisible);

  /*
  // Sensor Region (cuts 250 eV) for any particle within the DetBox logical volume
  detectorRegion = new G4Region(G4String("Detector"));
  logicalDetectorBox->SetRegion(detectorRegion);
  detectorRegion->AddRootLogicalVolume(logicalDetectorBox);
  */

  //  if(DetectorRegion) {G4cout << "\n###### DetectorRegion exits." << DetectorRegion << G4endl;}

  // Build the detector components within the DetBox
  if (verboseLevel > 1)
    G4cout << " iZip r " << Zip_Rout << " mm, z " << Zip_z << " mm" << G4endl;

  G4Tubs* solidZip = new G4Tubs("Zip_S", 0.0, Zip_Rout, Zip_z*0.5,  0, 360*deg);
  G4Material* zipMat = CDMSMaterialTable::GetMaterial("G4_Ge");
  G4LogicalVolume* logicalZip = new G4LogicalVolume(solidZip, zipMat, "Zip_L", 0,0,0);
  new G4PVPlacement(0, G4ThreeVector(0,0,0), logicalZip, "Zip_P", logicalDetectorBox, false, 0);

  // Visualization attributes
  G4VisAttributes* VisAttZip = new G4VisAttributes(G4Colour(128/255.,0/255.,0/255.));
  VisAttZip->SetForceSolid(true);
  logicalZip->SetVisAttributes(VisAttZip);

  //
  // Sensitive detectors
  //
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String detectorZipSDName = "Zip01";
  G4VSensitiveDetector* azipSD =
    SDman->FindSensitiveDetector(detectorZipSDName);

  if(!azipSD) {
    if (verboseLevel > 1)
      G4cout << " Creating sensitive detector " << detectorZipSDName << G4endl;

    azipSD = new CDMSZipSD(detectorZipSDName);
    SDman->AddNewDetector(azipSD);
  }

  logicalZip->SetSensitiveDetector(azipSD);

  return logicalDetectorBox;
}
