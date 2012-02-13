// $Id: CKGDetConstruction.cc,v 1.7 2011/07/22 21:09:17 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//  File:        CKGDetConstruction.cc                                //
//  Description: 100 kg configurable detector construction            //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        31 August 2010                                       //
//                                                                    //
//  20101210  M. Kelsey -- Follow CDMSVDetectorGeometry changes       //
//  20110629  M. Kelsey -- Initialize and use name arg in base class  //
//  20110722  M. Kelsey -- Discard matTable, using singleton instead. //
////////////////////////////////////////////////////////////////////////
 
#include "CDMSgeometry/detectors/CKGDetConstruction.hh"

#include "CDMSgeometry/detectors/CDMSMiniGeomMessenger.hh"
#include "CDMSgeometry/detectors/CKGZipSD.hh"
#include "CDMSg4base/CDMSMaterialTable.hh"
#include "G4Region.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"

 
CKGDetConstruction::CKGDetConstruction()
 :CDMSVDetectorGeometry("Detector")
{
  NTowers = 7;
  NTowerSides = 12;
  NStages = 7;
  ZipRad = 3.81*cm;
  ZipThick = 2.54*cm;
  ZipGap = 0.1*cm;

  ZipAxis1Len = 7.5*cm;
  ZipAxis2Len = 7.2*cm;
  HousingThickness = 0.4*cm;
  Zip_z = 1.0*cm;
  LidClearance = 0.2*cm;

  SpoolLen = 5.08*cm;
  SpoolThickness = 0.15875*cm;  // 1/16" 
  SpoolIR = 2.54*cm;
  CTubeIR = 2.54*cm;
  CTubeThick = 0.07112*cm;

  StageHeight = 2.54*cm;
  StageWallThick = 0.6*cm;
  StageGap = 0.2*cm;
  ZipClearance = 0.1*cm;

  SquidLen = 7.62*cm;
  SquidThick = 0.3*cm;
  SquidWidth = 1.5*cm;

  VesselExtraHeight = 10.0*cm;
  VesselDeltaHeight = 5.0*cm;

  // Wall thicknesses 1/8" and 3/16"
  for (G4int i = 0; i < 6; i++) VesselThick[i] = 0.3175*cm;
  VesselThick[4] = 0.47625*cm;
 
  Vessel0Rad = 35.0*cm;
  VesselDeltaRad = 5.0*cm;
  VesselGap = 15.0*cm;  // currently the vertical distance between
                        // vessels 4 and 5

  // Top and bottom lid thicknesses for outermost vessels
  LidThickness[0] = 0.9525*cm;
  LidThickness[1] = 0.3175*cm;
  LidThickness[2] = 0.9525*cm;

  StemRad = 2.0*cm;
  VacRad = 5.0*cm;

  PipeThick = 0.3*cm;
  PipeBaseLen = 100.0*cm;

  FETLen = 3.81*cm;

  StripThick = 0.3*cm;
  StripWidth = 2.1*cm;

  ExtraSpace = 1.0*cm;
  eps = 0.001*mm;        // tolerance to avoid coincident boundaries

  detectorRegion = 0;
  //  geomMessenger = new CDMSMiniGeomMessenger(this);
}

 
CKGDetConstruction::~CKGDetConstruction()
{
  //  delete geomMessenger;
}


G4LogicalVolume* CKGDetConstruction::BuildGeometry()
{
  // Clear out previously defined regions
  if (!detectorRegion) delete detectorRegion;

  // Detector Box
  G4double xDetBox = 4*m;
  G4double yDetBox = 2*m;
  G4double zDetBox = 2*m;
  G4Material* vacuum = CDMSMaterialTable::GetMaterial("G4_Galactic");
  G4Box* detectorBox = new G4Box("DetBox", xDetBox/2., yDetBox/2., zDetBox/2.);
  G4LogicalVolume* logicalDetectorBox = 
    new G4LogicalVolume(detectorBox, vacuum, GetName(), 0,0,0);

  G4VisAttributes* VisAttDetBox = new G4VisAttributes(G4Color(0, 0, 1));
  VisAttDetBox->SetForceSolid(false);
  logicalDetectorBox->SetVisAttributes(VisAttDetBox);
  logicalDetectorBox->SetVisAttributes(G4VisAttributes::Invisible);

  ///////////////////////////////////////////////////////////////////////////
  //                                                                       //
  //  Build zip tower volumes (vacuum)                                     //
  //                                                                       //
  ///////////////////////////////////////////////////////////////////////////

  G4int NZipsPerTower = NTowerSides;

  G4double housingIR = ZipRad + ZipClearance;
  G4double interTowerHalfDist = housingIR + HousingThickness + ExtraSpace;

  G4double zipHousingHeight = ZipThick + ZipGap;
  G4double zipTowerHeight = G4double(NZipsPerTower)*zipHousingHeight
                          + 2.*(HousingThickness + LidClearance);
  G4double supportTowerHeight = 
    SpoolLen + G4double(NStages)*StageHeight + 4.*StageGap + SquidLen + FETLen;

  G4double zEnds[2];
  zEnds[0] = zipTowerHeight/2.;
  zEnds[1] = -zipTowerHeight/2.;
  G4double rInner[2] = {0.0, 0.0};
  G4double rOuter[2];
  rOuter[0] = interTowerHalfDist - ExtraSpace/2.;
  rOuter[1] = interTowerHalfDist - ExtraSpace/2.;;
  G4Polyhedra* zipTowerShape = new G4Polyhedra("ZipTowerShape", 0, 360*deg,
                                               NTowerSides, 2, 
                                               zEnds, rInner, rOuter);

  std::vector<G4LogicalVolume*> zipTowerLogVolumes(NTowers);
  G4VisAttributes* towerVolAtt = new G4VisAttributes(G4Color(0.0,1.0,0.0));
  towerVolAtt->SetForceWireframe(true);

  for (G4int i = 0; i < NTowers; i++) {
    zipTowerLogVolumes[i] = 
      new G4LogicalVolume(zipTowerShape, vacuum, "LogicalZipTower", 0, 0, 0);
    //    zipTowerLogVolumes[i]->SetVisAttributes(towerVolAtt);
    zipTowerLogVolumes[i]->SetVisAttributes(G4VisAttributes::Invisible);
  }

  G4int stackOption = 0;
  // stackOption = 0 : circular pattern of zips with 1 central zip, tower
  //                 : rotation about z-axis by pi/Nsides
  //                 : valid for NTowers = 7, 9
  // stackOption = 1 : same as 0, but with no tower rotation
  //                 : vertex of each tower polygon pointing to circle center
  //                 : valid for NTowers = 9
  // stackOption = 2 : diamond pattern of zips, no tower rotation
  //                 : valid for NTowers = 9  
 
  std::vector<G4ThreeVector> towerPos(NTowers);
  G4double dd = 2.*interTowerHalfDist;

  G4double root2 = 1.4142136;
  if (stackOption == 2) {
    G4int j;
    for (G4int i = 0; i < NTowers/2; i++) {
      j = 2*i;
      towerPos[j].setX(root2*dd*std::cos(360*deg*j/(NTowers-1) ) );
      towerPos[j].setY(root2*dd*std::sin(360*deg*j/(NTowers-1) ) );
      towerPos[j].setZ(-zipTowerHeight/2.);
      j++;
      towerPos[j].setX(dd*std::cos(360*deg*j/(NTowers-1) ) );
      towerPos[j].setY(dd*std::sin(360*deg*j/(NTowers-1) ) );
      towerPos[j].setZ(-zipTowerHeight/2.);
    }
  } else {
    if (NTowers == 9) dd *= root2;
    for (G4int i = 0; i < NTowers - 1; i++) {
      towerPos[i].setX(dd*std::cos(360*deg*i/(NTowers-1) ) );
      towerPos[i].setY(dd*std::sin(360*deg*i/(NTowers-1) ) );
      towerPos[i].setZ(-zipTowerHeight/2.);
    }
  }

  towerPos[NTowers-1].setX(0.0);
  towerPos[NTowers-1].setY(0.0);
  towerPos[NTowers-1].setZ(-zipTowerHeight/2.);

  G4RotationMatrix* towerRotation = new G4RotationMatrix();
  if (stackOption != 1) towerRotation->rotateZ(180*deg/NTowerSides);

  G4Material* germanium = CDMSMaterialTable::GetMaterial("G4_Ge");
  G4Material* silicon = CDMSMaterialTable::GetMaterial("G4_Si");

  G4LogicalVolume* logicalHousing = BuildZipHousing(zipHousingHeight);

  G4LogicalVolume* logicalTopLid = BuildTowerLid(false);
  G4LogicalVolume* logicalBottomLid = BuildTowerLid(true);

  G4double rStrip;
  G4double xStrip;
  G4double yStrip;
  G4double zStrip;

  // Make support tower
  G4LogicalVolume* supportTowerLogVolume =
    BuildSupportTower(interTowerHalfDist - ExtraSpace/2. + 10*eps,
                      supportTowerHeight);
  supportTowerLogVolume->SetVisAttributes(G4VisAttributes::Invisible);

  //////////////////////////////////////////////////////////////////////
  //                                                                  //
  //  Loop to place elements inside the tower volumes                 //
  //                                                                  //
  //////////////////////////////////////////////////////////////////////

  G4double zPos;
  G4int zipCopyNum;

  G4LogicalVolume* logicalGeZip = BuildZip(germanium);
  G4LogicalVolume* logicalSiZip = BuildZip(silicon);

  // Set up sensitive detector and hit collection for zips
  G4String detectorZipSDName = "Zip";
  G4int collID = -1;
  CKGZipSD* aZipSD = 0;
  G4SDManager* SDman = G4SDManager::GetSDMpointer(); 
  collID = SDman->GetCollectionID(detectorZipSDName);
  if (collID < 0) {
    aZipSD = new CKGZipSD(detectorZipSDName, NTowers, NTowerSides);
    SDman->AddNewDetector(aZipSD);
  }

  logicalGeZip->SetSensitiveDetector(aZipSD);
  logicalSiZip->SetSensitiveDetector(aZipSD);

  G4LogicalVolume* logicalZip;

  for (G4int i = 0; i < NTowers; i++) {

    // Place bottom lid
    zPos = -(zipTowerHeight - HousingThickness - LidClearance)/2.;
    new G4PVPlacement(0, G4ThreeVector(0., 0., zPos), logicalBottomLid,
                      "PhysicalBottomLid", zipTowerLogVolumes[i], false, i);

    // Build zips and housings
    zPos += (HousingThickness + LidClearance)/2. + zipHousingHeight/2.;
    for (G4int j = 0; j < NZipsPerTower; j++) {
      logicalZip = logicalGeZip;
      // Can use 2-D map to decide which zip is placed
      // if (map[i][j] != 0) logicalZip = logicalSiZip;
      zipCopyNum = i*NZipsPerTower + j;
      new G4PVPlacement(0, G4ThreeVector(0., 0., zPos), logicalZip,
                        "PhysicalZip", zipTowerLogVolumes[i], false,
                        zipCopyNum);
      new G4PVPlacement(0, G4ThreeVector(0., 0., zPos), logicalHousing,
                        "PhysicalHousing", zipTowerLogVolumes[i], false,
                        zipCopyNum);
      zPos += zipHousingHeight;
    }

    // Place top lid
    zPos += (HousingThickness + LidClearance)/2. - zipHousingHeight/2.;
    new G4PVPlacement(0, G4ThreeVector(0., 0., zPos), logicalTopLid,
                      "PhysicalTopLid", zipTowerLogVolumes[i], false, i);

    // Put side strips in tower
    G4double stripLen;
    G4RotationMatrix stripRot;
    stripRot.rotateZ(180*deg/NTowerSides);
    G4LogicalVolume* logicalStrip = 0;
    rStrip = housingIR + HousingThickness + StripThick/2. + eps; 

    for (G4int j = 0; j < NTowerSides; j++) {
      stripLen = G4double(j+1)*zipHousingHeight + HousingThickness + LidClearance;
      zStrip = (zipTowerHeight - stripLen)/2.;
      yStrip = rStrip*std::sin(360*deg*(j + 0.5)/NTowerSides);
      xStrip = rStrip*std::cos(360*deg*(j + 0.5)/NTowerSides);

      logicalStrip = BuildSideStrip(stripLen);
      new G4PVPlacement(G4Transform3D(stripRot, G4ThreeVector(xStrip, yStrip, zStrip)),
                        logicalStrip, "PhysicalStripLower",
                        zipTowerLogVolumes[i], false, i*NTowerSides+j);
      stripRot.rotateZ(360*deg/NTowerSides);
    }

    // Put zip tower in lab volume
    new G4PVPlacement(towerRotation, towerPos[i], zipTowerLogVolumes[i],
                      "PhysicalTower", logicalDetectorBox, false, i);

    // Place support tower on top of zip tower
    towerPos[i].setZ(supportTowerHeight/2. + eps);
    new G4PVPlacement(towerRotation, towerPos[i], supportTowerLogVolume,
                      "PhysicalSupportTower", logicalDetectorBox, false, i);
  }

  //////////////////////////////////////////////////////////////////////
  //                                                                  //
  //  Build and place icebox vessels                                  //
  //                                                                  //
  //////////////////////////////////////////////////////////////////////

  G4Material* aluminum = CDMSMaterialTable::GetMaterial("G4_Al");

  G4LogicalVolume* logicalTopOfVessel = 0;
  G4LogicalVolume* logicalBottom = 0;
  G4Tubs* bottomDisk = 0;

  G4LogicalVolume* logicalRestOfVessel = 0;
  G4Tubs* fridgeTube = 0;

  G4LogicalVolume* logicalFridgePipe = 0;
  G4double vesselRad = Vessel0Rad;
  G4double vesselHeight;
  G4VisAttributes* vesselAtt = 0;
  G4double colorFac = 0.0;
  G4double zVessel;
  G4double zVesselTop[6];
  zVesselTop[0] = SpoolLen + 2.*StageHeight;
  zVesselTop[1] = zVesselTop[0] + 2.*StageHeight + StageGap;
  zVesselTop[2] = zVesselTop[1] + SquidLen + StageHeight;
  zVesselTop[3] = zVesselTop[2] + StageHeight + StageGap
                + FETLen + VesselGap;
  zVesselTop[4] = zVesselTop[3] + VesselDeltaHeight;
  zVesselTop[5] = zVesselTop[4] + VesselDeltaHeight;

  G4double stemRad[6];
  stemRad[0] = StemRad;
  stemRad[1] = 1.5*StemRad;
  stemRad[2] = 2.0*StemRad;
  stemRad[3] = 3.0*StemRad;
  stemRad[4] = 4.0*StemRad;
  stemRad[5] = 5.0*StemRad;

  G4double zStemHole[6];
  zStemHole[0] = 1.0*cm + StemRad;
  zStemHole[1] = zStemHole[0] + VesselDeltaHeight;
  zStemHole[2] = zStemHole[1] + VesselDeltaHeight;
  zStemHole[3] = zStemHole[2] + VesselDeltaHeight;
  zStemHole[4] = zStemHole[3] + VesselDeltaHeight;
  zStemHole[5] = zStemHole[4] + VesselDeltaHeight;

  G4double vacRad[3];
  vacRad[0] = VacRad; 
  vacRad[1] = vacRad[0] + 1.0*cm;
  vacRad[2] = vacRad[1] + 1.0*cm;

  G4double zVacHole[3];
  zVacHole[0] = 2.0*cm + VacRad;
  zVacHole[1] = zVacHole[0] + VesselDeltaHeight;
  zVacHole[2] = zVacHole[1] + VesselDeltaHeight;

  G4RotationMatrix* pipeRot = new G4RotationMatrix();
  pipeRot->rotateY(90*deg);
  G4double xPipe;
  G4double zPipe;
  G4double zBottomLid;

  // First do 3 innermost vacuum vessels

  G4Material* copper = CDMSMaterialTable::GetMaterial("G4_Cu");
  for (G4int i = 0; i < 3; i++) {

    // Make the lids of the vacuum vessels

    logicalTopOfVessel = BuildDiskWithHoles(vesselRad, VesselThick[i], 
                                            towerPos, towerRotation);
    vesselAtt = new G4VisAttributes(G4Color(colorFac, colorFac, 1.0));
    vesselAtt->SetForceSolid(true);
    logicalTopOfVessel->SetVisAttributes(vesselAtt);

    new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, zVesselTop[i]), logicalTopOfVessel,
                    "PhysicalTopOfVessel", logicalDetectorBox, false, i);

    // Now make the sides and bottom

    vesselHeight = VesselExtraHeight + zipTowerHeight + zVesselTop[i]
                 + VesselDeltaHeight*G4double(i);
    logicalRestOfVessel = BuildVessel(vesselRad, vesselHeight, VesselThick[i],
                                      stemRad[i], zStemHole[i]);
    logicalRestOfVessel->SetVisAttributes(vesselAtt);
    zVessel =  zVesselTop[i] - (vesselHeight + VesselThick[i])/2.;
    new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, zVessel), logicalRestOfVessel,
                      "PhysicalRestOfVessel", logicalDetectorBox, false, i);
                
    bottomDisk = new G4Tubs("BottomDisk", 0., vesselRad - VesselThick[i],
                            VesselThick[i]/2., 0, 360*deg); 
    logicalBottom = new G4LogicalVolume(bottomDisk, copper, "LogicalBottomLid",
                                        0, 0, 0);
    logicalBottom->SetVisAttributes(vesselAtt);

    zBottomLid = zVessel - vesselHeight/2. + VesselThick[i]/2.;
    new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, zVessel - vesselHeight/2.), logicalBottom,
                      "PhysicalBottomLid", logicalDetectorBox, false, i);

    // Build the fridge pipes
    fridgeTube = new G4Tubs("FridgeTube", stemRad[i] - PipeThick, stemRad[i],
                            (PipeBaseLen - vesselRad)/2., 0, 360*deg);
    logicalFridgePipe =
      new G4LogicalVolume(fridgeTube, copper, "LogicalFridgePipe", 0, 0, 0);
    logicalFridgePipe->SetVisAttributes(vesselAtt);
    xPipe = -vesselRad - (PipeBaseLen - vesselRad)/2. + StemRad;
    zPipe = zVessel - vesselHeight/2. + zStemHole[i];
    new G4PVPlacement(pipeRot, G4ThreeVector(xPipe, 0.0, zPipe), logicalFridgePipe,
                      "PhysicalFridgePipe", logicalDetectorBox, false, i);

    vesselRad += VesselDeltaRad;
    colorFac += 0.3;  
  }

  // Now build 3 outer vacuum vessels

  colorFac = 0.0;
  G4Tubs* thickDisk = 0;
  G4Tubs* vacTube = 0;
  G4LogicalVolume* logicalThickLid = 0;
  G4LogicalVolume* logicalVacPipe = 0;

  for (G4int i = 0; i < 3; i++) {
    vesselRad += VesselDeltaRad;
    vesselHeight = VesselExtraHeight + zipTowerHeight + zVesselTop[i+3]
                 + VesselDeltaHeight*G4double(i+3);

    // Make the lids of the vacuum vessels
    thickDisk = new G4Tubs("ThickDisk", 0.0, vesselRad,
                           LidThickness[i]/2., 0, 360*deg);
    logicalThickLid =
      new G4LogicalVolume(thickDisk, copper, "LogicalThickLid", 0, 0, 0);

    vesselAtt = new G4VisAttributes(G4Color(1.0, 1.0-colorFac, 1.0-colorFac));
    vesselAtt->SetForceSolid(true);
    logicalThickLid->SetVisAttributes(vesselAtt);

    new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, zVesselTop[i+3] - LidThickness[i]/2.),
                      logicalThickLid, "PhysicalThickLid", logicalDetectorBox, false, i);

    zBottomLid =  zVesselTop[i+3] - vesselHeight - 3.*LidThickness[i]/2.;
    new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, zBottomLid),
                      logicalThickLid, "PhysicalThickLid", logicalDetectorBox,
                      false, i+3);

    // Now make the sides
    logicalRestOfVessel =
      BuildVessel(vesselRad, vesselHeight, VesselThick[i+3],
                  stemRad[i+3], zStemHole[i+3], vacRad[i], zVacHole[i]);
    logicalRestOfVessel->SetVisAttributes(vesselAtt);

    zVessel = zVesselTop[i+3] - vesselHeight/2. - LidThickness[i];
    new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, zVessel), logicalRestOfVessel,
                      "PhysicalRestOfOuterVessel", logicalDetectorBox,
                      false, i);

    // Build fridge pipes
    fridgeTube = new G4Tubs("FridgeTube", stemRad[i+3] - PipeThick,
                            stemRad[i+3],
                            (PipeBaseLen - vesselRad)/2., 0, 360*deg);

    logicalFridgePipe =
      new G4LogicalVolume(fridgeTube, copper, "LogicalFridgePipe", 0, 0, 0);
    logicalFridgePipe->SetVisAttributes(vesselAtt);

    xPipe = -vesselRad - (PipeBaseLen - vesselRad)/2. + StemRad;
    zPipe = zVessel - vesselHeight/2. + zStemHole[i+3];
    new G4PVPlacement(pipeRot, G4ThreeVector(xPipe, 0.0, zPipe), logicalFridgePipe,
                      "PhysicalFridgePipe", logicalDetectorBox, false, i+3);

    colorFac += 0.3;
  }

  // Add the vacuum pipe
  vacTube = new G4Tubs("VacTube", vacRad[2] - PipeThick, vacRad[2],
                       (PipeBaseLen - vesselRad)/2., 0, 360*deg);

  logicalVacPipe = new G4LogicalVolume(vacTube, aluminum, "LogicalVacPipe",
                                       0, 0, 0);
  logicalVacPipe->SetVisAttributes(vesselAtt);

  xPipe = vesselRad + (PipeBaseLen - vesselRad)/2. - StemRad;
  zPipe = zVessel + vesselHeight/2. - zVacHole[2];
  new G4PVPlacement(pipeRot, G4ThreeVector(xPipe, 0.0, zPipe), logicalVacPipe,
                      "PhysicalVacPipe", logicalDetectorBox, false, 0);

// Sensor Region (cuts 250 eV) for any particle within the DetBox logical 
// volume
//  detectorRegion = new G4Region(G4String("Detector"));
//  logicalDetectorBox->SetRegion(detectorRegion);
//  detectorRegion->AddRootLogicalVolume(logicalDetectorBox);
//  if(DetectorRegion) {G4cout << "\n###### DetectorRegion exits." 
//     << DetectorRegion << G4endl;}

  return logicalDetectorBox;
}


G4LogicalVolume* CKGDetConstruction::BuildZip(G4Material* mat)
{
  G4Tubs* zipDisk =
    new G4Tubs("ZipDisk", 0, ZipRad, ZipThick/2., 0, 360*deg);
  G4Box* zipCutBox = 
    new G4Box("ZipCutBox", ZipAxis1Len/2., ZipAxis2Len/2., ZipThick/2. + eps);
  G4IntersectionSolid* zipShape =
    new G4IntersectionSolid("ZipShape", zipDisk, zipCutBox);
  G4LogicalVolume* logicalZip =
    new G4LogicalVolume(zipShape, mat, "LogicalZip", 0, 0, 0);

  G4VisAttributes* zipAtt;
  if (mat->GetName() == "G4_Si") {
    zipAtt = new G4VisAttributes(G4Color(0.0,0.0,1.0));
  } else {
    zipAtt = new G4VisAttributes(G4Color(1.0,0.0,0.0));
  }

  zipAtt->SetForceSolid(true);
  logicalZip->SetVisAttributes(zipAtt);

  return logicalZip;
}


G4LogicalVolume* CKGDetConstruction::BuildZipHousing(G4double length)
{
  G4double zEnds[2];
  length -= eps;
  zEnds[0] = length/2.;
  zEnds[1] = -length/2.;

  G4double zipHousingRinner = ZipRad + ZipClearance;
  G4double zipHousingRouter = zipHousingRinner + HousingThickness;
  G4double rInner[2];
  rInner[0] = zipHousingRinner;
  rInner[1] = zipHousingRinner;
  G4double rOuter[2];
  rOuter[0] = zipHousingRouter;
  rOuter[1] = zipHousingRouter;
  G4Polyhedra* housingShape = new G4Polyhedra("HousingShape", 0, 360*deg, NTowerSides, 2, 
                                              zEnds, rInner, rOuter);
  G4Material* copper = CDMSMaterialTable::GetMaterial("G4_Cu");
  G4LogicalVolume* logicalHousing =
     new G4LogicalVolume(housingShape, copper, "LogicalHousing", 0,0,0);
  G4VisAttributes* housingAtt = new G4VisAttributes(G4Color(1.0,0.6,0.4));
  housingAtt->SetForceSolid(true);
  logicalHousing->SetVisAttributes(housingAtt);

  return logicalHousing;
}


G4LogicalVolume*
CKGDetConstruction::BuildSupportTower(G4double radius, G4double length)
{
  G4Material* vacuum = CDMSMaterialTable::GetMaterial("G4_Galactic");

  std::vector<G4double> StageGaps(NStages);
  StageGaps[0] = StageGap;
  StageGaps[1] = StageGap;
  StageGaps[2] = StageGap;
  StageGaps[3] = SquidLen;
  StageGaps[4] = StageGap;
  StageGaps[5] = FETLen;
  StageGaps[6] = 0.;
  G4double StageGapSum = 0.;
  for (G4int i = 0; i < NStages; i++) StageGapSum += StageGaps[i];

  G4double cTubeLen = G4double(NStages)*StageHeight + StageGapSum;

  G4double zEnds[2];
  zEnds[0] = -length/2.; 
  zEnds[1] = length/2.;

  G4double rInner[2]; 
  rInner[0] = 0.0;
  rInner[1] = 0.0;

  G4double rOuter[2];
  rOuter[0] = radius;
  rOuter[1] = radius;

  G4Polyhedra* supportTowerShape = 
    new G4Polyhedra("SupportTowerShape", 0, 360*deg, NTowerSides, 2, 
                    zEnds, rInner, rOuter);

  G4LogicalVolume* supportTowerLogVolume = 
    new G4LogicalVolume(supportTowerShape, vacuum, "LogicalSupportTower",
                        0, 0, 0);
  
  // Copper spool
  G4Tubs* spoolTube = new G4Tubs("SpoolTube", SpoolIR, SpoolIR+SpoolThickness,
                                 SpoolLen/2., 0, 360*deg);
  G4Material* copper = CDMSMaterialTable::GetMaterial("G4_Cu");
  G4LogicalVolume* logicalSpool =
    new G4LogicalVolume(spoolTube, copper, "LogicalSpoolTube", 0, 0, 0);
  G4VisAttributes* spoolAtt = new G4VisAttributes(G4Color(1.0,0.6,0.5));
  spoolAtt->SetForceSolid(true);
  logicalSpool->SetVisAttributes(spoolAtt);

  //  G4double zPos = (SpoolLen - supportTowerHeight)/2.;
  G4double zPos = (SpoolLen - length)/2.;
  new G4PVPlacement(0, G4ThreeVector(0., 0., zPos), logicalSpool,
                    "PhysicalSpoolTube", supportTowerLogVolume, false, 0);

  //
  // Graphite support tube
  //
  G4Tubs* graphiteTube = new G4Tubs("GraphiteTube", CTubeIR,
                                    CTubeIR+CTubeThick,
                                    cTubeLen/2., 0, 360*deg);
  G4Material* carbon = CDMSMaterialTable::GetMaterial("G4_C");
  G4LogicalVolume* logicalGraphTube =
    new G4LogicalVolume(graphiteTube, carbon, "LogicalGraphTube", 0, 0, 0);
  G4VisAttributes* graphTubeAtt = new G4VisAttributes(G4Color(0.5,0.5,0.5));
  graphTubeAtt->SetForceSolid(true);
  logicalGraphTube->SetVisAttributes(graphTubeAtt);

  zPos = (length - cTubeLen)/2.;
  new G4PVPlacement(0, G4ThreeVector(0., 0., zPos), logicalGraphTube,
                    "PhysicalGraphTube", supportTowerLogVolume, false, 0);

  //
  // Stages separating the icebox layers
  //
  G4LogicalVolume* logicalStage = BuildStage();
  zPos = -length/2 + SpoolLen + StageHeight/2.;
  for (G4int i = 0; i < NStages; i++) {
    new G4PVPlacement(0, G4ThreeVector(0., 0., zPos), logicalStage,
                      "PhysicalStage", supportTowerLogVolume, false, i);
    zPos += StageHeight + StageGaps[i];
  }

  //
  // Put in the SQUID and FET cards
  //
  G4Box* squidCard =
    new G4Box("SquidCard", SquidThick/2., SquidWidth/2., SquidLen/2.);
  G4LogicalVolume* logicalSquidCard =
    new G4LogicalVolume(squidCard, carbon, "LogicalSquidCard", 0, 0, 0);
  G4VisAttributes* squidCardAtt = new G4VisAttributes(G4Color(0.4,1.0,0.4));
  squidCardAtt->SetForceSolid(true);
  logicalSquidCard->SetVisAttributes(squidCardAtt);

  G4Box* fetCard =
    new G4Box("FETCard", SquidThick/2., SquidWidth/2., FETLen/2.);
  G4LogicalVolume* logicalFETCard =
    new G4LogicalVolume(fetCard, carbon, "LogicalFETCard", 0, 0, 0);
  logicalFETCard->SetVisAttributes(squidCardAtt);

  G4RotationMatrix squidRot;
  squidRot.rotateZ(180*deg/NTowerSides);
  G4double rSquid = ZipRad + ZipClearance + StageWallThick - SquidThick/2.;
  G4double zSquid = -length/2. + SpoolLen + 4.*StageHeight + 3.*StageGap 
                    + SquidLen/2.; 
  G4double xSquid; 
  G4double ySquid;
  G4double zFET = -length/2. + SpoolLen + 6.*StageHeight + 4.*StageGap
                  + SquidLen + FETLen/2.;

  for (G4int i = 0; i < NTowerSides; i++) {
    xSquid = rSquid*std::cos(360*deg*(i + 0.5)/NTowerSides);
    ySquid = rSquid*std::sin(360*deg*(i + 0.5)/NTowerSides);
    new G4PVPlacement(G4Transform3D(squidRot, G4ThreeVector(xSquid, ySquid, zSquid)),
                      logicalSquidCard, "PhysicalSquidCard",
                      supportTowerLogVolume, false, i);
    new G4PVPlacement(G4Transform3D(squidRot, G4ThreeVector(xSquid, ySquid, zFET)),
                      logicalFETCard, "PhysicalFETCard",
                      supportTowerLogVolume, false, i);

    squidRot.rotateZ(360*deg/NTowerSides);
  }

  G4double stripLenLo = 4.*(StageHeight + StageGap) + SpoolLen;
  G4RotationMatrix stripRot2;
  stripRot2.rotateZ(180*deg/NTowerSides);
  G4LogicalVolume* logicalStripLo = BuildSideStrip(stripLenLo);
  G4double zStripLo = (stripLenLo - length)/2.;

  G4double stripLenMid = 2.*(StageHeight + StageGap);
  G4LogicalVolume* logicalStripMid = BuildSideStrip(stripLenMid);
  G4double zStripMid = length/2. - StageHeight - FETLen - stripLenMid/2.; 

  G4double stripLenHi = StageHeight;
  G4LogicalVolume* logicalStripHi = BuildSideStrip(stripLenHi);
  G4double zStripHi = (length - stripLenHi)/2.;
  G4double rStrip = ZipRad + ZipClearance + StageWallThick + StripThick/2. + eps;
  G4double xStrip;
  G4double yStrip;
  for (G4int i = 0; i < NTowerSides; i++) {
    yStrip = rStrip*std::sin(360*deg*(i + 0.5)/NTowerSides);
    xStrip = rStrip*std::cos(360*deg*(i + 0.5)/NTowerSides);
    new G4PVPlacement(G4Transform3D(stripRot2, G4ThreeVector(xStrip, yStrip, zStripLo)),
                      logicalStripLo, "PhysicalStripLower",
                      supportTowerLogVolume, false, i);
    new G4PVPlacement(G4Transform3D(stripRot2, G4ThreeVector(xStrip, yStrip, zStripMid)),
                      logicalStripMid, "PhysicalStripMiddle",
                      supportTowerLogVolume, false, i);
    new G4PVPlacement(G4Transform3D(stripRot2, G4ThreeVector(xStrip, yStrip, zStripHi)),
                      logicalStripHi, "PhysicalStripUpper",
                      supportTowerLogVolume, false, i);

    stripRot2.rotateZ(360*deg/NTowerSides);
  }

  return supportTowerLogVolume;
}


G4LogicalVolume* CKGDetConstruction::BuildTowerLid(G4bool invert)
{
  G4double zPlanes[4];
  zPlanes[0] = -(HousingThickness + LidClearance)/2.;
  if (invert) {
    zPlanes[1] = (HousingThickness - LidClearance)/2. - eps;
    zPlanes[2] = (HousingThickness - LidClearance)/2.;
  } else {
    zPlanes[1] = -(HousingThickness - LidClearance)/2.;
    zPlanes[2] = -(HousingThickness - LidClearance)/2. + eps;
  }
  zPlanes[3] = -zPlanes[0];

  G4double rInner[4];
  if(invert) {
    rInner[0] = 0.;
    rInner[1] = 0.;
    rInner[2] = ZipRad + ZipClearance;
    rInner[3] = rInner[2];
  } else {
    rInner[0] = ZipRad + ZipClearance;
    rInner[1] = rInner[0];
    rInner[2] = 0.;
    rInner[3] = 0.;
  }

  G4double rOuter[4];
  rOuter[0] = ZipRad + ZipClearance + HousingThickness;
  rOuter[1] = rOuter[0];
  rOuter[2] = rOuter[0];
  rOuter[3] = rOuter[0];

  G4Polyhedra* lidShape = new G4Polyhedra("TowerLidShape", 0, 360*deg, NTowerSides, 4,
                                          zPlanes, rInner, rOuter);
  G4Material* copper = CDMSMaterialTable::GetMaterial("G4_Cu");
  G4LogicalVolume* logicalLid =
    new G4LogicalVolume(lidShape, copper, "LogicalTowerLid", 0,0,0);
  G4VisAttributes* lidAtt = new G4VisAttributes(G4Color(1.0,0.6,0.4));
  lidAtt->SetForceSolid(true);
  logicalLid->SetVisAttributes(lidAtt);

  return logicalLid;
}


G4LogicalVolume* CKGDetConstruction::BuildStage()
{
  G4double zPlanes[2];
  zPlanes[0] = -StageHeight/2.;
  zPlanes[1] = StageHeight/2.;

  G4double rInner[2];
  rInner[0] = ZipRad + ZipClearance;
  rInner[1] = rInner[0];

  G4double rOuter[2];
  rOuter[0] = rInner[0] + StageWallThick;
  rOuter[1] = rOuter[0];

  G4Polyhedra* stageShape =
    new G4Polyhedra("StageShape", 0, 360*deg, NTowerSides, 2, zPlanes, rInner, rOuter);

  G4Material* copper = CDMSMaterialTable::GetMaterial("G4_Cu");
  G4LogicalVolume* logicalStage =
    new G4LogicalVolume(stageShape, copper, "LogicalStage", 0,0,0);

  G4VisAttributes* stageAtt = new G4VisAttributes(G4Color(1.0,0.6,0.4));
  stageAtt->SetForceSolid(true);
  logicalStage->SetVisAttributes(stageAtt);

  return logicalStage;
}


// NOTE:  Vector passed by value because function changes contents
G4LogicalVolume*
CKGDetConstruction::BuildDiskWithHoles(G4double diskRad, G4double diskThick,
                                       std::vector<G4ThreeVector> holePos,
                                       G4RotationMatrix* rotation)
{
  G4double zPlanes[2];
  zPlanes[0] = -diskThick;
  zPlanes[1] = diskThick;

  G4double rInner[2];
  rInner[0] = 0.;
  rInner[1] = 0.;

  G4double rOuter[2];
  // rOuter[0] = (ZipRad + ZipClearance)/std::cos(180*deg/NTowerSides) +
  //             HousingThickness + ExtraSpace/2. + eps;
  rOuter[0] = ZipRad + ZipClearance +
              HousingThickness + ExtraSpace/2. + eps;

  rOuter[1] = rOuter[0];

  G4Polyhedra* hole =
    new G4Polyhedra("Hole", 0, 360*deg, NTowerSides, 2, zPlanes, rInner, rOuter);

  G4Tubs* disk =
    new G4Tubs("Disk", 0.0, diskRad, diskThick/2., 0, 360*deg);

  std::vector<G4VSolid*> subSolidPtrs;
  subSolidPtrs.push_back(disk);
  for (G4int i = 0; i < NTowers; i++) subSolidPtrs.push_back(0);

  for (G4int i = 0; i < NTowers; i++) {
    holePos[i].setZ(0.0);
    subSolidPtrs[i+1] =
      new G4SubtractionSolid("DiskWithHoles", subSolidPtrs[i], hole, rotation,
                             holePos[i]);
  }

  G4Material* copper = CDMSMaterialTable::GetMaterial("G4_Cu");
  G4LogicalVolume* logicalDiskWithHoles =
    new G4LogicalVolume(subSolidPtrs[NTowers], copper, "LogDiskWithHoles", 0, 0, 0);

  return logicalDiskWithHoles;
}


G4LogicalVolume*
CKGDetConstruction::BuildVessel(G4double tubeRad, G4double tubeHeight,
                                G4double tubeThick,
                                G4double holeRad1, G4double zHole1,
                                G4double holeRad2, G4double zHole2)
{
  G4Tubs* tube = new G4Tubs("OuterVesselTube", tubeRad - tubeThick,
                            tubeRad, tubeHeight/2., 0, 360*deg);
  G4Tubs* hole1 = new G4Tubs("OuterStemHole", 0.0, holeRad1 + eps, holeRad1,
                             0, 360*deg);
  G4ThreeVector holePos1(-tubeRad, 0.0, -tubeHeight/2. + zHole1);
  G4RotationMatrix* rotation = new G4RotationMatrix();
  rotation->rotateY(90*deg);

  G4VSolid* tubeWithHoles =
    new G4SubtractionSolid("TubeWithHole", tube, hole1, rotation, holePos1);

  if (holeRad2 > 0.0) {
    G4Tubs* hole2 = new G4Tubs("OuterVacuumHole", 0.0, holeRad2 + eps, holeRad2, 0, 360*deg);
    G4ThreeVector holePos2(tubeRad, 0.0, tubeHeight/2. - zHole2);
    tubeWithHoles =
      new G4SubtractionSolid("TubeWithTwoHoles", tubeWithHoles, hole2, rotation, holePos2);
  }

  G4Material* copper = CDMSMaterialTable::GetMaterial("G4_Cu");
  G4LogicalVolume* logicalVessel =
    new G4LogicalVolume(tubeWithHoles, copper, "LogicalVessel", 0, 0, 0);
                                                            
  return logicalVessel;
}


G4LogicalVolume*
CKGDetConstruction::BuildSideStrip(G4double stripLen)
{
  G4double xStrip = StripThick;
  G4double yStrip = StripWidth;
  G4double zStrip = stripLen;
  G4Box* stripBox = new G4Box("StripBox", xStrip/2., yStrip/2., zStrip/2.);

  G4Material* carbon = CDMSMaterialTable::GetMaterial("G4_C");
  G4LogicalVolume* logicalStrip = 
    new G4LogicalVolume(stripBox, carbon, "LogicalStrip", 0, 0, 0);

  G4VisAttributes* stripAtt = new G4VisAttributes(G4Color(0.5,0.0,1.0));
  stripAtt->SetForceSolid(true);
  logicalStrip->SetVisAttributes(stripAtt);

  return logicalStrip;
}


void CKGDetConstruction::SetZipRad(G4double val)
{
//  Zip_Rout = val;
}

void CKGDetConstruction::SetZipLen(G4double val)
{
//  Zip_z = val;
}

void CKGDetConstruction::SetDetBoxShimThick(G4double val)
{
//  DetBoxShim = val;
}
