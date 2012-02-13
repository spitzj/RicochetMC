////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCTowerConstruction.hh                              //
//                                                                    //
//  Description: Construction of ZIP Tower with readout strings       //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Adapted from: Michael Kelsey (SLAC)                  //
//  Date:        14 January 2012                                      //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "RMCgeometry/detectors/RMCTowerConstruction.hh"
#include "RMCgeometry/detectors/RMCZipConstruction.hh"
#include "RMCgeometry/interface/RMCGeometryManager.hh"
#include "RMCgeometry/interface/RMCTowerMessenger.hh"
#include "RMCg4base/RMCMaterialTable.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4Polyhedra.hh"
#include "G4PVPlacement.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include <cmath>

RMCTowerConstruction::RMCTowerConstruction()
  : RMCVDetectorGeometry("ZipTower"), zipBuilder(0),
    NTowerSides(12), NZipsPerTower(12),  LidClearance(2*mm), ExtraSpace(1*cm),
    StripThick(3*mm), StripWidth(21*mm),
    housingRinner(0.), housingRouter(0.),
    housingThickness(0.), housingHeight(0.),
    towerRadius(0.), towerClearanceR(0.), towerHeight(0.), stripRadius(0.),
    messenger(new RMCTowerMessenger(this)) {
  zipBuilder = RMCGeometryManager::Instance()->GetZip();
}

RMCTowerConstruction::~RMCTowerConstruction() {
  delete zipBuilder; zipBuilder=0;
  delete messenger; messenger=0;
}


void RMCTowerConstruction::SetVerboseLevel(G4int verbose) {
  RMCVDetectorGeometry::SetVerboseLevel(verbose);
  if (zipBuilder) zipBuilder->SetVerboseLevel(verbose);
}


void RMCTowerConstruction::FillExtraParameters() {
  if (!zipBuilder) return;		// If no ZIPs, no tower can be built

  SetVerboseLevel(verboseLevel);	// Ensure that children have verbosity

  if (verboseLevel>1)
    G4cout << "RMCTowerConstruction::FillExtraParameters()" << G4endl;

  // Sanity check, number of ZIPs must be <= number of housing sides
  if (NZipsPerTower > NTowerSides) {
    G4cerr << "Tower has too many ZIPs for readout (max " << NZipsPerTower
	   << ").  Reducing count" << G4endl;
    NZipsPerTower = NTowerSides;
  }

  // Push tower configuration down to ZIP to ensure consistency
  zipBuilder->BuildZipWithHousing(false);	// Housing will be handled here
  zipBuilder->SetHousingSides(NTowerSides);
  zipBuilder->FillExtraParameters();

  // Get ZIP spacing and housing parameters
  zipDeltaHeight = zipBuilder->GetZipThick() + zipBuilder->GetZipClearanceZ();

  housingRinner = zipBuilder->GetHousingRinner();
  housingRouter = zipBuilder->GetHousingRouter();

  housingThickness = zipBuilder->GetHousingThickness();
  housingHeight = NZipsPerTower*zipBuilder->GetHousingHeight();

  stripRadius = housingRouter + StripThick/2. + tolerance;
  towerRadius = stripRadius + StripThick/2.;
  towerHeight = housingHeight + 2.*(housingThickness + LidClearance);

  towerClearanceR = towerRadius + ExtraSpace;
}


G4LogicalVolume* RMCTowerConstruction::BuildGeometry() {
  if (verboseLevel)
    G4cout << "RMCTowerConstruction::BuildGeometry()" << G4endl;

  FillExtraParameters();

  if (verboseLevel>1) PrintParameters(G4cout);

  // Empty containment volume, without inter-tower clearance
  G4double zEnds[2]  = { -towerHeight/2., towerHeight/2. };
  G4double rInner[2] = { 0.0, 0.0 };
  G4double rOuter[2] = { towerRadius, towerRadius };

  G4Polyhedra* zipTowerShape =
    new G4Polyhedra(GetName(), 0, 360*deg, NTowerSides,
		    2, zEnds, rInner, rOuter);

  G4VisAttributes* towerVolAtt = new G4VisAttributes(G4Color(0.0,1.0,0.0));
  towerVolAtt->SetForceWireframe(true);

  G4Material* vacuum = RMCMaterialTable::GetMaterial("G4_Galactic");
  G4LogicalVolume* zipTowerLogVolume = 
    new G4LogicalVolume(zipTowerShape, vacuum, GetName(), 0, 0, 0);
  zipTowerLogVolume->SetVisAttributes(G4VisAttributes::Invisible);

  // Start with the housing for the stack of ZIPs
  new G4PVPlacement(G4Transform3D(), BuildHousing(), "PhysicalZipHousing",
		    zipTowerLogVolume, false, copyNumber);

  // Build tower contents from the bottom up; bottom lid first
  if (verboseLevel>1) G4cout << " attaching bottom lid" << G4endl;

  G4double zPos = zEnds[0] + (housingThickness + LidClearance)/2.;

  new G4PVPlacement(0, G4ThreeVector(0.,0.,zPos), BuildLid(lowerLid),
		    "PhysicalBottomLid", zipTowerLogVolume, false, copyNumber);

  // Position ZIPs and readout strips in tower, bottom up
  zPos += (housingThickness + LidClearance + zipDeltaHeight)/2.;

  if (verboseLevel>1)
    G4cout << " stacking " << NZipsPerTower << " ZIPs and housings" << G4endl;

  for (int izip=0; izip<NZipsPerTower; izip++) {
    G4int zipIndex = copyNumber*NZipsPerTower + izip;

    G4LogicalVolume* theZip = zipBuilder->BuildGeometryCopy(zipIndex);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,zPos), theZip,
		      zipBuilder->GetName(),
		      zipTowerLogVolume, false, zipIndex);

    AttachSideStrip(izip, zipTowerLogVolume);

    zPos += zipDeltaHeight;
  }

  // Attach the top lid last
  if (verboseLevel>1) G4cout << " attaching top lid" << G4endl;

  zPos += (housingThickness + LidClearance - zipDeltaHeight)/2.;
  new G4PVPlacement(0, G4ThreeVector(0., 0., zPos), BuildLid(upperLid),
		    "PhysicalTopLid", zipTowerLogVolume, false, copyNumber);

  if (verboseLevel)
    G4cout << "Tower mass " << zipTowerLogVolume->GetMass()/kg << " kg"
	   << G4endl;

  return zipTowerLogVolume;
}

G4LogicalVolume* RMCTowerConstruction::BuildHousing() {
  if (verboseLevel>1) G4cout << "RMCZipConstruction::BuildHousing()" << G4endl;

  G4Material* mat = RMCMaterialTable::GetMaterial(zipBuilder->GetHousingMaterial());
  if (!mat) return 0;		// Can't build from non-existent material

  G4double fullHeight = housingHeight + housingThickness;	// w/lid steps

  G4double zEnds[2]  = { -fullHeight/2., fullHeight/2. };
  G4double rInner[2] = { housingRinner, housingRinner };
  G4double rOuter[2] = { housingRouter, housingRouter };

  G4VSolid* housingShape =
    new G4Polyhedra("HousingShape", 0, 360*deg, NTowerSides, 2, 
		    zEnds, rInner, rOuter);

  G4LogicalVolume* logicalHousing =
    new G4LogicalVolume(housingShape, mat, "LogicalHousing", 0,0,0);
  G4VisAttributes* housingAtt = new G4VisAttributes(G4Colour(1.0,0.6,0.4));
  housingAtt->SetForceSolid(false);
  logicalHousing->SetVisAttributes(housingAtt);

  return logicalHousing;
}

G4LogicalVolume* RMCTowerConstruction::BuildSideStrip(G4int sideIndex) {
  if (verboseLevel>1)
    G4cout << "RMCTowerConstruction::BuildSideStrip " << sideIndex << G4endl;

  G4double xStrip = StripThick;
  G4double yStrip = StripWidth;
  G4double zStrip = GetSideStripLength(sideIndex);
  G4Box* stripBox = new G4Box("StripBox", xStrip/2., yStrip/2., zStrip/2.);

  G4Material* carbon = RMCMaterialTable::GetMaterial("G4_C");
  G4LogicalVolume* logicalStrip = 
    new G4LogicalVolume(stripBox, carbon, "LogicalStrip", 0, 0, 0);

  G4VisAttributes* stripAtt = new G4VisAttributes(G4Color(0.5,0.0,1.0));
  stripAtt->SetForceSolid(true);
  logicalStrip->SetVisAttributes(stripAtt);

  return logicalStrip;
}

void 
RMCTowerConstruction::AttachSideStrip(G4int sideIndex, G4LogicalVolume* theTower) {
  if (verboseLevel>1)
    G4cout << "RMCTowerConstruction::AttachSideStrip " << sideIndex << G4endl;

  G4double stripLen = GetSideStripLength(sideIndex);
  G4double stripPhi = (360*sideIndex+180)*deg / NTowerSides;

  G4ThreeVector stripPos;
  stripPos.setRhoPhiZ(stripRadius, stripPhi, (towerHeight-stripLen)/2.);
  
  G4RotationMatrix stripRot;
  stripRot.rotateZ(stripPhi);
  
  new G4PVPlacement(G4Transform3D(stripRot,stripPos), BuildSideStrip(sideIndex),
		    "PhysicalStripLower", theTower, false, sideIndex);

  return;
}

G4double RMCTowerConstruction::GetSideStripLength(G4int sideIndex) {
  if (verboseLevel>2) {
    G4cout << "RMCTowerConstruction::GetSideStripLength " << sideIndex
	   << G4endl;
  }

  return (sideIndex+1)*zipDeltaHeight + housingThickness + LidClearance;
}

G4LogicalVolume* RMCTowerConstruction::BuildLid(TowerLidSide side) {
  if (verboseLevel>1) 
    G4cout << "RMCTowerConstruction::BuildLid " << side << G4endl;

  G4bool invert = (side == lowerLid);	// Bottom lid is upside down version of top

  // Lid fits into housing opening, with an outer lip which goes over end

  G4double zPlanes[4];
  G4double zstep = (housingThickness - LidClearance)/2.;

  zPlanes[0] = -(housingThickness + LidClearance)/2.;
  if (invert) {
    zPlanes[1] = zstep - tolerance;
    zPlanes[2] = zstep;
  } else {
    zPlanes[1] = -zstep;
    zPlanes[2] = -zstep + tolerance;
  }
  zPlanes[3] = -zPlanes[0];

  G4double rInner[4] = { 0., 0., 0., 0.};
  G4double rOuter[4];
  if (invert) {
    rOuter[0] = rOuter[1] = housingRouter;
    rOuter[2] = rOuter[3] = housingRinner;
  } else {
    rOuter[0] = rOuter[1] = housingRinner;
    rOuter[2] = rOuter[3] = housingRouter;
  }


  G4Polyhedra* lidShape =
    new G4Polyhedra("TowerLidShape", 0, 360*deg, NTowerSides,
		    4, zPlanes, rInner, rOuter);
  G4Material* copper = RMCMaterialTable::GetMaterial("G4_Cu");
  G4LogicalVolume* logicalLid =
    new G4LogicalVolume(lidShape, copper, "LogicalTowerLid", 0,0,0);
  G4VisAttributes* lidAtt = new G4VisAttributes(G4Color(0.5,0.3,0.2));
  lidAtt->SetForceSolid(true);
  logicalLid->SetVisAttributes(lidAtt);

  return logicalLid;
}


void RMCTowerConstruction::PrintParameters(std::ostream& os) const {
  os << "RMCTowerConstruction parameters\n NTowerSides " << NTowerSides
     << "\n LidClearance " << LidClearance << " mm"
     << "\n ExtraSpace " << ExtraSpace << " mm"
     << "\n StripThick " << StripThick << " mm"
     << "\n StripWidth " << StripWidth << " mm"
     << "\n NZipsPerTower " << NZipsPerTower
     << "\n zipDeltaHeight " << zipDeltaHeight << " mm"
     << "\n housingRinner " << housingRinner << " mm"
     << "\n housingRouter " << housingRouter << " mm"
     << "\n housingThickness " << housingThickness << " mm"
     << "\n housingHeight " << housingHeight << " mm"
     << "\n towerClearanceR " << towerClearanceR << " mm"
     << "\n towerHeight " << towerHeight << " mm"
     << "\n stripRadius " << stripRadius << " mm"
     << std::endl;
}
