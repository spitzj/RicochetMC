////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCTowerSupportConstruction.hh                      //
//                                                                    //
//  Description: Construction of tower support structures in cryostat //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        8 December 2010                                      //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "RMCgeometry/detectors/RMCTowerSupportConstruction.hh"
#include "RMCgeometry/detectors/RMCTowerConstruction.hh"
#include "RMCgeometry/interface/RMCGeometryManager.hh"
#include "RMCgeometry/interface/RMCTowerSupportMessenger.hh"
#include "RMCg4base/RMCMaterialTable.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4Polyhedra.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include <cmath>
#include <iostream>
#include <sstream>

RMCTowerSupportConstruction::RMCTowerSupportConstruction()
  : RMCVDetectorGeometry("TowerSupport"), towerBuilder(0), NStages(6),
    StageHeight(2.54*cm), StageWallThick(6*mm), StageGap(2*mm),
    VesselMountStage(0),
    SpoolLength(5.08*cm), SpoolIR(2.54*cm), SpoolThickness(1.27*cm),
    CTubeIR(2.54*cm), CTubeThick(0.7112*mm),
    SquidLength(7.62*cm), SquidThick(3*mm), SquidWidth(15*mm), 
    SquidStage(2), SquidStageIndex(2),
    FETLength(3.81*cm), FETStage(4), FETStageIndex(4),
    NTowerSides(0), housingRadius(0.), stripThick(0.), stripWidth(0.),
    supportHeight(0.), supportRadius(0.), CTubeLength(0.), StageGapSum(0.),
    cardRadius(0.), stripRadius(0.), squidHeight(0.), fetHeight(0.),
    stripLoLen(0.), stripLoHeight(0.), stripMidLen(0.), stripMidHeight(0.),
    stripHiLen(0.), stripHiHeight(0.),
    messenger(new RMCTowerSupportMessenger(this)) {}


void RMCTowerSupportConstruction::FillExtraParameters() {
  if (verboseLevel>1) 
    G4cout << "RMCTowerSupportConstruction::FillExtraParameters()" << G4endl;

  CopyTowerParameters();
  FillSupportParameters();
  FillReadoutParameters();
}

void RMCTowerSupportConstruction::CopyTowerParameters() {
  if (verboseLevel>1) 
    G4cout << "RMCTowerSupportConstruction::CopyTowerParameters()" << G4endl;

  if (!towerBuilder) {
    towerBuilder = RMCGeometryManager::Instance()->GetTower();
    towerBuilder->SetVerboseLevel(verboseLevel);
  }

  towerBuilder->FillExtraParameters();		// Make values up to date

  NTowerSides      = towerBuilder->GetNTowerSides();
  housingRadius    = towerBuilder->GetHousingRadius();
  stripThick       = towerBuilder->GetStripThick();
  stripWidth       = towerBuilder->GetStripWidth();
}

void RMCTowerSupportConstruction::FillSupportParameters() {
  if (verboseLevel>1) 
    G4cout << "RMCTowerSupportConstruction::FillSupportParameters()" << G4endl;

  G4int i;	// Generic loop variable, used several times below

  StageGaps.resize(NStages,0.);
  for (i=0; i<NStages-1; i++) StageGaps[i] = StageGap;

  // Stage indices below zero count from outside inward
  FETStageIndex = (FETStage>=0) ? FETStage : NStages-1 + FETStage;
  SquidStageIndex = (SquidStage>=0) ? SquidStage : NStages-1 + SquidStage;

  StageGaps[FETStageIndex]   = FETLength;
  StageGaps[SquidStageIndex] = SquidLength;

  StageGapSum = 0.;
  for (i=0; i<NStages; i++) StageGapSum += StageGaps[i];

  CTubeLength   = NStages*StageHeight + StageGapSum;
  supportHeight = SpoolLength + CTubeLength;
  supportRadius = housingRadius+StageWallThick+SquidThick/2.+stripThick+tolerance;

  // Z=0 is the mounting surface, at the innermost cryostat's lid
  supportOffset = 0.;
  for (i=0; i<VesselMountStage; i++) {
    supportOffset += (StageHeight + StageGaps[i]);
  }

  supportBottom = -SpoolLength - supportOffset;
  supportTop    =  CTubeLength - supportOffset;
}


void RMCTowerSupportConstruction::FillReadoutParameters() {
  if (verboseLevel>1) 
    G4cout << "RMCTowerSupportConstruction::FillReadoutParameters()" << G4endl;

  cardRadius  = housingRadius + StageWallThick - SquidThick/2.;
  stripRadius = cardRadius + SquidThick + tolerance;

  // Sanity check: stripRadius + thickness should be less than supportRadius
  if (verboseLevel > 3) {
    G4cout << " outside of strips at radius " << stripRadius+stripThick
	   << "\n mother volume outside radius " << supportRadius << G4endl;
  }

  // FIXME:  Assumes that FET stage is always higher (warmer) than SQUID

  // Kapton readout strip from top of support down to FET (one stage)
  stripHiLen = (NStages-1 - FETStageIndex)*StageHeight;
  stripHiHeight = supportTop - stripHiLen/2.;

  // FET card position at (N-1)th stage gap
  fetHeight = (stripHiHeight-stripHiLen/2.) - FETLength/2.;

  // Kapton readout strip from FET down to SQUID (one stage)
  stripMidLen = (FETStageIndex-SquidStageIndex)*(StageHeight+StageGap) - StageGap;
  stripMidHeight = (fetHeight-FETLength/2.) - stripMidLen/2.;
  
  // SQUID readout at (N-3)th stage gap
  squidHeight = (stripMidHeight-stripMidLen/2.) - SquidLength/2.;

  // Kapton readout strip from SQUID down to bottom of support
  stripLoLen = SquidStageIndex*(StageHeight+StageGap) + StageHeight + SpoolLength;
  stripLoHeight = supportBottom + stripLoLen/2.;

  // Sanity check -- bottom strip should be same height below SQUID as above end
  if (verboseLevel > 3) {
    G4double altStripLoHeight = (squidHeight-SquidLength/2.) - stripLoLen/2.;
    G4cout << " stripLoHeight computed from supportBottom: " << stripLoHeight
	   << "\n stripLoHeight computed from squidHeight:   " << altStripLoHeight
	   << "\n difference: " << stripLoHeight-altStripLoHeight << G4endl;
  }
}


G4LogicalVolume* RMCTowerSupportConstruction::BuildGeometry() {
  if (verboseLevel)
    G4cout << "RMCTowerSupportConstruction::BuildGeometry()" << G4endl;

  FillExtraParameters();

  if (verboseLevel > 1) PrintParameters(G4cout);

  G4double zEnds[2] = { supportBottom, supportTop };
  G4double rInner[2] = { 0., 0. };
  G4double rOuter[2] = { supportRadius, supportRadius };

  G4Polyhedra* supportTowerShape = 
    new G4Polyhedra(GetName(), 0, 360*deg, NTowerSides, 2, 
                    zEnds, rInner, rOuter);

  G4Material* vacuum = RMCMaterialTable::GetMaterial("G4_Galactic");

  G4LogicalVolume* supportTowerLogVolume = 
    new G4LogicalVolume(supportTowerShape, vacuum, GetName(),
                        0, 0, 0);
  supportTowerLogVolume->SetVisAttributes(G4VisAttributes::Invisible);

  // Fill support tower with structures
  AttachCopperSpool(supportTowerLogVolume);
  AttachSupportTube(supportTowerLogVolume);

  for (G4int istage=0; istage<NStages; istage++)
    AttachSupportStage(istage, supportTowerLogVolume);

  for (G4int istrip=0; istrip<NTowerSides; istrip++)
    AttachReadout(istrip, supportTowerLogVolume);

  if (verboseLevel)
    G4cout << "Tower support mass " << supportTowerLogVolume->GetMass()/kg
	   << " kg" << G4endl;

  return supportTowerLogVolume;
}


void 
RMCTowerSupportConstruction::AttachCopperSpool(G4LogicalVolume* theSupport) {
  if (verboseLevel>1)
    G4cout << "RMCTowerSupportConstruction::AttachCopperSpool()" << G4endl;

  G4Material* copper = RMCMaterialTable::GetMaterial("G4_Cu");

  G4Tubs* spoolTube = new G4Tubs("SpoolTube", SpoolIR, SpoolIR+SpoolThickness,
                                 SpoolLength/2., 0, 360*deg);

  G4LogicalVolume* logicalSpool =
    new G4LogicalVolume(spoolTube, copper, "LogicalSpoolTube", 0, 0, 0);

  G4VisAttributes* spoolAtt = new G4VisAttributes(G4Color(1.0,0.6,0.5));
  spoolAtt->SetForceSolid(true);
  logicalSpool->SetVisAttributes(spoolAtt);

  G4double zPos = supportBottom + SpoolLength/2.;	// Lower part of support

  if (verboseLevel>2)
    G4cout << " Placing " << SpoolLength << " (l) x " << SpoolIR+SpoolThickness
	   << " (r) mm spool at Z = " << zPos << G4endl;

  new G4PVPlacement(0, G4ThreeVector(0., 0., zPos), logicalSpool,
		    "PhysicalSpoolTube", theSupport, false, 0);
}

void 
RMCTowerSupportConstruction::AttachSupportTube(G4LogicalVolume* theSupport) {
  if (verboseLevel>1)
    G4cout << "RMCTowerSupportConstruction::AttachSupportTube()" << G4endl;

  G4Material* carbon = RMCMaterialTable::GetMaterial("G4_C");

  G4Tubs* graphiteTube = new G4Tubs("GraphiteTube", CTubeIR,
                                    CTubeIR+CTubeThick,
                                    CTubeLength/2., 0, 360*deg);

  G4LogicalVolume* logicalGraphTube =
    new G4LogicalVolume(graphiteTube, carbon, "LogicalGraphTube", 0, 0, 0);

  G4VisAttributes* graphTubeAtt = new G4VisAttributes(G4Color(0.5,0.5,0.5));
  graphTubeAtt->SetForceSolid(true);
  logicalGraphTube->SetVisAttributes(graphTubeAtt);

  G4double zPos = supportTop - CTubeLength/2.;	// Upper part of support

  if (verboseLevel>2)
    G4cout << " Placing " << CTubeLength << " (l) x " << CTubeIR+CTubeThick
	   << " (r) mm CF tube at Z = " << zPos << G4endl;

  new G4PVPlacement(0, G4ThreeVector(0., 0., zPos), logicalGraphTube,
		    "PhysicalGraphTube", theSupport, false, 0);
}

void
RMCTowerSupportConstruction::AttachSupportStage(G4int stage, 
						  G4LogicalVolume* theSupport) {
  if (verboseLevel>1)
    G4cout << "RMCTowerSupportConstruction::AttachSupportStage " << stage << G4endl;

  G4double stageRinner = housingRadius;
  G4double stageRouter = housingRadius + StageWallThick;

  G4double zPlanes[2] = { -StageHeight/2., StageHeight/2. };
  G4double rInner[2]  = { stageRinner, stageRinner };
  G4double rOuter[2]  = { stageRouter, stageRouter };

  G4Polyhedra* stageShape =
    new G4Polyhedra("StageShape", 0, 360*deg, NTowerSides, 2, zPlanes, rInner, rOuter);

  G4Material* copper = RMCMaterialTable::GetMaterial("G4_Cu");
  G4LogicalVolume* logicalStage =
    new G4LogicalVolume(stageShape, copper, "LogicalStage", 0,0,0);

  G4VisAttributes* stageAtt = new G4VisAttributes(G4Color(1.0,0.6,0.4));
  stageAtt->SetForceSolid(true);
  logicalStage->SetVisAttributes(stageAtt);

  // NOTE:  Coordinate system has support stages starting at Z=0
  G4double zPos = -supportOffset + (stage+0.5)*StageHeight;
  for (G4int i=0; i<stage; i++) zPos += StageGaps[i];

  std::ostringstream name;
  name << "PhysicalStage" << stage;

  if (verboseLevel>2)
    G4cout << " Placing copper stage disk at Z = " << zPos << G4endl;

  new G4PVPlacement(0, G4ThreeVector(0., 0., zPos), logicalStage,
		    name.str().c_str(), theSupport, false, stage);
}


// FIXME:  These must be directly placed to avoid volume overlaps with
//	   stage walls and vessel lids
void
RMCTowerSupportConstruction::AttachReadout(G4int side, 
					     G4LogicalVolume* theSupport) {
  if (verboseLevel>1)
    G4cout << "RMCTowerSupportConstruction::AttachReadout " << side << G4endl;

  G4double readoutPhi = (360*side+180)*deg / NTowerSides;
  G4RotationMatrix readoutRot;
  readoutRot.rotateZ(readoutPhi);

  G4ThreeVector pos;

  // SQUID and FET readout cards
  pos.setRhoPhiZ(cardRadius, readoutPhi, squidHeight);
  new G4PVPlacement(G4Transform3D(readoutRot, pos), BuildSquidCard(),
		    "PhysicalSquidCard", theSupport, false, side);

  pos.setRhoPhiZ(cardRadius, readoutPhi, fetHeight);
  new G4PVPlacement(G4Transform3D(readoutRot, pos), BuildFETCard(),
		    "PhysicalFETCard", theSupport, false, side);

  // Three kapton readout strips, from tower to SQUID, SQUID to FET, and FET out
  pos.setRhoPhiZ(stripRadius, readoutPhi, stripLoHeight);
  new G4PVPlacement(G4Transform3D(readoutRot, pos), BuildSideStrip(stripLoLen),
		    "PhysicalStripLower", theSupport, false, side);

  pos.setRhoPhiZ(stripRadius, readoutPhi, stripMidHeight);
  new G4PVPlacement(G4Transform3D(readoutRot, pos), BuildSideStrip(stripMidLen),
		    "PhysicalStripMiddle", theSupport, false, side);

  pos.setRhoPhiZ(stripRadius, readoutPhi, stripHiHeight);
  new G4PVPlacement(G4Transform3D(readoutRot, pos), BuildSideStrip(stripHiLen),
		    "PhysicalStripUpper", theSupport, false, side);
}

G4LogicalVolume* RMCTowerSupportConstruction::BuildSquidCard() {
  if (verboseLevel>1)
    G4cout << "RMCTowerSupportConstruction::BuildSquidCard()" << G4endl;

  G4Material* carbon = RMCMaterialTable::GetMaterial("G4_C");

  G4Box* squidCard =
    new G4Box("SquidCard", SquidThick/2., SquidWidth/2., SquidLength/2.);
  G4LogicalVolume* logicalSquidCard =
    new G4LogicalVolume(squidCard, carbon, "LogicalSquidCard", 0, 0, 0);
  G4VisAttributes* squidCardAtt = new G4VisAttributes(G4Color(0.4,1.0,0.4));
  squidCardAtt->SetForceSolid(true);
  logicalSquidCard->SetVisAttributes(squidCardAtt);

  return logicalSquidCard;
}

G4LogicalVolume* RMCTowerSupportConstruction::BuildFETCard() {
  if (verboseLevel>1)
    G4cout << "RMCTowerSupportConstruction::BuildFETCard()" << G4endl;

  G4Material* carbon = RMCMaterialTable::GetMaterial("G4_C");

  G4Box* fetCard =
    new G4Box("FETCard", SquidThick/2., SquidWidth/2., FETLength/2.);
  G4LogicalVolume* logicalFETCard =
    new G4LogicalVolume(fetCard, carbon, "LogicalFETCard", 0, 0, 0);
  G4VisAttributes* fetCardAtt = new G4VisAttributes(G4Color(0.4,1.0,0.4));
  fetCardAtt->SetForceSolid(true);
  logicalFETCard->SetVisAttributes(fetCardAtt);

  return logicalFETCard;
}

G4LogicalVolume* 
RMCTowerSupportConstruction::BuildSideStrip(G4double stripLen) {
  if (verboseLevel>1) {
    G4cout << "RMCTowerSupportConstruction::BuildSideStrip " << stripLen << " mm"
	   << G4endl;
  }

  G4Material* carbon = RMCMaterialTable::GetMaterial("G4_C");

  G4double xStrip = stripThick;
  G4double yStrip = stripWidth;
  G4double zStrip = stripLen;
  G4Box* stripBox = new G4Box("StripBox", xStrip/2., yStrip/2., zStrip/2.);

  G4LogicalVolume* logicalStrip = 
    new G4LogicalVolume(stripBox, carbon, "LogicalStrip", 0, 0, 0);

  G4VisAttributes* stripAtt = new G4VisAttributes(G4Color(0.5,0.0,1.0));
  stripAtt->SetForceSolid(true);
  logicalStrip->SetVisAttributes(stripAtt);

  return logicalStrip;
}


void RMCTowerSupportConstruction::PrintParameters(std::ostream& os) const {
  os << "RMCTowerSupportConstruction parameters\n NStages " << NStages
     << "\n StageHeight " << StageHeight << " mm"
     << "\n StageWallThick " << StageWallThick << " mm"
     << "\n StageGap " << StageGap << " mm"
     << "\n VesselMountStage " << VesselMountStage
     << "\n SpoolLength " << SpoolLength << " mm"
     << "\n SpoolIR " << SpoolIR << " mm"
     << "\n SpoolThickness " << SpoolThickness << " mm"
     << "\n CTubeIR " << CTubeIR << " mm"
     << "\n CTubeThick " << CTubeThick << " mm"
     << "\n CTubeLength " << CTubeLength << " mm"
     << "\n SquidLength " << SquidLength << " mm"
     << "\n SquidThick " << SquidThick << " mm"
     << "\n SquidWidth " << SquidWidth << " mm"
     << "\n SquidStage " << SquidStage << " (" << SquidStageIndex << ")"
     << "\n FETLength " << FETLength << " mm"
     << "\n FETStage " << FETStage << " (" << FETStageIndex << ")"
     << "\n supportHeight " << supportHeight << " mm"
     << "\n supportRadius " << supportRadius << " mm"
     << "\n supportOffset " << supportOffset << " mm"
     << "\n supportBottom " << supportBottom << " mm"
     << "\n supportTop " << supportTop << " mm"
     << "\n StageGaps ";
  PrintVectorParameter(os, StageGaps);

  os << "\n StageGapSum " << StageGapSum << " mm"
     << "\n cardRadius " << cardRadius << " mm"
     << "\n stripRadius " << stripRadius << " mm"
     << "\n fetHeight " << fetHeight << " mm"
     << "\n squidHeight " << squidHeight << " mm"
     << "\n stripLoLen " << stripLoLen << " mm"
     << "\n stripLoHeight " << stripLoHeight << " mm"
     << "\n stripMidLen " << stripMidLen << " mm"
     << "\n stripMidHeight " << stripMidHeight << " mm"
     << "\n stripHiLen " << stripHiLen << " mm"
     << "\n stripHiHeight " << stripHiHeight << " mm"
     << std::endl;
}
