// $Id: CDMSVesselConstruction.cc,v 1.29 2011/07/22 21:09:17 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSVesselConstruction.cc                            //
//                                                                    //
//  Description: Construction of vacuum cryostat with tower supports  //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        16 November 2010                                     //
//                                                                    //
//  20101130  M. Kelsey -- Add reporting.                             //
//  20101208  M. Kelsey -- Move tower support structure to new class  //
//  20101209  M. Kelsey -- Pass verbosity to children, allow for no   //
//		towers (empty cryostat), fix overlapping volume bugs, //
//		build mother volume as envelope of vessel and pipes   //
//  20101210  M. Kelsey -- Position tower support at inner vessel lid //
//		Fix calculations of vessel Z coordinates (towers @ 0) //
//  20101211  M. Kelsey -- Drop unused parameters, use new pattern    //
//		offset interface.  Do not subtract holes for towers.  //
//  20101215  M. Kelsey -- Add number of vacuum vessels, multi pipes; //
//		decouple tower-support stages from cryostat count.    //
//		Subtract chunk from end of envelope to help shielding //
//  20101218  M. Kelsey -- Convert "extra param" arrays to vectors.   //
//  20110105  M. Kelsey -- Add function for cylindrical outer radius  //
//  20110105  M. Kelsey -- Compute offset of Z center from origin;    //
//		replace use of "[outer]" with STL .back(); replace    //
//		PipeBaseLen with PipeOutsideLen, and adjust meaning;  //
//		add offsets from ends for cryo and vacuum pipes.      //
//  20110107  M. Kelsey -- Check nTowers>0 before building supports.  //
//  20110215  M. Kelsey -- Report mass at end of construction.        //
//  20110629  M. Kelsey -- Initialize and use name arg in base class  //
//  20110722  M. Kelsey -- Discard matTable, using singleton instead. //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSgeometry/detectors/CDMSVesselConstruction.hh"
#include "CDMSgeometry/detectors/CDMSTowerConstruction.hh"
#include "CDMSgeometry/detectors/CDMSTowerSupportConstruction.hh"
#include "CDMSgeometry/detectors/CDMSTowerPattern.hh"
#include "CDMSgeometry/interface/CDMSGeometryManager.hh"
#include "CDMSgeometry/interface/CDMSVesselMessenger.hh"
#include "CDMSg4base/CDMSMaterialTable.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4Polyhedra.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4Transform3D.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4VisAttributes.hh"
#include <cmath>
#include <iostream>


CDMSVesselConstruction::CDMSVesselConstruction()
  : CDMSVDetectorGeometry("Cryostat"),
    theManager(CDMSGeometryManager::Instance()),
    towerPattern(0), towerBuilder(0), supportBuilder(0),
    NStages(6), NVacuum(3), NTowers(0), Vessel0Rad(35*cm), VesselDeltaRad(5*cm),
    VesselDeltaHeight(5*cm), VesselExtraHeight(10*cm), VesselGap(15*cm),
    StemInnerRad(2*cm), VacInnerRad(5*cm), StemZInset(1*cm), VacZInset(2*cm),
    PipeThick(3*mm), PipeOutsideLen(40*cm),
    NTowerSides(0), ExtraSpace(1*cm), housingRadius(0.), housingThickness(0.),
    towerHeight(0.), towerRadius(0.), supportNStages(0),
    supportMountStage(0), supportSquidStage(0), supportFETStage(0),
    supportRadius(0.), supportBottom(0.), spoolLength(0.), stageHeight(0.),
    messenger(0) {
  // Vacuum vessel parameters are arrays, and must be set accordingly
  // Wall thicknesses 1/8" and 3/16"

  VesselThick.resize(6,0.);
  for (G4int i = 0; i < 6; i++) VesselThick[i] = 0.3175*cm;
  VesselThick[4] = 0.47625*cm;

  LidThickness.resize(3,0.);
  LidThickness[0] = 0.9525*cm;
  LidThickness[1] = 0.3175*cm;
  LidThickness[2] = 0.9525*cm;

  // Alway have tower layout pattern, even if zero
  towerPattern = theManager->GetTowerPattern();

  // Create Messenger after the parameter lists have been filled above
  messenger = new CDMSVesselMessenger(this);
}

CDMSVesselConstruction::~CDMSVesselConstruction() {
  delete messenger;
}


void CDMSVesselConstruction::SetVerboseLevel(G4int verbose) {
  CDMSVDetectorGeometry::SetVerboseLevel(verbose);
  if (towerBuilder)   towerBuilder->SetVerboseLevel(verbose);
  if (towerPattern)   towerPattern->SetVerboseLevel(verbose);
  if (supportBuilder) supportBuilder->SetVerboseLevel(verbose);
}
  

// Fill vector parameters (usually from commands)
void 
CDMSVesselConstruction::SetVesselThick(const std::vector<G4double>& value) {
  VesselThick = value;
  if (value.size() < (size_t)NStages) 
    for (G4int i=value.size(); i<NStages; i++) VesselThick.push_back(0.);
}

void  
CDMSVesselConstruction::SetLidThickness(const std::vector<G4double>& value) {
  LidThickness = value;
  if (value.size() < (size_t)NVacuum) 
    for (G4int i=value.size(); i<NVacuum; i++) LidThickness.push_back(0.);
}


// Fill computed and copied quantities

void CDMSVesselConstruction::FillExtraParameters() {
  SetVerboseLevel(verboseLevel);	// Ensure that children have verbosity

  //if (verboseLevel>1) 
    G4cout << "CDMSVesselConstruction::FillExtraParameters()" << G4endl;

  UpdateTowerPattern();	
  CopyTowerParameters();
  CopySupportParameters();
  FillVesselParameters();

  // Shift so that outermost volume is centered in external coordinates
  position.setZ(-zVessel.back());
}

void CDMSVesselConstruction::CopyTowerParameters() {
  if (!towerBuilder) return;

  if (verboseLevel>1) 
    G4cout << "CDMSVesselConstruction::CopyTowerParameters()" << G4endl;

  towerBuilder->FillExtraParameters();		// Make values up to date

  NTowerSides      = towerBuilder->GetNTowerSides();
  ExtraSpace       = towerBuilder->GetExtraSpace();
  housingRadius    = towerBuilder->GetHousingRadius();
  housingThickness = towerBuilder->GetHousingThickness();
  towerHeight      = towerBuilder->GetTowerHeight();
  towerRadius	   = towerBuilder->GetTowerRadius();
}

void CDMSVesselConstruction::CopySupportParameters() {
  if (!supportBuilder) {
    stageGaps.clear();				// No extra gaps
    stageGaps.push_back(0.);
    supportNStages = supportMountStage = 0;
    supportSquidStage = supportFETStage = 0;
    supportRadius = supportBottom = spoolLength = stageHeight = 0.;
    return;
  }

  if (verboseLevel>1) 
    G4cout << "CDMSVesselConstruction::CopySupportParameters()" << G4endl;

  supportBuilder->FillExtraParameters();

  supportNStages    = supportBuilder->GetNStages();
  supportMountStage = supportBuilder->GetVesselMountStage();
  supportSquidStage = supportBuilder->GetSquidStage();
  supportFETStage   = supportBuilder->GetFETStage();
  supportRadius     = supportBuilder->GetSupportRadius();
  supportBottom     = supportBuilder->GetSupportBottom();
  spoolLength       = supportBuilder->GetSpoolLength();
  stageHeight       = supportBuilder->GetStageHeight();
  stageGaps         = supportBuilder->GetStageGaps();
}

void CDMSVesselConstruction::FillVesselParameters() {
  if (verboseLevel>1) 
    G4cout << "CDMSVesselConstruction::FillVesselParameters()" << G4endl;

  // Pre-allocate all "arrays" to allow initialization by subscript
  vesselRadius.resize(NStages,0.);
  vesselHeight.resize(NStages,0.);
  zVesselTop.resize(NStages,0.);
  zVessel.resize(NStages,0.);
  zVesselBottom.resize(NStages,0.);
  pipeLen.resize(NStages,0.);
  stemRad.resize(NStages,0.);
  zStemHole.resize(NStages,0.);
  vacRad.resize(NVacuum,0.);
  zVacHole.resize(NVacuum,0.);
  vesselVisAtt.resize(NStages,0);
  
  ComputeVesselTopCoordinates();
  ComputeVesselDimensions();
  ComputePipeCoordinates();

  // Special:  Each vacuum vessel has particular color settings
  G4double color = 0.3;
  for (G4int i=0; i<NStages; i++) {
    if (i<NStages-NVacuum) {		// Inner vessels
      vesselVisAtt[i] = new G4VisAttributes(G4Color(i*color, i*color, 1.0));
    } else {				// Outer vessels
      vesselVisAtt[i] = new G4VisAttributes(G4Color(1.0, 1.0-i*color, 1.0-i*color));
    }
    vesselVisAtt[i]->SetForceSolid(true);
  }
}

void CDMSVesselConstruction::ComputeVesselTopCoordinates() {
  if (verboseLevel>1) 
    G4cout << "CDMSVesselConstruction::ComputeVesselTopCoordinates()" << G4endl;
  
  G4int i;	// Generic loop variable, used several times below
  
  // Cryostat vessels are nested with the innermost to have towers at Z=0
  // The tower support has its internal Z=0 at the lid of the innermost
  // The SQUID readout card is inside the last cryogen-only vessel
  // The FET readout card is inside the first (inner) vacuum-pumped vessel
  
  G4int FETVessel   = NStages - NVacuum;
  G4int SquidVessel = FETVessel - 1;

  zVesselTop[0] = -supportBottom + towerHeight/2.;

  // Number of support stages and vessels up to SQUID card
  G4int deltaStages  = supportSquidStage; //*** - supportMountStage;
  G4int deltaVessels = SquidVessel - 1;
  G4int stagesPerVessel = deltaStages / deltaVessels;

  if (verboseLevel > 2) {
    G4cout << " computing zVesselTop for vessels 1 to " << SquidVessel-1
	   << "\n deltaStages " << deltaStages 
	   << " deltaVessels " << deltaVessels
	   << " stagesPerVessel " << stagesPerVessel << G4endl;
  }

  for (i=1; i<SquidVessel; deltaStages-=stagesPerVessel, i++) {
    zVesselTop[i] = zVesselTop[i-1];
    G4int jOffset = (i-1)*stagesPerVessel; //*** + supportMountStage;

    if (verboseLevel > 2) 
      G4cout << " computing zVesselTop[" << i << "] with jOffset "
	     << jOffset << G4endl;

    for (G4int j=0; j<stagesPerVessel; j++) {
      zVesselTop[i] += stageHeight + stageGaps[j+jOffset];
    }
  }

  // SQUID card vessel may include support stages
  G4int nextStage = supportSquidStage;
  if (deltaStages > 0) nextStage -= deltaStages;

  if (verboseLevel > 2) {
    G4cout << " computing SQUID zVesselTop[" << SquidVessel << "]"
	   << " with stages " << nextStage << " to " << supportSquidStage
	   << G4endl;
  }

  G4cout << stagesPerVessel << G4endl;
  G4cout << stageHeight << G4endl;;
  G4cout << stageGaps[0] << G4endl;
  G4cout << stageGaps[1] << G4endl;
  zVesselTop[SquidVessel] = zVesselTop[SquidVessel-1];
  for (i=nextStage; i<=supportSquidStage; i++) {
    zVesselTop[SquidVessel] += stageHeight + stageGaps[i];
  }

  // FET card vessel may include support stages
  nextStage = supportSquidStage + 1;

  if (verboseLevel > 2) {
    G4cout << " computing FET zVesselTop[" << FETVessel << "]"
	   << " with stages " << nextStage << " to " << supportNStages-1
	   << G4endl;
  }

  zVesselTop[FETVessel] = zVesselTop[FETVessel-1] + VesselGap;
  for (i=nextStage; i<supportNStages; i++) {
    zVesselTop[FETVessel] += stageHeight + stageGaps[i];
  }

  // Outer vessels are empty
  for (i=FETVessel+1; i<NStages; i++) {
    if (verboseLevel > 2)
      G4cout << " computing zVesselTop[" << i << "]" << G4endl;
    zVesselTop[i] = zVesselTop[i-1] + VesselDeltaHeight;
  }
}

void CDMSVesselConstruction::ComputeVesselDimensions() {
  if (verboseLevel>1) 
    G4cout << "CDMSVesselConstruction::ComputeVesselDimensions()" << G4endl;

  // Other vessel dimensions are relative to top
  for (G4int i=0; i<NStages; i++) {
    vesselRadius[i] = Vessel0Rad + i*VesselDeltaRad;

    vesselHeight[i] = (VesselExtraHeight + towerHeight/2. + i*VesselDeltaHeight
		       + zVesselTop[i]);	// Top[0] includes half-tower

    G4double halfFull = (vesselHeight[i] + VesselThick[i])/2.;
    zVessel[i] = zVesselTop[i] - halfFull;
    zVesselBottom[i] = zVessel[i] - halfFull;
  }
}

void CDMSVesselConstruction::ComputePipeCoordinates() {
  if (verboseLevel>1) 
    G4cout << "CDMSVesselConstruction::ComputeVesselDimensions()" << G4endl;

  // Nested cryogen pipes -- FIXME:  This is hardwired for six stages!
  stemRad[0] = StemInnerRad;
  stemRad[1] = 1.5*StemInnerRad;
  stemRad[2] = 2.0*StemInnerRad;
  stemRad[3] = 3.0*StemInnerRad;
  stemRad[4] = 4.0*StemInnerRad;
  stemRad[5] = 5.0*StemInnerRad;

  // Cryogen pipes are near the bottom of vessels, vacuum pipes near the top
  G4double zoffset = 0.;

  G4int i;
  for (i=0; i<NStages; i++) {
    zoffset = StemInnerRad + StemZInset + i*VesselDeltaHeight;
    zStemHole[i] = zVesselBottom[i] + zoffset;
    pipeLen[i]   = PipeOutsideLen + (vesselRadius.back()-vesselRadius[i]);
  }

  for (i=0; i<NVacuum; i++) {
    G4int stage = NStages-NVacuum + i;
    zoffset = VacInnerRad + VacZInset + i*VesselDeltaHeight;
    zVacHole[i] = zVesselTop[stage] - zoffset;
    vacRad[i]   = VacInnerRad + i*(1.0*cm);
  }
}


void CDMSVesselConstruction::UpdateTowerPattern() {
  if (!towerPattern) return;		// No pattern means no towers to use

  if (verboseLevel>1) 
    G4cout << "CDMSVesselConstruction::UpdateTowerPattern()" << G4endl;

  NTowers = towerPattern->GetNumberOfTowers();
  if (NTowers < 1) {			// No towers to be built
    towerBuilder = 0;
    supportBuilder = 0;
    return;
  }

  // Get builders for towers and supports
  if (!towerBuilder)   towerBuilder   = theManager->GetTower();
  if (!supportBuilder) supportBuilder = theManager->GetTowerSupport();

  if (towerBuilder) towerBuilder->FillExtraParameters();
  towerPattern->GeneratePattern();
}


G4LogicalVolume* CDMSVesselConstruction::BuildGeometry() {
  if (verboseLevel)
    G4cout << "CDMSVesselConstruction::BuildGeometry()" << G4endl;

  FillExtraParameters();

  if (verboseLevel > 1) PrintParameters(G4cout);

  // Psuedo "world volume" in which to place all detector structures
  G4LogicalVolume* logicalIcebox = BuildEnvelope();

  // Now start placing everything (mostly at the origin)
  for (G4int i=0; i<NStages; i++) PlaceVesselStage(i, logicalIcebox);

  if (verboseLevel)
    G4cout << "Cryostat mass " << logicalIcebox->GetMass()/kg << " kg"
	   << " without tower supports" << G4endl;

  PlaceTowerSupports(logicalIcebox);

  if (verboseLevel)
    G4cout << "Cryostat mass " << logicalIcebox->GetMass()/kg << " kg"
	   << G4endl;

  return logicalIcebox;
}

G4LogicalVolume* CDMSVesselConstruction::BuildEnvelope() {
  if (verboseLevel > 1)
    G4cout << "CDMSVesselConstruction::BuildEnvelope()" << G4endl;

  G4double clearance = VesselThick.back() + tolerance;
  G4double envBelow  = -zVesselBottom.back() + clearance;
  G4double envAbove  = zVesselTop.back() + clearance;
  G4double envRadius = vesselRadius.back() + clearance;

  G4double envHalfHeight = (envBelow>envAbove) ? envBelow : envAbove;

  if (verboseLevel > 2) {
    G4cout << " envBelow " << envBelow << " mm"
	   << "\n envAbove " << envAbove << " mm"
	   << "\n envRadius " << envRadius << " mm"
	   << "\n envHalfHeight " << envHalfHeight << " mm"
	   << G4endl;
  }

  G4VSolid* envelope = new G4Tubs("VesselEnvelope", 0., envRadius,
				  envHalfHeight, 0., 360*deg);

  // Get Z coords right by removing "extra" cylinder from one end of envelope
  G4double excessHeight = std::abs(envAbove-envBelow);

  if (excessHeight > 0.) {
    G4double excessPosition = 0.;
    if (envBelow>envAbove) excessPosition = envAbove + excessHeight/2.;
    else excessPosition = -envBelow - excessHeight/2.;
    // NOTE:  Position computed above is below bottom of envelope, but works!

    excessHeight += clearance;
    excessPosition += (excessPosition>0) ? clearance/2. : -clearance/2.;

    if (verboseLevel > 2) {
      G4cout << " Removing extra cylindrical piece at "
	     << (excessPosition>=0.?"top":"bottom") << " of shape"
	     << "\n excessHeight " << excessHeight << " mm"
	     << "\n excessPosition " << excessPosition << " mm"
	     << G4endl;
    }

    G4Tubs* excessVessel = new G4Tubs("VesselExcess", 0., envRadius+clearance,
				      excessHeight/2., 0., 360*deg);
    envelope =
      new G4SubtractionSolid("VesselEnvelope", envelope, excessVessel,
			     0, G4ThreeVector(0.,0.,excessPosition));
  }

  // Combine cylindrical vessel with pipe stems
  if (verboseLevel > 2)
    G4cout << " Joining vacuum and cryo pipes to vessel shape" << G4endl;

  G4double pipeLength = PipeOutsideLen + tolerance;

  G4Tubs* vacTube = new G4Tubs("VacEnvelope", 0., vacRad.back()+tolerance,
			       pipeLength/2., 0, 360*deg);
  G4Tubs* cryoTube = new G4Tubs("CryoEnvelope", 0., stemRad.back()+tolerance,
				pipeLength/2., 0, 360*deg);

  G4int outer = NStages-1;
  envelope = new G4UnionSolid("VesselPlusVac", envelope, vacTube,
			      GetVacuumPipePosition(outer));

  envelope = new G4UnionSolid(GetName(), envelope, cryoTube,
			      GetCryoPipePosition(outer));

  // Create final mother volume for construction and export
  G4Material* vacuum = CDMSMaterialTable::GetMaterial("G4_Galactic");
  G4LogicalVolume* logicalIcebox =
    new G4LogicalVolume(envelope, vacuum, GetName(), 0, 0, 0);

  return logicalIcebox;
}


// Remove cylindrical spaces for each tower from mother's physical volume

G4VSolid* CDMSVesselConstruction::SubtractTowerPositions(G4VSolid* icebox) {
  if (!towerPattern) return icebox;	// Avoid unnecessary work

  if (verboseLevel>1)
    G4cout << "CDMSVesselConstruction::SubtractTowerPositions()" << G4endl;

  G4double zEnds[2] = { -towerHeight/2., towerHeight/2. };
  G4double rInner[2] = { 0., 0. };
  G4double rOuter[2] = { towerRadius, towerRadius };
  
  G4Polyhedra* towerShape = 
    new G4Polyhedra("TowerShape", 0, 360*deg, NTowerSides, 2, zEnds,
		    rInner, rOuter);

  towerPattern->SetCenter();			// Towers at origin
  for (G4int i=0; i<towerPattern->GetNumberOfTowers(); i++) {
    icebox = new G4SubtractionSolid(GetName(), icebox, towerShape,
				    towerPattern->GetTransform3D(i));
  }

  return icebox;
}


void CDMSVesselConstruction::PlaceTowerSupports(G4LogicalVolume* world) {
  // Bare cryostat won't have any towers, so no supports are required
  G4int nTowers = towerPattern ? towerPattern->GetNumberOfTowers() : 0;
  if (nTowers < 1) return;

  if (verboseLevel>1)
    G4cout << "CDMSVesselConstruction::PlaceTowerSupports()" << G4endl;

  // Reference plane (spool/CF boundary) at lid of innermost cryostat
  towerPattern->SetCenter(0.,0.,zVesselTop[0]);

  // Make single copy of tower support struture, and place it everywhere
  G4LogicalVolume* theSupport = supportBuilder->BuildGeometry();

  for (G4int i=0; i<nTowers; i++) {
    if (verboseLevel>2) {
      G4cout << " placing support #" << i << " of " << nTowers << " at "
	     << towerPattern->GetPosition(i) << " mm" << G4endl;
    }

    new G4PVPlacement(towerPattern->GetTransform3D(i), theSupport,
		      "PhysicalSupportTower", world, false, i);
  }
}

// NOTE:  Components are directly placed at top level to avoid volume overlaps
void 
CDMSVesselConstruction::PlaceVesselStage(G4int stage, G4LogicalVolume* world) {
  if (verboseLevel>1)
    G4cout << "CDMSVesselConstruction::PlaceVesselStage " << stage << G4endl;

  G4ThreeVector pos;

  pos.setZ(zVesselTop[stage]);
  new G4PVPlacement(0, pos, BuildVesselTop(stage), "PhysicalTopOfVessel",
		    world, false, stage);
  
  pos.setZ(zVessel[stage]);
  new G4PVPlacement(0, pos, BuildVesselSide(stage), "PhysicalSideOfVessel",
		    world, false, stage);

  pos.setZ(zVesselBottom[stage]);
  new G4PVPlacement(0, pos, BuildVesselBottom(stage), "PhysicalBottomOfVessel",
		    world, false, stage);

  new G4PVPlacement(GetCryoPipePosition(stage), BuildCryoPipe(stage),
		    "PhysicalFridgePipe", world, false, stage);

  G4LogicalVolume* vacuumPipe = BuildVacuumPipe(stage);		// May be null!
  if (vacuumPipe) new G4PVPlacement(GetVacuumPipePosition(stage), vacuumPipe,
				    "PhysicalVacuumPipe", world, false, stage);
}


G4LogicalVolume* CDMSVesselConstruction::BuildVesselTop(G4int stage) {
  if (verboseLevel>1)
    G4cout << "CDMSVesselConstruction::BuildVesselTop " << stage << G4endl;

  G4double zPlanes[2] = { -VesselThick[stage]/2., VesselThick[stage]/2. };
  G4double rInner[2] = { 0., 0. };
  G4double rOuter[2] = { supportRadius+tolerance, supportRadius+tolerance };

  // Start with a full disk, and take away the holes
  G4VSolid* solidTop = new G4Tubs("Disk", 0.0, vesselRadius[stage],
				  VesselThick[stage]/2., 0, 360*deg);

  // Bare cryostat won't have any towers, so no holes are required
  G4int nTowers = towerPattern ? towerPattern->GetNumberOfTowers() : 0;

  // Only inner cryostats need holes for tower supports
  if (nTowers > 0 && stage < (NStages-NVacuum)) {
    if (verboseLevel>2)
      G4cout << " removing " << nTowers << " holes for tower supports"
	     << G4endl;

    G4Polyhedra* hole = new G4Polyhedra("Hole", 0, 360*deg, NTowerSides,
					2, zPlanes, rInner, rOuter);

    towerPattern->SetCenter();
    for (G4int i=0; i<nTowers; i++) {
      solidTop = new G4SubtractionSolid("DiskWithHoles", solidTop, hole,
					towerPattern->GetTransform3D(i));
    }
  }

  // Now create the actual lid with perforations
  G4Material* copper = CDMSMaterialTable::GetMaterial("G4_Cu");
  G4LogicalVolume* logicalTop =
    new G4LogicalVolume(solidTop, copper, "LogDiskWithHoles", 0, 0, 0);
  logicalTop->SetVisAttributes(vesselVisAtt[stage]);
  
  return logicalTop;
}

G4LogicalVolume* CDMSVesselConstruction::BuildVesselSide(G4int stage) {
  if (verboseLevel>1)
    G4cout << "CDMSVesselConstruction::BuildVesselSide " << stage << G4endl;

  G4double vesselIR = vesselRadius[stage] - VesselThick[stage];
  G4Tubs* tube = new G4Tubs("OuterVesselTube", vesselIR, vesselRadius[stage],
			    vesselHeight[stage]/2., 0, 360*deg);

  G4RotationMatrix* pipeRot = new G4RotationMatrix;
  pipeRot->rotateY(90*deg);

  // All cryostat vessels have cryogen pipe on the left side
  if (verboseLevel>2) G4cout << " opening hole for cryo pipes" << G4endl;

  G4Tubs* hole1 = new G4Tubs("OuterStemHole", 0.0, stemRad[stage]+tolerance,
			     stemRad[stage], 0, 360*deg);

  // Z position here is relative to center of vessel, not absolute coordinate
  G4ThreeVector holePos1(-vesselRadius[stage], 0.0,
			 zStemHole[stage]-zVessel[stage]);

  G4VSolid* tubeWithHoles =
    new G4SubtractionSolid("TubeWithHole", tube, hole1, pipeRot, holePos1);

  // Outer vessels have a vacuum pipe on the right side
  if (stage >= (NStages-NVacuum)) {
    if (verboseLevel>2) G4cout << " opening hole for vacuum pipes" << G4endl;

    G4int vacStage = stage - NVacuum;
    G4Tubs* hole2 =
      new G4Tubs("OuterVacuumHole", 0.0, vacRad[vacStage]+tolerance,
		 vacRad[vacStage], 0, 360*deg);

    // Z position here is relative to center of vessel, not absolute coordinate
    G4ThreeVector holePos2(vesselRadius[stage], 0.0,
			   zVacHole[vacStage]-zVessel[stage]);

    tubeWithHoles =
      new G4SubtractionSolid("TubeWithTwoHoles", tubeWithHoles, hole2, pipeRot, holePos2);
  }

  G4Material* copper = CDMSMaterialTable::GetMaterial("G4_Cu");
  G4LogicalVolume* logicalSide =
    new G4LogicalVolume(tubeWithHoles, copper, "LogicalVessel", 0, 0, 0);
                                                            
  logicalSide->SetVisAttributes(vesselVisAtt[stage]);

  return logicalSide;
}

G4LogicalVolume* CDMSVesselConstruction::BuildVesselBottom(G4int stage) {
  if (verboseLevel>1)
    G4cout << "CDMSVesselConstruction::BuildVesselBottom " << stage << G4endl;

  G4Tubs* bottomDisk =
    new G4Tubs("BottomDisk", 0., vesselRadius[stage]-VesselThick[stage],
	       VesselThick[stage]/2., 0, 360*deg); 

  G4Material* copper = CDMSMaterialTable::GetMaterial("G4_Cu");
  G4LogicalVolume* logicalBottom = 
    new G4LogicalVolume(bottomDisk, copper, "LogicalBottomLid", 0, 0, 0);
  logicalBottom->SetVisAttributes(vesselVisAtt[stage]);
  
  return logicalBottom;
}
 
G4LogicalVolume* CDMSVesselConstruction::BuildCryoPipe(G4int stage) {
  if (verboseLevel>1)
    G4cout << "CDMSVesselConstruction::BuildCryoPipe " << stage << G4endl;

  G4double stemIR = stemRad[stage]-PipeThick;
  G4Tubs* cryoTube = new G4Tubs("CryoTube", stemIR, stemRad[stage],
				pipeLen[stage]/2., 0, 360*deg);
  
  G4Material* copper = CDMSMaterialTable::GetMaterial("G4_Cu");
  G4LogicalVolume* logicalCryoPipe =
    new G4LogicalVolume(cryoTube, copper, "LogicalCryoPipe", 0, 0, 0);
  logicalCryoPipe->SetVisAttributes(vesselVisAtt[stage]);
  
  return logicalCryoPipe;
}

G4LogicalVolume* CDMSVesselConstruction::BuildVacuumPipe(G4int stage) {
  G4int vacStage = stage - NVacuum;		// Outer vessels are vacuum
  if (vacStage < 0) return 0;

  if (verboseLevel>1)
    G4cout << "CDMSVesselConstruction::BuildVacuumPipe " << stage << G4endl;

  G4double vacIR = vacRad[vacStage]-PipeThick;
  G4Tubs* vacTube = new G4Tubs("VacuumTube", vacIR, vacRad[vacStage],
				pipeLen[stage]/2., 0, 360*deg);

  G4Material* aluminum = CDMSMaterialTable::GetMaterial("G4_Al");
  G4LogicalVolume* logicalVacPipe =
    new G4LogicalVolume(vacTube, aluminum, "LogicalVacPipe", 0, 0, 0);
  logicalVacPipe->SetVisAttributes(vesselVisAtt[stage]);

  return logicalVacPipe;
}


// Compute position and rotation for pipe positions (also used for envelope)

G4Transform3D CDMSVesselConstruction::GetCryoPipePosition(G4int stage) {
  G4RotationMatrix pipeRot;
  pipeRot.rotateY(90*deg);

  G4double rPipe = pipeLen[stage]/2. + vesselRadius[stage] - VesselThick[stage];
  return G4Transform3D(pipeRot, G4ThreeVector(-rPipe, 0., zStemHole[stage]));
}

G4Transform3D CDMSVesselConstruction::GetVacuumPipePosition(G4int stage) {
  G4int vacStage = stage - NVacuum;		// Outer vessels are vacuum
  if (vacStage < 0) return G4Transform3D::Identity;

  G4RotationMatrix pipeRot;
  pipeRot.rotateY(90*deg);

  G4double rPipe = pipeLen[stage]/2. + vesselRadius[stage] - VesselThick[stage];
  return G4Transform3D(pipeRot, G4ThreeVector(rPipe, 0., zVacHole[vacStage]));
}


G4double CDMSVesselConstruction::GetRadius() const {
  return (vesselRadius.back() + VesselThick.back());
}

G4double CDMSVesselConstruction::GetLength() const {
  // These are the (absolute) half-heights in each direction
  G4double heightTop    = zVesselTop.back() + VesselThick.back();
  G4double heightBottom = -zVesselBottom.back() + VesselThick.back();

  return (heightTop+heightBottom);
}


void CDMSVesselConstruction::PrintParameters(std::ostream& os) const {
  os << "CDMSVesselConstruction parameters\n NStages " << NStages
     << "\n NVacuum " << NVacuum
     << "\n VesselThick";
  PrintVectorParameter(os, VesselThick);

  os << "\n Vessel0Rad " << Vessel0Rad << " mm"
     << "\n VesselDeltaRad " << VesselDeltaRad << " mm"
     << "\n VesselDeltaHeight " << VesselDeltaHeight << " mm"
     << "\n VesselExtraHeight " << VesselExtraHeight << " mm"
     << "\n VesselGap " << VesselGap << " mm"
     << "\n LidThickness";
  PrintVectorParameter(os, LidThickness);

  os << "\n StemInnerRad " << StemInnerRad << " mm"
     << "\n VacInnerRad " << VacInnerRad << " mm"
     << "\n PipeThick " << PipeThick << " mm"
     << "\n PipeOutsideLen " << PipeOutsideLen << " mm"
     << "\n vesselRadius";
  PrintVectorParameter(os, vesselRadius);

  os << "\n vesselHeight";
  PrintVectorParameter(os, vesselHeight);

  os << "\n zVesselTop";
  PrintVectorParameter(os, zVesselTop);

  os << "\n zVessel";
  PrintVectorParameter(os, zVessel);

  os << "\n zVesselBottom";
  PrintVectorParameter(os, zVesselBottom);

  os << "\n stemRad";
  PrintVectorParameter(os, stemRad);

  os << "\n zStemHole";
  PrintVectorParameter(os, zStemHole);

  os << "\n vacRad";
  PrintVectorParameter(os, vacRad);

  os << "\n zVacHole";
  PrintVectorParameter(os, zVacHole);

  os << std::endl;
}
