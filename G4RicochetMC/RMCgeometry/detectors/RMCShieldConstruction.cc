////////////////////////////////////////////////////////////////////////
//  File:        RMCShieldConstruction.cc                            //
//                                                                    //
//  Description: RMC shield and veto arrangement outside of icebox   //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        3 December 2010                                      //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "RMCgeometry/detectors/RMCShieldConstruction.hh"
#include "RMCgeometry/detectors/RMCVesselConstruction.hh"
#include "RMCgeometry/interface/RMCShieldMessenger.hh"
#include "RMCgeometry/interface/RMCGeometryManager.hh"
#include "RMCg4base/RMCMaterialTable.hh"
#include "G4Box.hh"
#include "G4Color.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4Point3D.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4Polyhedra.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4VSolid.hh"
#include "G4VisAttributes.hh"
#include "G4VSensitiveDetector.hh"
#include "G4SDManager.hh"
#include "RMCVetoSD.hh"
#include <algorithm>
#include <cmath>
#include <iostream>


// Constructor and destructor for shielding layer

RMCShieldConstruction::LayerData::
LayerData(const char* name, G4double gap, G4double side, G4double top,
	  G4double bottom, const char* mat)
  : Gap(gap), SideThick(side), TopThick(top), BottomThick(bottom),
    Radius(0.), Length(0.), Mid(0.), Z(0.), Name(name), matName(mat),
    visAtt(0) {}

RMCShieldConstruction::LayerData::~LayerData() {
  delete visAtt; visAtt=0;
}


// Fill parameters for single layer of shielding (updates input parameters)
// This works for the outer shielding, where we must work from the inner
// layers out, relative to the outer icebox boundary.

void 
RMCShieldConstruction::LayerData::FillOuterShield(G4double& radius, G4double& length,
					G4double& z) {
  Radius = (radius += Gap + SideThick);
  Length = (length += 2.*Gap + TopThick + BottomThick);
  Mid    = (BottomThick - TopThick)/2.;
  Z      = (z -= Mid);
}


// Fill parameters for single layer of shielding (updates input parameters)
// This works for the inner shielding, where we must work from the outer
// layers in, relative to the inner icebox boundary.

void 
RMCShieldConstruction::LayerData::FillInnerShield(G4double& radius, G4double& length,
					G4double& z) {
  Radius = (radius -= (Gap + SideThick));
  Length = (length -= (2.*Gap + TopThick + BottomThick));
  Mid    = (BottomThick - TopThick)/2.;
  Z      = (z += Mid);
}


// Constructor and destructor

RMCShieldConstruction::RMCShieldConstruction()
  : RMCVDetectorGeometry("VetoShield"), icebox(new RMCVesselConstruction()), DetectorName("Veto"),
    iceboxRad(61*cm), iceboxLength(136*cm), iceboxZ(0*cm),
    fridgeStemRad(11.0*cm), vacStemRad(8.0*cm),
    zFridgeStem(-40*cm), zVacStem(40*cm),
    mother("VetoShield",0.,1.4*m,1.5*m,1.5*m,"G4_AIR"),
    muMetal("MuMetalShield",1.0*cm,0.0381*cm,0.0381*cm,0.0381*cm,"MuMetal"),  
    innerPoly("InnerPoly",0.4826*cm,1.0*cm,1.0*cm,1.0*cm,"G4_POLYETHYLENE"),
    innerLead("InnerLead",0.1*cm,1.0*cm,1.0*cm,1.0*cm,"G4_Pb"),
    outerLead("OuterLead",0.47625*cm,15.0*cm,17.78*cm,17.78*cm,"G4_Pb"),
    outerPoly("OuterPoly",0.79375*cm,25.0*cm,40.64*cm,40.64*cm,"G4_POLYETHYLENE"),
    outerShieldRad(0.0), outerShieldLength(0.0), outerShieldZ(0.0),
    scintThick(5.08*cm), NTopPanelsX(3), NTopPanelsY(2), NTopPanels(6),
    overhang(20.0*cm), overlapX(5.08*cm), overlapY(10.16*cm),
    topVetoSupportHeight(5.72*cm), bottomVetoSupportHeight(8.4074*cm),
    topVetoLength(0.), topPanelLength(0.), topPanelWidth(0.),
    NSidePanels(8), sidePanelClearance(5.08*cm), sidePanelOverlap(5.08*cm),
    sideCornerOverlap(5.08*cm), sideRingRadius1(0.), sideRingRadius2(0.),
    sidePanelWidth(0.), sidePanelShift(0.), midPanelShift(0.),
    loPanelBottom(0.), loPanelTop(0.), loPanelHeight(0.), loPanelZ(0.),
    midPanelHeight(0.), midPanelZ(0.),
    hiPanelBottom(0.), hiPanelTop(0.), hiPanelHeight(0.), hiPanelZ(0.),
    iVacPanel(0), iFridgePanel(0),
    messenger(new RMCShieldMessenger(this)) {}

RMCShieldConstruction::~RMCShieldConstruction() {
  delete messenger; messenger=0;
}


void RMCShieldConstruction::SetVerboseLevel(G4int verbose) {
  RMCVDetectorGeometry::SetVerboseLevel(verbose);
  if (icebox) icebox->SetVerboseLevel(verbose);
}


// Specify that shielding surrounds a cryostat (fetched from manager)

void RMCShieldConstruction::UseIcebox() {
  icebox = RMCGeometryManager::Instance()->GetVessel();
}


// Check here for valid polygons, report error
void RMCShieldConstruction::SetNSidePanels(G4int value) {
  if (value > 0 && value < 4) {
    value = 3;
    G4cerr << " RMCShieldConstruction: side panels must form a polygon."
	   << "  N set to 3" << G4endl;
  }

  NSidePanels = value;
}


void RMCShieldConstruction::FillExtraParameters() {
  SetVerboseLevel(verboseLevel);	// Ensure that children have verbosity

  if (verboseLevel)
    G4cout << "RMCShieldConstruction::FillExtraParameters()" << G4endl;

  FillIceboxParameters();
  FillShieldParameters();
  FillVetoParameters();
  FillMotherParameters();
}

void RMCShieldConstruction::FillIceboxParameters() {
  if (!icebox) return;

  if (verboseLevel>1)
    G4cout << "RMCShieldConstruction::FillIceboxParameters()" << G4endl;

  icebox->FillExtraParameters();
  iceboxRad    = icebox->GetRadius();
  iceboxLength = icebox->GetLength();
  iceboxZ      = icebox->GetPosition().z();

  // Include clearance to avoid 
  fridgeStemRad = icebox->GetStemRadList().back() + 2.*tolerance;
  vacStemRad    = icebox->GetVacRadList().back() + 2.*tolerance;
  
  zFridgeStem   = icebox->GetZStemHoleList().back() + iceboxZ;
  zVacStem      = icebox->GetZVacHoleList().back() + iceboxZ;
}

void RMCShieldConstruction::FillShieldParameters() {
  if (verboseLevel>1)
    G4cout << "RMCShieldConstruction::FillShieldParameters()" << G4endl;

  // All of the sheilding layers' parameters are computed incrementally
  G4double outerShieldRadius = iceboxRad;
  G4double innerShieldRadius = icebox->GetVessel0Rad();
  G4double outerShieldLen = iceboxLength;
  const std::vector<G4double>& heights = icebox->GetVesselHeight();
  G4double innerShieldLen = heights[0];
  G4double outerShieldZPos   = 0.;
  const std::vector<G4double>& zpositions = icebox->GetVesselZPos();
  G4double innerShieldZPos   = zpositions[0] - zpositions[5];

  // Fill the inner layers incrementally, working from outermost
  // to innermost
  innerPoly.FillInnerShield(innerShieldRadius, innerShieldLen, innerShieldZPos);
  innerPoly.visAtt = new G4VisAttributes(G4Color(0.5, 1.0, 0.5));
  innerPoly.visAtt->SetForceSolid(true);

  innerLead.FillInnerShield(innerShieldRadius, innerShieldLen, innerShieldZPos);
  innerLead.visAtt = new G4VisAttributes(G4Color(0.5, 0.5, 0.5));
  innerLead.visAtt->SetForceSolid(true);

  muMetal.FillInnerShield(innerShieldRadius, innerShieldLen, innerShieldZPos);
  muMetal.visAtt = new G4VisAttributes(G4Color(1, 1, 0));
  muMetal.visAtt->SetForceSolid(true);

  // Fill the outer layers incrementally, working from innermost
  // to outermost
  outerPoly.FillOuterShield(outerShieldRadius, outerShieldLen, outerShieldZPos);
  outerPoly.visAtt = new G4VisAttributes(G4Color(0.5, 1.0, 0.5));
  outerPoly.visAtt->SetForceSolid(true);

  outerLead.FillOuterShield(outerShieldRadius, outerShieldLen, outerShieldZPos);
  outerLead.visAtt = new G4VisAttributes(G4Color(1, 1, 0));
  outerLead.visAtt->SetForceSolid(true);

  // Copy outermost sheild layer data for use by veto ayer
  outerShieldRad    = outerShieldRadius;
  outerShieldLength = outerShieldLen;
  outerShieldZ      = outerShieldZPos;
}

void RMCShieldConstruction::FillVetoParameters() {
  if (verboseLevel>1)
    G4cout << "RMCShieldConstruction::FillVetoParameters()" << G4endl;

  FillEndcapVetoParameters();
  FillBarrelVetoParameters();
}

void RMCShieldConstruction::FillEndcapVetoParameters() {
  NTopPanels = NTopPanelsX*NTopPanelsY;
  if (0 == NTopPanels) return;			// User suppressed panels

  if (verboseLevel>1)
    G4cout << "RMCShieldConstruction::FillEndcapVetoParameters()" << G4endl;

  // Endcap veto is built as a square centered on the cylindrical shield
  topVetoLength = 2.*(outerShieldRad + overhang);
  topPanelWidth  = (topVetoLength - overlapX) / NTopPanelsX + overlapX;
  topPanelLength = (topVetoLength - overlapY) / NTopPanelsY + overlapY;

  panelX.resize(NTopPanels,0.);		// Preallocate for assignment by index
  panelY.resize(NTopPanels,0.);
  panelZtop.resize(NTopPanels,0.);
  panelZbot.resize(NTopPanels,0.);

  // Panels overlap in both X and Y, with staggered Z offsets
  G4double x0  = (-topVetoLength + topPanelWidth)/2.;
  G4double y0  = (-topVetoLength + topPanelLength)/2.;
  G4double z0t = outerShieldZ + outerShieldLength/2. + topVetoSupportHeight;
  G4double z0b = outerShieldZ - outerShieldLength/2. - bottomVetoSupportHeight;

  G4double deltaX = topPanelWidth - overlapX;
  G4double deltaY = topPanelLength - overlapY;
  G4double deltaZ = scintThick + tolerance;

  G4int ix, iy;
  for (G4int i=0; i<NTopPanels; i++) {
    ix = i % NTopPanelsX;		// X index counts "first", Y "second"
    iy = i / NTopPanelsX;
    G4int zStep = 2*(ix%2) + iy%2;	// Staggered in double layers

    panelX[i] = x0 + ix*deltaX;
    panelY[i] = y0 + iy*deltaY;
    panelZtop[i] = z0t + zStep*deltaZ;
    panelZbot[i] = z0b - zStep*deltaZ;
  }
}

void RMCShieldConstruction::FillBarrelVetoParameters() {
  if (0 == NSidePanels) return;		// User suppressed panels

  if (verboseLevel>1)
    G4cout << "RMCShieldConstruction::FillBarrelVetoParameters()" << G4endl;

  // Panels are placed in regular polygon with sides shifted 
  const G4double deltaPhi = 360*deg / NSidePanels;
  const G4double sinPhi = sin(deltaPhi);
  const G4double tanHalf = tan(deltaPhi/2.);

  // Middle panels are placed outside top and bottom panels
  sideRingRadius1 = outerShieldRad + scintThick/2. + sidePanelClearance;
  sideRingRadius2 = sideRingRadius1 + scintThick + tolerance;

  sidePanelWidth = 2.*sideRingRadius1*tanHalf + sideCornerOverlap/sinPhi;
  sidePanelShift = (sideCornerOverlap/2. + scintThick) / sinPhi;
  midPanelShift  = sidePanelShift + scintThick/sinPhi;

  // Bottom panel runs from lower endcap veto to just below cryogen pipe
  loPanelBottom = panelZbot[0] + scintThick/2. + tolerance;
  loPanelTop    = zFridgeStem - fridgeStemRad;

  loPanelHeight = loPanelTop - loPanelBottom;
  loPanelZ      = (loPanelTop + loPanelBottom)/2.;

  // Top panel runs from upper endcap veto to just above vacuum pipe
  hiPanelBottom = zVacStem + vacStemRad;
  hiPanelTop    = panelZtop[0] - scintThick/2. - tolerance;

  hiPanelHeight = hiPanelTop - hiPanelBottom;
  hiPanelZ      = (hiPanelTop + hiPanelBottom)/2.;

  // Middle panel overlaps the top and bottom panels (so pipes go through it)
  midPanelHeight = hiPanelBottom - loPanelTop + 2.*sidePanelOverlap;
  midPanelZ      = (hiPanelBottom + loPanelTop)/2.;

  // Calculate placements
  topTransform.resize(NSidePanels);	// Preallocate for assignment by index
  middleTransform.resize(NSidePanels);
  bottomTransform.resize(NSidePanels);

  G4RotationMatrix rot;			// Buffers to reduce memory churn
  G4ThreeVector base, shift;

  for (G4int i = 0; i < NSidePanels; i++) {
    base.setRhoPhiZ(sideRingRadius1, i*deltaPhi, hiPanelZ);
    shift.setRhoPhiZ(sidePanelShift, i*deltaPhi-90*deg, 0.);
    topTransform[i] = G4Transform3D(rot, base+shift);

    base.setZ(loPanelZ);		// Same (x,y) position as top panel
    bottomTransform[i] = G4Transform3D(rot, base+shift);

    base.setRhoPhiZ(sideRingRadius2, i*deltaPhi, midPanelZ);
    shift.setRhoPhiZ(midPanelShift, i*deltaPhi-90*deg, 0.);
    middleTransform[i] = G4Transform3D(rot, base+shift);

    rot.rotateZ(deltaPhi);		// Step orientation around cylinder
  }

  // Pipe stems are along +x (vacuum) and -x (cryogen)
  iVacPanel = 0;
  iFridgePanel = NSidePanels / 2;	// FIXME:  Assumes even number of panels
}

void RMCShieldConstruction::FillMotherParameters() {
  if (verboseLevel>1)
    G4cout << "RMCShieldConstruction::FillMotherParameters()" << G4endl;

  // Start with dimensions of shielding
  mother.Radius = outerShieldRad + tolerance;
  mother.Length = outerShieldLength + 2*std::abs(outerShieldZ) + tolerance;

  // FIXME:  Would like asymmetric mother, but without recomputing everything
  mother.Z = 0.;

  // Use panel transforms to determine coordinate of outermost corner
  if (NSidePanels > 0) {
    // NOTE: Corner point is -Y because shift is "downward" (-phi)
    G4double maxX=0., maxY=0.;
    G4Point3D base(scintThick/2., -sidePanelWidth/2., 0.), ipos;
    for (int i=0; i<NSidePanels; i++) {
      ipos = middleTransform[i]*base;
      maxX = std::max(maxX, std::abs(ipos.x()));
      maxY = std::max(maxY, std::abs(ipos.y()));
    }

    mother.Radius = std::max(std::max(maxX,maxY), topVetoLength/2.) + tolerance;
  }

  // Outside lengths are the outer surface of the outermost veto panels
  if (NTopPanels > 0) {
    G4double topZmax = *std::max_element(panelZtop.begin(), panelZtop.end());
    G4double bottomZmax = *std::min_element(panelZbot.begin(), panelZbot.end());
    
    mother.Length = 2.*std::max(topZmax, -bottomZmax) + scintThick + tolerance;
  }

  if (icebox) icebox->SetPipeOutsideLen(mother.Radius - iceboxRad);

  mother.visAtt = new G4VisAttributes(G4Color(0.5, 0.5, 1));
  mother.visAtt->SetForceWireframe(true);
}


G4LogicalVolume* RMCShieldConstruction::BuildGeometry() {
  if (verboseLevel)
    G4cout << "RMCShieldConstruction::BuildGeometry()" << G4endl;

  FillExtraParameters();
  if (verboseLevel > 1) PrintParameters(G4cout);

  // Build mother volume for shield and veto
  G4LogicalVolume* mother = CreateMotherVolume();

  PlaceShielding(mother);
  PlaceVeto(mother);

  if (verboseLevel)
    G4cout << "Shielding and veto mass " << mother->GetMass()/kg << " kg"
	   << G4endl;

  return mother;
}


G4LogicalVolume* RMCShieldConstruction::CreateMotherVolume() {
  if (verboseLevel)
    G4cout << "RMCShieldConstruction::CreateMotherVolume()" << G4endl;

  // Make mother volume a rectangular box, to enclose scintillator lids
  G4Box* box = new G4Box(GetName(), mother.Radius, mother.Radius,
			 mother.Length/2.);
  
  G4LogicalVolume* logicalShieldVolume = 
    new G4LogicalVolume(box, RMCMaterialTable::GetMaterial(mother.matName),
                        mother.Name, 0, 0, 0);
  logicalShieldVolume->SetVisAttributes(mother.visAtt);

  return logicalShieldVolume;
}


void RMCShieldConstruction::PlaceShielding(G4LogicalVolume* mother)
{
  if (verboseLevel)
    G4cout << "RMCShieldConstruction::PlaceShielding()" << G4endl;

  // Get the mother volume for the inner shielding.  This should just
  // be the innermost vessel of the icebox.
  //G4LogicalVolume* innerShieldMotherVolume = icebox->

  PlaceShieldLayer(muMetal, mother, true);
  PlaceShieldLayer(innerPoly, mother, true);
  PlaceShieldLayer(innerLead, mother, true);
  PlaceShieldLayer(outerLead, mother, false);
  PlaceShieldLayer(outerPoly, mother, false);
}

void RMCShieldConstruction::PlaceShieldLayer(const LayerData& layer,
					      G4LogicalVolume* mother, G4bool innerType) {
  if (verboseLevel > 1)
    G4cout << "RMCShieldConstruction::PlaceShieldLayer " << layer.Name
	   << G4endl;

  G4LogicalVolume* logicalLayer = BuildShieldLayer(layer, innerType);
  new G4PVPlacement(0, G4ThreeVector(0.,0.,layer.Z), logicalLayer, layer.Name,
		    mother, false, 0);
}

// Shield layers consist of a hollow cylinder (possibly offset vertically)
// with two side penetrations

G4LogicalVolume*
RMCShieldConstruction::BuildShieldLayer(const LayerData& layer, G4bool innerType) {
  if (verboseLevel > 1)
    G4cout << "RMCShieldConstruction::BuildShieldLayer " << layer.Name
	   << G4endl;
  
  // Start with a hollow cylinder (hole may be offset)
  G4double coreRad = layer.Radius-layer.SideThick;
  G4double coreLen = layer.Length - layer.TopThick - layer.BottomThick;
  G4ThreeVector pos(0.0, 0.0, layer.Mid);

  G4Tubs* outer = new G4Tubs("Tube1", 0.0, layer.Radius, layer.Length/2., 0, 360*deg);
  G4Tubs* core  = new G4Tubs("Tube2", 0.0, coreRad, coreLen/2., 0, 360*deg);
  G4VSolid* shell = new G4SubtractionSolid("Shell", outer, core, 0, pos);

  // Remove the cryogen and vacuum pipe penetrations
  G4double holeX = layer.Radius - layer.SideThick/2.;	// Middle of material

  // Extra length on hole accounts for curvature of shield layer
  if(innerType == false)
  {
    shell = MakePipeHoleInSolid(shell, fridgeStemRad, layer.SideThick+1.5*cm,
				-holeX, 0., zFridgeStem-layer.Z);

    shell = MakePipeHoleInSolid(shell, vacStemRad, layer.SideThick+1.5*cm,
				holeX, 0., zVacStem-layer.Z);
  }

  // Return completed structure
  G4LogicalVolume* logicalShieldLayer =
    new G4LogicalVolume(shell, RMCMaterialTable::GetMaterial(layer.matName),
			("Logical"+layer.Name).c_str(), 0, 0, 0);
  logicalShieldLayer->SetVisAttributes(layer.visAtt);

  if (verboseLevel)
    G4cout << "Shielding " << layer.Name << " mass "
	   << logicalShieldLayer->GetMass()/kg << " kg" << G4endl;

  return logicalShieldLayer;
}


// Add pipe penetration (along Y axis) through solid at specified location

G4VSolid*
RMCShieldConstruction::MakePipeHoleInSolid(G4VSolid* solid,
					    G4double radius, G4double thick,
					    const G4ThreeVector& hole) {
  if (verboseLevel > 1) {
    G4cout << "RMCShieldConstruction::MakePipeHoleInSolid"
	   << " r " << radius << " l " << thick << " at " << hole << G4endl;
  }

  // Both pipe penetrations are horizontal, along X axis
  static G4RotationMatrix rotation(CLHEP::HepRotationY(90*deg));

  G4Tubs* pipe = new G4Tubs("PipeHole", 0., radius, thick, 0, 360*deg);

  return new G4SubtractionSolid(solid->GetName(), solid, pipe, &rotation, hole);
}


// Veto scintillator panels surround outermost shielding layer

void RMCShieldConstruction::PlaceVeto(G4LogicalVolume* mother)
{
  if (verboseLevel) G4cout << "RMCShieldConstruction::PlaceVeto()" << G4endl;

  // Users may suppress veto panels by setting counts to zero
  if (NTopPanels > 0)  PlaceEndcapVeto(mother);
  if (NSidePanels > 0) PlaceBarrelVeto(mother);
}

//  Top and bottom scintillator panels (bottom array is mirror of top)

void RMCShieldConstruction::PlaceEndcapVeto(G4LogicalVolume* mother)
{
  if (0 == NTopPanels) return;		// User suppressed veto

  if (verboseLevel>1)
    G4cout << "RMCShieldConstruction::PlaceEndcapVeto()" << G4endl;

  G4Material* scint = RMCMaterialTable::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  G4Box* topPanel = new G4Box("TopPanel", topPanelWidth/2.,
                              topPanelLength/2., scintThick/2.);

  G4LogicalVolume* logTopPanel =
    new G4LogicalVolume(topPanel, scint, "LogTopScintPanel");

  G4VisAttributes* VisAttScint = new G4VisAttributes(G4Color(1.0, 0.0, 1.0));
  VisAttScint->SetForceSolid(true);
  logTopPanel->SetVisAttributes(VisAttScint);

  logTopPanel->SetSensitiveDetector(BuildSensitiveDetector());

  G4ThreeVector pos;		// Buffer to reduce memory churn
  for (G4int i=0; i < NTopPanels; i++) {
    pos.set(panelX[i], panelY[i], panelZtop[i]);
    new G4PVPlacement(0,pos,logTopPanel,"PhysTopScintPanel",mother,false,copyNumber);

    pos.setZ(panelZbot[i]);
    new G4PVPlacement(0,pos,logTopPanel,"PhysBottomScintPanel",mother,false,
		      copyNumber+NTopPanels);
    copyNumber++;
  }
  copyNumber += NTopPanels; 
}


void RMCShieldConstruction::PlaceBarrelVeto(G4LogicalVolume* mother)
{
  if (0 == NSidePanels) return;		// User suppressed veto

  if (verboseLevel>1)
    G4cout << "RMCShieldConstruction::PlaceBarrelVeto()" << G4endl;

  PlaceTopSideVeto(mother);
  PlaceMiddleSideVeto(mother);
  PlaceBottomSideVeto(mother);
}

void RMCShieldConstruction::PlaceTopSideVeto(G4LogicalVolume* mother)
{
  if (verboseLevel>1)
    G4cout << "RMCShieldConstruction::PlaceTopSideVeto()" << G4endl;

  G4Material* scint = RMCMaterialTable::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  // Build panels oriented vertically in Y-Z plane
  G4Box* topSidePanel = new G4Box("TopSidePanel",  scintThick/2.,
				  sidePanelWidth/2., hiPanelHeight/2.);

  G4LogicalVolume* logTopSidePanel =
    new G4LogicalVolume(topSidePanel, scint, "LogTopSideScint");
  G4VisAttributes* VisAttSideScint = new G4VisAttributes(G4Color(1., 0.5, 0.5));
  VisAttSideScint->SetForceSolid(true);
  logTopSidePanel->SetVisAttributes(VisAttSideScint);

  logTopSidePanel->SetSensitiveDetector(BuildSensitiveDetector());

  for (G4int i=0; i<NSidePanels; i++) {
    new G4PVPlacement(topTransform[i], logTopSidePanel, "PhysTopSideScint",
                      mother, false, copyNumber);
    copyNumber++;
  }
}

void RMCShieldConstruction::PlaceMiddleSideVeto(G4LogicalVolume* mother)
{
  if (verboseLevel>1)
    G4cout << "RMCShieldConstruction::PlaceMiddleSideVeto()" << G4endl;

  G4Material* scint = RMCMaterialTable::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  // Build panels oriented vertically in Y-Z plane
  G4Box* middleSidePanel = new G4Box("MiddleSidePanel", scintThick/2.,
				     sidePanelWidth/2., midPanelHeight/2.);

  // Basic panel can be placed everywhere with one pointer
  G4LogicalVolume* logMidSidePanel =
    new G4LogicalVolume(middleSidePanel, scint, "LogMidSideScint");

  G4VisAttributes* VisAttMidScint = new G4VisAttributes(G4Color(1.0, 1.0, 0.));
  VisAttMidScint->SetForceSolid(true);
  logMidSidePanel->SetVisAttributes(VisAttMidScint);

  logMidSidePanel->SetSensitiveDetector(BuildSensitiveDetector());

  // Middle panels intersect the cryostat plumbing, will need to make holes
  G4double holeOffset;

  for (G4int i = 0; i<NSidePanels; i++) {
    G4LogicalVolume* placementPanel = logMidSidePanel;

    if (PanelIntersectsVacStem(i, holeOffset)) {
      placementPanel = MakeVacHoleInPanel(middleSidePanel, i, holeOffset);
      placementPanel->SetVisAttributes(VisAttMidScint);
    } else if (PanelIntersectsFridgeStem(i, holeOffset)) {
      placementPanel = MakeFridgeHoleInPanel(middleSidePanel, i, holeOffset);
      placementPanel->SetVisAttributes(VisAttMidScint);
    }

    new G4PVPlacement(middleTransform[i], placementPanel,
		      "PhysMiddleSideScint", mother, false, copyNumber);
    copyNumber++;
  }
}

void RMCShieldConstruction::PlaceBottomSideVeto(G4LogicalVolume* mother)
{
  if (verboseLevel>1)
    G4cout << "RMCShieldConstruction::PlaceBottomSideVeto()" << G4endl;

  G4Material* scint = RMCMaterialTable::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  // Build panels oriented vertically in Y-Z plane
  G4Box* bottomSidePanel = new G4Box("BottomSidePanel", scintThick/2.,
				     sidePanelWidth/2., loPanelHeight/2.);

  G4LogicalVolume* logBottomSidePanel =
    new G4LogicalVolume(bottomSidePanel, scint, "LogBottomSideScint");
  G4VisAttributes* VisAttSideScint = new G4VisAttributes(G4Color(1., 0.5, 0.5));
  VisAttSideScint->SetForceSolid(true);
  logBottomSidePanel->SetVisAttributes(VisAttSideScint);

  logBottomSidePanel->SetSensitiveDetector(BuildSensitiveDetector());

  for (G4int i = 0; i < NSidePanels; i++) {
    new G4PVPlacement(bottomTransform[i], logBottomSidePanel,
                      "PhysBottomSideScint", mother, false, copyNumber);
    copyNumber++;
  }
}


// Determine whether final panel position contacts cryostat plumbing

void RMCShieldConstruction::GetCornersOfPanel(G4int index, G4Point3D& right,
					       G4Point3D& left) {
  static const G4Point3D rightCorner(-scintThick/2., -sidePanelWidth/2., 0.);
  static const G4Point3D leftCorner(-scintThick/2., sidePanelWidth/2., 0.);

  // Apply positioning transformation to near corners, and return
  right = middleTransform[index] * rightCorner;
  left  = middleTransform[index] * leftCorner;

  if (verboseLevel > 2) {
    G4cout << "RMCShieldConstruction::GetCornersOfPanel " << index
	   << "\n right = " << right << " left = " << left << G4endl;
  }
}

G4bool RMCShieldConstruction::PanelIntersectsVacStem(G4int index,
						      G4double& panelY) {
  G4Point3D rightCorner, leftCorner;
  GetCornersOfPanel(index, rightCorner, leftCorner);

  // Vacuum stem is on positive X side
  // Intersects if right corner below +Y edge AND left corner above -Y edge
  G4bool onSide  = (rightCorner.x() > 0. || leftCorner.x() > 0.);
  G4bool hitPipe = (rightCorner.y() < vacStemRad && 
		    leftCorner.y() > -vacStemRad);

  if (verboseLevel > 2) {
    G4cout << " checking panel " << index << " for vacStem r "
	   << vacStemRad << " mm at x > 0" << G4endl
	   << " right.x() " << rightCorner.x() << " left.x() " << leftCorner.x()
	   << " => on side? " << onSide << G4endl
	   << " right.y() " << rightCorner.y() << " left.y() " << leftCorner.y()
	   << " => hit pipe? " << hitPipe << G4endl;
  }

  // Find distance of stem in X-Y plane from center of panel
  if (onSide && hitPipe) panelY = DistanceToPipe(rightCorner, leftCorner);

  return (onSide && hitPipe);
}

G4bool RMCShieldConstruction::PanelIntersectsFridgeStem(G4int index,
							 G4double& panelY) {
  G4Point3D rightCorner, leftCorner;
  GetCornersOfPanel(index, rightCorner, leftCorner);

  // Fridge stem is on negative X side
  // Intersects if left corner below +Y edge AND right corner above -Y edge
  G4bool onSide  = (rightCorner.x() < 0. || leftCorner.x() < 0.);
  G4bool hitPipe = (rightCorner.y() > -fridgeStemRad && 
		    leftCorner.y() < fridgeStemRad);

  if (verboseLevel > 2) {
    G4cout << " checking panel " << index << " for fridgeStem r "
	   << fridgeStemRad << " mm at x < 0" << G4endl
	   << " right.x() " << rightCorner.x() << " left.x() " << leftCorner.x()
	   << " => on side? " << onSide << G4endl
	   << " right.y() " << rightCorner.y() << " left.y() " << leftCorner.y()
	   << " => hit pipe? " << hitPipe << G4endl;
  }
    
  // Find distance of stem in X-Y plane from center of panel
  if (onSide && hitPipe) panelY = DistanceToPipe(rightCorner, leftCorner);

  return (onSide && hitPipe);
}

// Compute signed distance to pipe stem (lies on X axis) from panel midpoint

G4double RMCShieldConstruction::DistanceToPipe(const G4Point3D& right,
						const G4Point3D& left) {
  if (verboseLevel > 2) 
    G4cout << "RMCShieldConstruction::DistanceToPipe" << G4endl;

  // Panel may be aligned with Y, but never with X)
  G4double invslope = (left.x()-right.x()) / (left.y()-right.y());
  G4double XatXaxis = right.x() - right.y()*invslope;
  
  G4double Xmid = (right.x()+left.x()) / 2.;
  G4double Ymid = (right.y()+left.y()) / 2.;

  // Distance to X axis along front face of panel [see GetsCornersOfPanel()]
  G4double Ydist = std::sqrt((Xmid-XatXaxis)*(Xmid-XatXaxis) + Ymid*Ymid);

  if (verboseLevel > 2) {
    G4cout << " midpoint x " << Xmid << " y " << Ymid
	   << " invslope " << invslope << " XatXaxis " << XatXaxis
	   << "\n distance from midpoint on face to X axis " << Ydist << G4endl;
  }

  // Sign uses position and orientation of panel (right -ve, left +ve)
  Ydist *= (Ymid>0.?-1.:1.) * (left.y()>right.y()?1.:-1);

  // Correct signed distance along midplane, to account for tilt
  Ydist += scintThick/2. * invslope;

  if (verboseLevel > 2)
    G4cout << " signed distance from center in midplane " << Ydist << G4endl;

  return Ydist;
}


// Add properly oriented hole to veto panel before placement

G4LogicalVolume* 
RMCShieldConstruction::MakeVacHoleInPanel(G4VSolid* panel, G4int index,
					   G4double holeY) {
  if (verboseLevel > 1)
    G4cout << "RMCShieldConstruction::MakeVacHoleInPanel " << index << G4endl;

  return MakePipeHoleInPanel(panel, index, vacStemRad, holeY, zVacStem);
}

G4LogicalVolume* 
RMCShieldConstruction::MakeFridgeHoleInPanel(G4VSolid* panel, G4int index,
					      G4double holeY) {
  if (verboseLevel > 1)
    G4cout << "RMCShieldConstruction::MakeFridgeHoleInPanel " << index
	   << G4endl;

  return MakePipeHoleInPanel(panel, index, fridgeStemRad, holeY, zFridgeStem);
}

G4LogicalVolume* 
RMCShieldConstruction::MakePipeHoleInPanel(G4VSolid* panel, G4int index,
					    G4double radius, G4double holeY,
					    G4double holeZ) {
  if (verboseLevel > 2) {
    G4cout << "RMCShieldConstruction::MakePipeHoleInPanel " << index 
	   << " r " << radius << " y " << holeY << " z(abs) " << holeZ
	   << " mm" << G4endl;
  }

  // Rotations are around Z axis, given by angle from rotated X axis to base
  G4double phi = middleTransform[index].getRotation().phiX();

  // Penetration length needs to account for rotation angle
  G4double thick = (scintThick/std::abs(cos(phi)) + 2.*radius*std::abs(tan(phi))
		    + 2.*tolerance);

  // Pipe penetration must be rotated to horizontal before placement
  G4RotationMatrix holerot;
  holerot.rotateZ(phi).rotateY(90*deg);	// CLHEP rotates coord axes, not object

  G4ThreeVector holepos(0., holeY, holeZ-midPanelZ);

  if (verboseLevel > 2) {
    G4cout << " Placing hole " << thick << " mm long, angled " << phi/deg
	   << " deg at " << holepos << G4endl;
  }

  // Add penetration to input panel, and wrap in logical volume
  G4Tubs* pipe = new G4Tubs("PipeHole", 0., radius, thick/2., 0, 360*deg);

  G4VSolid* panelWithHole =
    new G4SubtractionSolid(panel->GetName(), panel, pipe, &holerot, holepos);

  G4Material* scint = RMCMaterialTable::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  return new G4LogicalVolume(panelWithHole, scint, "LogPanelWithPipeHole");
}


G4VSensitiveDetector* RMCShieldConstruction::BuildSensitiveDetector()
{
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4VSensitiveDetector* vetoSD =
    SDman->FindSensitiveDetector(DetectorName, false);
  if (!vetoSD) {
    if (verboseLevel > 1)
      G4cout << " Creating sensitive detector " << DetectorName << G4endl;

    vetoSD = new RMCVetoSD(DetectorName);
    SDman->AddNewDetector(vetoSD);
  }

  return vetoSD;
}


// Report values of all configuration parameters

void RMCShieldConstruction::PrintParameters(std::ostream& os) const {
  os << "RMCShieldConstruction parameters"
     << "\n iceboxRad " << iceboxRad << " mm"
     << "\n iceboxLength " << iceboxLength << " mm"
     << "\n iceboxZ " << iceboxZ << " mm"
     << "\n fridgeStemRad " << fridgeStemRad << " mm"
     << "\n vacStemRad " << vacStemRad << " mm"
     << "\n zFridgeStem " << zFridgeStem << " mm"
     << "\n zVacStem " << zVacStem << " mm"
     << "\n mother " << mother
     << "\n muMetal " << muMetal
     << "\n innerPoly " << innerPoly
     << "\n innerLead " << innerLead
     << "\n outerLead " << outerLead
     << "\n outerPoly " << outerPoly

     << "\n outerShieldRad " << outerShieldRad << " mm"
     << "\n outerShieldLength " << outerShieldLength << " mm"
     << "\n outerShieldZ " << outerShieldZ << " mm"
     << "\n scintThick " << scintThick << " mm"

     << "\n NTopPanelsX " << NTopPanelsX << " NTopPanelsY " << NTopPanelsY
     << "\n overhang " << overhang << " mm"
     << "\n overlapX " << overlapX << " mm"
     << "\n overlapY " << overlapY << " mm"
     << "\n topVetoSupportHeight " << topVetoSupportHeight << " mm"
     << "\n bottomVetoSupportHeight " << bottomVetoSupportHeight << " mm"
     << "\n topVetoLength " << topVetoLength << " mm"
     << "\n topPanelLength " << topPanelLength << " mm"
     << "\n topPanelWidth " << topPanelWidth << " mm";

  os << "\n panelX[" << panelX.size() << "] "; PrintVectorParameter(os, panelX);
  os << "\n panelY[" << panelX.size() << "] "; PrintVectorParameter(os, panelY);
  os << "\n panelZtop[" << panelX.size() << "] "; PrintVectorParameter(os, panelZtop);
  os << "\n panelZbot[" << panelX.size() << "] "; PrintVectorParameter(os, panelZbot);

  os << "\n NSidePanels " << NSidePanels
     << "\n sidePanelClearance " << sidePanelClearance << " mm"
     << "\n sidePanelOverlap " << sidePanelOverlap << " mm"
     << "\n sideCornerOverlap " << sideCornerOverlap << " mm"
     << "\n sideRingRadius1 " << sideRingRadius1 << " mm"
     << "\n sideRingRadius2 " << sideRingRadius2 << " mm"
     << "\n sidePanelWidth " << sidePanelWidth << " mm"
     << "\n sidePanelShift " << sidePanelShift << " mm"
     << "\n midPanelShift " << midPanelShift << " mm"
     << "\n loPanelBottom " << loPanelBottom << " loPanelTop " << loPanelTop
     << "\n loPanelHeight " << loPanelHeight << " mm"
     << "\n loPanelZ " << loPanelZ << " mm"
     << "\n midPanelHeight " << midPanelHeight << " mm"
     << "\n midPanelZ " << midPanelZ << " mm"
     << "\n hiPanelBottom " << hiPanelBottom << " hiPanelTop " << hiPanelTop
     << "\n hiPanelHeight " << hiPanelHeight << " mm"
     << "\n hiPanelZ " << hiPanelZ << " mm"

     << std::endl;
}

// Print out parameters of single shielding layer

std::ostream& 
operator<<(std::ostream& os, const RMCShieldConstruction::LayerData& layer) {
  os << layer.Name << "(" << layer.matName << ")"
     << "\n  gap " << layer.Gap << " side " << layer.SideThick
     << " top " << layer.TopThick  << " bottom " << layer.BottomThick
     << " r " << layer.Radius << " l " << layer.Length << " mid " << layer.Mid
     << " z " << layer.Z << " [mm]";

  return os;
}
