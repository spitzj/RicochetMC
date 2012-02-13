////////////////////////////////////////////////////////////////////////
// $Id: Am241SourceHolder.cc,v 1.16 2011/07/22 21:09:17 kelsey Exp $
//  File:        Am241SourceHolder.hh                                 //
//  Description: Support plate and canisters for Am241 sources        //
//                                                                    //
//  This geometry defines the structure for the Am-241 source used    //
//  with single-ZIP test facilities.  It also directly instantiates   //
//  an instance of CDMSsources/Am241Source configured for each of     //
//  the canisters on the plate.                                       //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        20 April 2011                                        //
//                                                                    //
//  20110422  M. Kelsey -- Add names for each active source (Region), //
//		fix holder positions (mirrored vs. actual install),   //
//		add tiny gaps between source canister components.     //
//  20110426  M. Kelsey -- Make real source, CDMSVSourceConstruction, //
//		add event-generator functions, source activity (uCi). //
//  20110427  M. Kelsey -- Drop use of G4Region, use Messenger, fix   //
//		fix bug in throwing source atom position.  Option to  //
//		use hard-coded gamma lines instead of Am-241 decay.   //
//  20110428  M. Kelsey -- Use CLHEP RandExponential.h directly (not  //
//		included in GEANT4 mappings).  Carve out recess in    //
//		source puck so that active surface is exposed.        //
//  20110429  M. Kelsey -- Redesign canister as hollow cylinder       //
//  20110502  M. Kelsey -- Cast std::floor() to G4int for warnings.   //
//  20110617  M. Kelsey -- Allow zero thickness to leave out pieces   //
//  20110623  M. Kelsey -- Add configurable foil material             //
//  20110624  M. Kelsey -- Move active-position Z back behind surface //
//  20110629  M. Kelsey -- Initialize and use name arg in base class  //
//  20110722  M. Kelsey -- Discard matTable, using singleton instead. //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSgeometry/detectors/Am241SourceHolder.hh"

#include "CDMSgeometry/detectors/CDMSZipConstruction.hh"
#include "CDMSgeometry/interface/Am241HolderMessenger.hh"
#include "CDMSgeometry/interface/CDMSGeometryManager.hh"
#include "CDMSsources/Am241Lines.hh"
#include "CDMSg4base/CDMSMaterialTable.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Polyhedra.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "Randomize.hh"
#include "CLHEP/Random/RandExponential.h"
#include <cmath>


// Static class member, defined for convenience (since G4ThreeVector doesn't)

const G4ThreeVector Am241SourceHolder::origin;


// Constructor and destructor

Am241SourceHolder::Am241SourceHolder()
  : CDMSVSourceConstruction("Am241Source"), zipBuilder(0), drawSolid(true),
    gammaLines(0), PlateMaterial("G4_Cu"),
    PlateSides(6), PlateRadius(101.35*mm), PlateThickness(1.09*mm),
    CanHeight(11.1*mm), CanThickness(1.6*mm), SourcePuckHeight(6.35*mm),
    LeadDiskThickness(1.6*mm), AlphaFoilThickness(25.4*um), 
    AlphaFoilMaterial("G4_Al"),
    LeadHoleRadius(0.2032*mm), PlateHoleRadius(1*mm),
    NumberOfCans(4), SourceName(4,""), SourceActivity(4,0.),
    CanInnerRadius(4, 0.), SourcePuckRadius(4,0.), ActiveRadius(4,0.),
    ActiveHeight(4,3.175*mm), CanPositionR(4,0.),
    CanPositionPhi(4,0.), HolePositionR(4,0.), HolePositionPhi(4,0.),
    MotherZmin(0.), MotherZmax(0.), sourceAm241(0), 
    ActiveSurface(4,0.), ActiveVisible(4,0.), decayTime(4,0.),
    messenger(new Am241HolderMessenger(this)) {
  // All sources are the same; pointer is for convenience later
  sourceAm241 = G4ParticleTable::GetParticleTable()->GetIon(95,241,0.*keV);

  // Set default configuration for UMN test facility (G101)
  SourceName[0] = "Am241Center";
  SourceName[1] = "Am241LargeQ3";		// FIXME:  Want to use quads
  SourceName[2] = "Am241SmallQ4";
  SourceName[3] = "Am241SmallQ1";

  SourceActivity[0] = SourceActivity[1]	=	// All sources equal strength
    SourceActivity[2] = SourceActivity[3] = 1e-6*curie;

  SourcePuckRadius[0] = SourcePuckRadius[1] = 12.7*mm;		// 1" diam
  SourcePuckRadius[2] = SourcePuckRadius[3] = 4.7625*mm;	// 3/8" diam

  ActiveRadius[0] = ActiveRadius[1] = 2.5*mm;
  ActiveRadius[2] = ActiveRadius[3] = 1.5875*mm;		// 1/4" diam

  // Phi angles chosen because G4 puts corner of polygon at 0, not side
  CanPositionR[0] = 0.*mm;     CanPositionPhi[0] =   0.*deg;	// Center
  CanPositionR[1] = 30.200*mm; CanPositionPhi[1] = 225.*deg;
  CanPositionR[2] = 44.323*mm; CanPositionPhi[2] = 315.*deg;
  CanPositionR[3] = 39.345*mm; CanPositionPhi[3] =  45.*deg;

  HolePositionR[2] = 1.24*mm; HolePositionPhi[2] = 315.*deg;	// Offset hole

  UseZipBuilder(true);		// Default is to always use ZIP housing
}

Am241SourceHolder::~Am241SourceHolder() {
  delete messenger; messenger=0;
  delete gammaLines; gammaLines=0;
}


// Configuration functions

void Am241SourceHolder::SetVerboseLevel(G4int verbose) {
  CDMSVDetectorGeometry::SetVerboseLevel(verbose);
  if (zipBuilder) zipBuilder->SetVerboseLevel(verbose);
  if (gammaLines) gammaLines->SetVerboseLevel(verbose);
}

void Am241SourceHolder::UseZipBuilder(G4bool useZip) {
  if (verboseLevel)
    G4cout << "Am241SourceHolder::UseZipBuilder " << useZip << G4endl;

  zipBuilder = useZip ? CDMSGeometryManager::Instance()->GetZip() : 0;
}

void Am241SourceHolder::UseCDMSGammas(G4bool useGammas) {
  if (verboseLevel)
    G4cout << "Am241SourceHolder::UseCDMSGammas " << useGammas << G4endl;

  if (useGammas && gammaLines) return;	// Already configured as requested
  else delete gammaLines;		// Avoid memory leaks

  gammaLines = useGammas ? new Am241Lines(origin,verboseLevel) : 0;
}


// Generating events once geometry has been configured

void Am241SourceHolder::GeneratePrimaries(G4Event* event) {
  if (verboseLevel>1) G4cout << "Am241CanSource::GeneratePrimaries" << G4endl;

  InitializeGun();		// Configure ParticleGun for Am-241 decays

  G4int iCan = ChooseSourceCanister();		// Choose decay location
  GenerateSourceAtom(iCan);

  particleGun->GeneratePrimaryVertex(event);	// Generate sensible decay
}

void Am241SourceHolder::InitializeGun() {
  if (verboseLevel>1)
    G4cout << "Am241CanSourceHolder::InitializeGun " << G4endl;

  if (gammaLines) {		// For hard-wired gammas, do generation here
    SetParticleType(G4Gamma::Definition());
    gammaLines->SetDirection(GetDirection());
    G4LorentzVector theGamma = gammaLines->shoot();
    particleGun->SetParticleMomentum(theGamma.vect());
    if (verboseLevel>2) G4cout << " generated gamma " << theGamma << G4endl;
  } else {	// Am-241 in holder; RadioactiveDecay will make gammas, alphas
    SetParticleType(sourceAm241);
    particleGun->SetParticleEnergy(0.);
    particleGun->SetParticleMomentumDirection(GetDirection());
  }

  particleGun->SetParticleCharge(0.*eplus);
}

// Select canister for next decay

G4int Am241SourceHolder::ChooseSourceCanister() {
  if (verboseLevel>1)
    G4cout << "Am241SourceHolder::ChooseSourceCanister()" << G4endl;

  G4double randCan = NumberOfCans * G4UniformRand();

  return (G4int)std::floor(randCan);		// Truncate to get 0 to (N-1)
}

// Position and time (seconds) of decay

void Am241SourceHolder::GenerateSourceAtom(G4int iCan) {
  if (iCan < 0 || iCan >= NumberOfCans) return;		// Sanity check

  if (verboseLevel>1)
    G4cout << "Am241SourceHolder::GenerateSourceAtom " << iCan << G4endl;

  particleGun->SetParticlePosition(GenerateDecayPos(iCan));
  particleGun->SetParticleTime(GenerateDecayTime(iCan));
}

G4ThreeVector Am241SourceHolder::GenerateDecayPos(G4int iCan) {
  if (iCan < 0 || iCan >= NumberOfCans) return origin;	// Sanity check

  // Sample disk uniformly (requires Cartesian throws with accept/reject)
  G4double radius = LeadHoleRadius;

  G4double atomX, atomY, atomR;
  do {
    atomX = radius*(2.*G4UniformRand()-1.);	// Sample (-R,+R) square box
    atomY = radius*(2.*G4UniformRand()-1.);
    atomR = std::sqrt(atomX*atomX+atomY*atomY);
    if (verboseLevel>3) 
      G4cout << " trying atom @ (" << atomX << "," << atomY << ") mm"
	     << " R " << atomR << G4endl;
  } while (atomR > radius);

  G4ThreeVector atomPos = GetActivePos(iCan);
  atomPos.setX(atomPos.x() + atomX);	// Locate atom relative to center
  atomPos.setY(atomPos.y() + atomY);

  if (verboseLevel>2) 
    G4cout << " Am-241 decay @ " << atomPos << " mm" << G4endl;

  return atomPos;
}

// Use source activity and relative area to generate exponential times
// NOTE:  GEANT4 uses nanoseconds as basic time unit

G4double Am241SourceHolder::GenerateDecayTime(G4int iCan) {
  if (iCan < 0 || iCan >= NumberOfCans) return -1.;	// Sanity check

  G4double decayTau = (ActiveVisible[iCan]/ActiveSurface[iCan]
		       / SourceActivity[iCan]);
  if (verboseLevel>3) G4cout << " decayTau = " << decayTau << " ns" << G4endl;

  G4double decayInterval = CLHEP::RandExponential::shoot(decayTau);

  decayTime[iCan] += decayInterval;	// Increment timestamp for this source
  if (verboseLevel>2) 
    G4cout << " Am-241 decay @ " << decayTime[iCan] << " ns" << G4endl;

  return decayTime[iCan];
}


// Initialize parameters computed from other objects

void Am241SourceHolder::FillExtraParameters() {
  if (verboseLevel>1)
    G4cout << "Am241SourceHolder::FillExtraParameters()" << G4endl;

  SetVerboseLevel(verboseLevel);	// Ensure that children have verbosity

  if (zipBuilder) {		// Override local values with ZIP configuration
    if (verboseLevel>2)
      G4cout << " Configuring plate to match ZIP housing" << G4endl;

    zipBuilder->FillExtraParameters();	// Make sure everything is up to date

    PlateMaterial = zipBuilder->GetHousingMaterial();
    PlateSides    = zipBuilder->GetHousingSides();
    PlateRadius   = zipBuilder->GetHousingRinner() - tolerance;
  }

  // Compute actual height of canister with bottom lid plus stack
  CanHeight = CanThickness + SourcePuckHeight + 2*LeadDiskThickness
    + AlphaFoilThickness + 3*tolerance;

  CanInnerRadius.resize(NumberOfCans,0.);
  for (G4int i=0; i<NumberOfCans; i++) {
    CanInnerRadius[i] = SourcePuckRadius[i] + tolerance;
  }

  // Create locations for canisters, holes, and collimators (plate)
  CanPos.resize(NumberOfCans,origin);
  HolePos.resize(NumberOfCans,origin);
  CollimatorPos.resize(NumberOfCans,origin);
  SourcePos.resize(NumberOfCans,origin);

  ActiveSurface.resize(NumberOfCans,0.);	// Surface areas for decay times
  ActiveVisible.resize(NumberOfCans,0.);

  // Local coordinate system puts midplane of plate at Z=0
  G4double canZ = PlateThickness/2. + tolerance + CanHeight/2.;
  G4double sourceZ = (PlateThickness/2. + LeadDiskThickness
		      + AlphaFoilThickness + 3*tolerance);

  for (G4int i=0; i<NumberOfCans; i++) {
    CanPos[i].setRhoPhiZ(CanPositionR[i], CanPositionPhi[i], canZ);
    HolePos[i].setRhoPhiZ(HolePositionR[i], HolePositionPhi[i], 0.);

    (CollimatorPos[i]=CanPos[i]+HolePos[i]).setZ(0.);	// Plate coordinates

    // Point on active region aligned with collimator
    SourcePos[i] = CollimatorPos[i];
    SourcePos[i].setZ(sourceZ + (SourcePuckHeight+ActiveHeight[i])/2.);
    SourcePos[i] += GetPosition();

    // Compute surface areas of each source, to be used for event timestamps
    ActiveSurface[i] = 2.*pi*ActiveRadius[i]*ActiveRadius[i]*ActiveHeight[i];
    ActiveVisible[i] = pi*LeadHoleRadius*LeadHoleRadius;
  }

  // Compute (internal) Z coordinate range so plate can be positioned as lid
  MotherZmin = -PlateThickness/2.;	// Midline of plate is X-Y plane
  MotherZmax = PlateThickness/2. + CanHeight;

  // Brand new source, decay event times start at zero
  decayTime.resize(NumberOfCans, 0.*second);
}


// Construct copper end-plate with source canisters

G4LogicalVolume* Am241SourceHolder::BuildGeometry() {
  if (verboseLevel) G4cout << "Am241SourceHolder::BuildGeometry()" << G4endl;

  FillExtraParameters();

  if (verboseLevel>1) PrintParameters(G4cout);

  G4LogicalVolume* theSourceHolder = BuildMother();

  PlacePlate(theSourceHolder);
  for (int i=0; i<NumberOfCans; i++) PlaceCanister(i, theSourceHolder);

  return theSourceHolder;
}

// Empty containment volume, with Z origin in middle of base plate

G4LogicalVolume* Am241SourceHolder::BuildMother() const {
  if (verboseLevel>1) G4cout << "Am241SourceHolder::BuildMother()" << G4endl;

  G4double zEnds[2]  = { MotherZmin-tolerance, MotherZmax+tolerance };
  G4double rInner[2] = { 0.0, 0.0 };
  G4double rOuter[2] = { GetRadius()+tolerance, GetRadius()+tolerance };

  G4Polyhedra* holderShape =
    new G4Polyhedra(GetName(), 0, 360*deg, PlateSides,
		    2, zEnds, rInner, rOuter);

  G4Material* vacuum = CDMSMaterialTable::GetMaterial("G4_Galactic");
  G4LogicalVolume* LogicalSourceHolder = 
    new G4LogicalVolume(holderShape, vacuum, GetName(), 0, 0, 0);
  LogicalSourceHolder->SetVisAttributes(G4VisAttributes::Invisible);

  return LogicalSourceHolder;
}

void Am241SourceHolder::PlacePlate(G4LogicalVolume* mother) const {
  if (verboseLevel>1) G4cout << "Am241SourceHolder::PlacePlate()" << G4endl;

  // Mother volume arranged so that plate is centered at origin
  new G4PVPlacement(0, origin, BuildPlate(), "PhysicalBasePlate", mother,
		    false, 0);
}

G4LogicalVolume* Am241SourceHolder::BuildPlate() const {
  if (verboseLevel>1) G4cout << "Am241SourceHolder::BuildPlate()" << G4endl;

  // Start with simple copper plate, then drill collimator holes
  G4double zEnds[2]  = { -PlateThickness/2., PlateThickness/2. };
  G4double rInner[2] = { 0.0, 0.0 };
  G4double rOuter[2] = { PlateRadius, PlateRadius };

  G4VSolid* basePlate =
    new G4Polyhedra("basePlate", 0, 360*deg, PlateSides,
		    2, zEnds, rInner, rOuter);

  // All collimator holes are identical, but only if sources attached!
  if (NumberOfCans > 0) {
    if (verboseLevel>1)
      G4cout << " adding " << NumberOfCans << " collimator holes" << G4endl;

    G4Tubs* collimator = new G4Tubs("PlateHole", 0., PlateHoleRadius,
				    PlateThickness/2.+tolerance, 0., 360*deg);

    for (G4int i=0; i<NumberOfCans; i++) {
      if (verboseLevel>2)
	G4cout << " collimator " << i << " @ " << CollimatorPos[i] << G4endl;

      basePlate = new G4SubtractionSolid("basePlate", basePlate,
					 collimator, 0, CollimatorPos[i]);
    }
  }

  // Default is copper, but user may change it (must be same as ZIP housing)
  G4Material* copper = CDMSMaterialTable::GetMaterial(PlateMaterial);

  G4VisAttributes* plateAtt =
    SetSolidOrWireFrame(new G4VisAttributes(G4Colour(1.0,0.6,0.4)));

  G4LogicalVolume* LogicalBasePlate = 
    new G4LogicalVolume(basePlate, copper, "LogicalBasePlate", 0, 0, 0);
  LogicalBasePlate->SetVisAttributes(plateAtt);

  return LogicalBasePlate;
}

// Cylindrical can, open at bottom with source and shielding inside

void Am241SourceHolder::PlaceCanister(G4int i, G4LogicalVolume* mother) const {
  if (verboseLevel>1)
    G4cout << "Am241SourceHolder::PlaceCanister " << i << G4endl;

  // Mother volume arranged so that plate is centered at origin
  new G4PVPlacement(0, CanPos[i], BuildCanister(i),
		    "SourceCanister", mother, false, i);
}


// Canister contains lead shields, alpha foil, and source disk

G4LogicalVolume* Am241SourceHolder::BuildCanister(G4int i) const {
  if (verboseLevel>1)
    G4cout << "Am241SourceHolder::BuildCanister " << i << G4endl;

  // Create a mother volume for all components
  G4Tubs* canister = new G4Tubs("Canister", 0., CanInnerRadius[i]+CanThickness,
				CanHeight/2., 0., 360*deg);
  G4Material* vacuum = CDMSMaterialTable::GetMaterial("G4_Galactic");
  G4LogicalVolume* canMother = 
    new G4LogicalVolume(canister, vacuum, "CanisterMother", 0, 0, 0);
  canMother->SetVisAttributes(G4VisAttributes::Invisible);

  // Place can, source and shielding components within mother volume
  G4ThreeVector nextPos = origin;
  new G4PVPlacement(0, nextPos, BuildCopperCan(i), "Canister",
		    canMother, false, i);

  G4double nextZ = -CanHeight/2.;	// Start at bottom of can

  if (LeadDiskThickness > 0.) {
    nextPos.setZ(nextZ+=LeadDiskThickness/2.);
    if (verboseLevel>2) G4cout << " lead collimator @ " << nextPos << G4endl;
    new G4PVPlacement(0, nextPos, BuildLeadDisk(i,true), "LeadCollimator",
		      canMother, false, i);
  }

  if (AlphaFoilThickness > 0.) {
    nextPos.setZ(nextZ+=LeadDiskThickness/2.+tolerance+AlphaFoilThickness/2.);
    if (verboseLevel>2) G4cout << " alpha foil @ " << nextPos << G4endl;
    new G4PVPlacement(0, nextPos, BuildAlphaFoil(i), "AlphaFoil",
		      canMother, false, i);
  }

  nextPos.setZ(nextZ+=AlphaFoilThickness/2.+tolerance+SourcePuckHeight/2.);
  if (verboseLevel>2) G4cout << " source @ " << nextPos << G4endl;
  new G4PVPlacement(0, nextPos, BuildSource(i), "Am241Source",
		    canMother, false, i);

  if (LeadDiskThickness > 0.) {
    nextPos.setZ(nextZ+=SourcePuckHeight/2.+tolerance+LeadDiskThickness/2.);
    if (verboseLevel>2) G4cout << " lead shield @ " << nextPos << G4endl;
    new G4PVPlacement(0, nextPos, BuildLeadDisk(i,false), "LeadShield",
		      canMother, false, i);
  }

  // Canister completed, return
  return canMother;
}

// Copper canister surrounding source and shielding

G4LogicalVolume* Am241SourceHolder::BuildCopperCan(G4int i) const {
  if (verboseLevel>1)
    G4cout << "Am241SourceHolder::BuildCopperCan " << i << G4endl;

  G4VSolid* canister =
    new G4Tubs("Canister", 0., CanInnerRadius[i]+CanThickness,
	       CanHeight/2., 0., 360*deg);

  G4Tubs* canInterior = new G4Tubs("CanInterior", 0., CanInnerRadius[i],
				   (CanHeight-CanThickness)/2., 0., 360*deg);

  G4ThreeVector hollowPos(0., 0., -CanThickness/2.);
  canister = new G4SubtractionSolid("Canister", canister, canInterior, 
				    0, hollowPos);

  // Default is copper, but user may change it (must be same as ZIP housing)
  G4Material* copper = CDMSMaterialTable::GetMaterial(PlateMaterial);

  G4VisAttributes* canAtt = 
    SetSolidOrWireFrame(new G4VisAttributes(G4Colour(1.0,0.6,0.4)));

  G4LogicalVolume* canShell =
    new G4LogicalVolume(canister, copper, "CanisterShell", 0, 0, 0);
  canShell->SetVisAttributes(canAtt);

  return canShell;
}

// Plastic source puck from vendor; active source will be at center

G4LogicalVolume* Am241SourceHolder::BuildSource(G4int i) const {
  if (verboseLevel>1)
    G4cout << "Am241SourceHolder::BuildSource " << i << G4endl;

  G4VSolid* sourcePuck = new G4Tubs("SourcePuck", 0., SourcePuckRadius[i],
				    SourcePuckHeight/2., 0., 360*deg);

  // Plastic disk has recess on bottom side, so source is exposed
  G4double recessH = SourcePuckHeight/2. - ActiveHeight[i]/2.;
  G4double recessR = SourcePuckRadius[i] - 1.5875*mm;	// 1/16" rim
  G4Tubs* recess = new G4Tubs("puckRecess", 0., recessR, recessH/2.,
			      0., 360*deg);

  G4ThreeVector recessPos(0., 0., (recessH-SourcePuckHeight)/2.);
  sourcePuck = new G4SubtractionSolid("SourcePuck", sourcePuck, recess,
				      0, recessPos);

  // NOTE:  Physical source won't have any actual americium in it
  G4Material* acrylic = CDMSMaterialTable::GetMaterial("G4_POLYCARBONATE");

  G4VisAttributes* puckAtt = 
    SetSolidOrWireFrame(new G4VisAttributes(G4Colour(0.75,0.,0.)));

  G4LogicalVolume* logicalSourcePuck = 
    new G4LogicalVolume(sourcePuck, acrylic, "SourcePuck", 0, 0, 0);
  logicalSourcePuck->SetVisAttributes(puckAtt);

  // Overlay active source region for visualization; not used
  AddActiveSource(i, logicalSourcePuck);

  return logicalSourcePuck;
}

void Am241SourceHolder::AddActiveSource(G4int i, G4LogicalVolume* puck) const {
  if (verboseLevel>1)
    G4cout << "Am241SourceHolder::AddActiveSource " << i << G4endl;

  // Active source is at center of plastic puck
  G4LogicalVolume* logicalSource = BuildActiveSource(i);
  new G4PVPlacement(0, origin, logicalSource, "ActiveSource", puck, false, i);
}

G4LogicalVolume* Am241SourceHolder::BuildActiveSource(G4int i) const {
  if (verboseLevel>1)
    G4cout << "Am241SourceHolder::BuildActiveSource " << i << G4endl;

  G4Tubs* source = new G4Tubs("ActiveSource", 0., ActiveRadius[i],
			      ActiveHeight[i]/2., 0., 360*deg);

  // NOTE:  Physical source won't have any actual americium in it
  G4Material* acrylic = CDMSMaterialTable::GetMaterial("G4_POLYCARBONATE");

  G4VisAttributes* srcAtt =
    SetSolidOrWireFrame(new G4VisAttributes(G4Colour(0.95,0.,0.)));

  G4LogicalVolume* logicalSource = 
    new G4LogicalVolume(source, acrylic, "ActiveSource", 0, 0, 0);
  logicalSource->SetVisAttributes(srcAtt);

  return logicalSource;
}


// Shielding disk for alphas

G4LogicalVolume* Am241SourceHolder::BuildAlphaFoil(G4int i) const {
  if (AlphaFoilThickness <= 0.) return 0;	// Foil is not used

  if (verboseLevel>1)
    G4cout << "Am241SourceHolder::BuildAlphaFoil " << i << G4endl;

  G4VSolid* alphaFoil = new G4Tubs("alphaFoil", 0., SourcePuckRadius[i],
				  AlphaFoilThickness/2., 0., 360*deg);

  G4Material* foil = CDMSMaterialTable::GetMaterial(AlphaFoilMaterial);
  if (!foil) {
    G4cerr << "Cannot make alpha foil from " << AlphaFoilMaterial << G4endl;
    return 0;
  }

  G4VisAttributes* foilAtt =
    SetSolidOrWireFrame(new G4VisAttributes(G4Colour(0.56,0.64,0.72)));

  G4LogicalVolume* logicalAlphaFoil = 
    new G4LogicalVolume(alphaFoil, foil, "LogicalAlphaFoil", 0, 0, 0);
  logicalAlphaFoil->SetVisAttributes(foilAtt);

  return logicalAlphaFoil;
}

// Shielding disk, one above and one below (with collimator hole) source

G4LogicalVolume* 
Am241SourceHolder::BuildLeadDisk(G4int i, G4bool hole) const {
  if (LeadDiskThickness <= 0.) return 0;	// Collimator is not used

  if (verboseLevel>1) {
    G4cout << "Am241SourceHolder::BuildLeadDisk " << i << " with"
	   << (hole?"":"out") << " collimator hole" << G4endl;
  }

  G4VSolid* leadDisk = new G4Tubs("leadDisk", 0., SourcePuckRadius[i],
				  LeadDiskThickness/2., 0., 360*deg);

  if (hole && LeadHoleRadius > 0.) {
    if (verboseLevel>2)
      G4cout << " collimator @ " << HolePos[i] << " on disk" << G4endl;

    G4Tubs* collimator = new G4Tubs("diskHole", 0., LeadHoleRadius,
				    LeadDiskThickness/2.+tolerance,0.,360*deg);
    leadDisk = new G4SubtractionSolid("leadDisk", leadDisk, collimator,
				      0, HolePos[i]);
  }

  G4Material* lead = CDMSMaterialTable::GetMaterial("G4_Pb");

  G4VisAttributes* leadAtt = 
    SetSolidOrWireFrame(new G4VisAttributes(G4Colour(0.28,0.32,0.36)));

  G4LogicalVolume* logicalLeadDisk = 
    new G4LogicalVolume(leadDisk, lead, "LogicalLeadDisk", 0, 0, 0);
  logicalLeadDisk->SetVisAttributes(leadAtt);

  return logicalLeadDisk;
}


// Configure drawing attributes for opaque or clear objects
G4VisAttributes*
Am241SourceHolder::SetSolidOrWireFrame(G4VisAttributes* att) const {
  if (drawSolid) att->SetForceSolid(true);
  else att->SetForceWireframe(true);
  return att;
}


// Report configuration

void Am241SourceHolder::PrintParameters(std::ostream& os) const {
  os << "Am241SourceHolder parameters"
     << "\n MotherZmin " << MotherZmin << " MotherZmax " << MotherZmax
     << "\n PlateMaterial " << PlateMaterial
     << "\n PlateSides " << PlateSides
     << "\n PlateRadius " << PlateRadius << " mm"
     << "\n PlateThickness " << PlateThickness << " mm"
     << "\n CanHeight " << CanHeight << " mm"
     << "\n CanThickness " << CanThickness << " mm"
     << "\n SourcePuckHeight " << SourcePuckHeight << " mm"
     << "\n LeadDiskThickness " << LeadDiskThickness << " mm"
     << "\n AlphaFoilThickness " << AlphaFoilThickness << " mm"
     << "\n LeadHoleRadius " << LeadHoleRadius << " mm"
     << "\n PlateHoleRadius " << PlateHoleRadius << " mm"
     << "\n NumberOfCans " << NumberOfCans << " mm";

  os << "\n SourceNames";
  PrintVectorParameter(os, SourceName);

  os << "\n SourceActivity";
  PrintVectorParameter(os, SourceActivity); os << " /ns";

  os << "\n CanInnerRadius";
  PrintVectorParameter(os, CanInnerRadius); os << " mm";

  os << "\n SourcePuckRadius";
  PrintVectorParameter(os, SourcePuckRadius); os << " mm";

  os << "\n ActiveRadius";
  PrintVectorParameter(os, ActiveRadius); os << " mm";

  os << "\n ActiveHeight";
  PrintVectorParameter(os, ActiveHeight); os << " mm";

  os << "\n CanPositionR";
  PrintVectorParameter(os, CanPositionR); os << " mm";

  os << "\n CanPositionPhi";
  PrintVectorParameter(os, CanPositionPhi); os << " mm";

  os << "\n HolePositionR";
  PrintVectorParameter(os, HolePositionR); os << " mm";

  os << "\n HolePositionPhi";
  PrintVectorParameter(os, HolePositionPhi); os << " mm";

  os << "\n CanPos";
  PrintVectorParameter(os, CanPos); os << " mm";

  os << "\n HolePos";
  PrintVectorParameter(os, HolePos); os << " mm";

  os << "\n CollimatorPos";
  PrintVectorParameter(os, CollimatorPos); os << " mm";

  os << "\n SourcePos";
  PrintVectorParameter(os, SourcePos); os << " mm";

  os << G4endl;
}
