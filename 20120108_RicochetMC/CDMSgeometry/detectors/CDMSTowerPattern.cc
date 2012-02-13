// $Id: CDMSTowerPattern.cc,v 1.16 2011/05/03 06:13:52 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSTowerPattern.cc                                  //
//                                                                    //
//  Description: Generate coordinates for multiple-tower geometries   //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        5 November 2010                                      //
//                                                                    //
//  20101130  M. Kelsey -- Add verbosity reporting; cache tower ptr.  //
//  20101201  M. Kelsey -- Cache transformation matrices.             //
//  20101211  M. Kelsey -- Add origin offset for general geometries   //
//  20110121  M. Kelsey -- Use enumerator for layout codes, split     //
//		different pattern calculations into functions.  Add   //
//		new close-packed hexagonal layout for large Ntowers.  //
//  20110123  M. Kelsey -- Expand "hexgonal" layout to support offset //
//		pattern with three towers around symetry point.       //
//  20110125  M. Kelsey -- For heaxongal layout, sort positions by    //
//		radius before filling.                                //
//  20110126  M. Kelsey -- Redesign hexagon layouts do allow sorting  //
//		list of positions by radius; move SetCenter() to .cc  //
//  20110204  M. Kelesy -- Generalize layout code selection, add      //
//		utility to work out best-guess tower and ZIP layout   //
//  20110210  M. Kelsey -- Check each parameter for changes, and redo //
//		pattern calculations if necessary, even if usingTower //
//		is same as previous value.                            //
//  20110326  M. Kelsey -- Use Manager to get ZIP information. Change //
//		default numberOfTowers to zero.                       //
//  20110502  M. Kelsey -- Cast std::ceil() to G4int for warnings.    //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSgeometry/detectors/CDMSTowerPattern.hh"
#include "CDMSgeometry/detectors/CDMSTowerConstruction.hh"
#include "CDMSgeometry/detectors/CDMSZipConstruction.hh"
#include "CDMSgeometry/interface/CDMSGeometryManager.hh"
#include "CDMSgeometry/interface/CDMSTowerPatternMessenger.hh"
#include "G4Material.hh"


// Constructor and destructor

CDMSTowerPattern::CDMSTowerPattern(const CDMSTowerConstruction* theTower)
  : verboseLevel(0), layoutCode(RING), numberOfTowers(0), towerSides(12),
    towerRadius(75*mm), towerHeight(350*mm), towerSpace(170*mm), valid(false),
    usingTower(theTower), messenger(new CDMSTowerPatternMessenger(this)) {}

CDMSTowerPattern::~CDMSTowerPattern() {
  delete messenger;
}


// Mapping between layout codes and strings

const char* CDMSTowerPattern::layoutName(CDMSTowerPattern::Layout layout) {
  switch (layout) {
  case RING:      return "ring"; break;
  case RINGROT:   return "ringRot"; break;
  case GRID:      return "grid"; break;
  case HEX:       return "hex"; break;
  case HEXCENTER: return "central"; break;	// Internal use only!
  case HEXCORNER: return "triplet"; break;	// Internal use only!
  default: ;
  }

  G4cerr << "CDMSTowerPattern: ERROR unrecognized layout " << layout << G4endl;
  return 0;	// Not just a null string, but a job killer!
}

CDMSTowerPattern::Layout 
CDMSTowerPattern::layoutCodeFromName(const G4String& name) {
  if (name == "ring")    return RING;
  if (name == "ringRot") return RINGROT;
  if (name == "grid")    return GRID;
  if (name == "hex")     return HEX;

  G4cerr << "CDMSTowerPattern: ERROR invalid layout name " << name << G4endl;
  return UNKNOWN;
}


// Set parameters based on input structure and build

void CDMSTowerPattern::FillTowerParameters() {
  if (verboseLevel>1) 
    G4cout << "CDMSTowerPattern::FillTowerParameters()" << G4endl;

  if (numberOfTowers<1) {		// No towers to be built
    SetAndFlagChange(towerSides, 0);
    SetAndFlagChange(towerHeight, 0.);
    SetAndFlagChange(towerRadius, 0.);
    SetAndFlagChange(towerSpace, 0.);

    position.clear();
    rotation = G4RotationMatrix::IDENTITY;
    return;
  }

  usingTower = CDMSGeometryManager::Instance()->GetTower();

  SetAndFlagChange(towerSides, usingTower->GetNTowerSides());
  SetAndFlagChange(towerHeight, usingTower->GetTowerHeight());
  SetAndFlagChange(towerSpace, 2.*usingTower->GetTowerClearanceR());
  SetAndFlagChange(towerRadius,
		   usingTower->GetStripRadius() + usingTower->GetStripThick()/2.);
}


// Shift origin of coordinates for pattern

void CDMSTowerPattern::SetCenter(const G4ThreeVector& pos) {
  if (valid) {		// Avoid unnecessary effort by just shifting current
    if (verboseLevel>1)
      G4cout << "CDMSTowerPattern::SetCenter moving from " << center
	     << " to " << pos << G4endl;

    for (G4int i=0; i<numberOfTowers; i++) {
      position[i] -= center - pos;	// Move from old center to new (sign!)
      transform[i] = G4Transform3D(rotation, position[i]);
    }
  }

  center = pos;		// Replace reference position
}

void CDMSTowerPattern::SetCenter(G4double x, G4double y, G4double z) {
  SetCenter(G4ThreeVector(x,y,z));
}


// Select layout pattern (including centered or offset hex) based on tower count

CDMSTowerPattern::Layout 
CDMSTowerPattern::ChooseLayoutPattern(G4int nTowers) const {
  if (0 == nTowers) return UNKNOWN;
  if (9 == nTowers && layoutCode != HEX) return GRID;	// 3 x 3 Diamond array

  // Use triplet pattern for exact multiples of 3, centered otherwise
  return (nTowers%3 == 0) ? HEXCORNER : HEXCENTER;
}


// Build pattern based on current parameters

void CDMSTowerPattern::GeneratePattern() {
  if (valid) return;			// Already current, avoid extra work

  if (verboseLevel) G4cout << "CDMSTowerPattern::GeneratePattern()" << G4endl;

  FillTowerParameters();	// Update tower dimensions

  // Discard old information and preload vectors to fill by index
  static const G4ThreeVector origin;

  position.clear();
  position.resize(numberOfTowers, origin);

  transform.clear();
  transform.resize(numberOfTowers, G4Transform3D::Identity);

  if (numberOfTowers < 1) return;	// User requested empty detector

  rotation = G4RotationMatrix::IDENTITY;
  if (layoutCode != RING) rotation.rotateZ(180*deg/towerSides);

  switch (layoutCode) {
  case RING:
  case RINGROT: GenerateRingPositions(); break;
  case GRID:    GenerateGridPositions(); break;
  case HEXCENTER:
  case HEXCORNER:
  case HEX:	GenerateHexagonPositions(); break;
  default: ;
  }

  // Construct G4Transform3D's from rotation and positions, with origin offset
  for (G4int i=0; i<numberOfTowers; i++) {
    position[i] += center;
    transform[i] = G4Transform3D(rotation, position[i]);
  }

  if (verboseLevel > 1) Print();		// Report configuration

  valid = true;
}


void CDMSTowerPattern::GenerateRingPositions() {
  if (verboseLevel>1) G4cout << " calculating ring pattern" << G4endl;

  G4double towerAng = 360*deg / (numberOfTowers-1);
  G4double space = (numberOfTowers==9) ? std::sqrt(2.)*towerSpace : towerSpace;

  for (G4int i=0; i<numberOfTowers-1; i++) {
    position[i].setRhoPhiZ(space,i*towerAng,0.);
  }

  position[numberOfTowers-1].set(0.,0.,0.);	  // Final tower always center
}

void CDMSTowerPattern::GenerateGridPositions() {
  if (verboseLevel>1) G4cout << " calculating diamond pattern" << G4endl;

  G4double towerAng = 360*deg / (numberOfTowers-1);
  G4double r2space = std::sqrt(2.)*towerSpace;

  for (G4int i=0; i<numberOfTowers-1; i+=2) {
    position[i].setRhoPhiZ(r2space,i*towerAng,0.);
    position[i+1].setRhoPhiZ(towerSpace,(i+1)*towerAng,0.);
  }

  position[numberOfTowers-1].set(0.,0.,0.);	  // Final tower always center
}

void CDMSTowerPattern::GenerateHexagonPositions() {
  Layout hexType = ChooseLayoutPattern(numberOfTowers);
  if (hexType == GRID) hexType = HEXCORNER;	// Catch special case

  if (hexType != HEXCENTER && hexType != HEXCORNER) {	// Sanity check
    G4cerr << " GenerateHexagonPositions() got invalid layout " << hexType
	   << G4endl;
    return;
  }

  if (verboseLevel>1) 
    G4cout << " calculating " << layoutName(hexType) << " hexagonal pattern"
	   << G4endl;

  G4int ringSize, nDone, nLeft;		// Size of ring, interior done, residue
  G4int lastRing =
    GetHexRingParameters(numberOfTowers-1, hexType, ringSize, nDone, nLeft);

  std::vector<G4ThreeVector> slots;	// Buffer for tower positions in ring

  G4int iTower=0;		// Index of next tower to be placed
  for (G4int iRing=0; iRing<=lastRing; iRing++) {
    ringSize = GetRingPositions(iRing, hexType, slots);
    nLeft = numberOfTowers - iTower;

    if (verboseLevel>2) 
      G4cout << " ring " << iRing << " has " << ringSize << " slots,"
	     << " vector filled with " << slots.size() << G4endl;

    if (nLeft < ringSize) {		// For partial ring, sort by radius
      if (verboseLevel>2)
	G4cout << " sorting and using first " << nLeft << " slots" << G4endl;

      std::sort(slots.begin(), slots.end(), CompareVectorRadii);
    }

    for (G4int i=0; i<ringSize && i<nLeft; i++) {
      position[iTower] = slots[i];
      iTower++;
    }
  }
}

// Function to allow sorting three-vectors by distance from origin

G4bool CDMSTowerPattern::CompareVectorRadii(const G4ThreeVector& a,
					    const G4ThreeVector& b) {
  return (a.r() < b.r());
}

// Fill list with coordinates of all towers in specified ring

G4int CDMSTowerPattern::GetRingPositions(G4int iRing, Layout hexType,
					 std::vector<G4ThreeVector>& slots) {
  if (hexType != HEXCENTER && hexType != HEXCORNER) {	// Sanity check
    G4cerr << " GetRingPositions() got invalid layout " << hexType << G4endl;
    return 0;
  }

  if (verboseLevel>2) {
    G4cout << " filling tower positions for " << layoutName(hexType)
	   << " hex ring " << iRing << G4endl;
  }

  static const G4ThreeVector origin;
  slots.clear();

  // Handle special-case first -- ring 0 of HEXCENTER is a single tower
  if (iRing == 0 && hexType == HEXCENTER) {
    slots.resize(1, origin);
    return 1;
  }

  // Triplet has 6n+3 per ring; centered hexes have 6n in outer rings
  G4int nRing = 6*iRing + (hexType==HEXCORNER ? 3 : 0);
  slots.resize(nRing, origin);

  for (int iphi=0; iphi<nRing; iphi++) {
    if (hexType==HEXCENTER) PositionInHexCenter(iphi, iRing, slots[iphi]);
    else PositionInHexCorner(iphi, iRing, slots[iphi]);
  }

  return nRing;
}


// Towers positioned at corners of ring and along hexagonal "sides"

void
CDMSTowerPattern::PositionInHexCenter(G4int iphi, G4int iRing,
				      G4ThreeVector& pos) {
  // Each hexagonal ring has "iRing" slots per side

  G4int iside  = iphi / iRing;		// Which side of hexagon are we on?
  G4int ialong = iphi % iRing;		// Which position from corner are we?

  if (verboseLevel>2)
    G4cout << " placing tower in slot " << iphi << " side " << iside
	   << " offset " << ialong << G4endl;

  // For close-packed hexagons, easier to compute (x,y) directly than (rho,phi)
  G4double phiCorner = iside * 60*deg;		// Position angle to corner
  G4double ringRadius = iRing * towerSpace;
  G4double distAlong  = ialong * towerSpace;

  HexagonCoordinates(phiCorner, ringRadius, 120.*deg, distAlong, pos);
}

// Rings are "truncated triangles" with long sides and short sides

void
CDMSTowerPattern::PositionInHexCorner(G4int iphi, G4int iRing,
				      G4ThreeVector& pos) {
  const G4double cornerRad = 0.5*towerSpace/cos(30*deg);

  // Central core is triangle of towers, with reference corner on -Y axis
  if (iRing == 0) {
    G4double cornerDist = cornerRad;
    G4double phiCorner  = (120*iphi - 90) * deg;
    HexagonCoordinates(phiCorner, cornerDist, 0., 0., pos);
    return;
  }

  G4int nSide = 2*iRing + 1;		// Number of slots per "triangle side"

  G4int nlong = iRing + 2;		// Per long side, including corners
  G4int dlong = iRing;			// ... "empty slots" from "long" origin
  
  G4int nshort = nSide - nlong;		// Per short side, excluding corners
  G4int dshort = nlong;			// ... "empty slots" from "short" origin

  G4int longbase  = nlong+2*nshort+1;	// Length of "long side" triangle
  G4int shortbase = nshort+2*nlong-1;	// Length of "short side" triangle

  G4int iside  = iphi / nSide;		// Which "triangle side" to fill
  G4int ialong = iphi % nSide;		// True slot along "triangle side"

  if (verboseLevel>2)
    G4cout << " placing tower in slot " << iphi << " side " << iside
	   << " offset " << ialong << G4endl;

  G4double cornerDist, phiCorner, distAlong;
  if (ialong < nlong) {			// First slots filled along long side
    phiCorner  = (120*iside - 90) * deg;	// "Triangle" points downward
    cornerDist = longbase * cornerRad;
    distAlong  = (ialong+dlong) * towerSpace;
  } else {				// End of each "side" is short side
    phiCorner  = (120*iside - 30) * deg;	// "Triangle" points upward
    cornerDist = shortbase * cornerRad;
    distAlong  = (ialong-nlong+dshort) * towerSpace;
  }

  HexagonCoordinates(phiCorner, cornerDist, 150.*deg, distAlong, pos);
}


// Compute indexing parameters for hexagonal ring currently being filled

G4int 
CDMSTowerPattern::GetHexRingParameters(G4int i,
				       CDMSTowerPattern::Layout hexType, 
				       G4int& ringSize, G4int& filledSlots,
				       G4int& outerSlots) {
  //***  if (!good(i)) return -1;		// Sanity check

  // Find ring index currently being filled (c.f. "filledSlots" calculation)
  // ==> Condition is smallest index such that (iring+1)*(3*iring+nCore) >= i+1

  // Innermost set of towers is "core", indexed as "ring 0"
  G4int nCore = hexType==HEXCENTER ? 1 : hexType==HEXCORNER ? 3 : 0;

  // Solve quadratic eq., take integer just above solution
  G4int nc3 = nCore + 3;
  G4int iring = (G4int)std::ceil((std::sqrt(nc3*nc3 - 12*(nCore-i-1))-nc3) / 6.);

  ringSize = 6*iring + 3*(nCore/3);
  filledSlots = iring * (3*(iring-1) + nCore);
  
  // Determine whether current ring is partial (outermost)
  if (numberOfTowers <= filledSlots) outerSlots = 0;
  else outerSlots = numberOfTowers - filledSlots;

  if (verboseLevel>2) {
    G4cout << " ring " << iring << " has " << ringSize << " slots; "
	   << filledSlots << " completed; "  << outerSlots << " remaining"
	   << G4endl;
  }

  return iring;
}


// Compute coordinates of a hexagon identified by a (rho,phiC) "corner" and
// an offset along a different (phiS) direction

void
CDMSTowerPattern::HexagonCoordinates(G4double phiCorner, G4double ringRadius,
				     G4double phiSide, G4double distAlong,
				     G4ThreeVector& pos) {
  if (verboseLevel>2) {
    G4cout << " phiCorner " << phiCorner/deg << " ringRadius " << ringRadius
	   << " phiSide " << phiSide/deg << " distAlong " << distAlong
	   << G4endl;
  }

  G4double x = ringRadius*cos(phiCorner) + distAlong*cos(phiCorner+phiSide);
  G4double y = ringRadius*sin(phiCorner) + distAlong*sin(phiCorner+phiSide);
  
  if (verboseLevel>2) G4cout << " position " << x << " " << y << G4endl;

  pos.set(x, y, 0.);
}


// Return requested configuration data, regenerating if necessary

const G4RotationMatrix& CDMSTowerPattern::GetRotation() const {
  if (!valid) const_cast<CDMSTowerPattern*>(this)->GeneratePattern();

  return rotation;
}

const G4ThreeVector& CDMSTowerPattern::GetPosition(G4int i) const {
  if (!valid) const_cast<CDMSTowerPattern*>(this)->GeneratePattern();

  return good(i) ? position[i] : center;
}

const G4Transform3D& CDMSTowerPattern::GetTransform3D(G4int i) const {
  if (!valid) const_cast<CDMSTowerPattern*>(this)->GeneratePattern();

  return good(i) ? transform[i] : G4Transform3D::Identity;
}


// Dump current configuraiton for diagnostics

void CDMSTowerPattern::Print() const {
  G4cout << " Layout " << GetLayoutName()
	 << " for " << numberOfTowers << " towers"
	 << " centered at " << center << " mm"
	 << "\n sides " << towerSides << " radius " << towerRadius << " mm"
	 << " height " << towerHeight << " mm"
	 << "\n separation " << towerSpace << "mm"
	 << " rotation " << rotation.delta()/deg << " deg"
	 << G4endl;
  
  for (G4int i=0; i<numberOfTowers; i++) {
    G4cout << " Tower #" << i << " at rho " << position[i].rho()
	   << " mm, phi " << position[i].phi()/deg << " deg, z "
	   << position[i].z() << " mm" << G4endl;
  }
}


// Derive tower parameters with constraints (X are not implemented yet):
//
// 1) Total mass of Ge should be about 100 kg (allowing for integer ZIPs)
// 2) All towers should have same number of ZIPs (consistent construction)
// 3) Tower pattern should be close to circular/convex (no veto shadows)
// X) Maximize number of ZIPs per tower (adjcent surfaces veto contaminants)
// X) Minimize housing material (source of radiogenic neutrons)
//
// Returns FALSE if any of the constraints failed or if no solution

G4bool CDMSTowerPattern::Optimize(G4double activeMass) {
  if (verboseLevel)
    G4cout << "CDMSTowerPattern::Optimize" << activeMass/kg << " kg" << G4endl;

  if (!usingTower) return false;	// Cannot optimize without material
  FillTowerParameters();		// Ensure local parameters up to date

  // Get ZIP information
  const CDMSZipConstruction* theZip = CDMSGeometryManager::Instance()->GetZip();
  if (!theZip) return false;		// Cannot optimize without ZIP info

  // Get mass of single ZIP
  G4Material* zipMat = theZip->GetZipG4Material();
  G4double zipDens = zipMat ? zipMat->GetDensity() : 5.323*g/cm3;
  G4double zipR = theZip->GetRadius();
  G4double zipL = theZip->GetLength();
  G4double zipMass = pi*zipR*zipR * zipL * zipDens;

  // Ideal total number of ZIPs to deploy to reach mass target
  G4int nTotalZIPs = (G4int)std::ceil(activeMass / zipMass);

  if (verboseLevel > 1)
    G4cout << " Each ZIP is " << zipMass/kg << " kg; need " << nTotalZIPs
	   << " for target mass (" << nTotalZIPs*zipMass/kg << " kg" << G4endl;

  // Make initial guess using maximum length towers
  G4int maxZIPs = towerSides;		// Max of one ZIP per side

  Layout layout = layoutCode;
  G4int nZIPs   = maxZIPs;
  G4int nTowers = nTotalZIPs/nZIPs;
  if (nTowers*nZIPs < nTotalZIPs) nTowers++;

  if (verboseLevel > 1)
    G4cout << " Starting with " << nTowers << " towers of " << nZIPs
	   << " ZIPs each" << G4endl;

  // Adjust parameters until successful or impossible
  G4int  trials=0, maxTrials=20;

  G4bool finished = false;
  while (!finished && trials<maxTrials) {
    trials++;				// Number of adjustments attempted
    if (verboseLevel > 1) G4cout << " Starting trial " << trials << G4endl;

    // Adjust ZIPs and towers to approach mass target
    G4int nMissingZIPs = nTotalZIPs - nTowers*nZIPs;
    if (nMissingZIPs == 0) finished = true;
    if (nMissingZIPs != 0) {
      if (verboseLevel > 1)
	G4cout << " Need to " << (nMissingZIPs>0 ? "add " : "remove ")
	       << std::abs(nMissingZIPs) << " ZIPs" << G4endl;

      if (std::abs(nMissingZIPs) > nZIPs) nTowers += nMissingZIPs/nZIPs;
      else if (nMissingZIPs > 0 && nZIPs<maxZIPs) nZIPs++;
      else if (nMissingZIPs < 0) nZIPs--;
      
      if (verboseLevel > 1)
	G4cout << " now trying " << nTowers << " towers of " << nZIPs
	       << " ZIPs each" << G4endl;
    }

    layout = ChooseLayoutPattern(nTowers);	// See if good pattern

    if (GRID == layout || RING == layout) {	// Exact pattern possible
      finished = true;
      break;
    }

    if (HEXCENTER == layout || HEXCORNER == layout) {
      G4int ringSize, nDone, nLeft;	// Size of ring, interior done, residue
      GetHexRingParameters(nTowers-1, layout, ringSize, nDone, nLeft);

      nLeft = nTowers - nDone;			// Filled slots in outer ring
      G4int nEmpty = ringSize - nLeft;		// Empty slots in outer ring

      if (verboseLevel > 1)
	G4cout << " hex " << layoutName(layout) << " has " << nEmpty
	       << " missing positions" << G4endl;

      if (nEmpty >= 3) {	// Trade ZIPs for towers to fill one gap
	nZIPs--;
	nTowers++;
	finished = false;
      }
    }
  }	// while (!finished

  // See if "optimization" was successful or not
  G4bool goodResult = (layout != UNKNOWN) && (nTowers > 0) && (nZIPs > 0);

  if (goodResult) {	// Copy derived values for building detector
    G4cout << " CDMS Detector Optimization successful" << G4endl;

    if (verboseLevel>1) {
      G4cout << " Recommend " << GetLayoutName() << " layout with "
	     << nTowers << " towers of " << nZIPs << " ZIPs each"
	     << " (" << zipR*2. << " x " << zipL << " mm)" 
	     << "\n Target mass " << nTowers*nZIPs*zipMass/kg << " kg"
	     << " vs. goal " << activeMass/kg << " kg" << G4endl;
    }

    layoutCode = layout;
    numberOfTowers = nTowers;
    const_cast<CDMSTowerConstruction*>(usingTower)->SetNZipsPerTower(nZIPs);
  } else {		// Report failure, leave existing parameters alone
    G4cerr << " CDMS Detector Optimization failed" << G4endl;

    if (verboseLevel>1) {
      G4cerr << " After " << trials << " attempts, got " << GetLayoutName()
	     << " layout with " << nTowers << " towers of " << nZIPs
	     << " ZIPs each (" << zipR*2. << " x " << zipL << " mm)"
	     << "\n Target mass " << nTowers*nZIPs*zipMass/kg << " kg"
	     << " vs. goal " << activeMass/kg << " kg"
	     << "\n Please set layout pattern manually." << G4endl;
    }
  }

  return goodResult;
}
