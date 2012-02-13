#ifndef Am241SourceHolder_hh
#define Am241SourceHolder_hh 1
////////////////////////////////////////////////////////////////////////
// $Id: Am241SourceHolder.hh,v 1.8 2011/07/22 21:09:17 kelsey Exp $
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
//  20110422  M. Kelsey -- Add names for each active source (Region)  //
//  20110426  M. Kelsey -- Make real source, CDMSVSourceConstruction, //
//		add event-generator functions, source activity (uCi). //
//  20110427  M. Kelsey -- Add missing Set/GetActivity functions, add //
//		flag to select solid vs. wire-frame drawing (for vis) //
//		Option to use hard-coded gamma lines instead of       //
//		RadioactiveDecay database.                            //
//  20110429  M. Kelsey -- Redesign canister as hollow cylinder       //
//  20110623  M. Kelsey -- Add configurable foil material             //
//  20110722  M. Kelsey -- Remove materials table; use singleon in .cc//
//////////////////////////////////////////////////////////////////////// 

#include "CDMSg4base/CDMSVSourceConstruction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include <vector>

class Am241HolderMessenger;
class Am241Lines;
class CDMSZipConstruction;
class G4LogicalVolume;
class G4VSensitiveDetector;
class G4VSolid;
class G4VisAttributes;


class Am241SourceHolder : public CDMSVSourceConstruction {
public:
  Am241SourceHolder();
  virtual ~Am241SourceHolder();

  void GeneratePrimaries(G4Event*);		// Shoot Am-241 decays

  virtual G4double GetRadius() const { return PlateRadius; }
  virtual G4double GetLength() const { return PlateThickness+CanHeight; }

  virtual G4LogicalVolume* BuildGeometry();	// Construct for positioning

  virtual void FillExtraParameters();		// Store housing parameters

  virtual void PrintParameters(std::ostream& os) const;

  virtual void SetVerboseLevel(G4int verbose=0);

  void UseZipBuilder(G4bool useZip=true);	// Get dimensions from housing
  void UseCDMSGammas(G4bool useGammas=true);	// Use hard-coded gamma lines

  // Select opaque vs. wire-frame drawing, in order to see hits
  void DrawSolid(G4bool solid=true) { drawSolid = solid; }
  
  void SetPlateMaterial(const G4String& val) { PlateMaterial = val; }
  void SetPlateSides(const G4int& val) { PlateSides = val; }
  void SetPlateRadius(const G4double& val) { PlateRadius = val; }
  void SetPlateThickness(const G4double& val) { PlateThickness = val; }
  void SetCanHeight(const G4double& val) { CanHeight = val; }
  void SetCanThickness(const G4double& val) { CanThickness = val; }
  void SetSourcePuckHeight(const G4double& val) { SourcePuckHeight = val; }
  void SetLeadDiskThickness(const G4double& val) { LeadDiskThickness = val; }
  void SetAlphaFoilThickness(const G4double& val) { AlphaFoilThickness = val; }
  void SetAlphaFoilMaterial(const G4String& val) { AlphaFoilMaterial = val; }
  void SetLeadHoleRadius(const G4double& val) { LeadHoleRadius = val; }
  void SetPlateHoleRadius(const G4double& val) { PlateHoleRadius = val; }
  void SetNumberOfCans(const G4int& val) { NumberOfCans = val; }
  void SetSourceName(const std::vector<G4String>& val) { SourceName = val; }
  void SetSourceActivity(const std::vector<G4double>& val) { SourceActivity = val; }
  void SetCanInnerRadius(const std::vector<G4double>& val) { CanInnerRadius = val; }
  void SetSourcePuckRadius(const std::vector<G4double>& val) { SourcePuckRadius = val; }
  void SetActiveRadius(const std::vector<G4double>& val) { ActiveRadius = val; }
  void SetActiveHeight(const std::vector<G4double>& val) { ActiveHeight = val; }
  void SetCanPositionR(const std::vector<G4double>& val) { CanPositionR = val; }
  void SetCanPositionPhi(const std::vector<G4double>& val) { CanPositionPhi = val; }
  void SetHolePositionR(const std::vector<G4double>& val) { HolePositionR = val; }
  void SetHolePositionPhi(const std::vector<G4double>& val) { HolePositionPhi = val; }

  const G4String& GetPlateMaterial() const { return PlateMaterial; }
  G4int    GetPlateSides() const { return PlateSides; }
  G4double GetPlateRadius() const { return PlateRadius; }
  G4double GetPlateThickness() const { return PlateThickness; }
  G4double GetCanHeight() const { return CanHeight; }
  G4double GetCanThickness() const { return CanThickness; }
  G4double GetSourcePuckHeight() const { return SourcePuckHeight; }
  G4double GetLeadDiskThickness() const { return LeadDiskThickness; }
  G4double GetAlphaFoilThickness() const { return AlphaFoilThickness; }
  const G4String& GetAlphaFoilMaterial() const { return AlphaFoilMaterial; }
  G4double GetLeadHoleRadius() const { return LeadHoleRadius; }
  G4double GetPlateHoleRadius() const { return PlateHoleRadius; }
  G4int GetNumberOfCans() const { return NumberOfCans; }

  const std::vector<G4String>& GetSourceName() const { return SourceName; }
  const std::vector<G4double>& GetSourceActivity() const { return SourceActivity; }
  const std::vector<G4double>& GetCanInnerRadius() const { return CanInnerRadius; }
  const std::vector<G4double>& GetSourcePuckRadius() const { return SourcePuckRadius; }
  const std::vector<G4double>& GetActiveRadius() const { return ActiveRadius; }
  const std::vector<G4double>& GetActiveHeight() const { return ActiveHeight; }
  const std::vector<G4double>& GetCanPositionR() const { return CanPositionR; }
  const std::vector<G4double>& GetCanPositionPhi() const { return CanPositionPhi; }
  const std::vector<G4double>& GetHolePositionR() const { return HolePositionR; }
  const std::vector<G4double>& GetHolePositionPhi() const { return HolePositionPhi; }

  const G4String& GetSourceName(size_t i) const { return SourceName[i]; }
  G4double GetSourceActivity(size_t i) const { return SourceActivity[i]; }
  G4double GetCanInnerRadius(size_t i) const { return CanInnerRadius[i]; }
  G4double GetSourcePuckRadius(size_t i) const { return SourcePuckRadius[i]; }
  G4double GetActiveRadius(size_t i) const { return ActiveRadius[i]; }
  G4double GetActiveHeight(size_t i) const { return ActiveHeight[i]; }
  G4double GetCanPositionR(size_t i) const { return CanPositionR[i]; }
  G4double GetCanPositionPhi(size_t i) const { return CanPositionPhi[i]; }
  G4double GetHolePositionR(size_t i) const { return HolePositionR[i]; }
  G4double GetHolePositionPhi(size_t i) const { return HolePositionPhi[i]; }

  // Absolute location of active source, for use with ParticleGun
  const G4ThreeVector& GetActivePos(size_t i) const { return SourcePos[i]; }

private:
  // Support functions for generating useful events
  void InitializeGun();
  G4int ChooseSourceCanister();			// Which source decays next?
  void GenerateSourceAtom(G4int iCan);		// Choose point in active area
  G4ThreeVector GenerateDecayPos(G4int iCan);
  G4double GenerateDecayTime(G4int iCan);

  // Support functions for building composite geometry
  G4LogicalVolume* BuildMother() const;
  G4LogicalVolume* BuildPlate() const;
  G4LogicalVolume* BuildCanister(G4int i) const;
  G4LogicalVolume* BuildCopperCan(G4int i) const;
  G4LogicalVolume* BuildSource(G4int i) const;
  G4LogicalVolume* BuildActiveSource(G4int i) const;
  G4LogicalVolume* BuildAlphaFoil(G4int i) const;
  G4LogicalVolume* BuildLeadDisk(G4int i, G4bool bottom=true) const;

  void PlacePlate(G4LogicalVolume* mother) const;
  void PlaceCanister(G4int i, G4LogicalVolume* mother) const;
  void AddActiveSource(G4int i, G4LogicalVolume* puck) const;

  G4VisAttributes* SetSolidOrWireFrame(G4VisAttributes* att) const;

  // Three-vector canister and collimator positions for building
  const G4ThreeVector& GetCanPosition(size_t i) const { return CanPos[i]; }
  const G4ThreeVector& GetHolePosition(size_t i) const { return HolePos[i]; }

  // Absolute (relative to base plate) location of collimator holes
  G4double GetCollimatorR(size_t i) const { return CollimatorPos[i].rho(); }
  G4double GetCollimatorPhi(size_t i) const { return CollimatorPos[i].phi(); }
  const G4ThreeVector& GetCollimatorPos(size_t i) const {
    return CollimatorPos[i];
  }

  CDMSZipConstruction* zipBuilder;	// For dimensions and material

  G4bool drawSolid;			// Opaque vs. wire-frame drawing
  Am241Lines* gammaLines;		// Hard-coded spectrum; null uses RadDcy

  // Mounting plate acts as a (recessed) lid for the ZIP housing
  G4String PlateMaterial;		// Same as ZIP housing material
  G4int    PlateSides;			// Same as ZIP housing sides
  G4double PlateRadius;
  G4double PlateThickness;

  // All source canisters are identical height and structure
  G4double CanHeight;
  G4double CanThickness;
  G4double SourcePuckHeight;
  G4double LeadDiskThickness;
  G4double AlphaFoilThickness;
  G4String AlphaFoilMaterial;		// Aluminum, copper, or Kapton

  G4double LeadHoleRadius;
  G4double PlateHoleRadius;

  G4int NumberOfCans;
  std::vector<G4String> SourceName;	// Used to assign G4Region
  std::vector<G4double> SourceActivity;	// Strength (curies) for event rate

  std::vector<G4double> CanInnerRadius;	// Computed from source + tolerance
  std::vector<G4double> SourcePuckRadius;
  std::vector<G4double> ActiveRadius;
  std::vector<G4double> ActiveHeight;
  std::vector<G4double> CanPositionR;
  std::vector<G4double> CanPositionPhi;
  std::vector<G4double> HolePositionR;	// Zero means centered on can
  std::vector<G4double> HolePositionPhi;

  std::vector<G4ThreeVector> CanPos;
  std::vector<G4ThreeVector> HolePos;
  std::vector<G4ThreeVector> CollimatorPos;	// Absolute  == canPos+holePos
  std::vector<G4ThreeVector> SourcePos;		// Center of active source
  
  G4double MotherZmin;			// Vertical range for coordinates
  G4double MotherZmax;

  G4ParticleDefinition* sourceAm241;	// Avoid repeated function calls
  std::vector<G4double> ActiveSurface;	// Total surface of each source
  std::vector<G4double> ActiveVisible;	// Surface of source at collimator
  std::vector<G4double> decayTime;	// Timestamp of last decay per can

  Am241HolderMessenger* messenger;

  static const G4ThreeVector origin;	// For convenience, no memory churn
};

#endif	/* Am241SourceHolder_hh */
