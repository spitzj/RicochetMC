#ifndef CDMSTowerPattern_hh
#define CDMSTowerPattern_hh 1
// $Id: CDMSTowerPattern.hh,v 1.12 2011/03/29 04:50:13 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//  File:        CDMSTowerPattern.hh                                  //
//                                                                    //
//  Description: Generate coordinates for multiple-tower geometries   //
//                                                                    //
// layout = 0 : circular pattern of zips with 1 central zip, tower    //
//  "RING"    : sides face center of pattern.                         //
//            : valid for NTowers = 7, 9                              //
// layout = 1 : same as 0, but with tower vertices facing center.     //
//  "RINGROT" : valid for NTowers = 9                                 //
// layout = 2 : diamond pattern of zips.                              //
//  "GRID"    : valid for NTowers = 9                                 //
// layout = 3 : close packed hexagonal pattern, equivalent to 0 for   //
//  "HEX"     : NTowers=7; valid for any NTowers >= 7                 //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        5 November 2010                                      //
//                                                                    //
//  20101130  M. Kelsey -- Add cache of tower pointer to reduce work. //
//  20101201  M. Kelsey -- Add cache of G4Transform3Ds, return refs.  //
//  20101211  M. Kelsey -- Add origin offset for general geometries   //
//  20110121  M. Kelsey -- Add "hexagonal" layout with multiple rings //
//		and define enum for layouts instead of magic numbers, //
//		and provide enum<->string translations.  Separate     //
//		different pattern calculations into functions.        //
//  20110123  M. Kelsey -- Expand "hexgonal" layout to support offset //
//		pattern with three towers around symetry point.       //
//  20110126  M. Kelsey -- Redesign hexagon layouts do allow sorting  //
//		list of positions by radius; move SetCenter() to .cc  //
//  20110204  M. Kelsey -- Pass layout type to GetRingParameters()    //
//  		Generalize layout code selection, add Optimize(mass)  //
//		utility to work out best-guess tower and ZIP layout   //
//  20110210  M. Kelsey -- Add general function to handle changing    //
//		parameters; split parameter setting from UseTower()   //
//////////////////////////////////////////////////////////////////////// 

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include <vector>

class CDMSTowerConstruction;
class CDMSTowerPatternMessenger;


class CDMSTowerPattern {
public:
  CDMSTowerPattern(const CDMSTowerConstruction* theTower=0);
  virtual ~CDMSTowerPattern();

  void SetVerboseLevel(G4int verbose) { verboseLevel = verbose; }

  // layout = 0 : circular pattern of zips with 1 central zip, tower
  //  "RING"    : sides face center of pattern.
  //            : valid for NTowers = 7, 9
  // layout = 1 : same as 0, but with tower vertices facing center.
  //  "RINGROT" : valid for NTowers = 9
  // layout = 2 : diamond pattern of zips.
  //  "GRID"    : valid for NTowers = 9
  // layout = 3 : close packed hexagonal pattern, equivalent to 0 for
  //  "HEX"     : NTowers=7; valid for any NTowers >= 7
  //            : NOTE: HEXCENTER puts one tower at central point,
  //		:	HEXCORNER puts three towers around central point
  enum Layout { RING, RINGROT, GRID, HEX, HEXCENTER, HEXCORNER, UNKNOWN=-1 };

  // Compute "best guess" tower and ZIP counts and pattern for given mass
  G4bool Optimize(G4double activeMass=100.*kg);

  void GeneratePattern();	// Build or replace position vectors

  void FillTowerParameters();		// Replace parameters using tower data

  void SetLayout(Layout layout)        { SetAndFlagChange(layoutCode, layout); }
  void SetLayout(const G4String& name) { SetLayout(layoutCodeFromName(name)); }

  void SetNumberOfTowers(G4int ntower) { SetAndFlagChange(numberOfTowers, ntower); }
  void SetTowerSides(G4int nsides)     { SetAndFlagChange(towerSides, nsides); }
  void SetTowerRadius(G4double radius) { SetAndFlagChange(towerRadius,radius); }
  void SetTowerHeight(G4double height) { SetAndFlagChange(towerHeight,height); }
  void SetTowerSpace(G4double space)   { SetAndFlagChange(towerSpace, space); }

  void SetCenter(const G4ThreeVector& pos);
  void SetCenter(G4double x=0., G4double y=0., G4double z=0.);

  Layout GetLayout() const        { return layoutCode; }
  const char* GetLayoutName() const { return layoutName(layoutCode); }
 
  G4int GetNumberOfTowers() const { return numberOfTowers; }
  G4int GetTowerSides() const     { return towerSides; }
  G4double GetTowerRadius() const { return towerRadius; }
  G4double GetTowerHeight() const { return towerHeight; }
  G4double GetTowerSpace() const  { return towerSpace; }
  const G4ThreeVector& GetCenter() const { return center; }

  const G4RotationMatrix& GetRotation() const;
  const G4ThreeVector& GetPosition(G4int i) const;
  const G4Transform3D& GetTransform3D(G4int i) const;

  void Print() const;			// Report current configuration

protected:
  // Translate between enumerator layout codes and labels
  static const char* layoutName(Layout layout);
  static Layout layoutCodeFromName(const G4String& name);

  // Special function to copy value into data member and flag if different
  template <typename T>
  void SetAndFlagChange(T& param, const T& value) {
    if (value != param) reset(); param = value;
  }

  void reset() { valid = false; }	// Trigger to regenerate positions

  bool good(G4int i) const { return (i>=0 && i<numberOfTowers); }

  // Select "best" layout pattern based on number of towers
  Layout ChooseLayoutPattern(G4int nTowers) const;

  void GenerateRingPositions();		// Fill position vectors for patterns
  void GenerateGridPositions();
  void GenerateHexagonPositions();

  // Utility functions to compute hexagonal close-packing parameters
  
  // Get parameters for ring currently being filled (nCore=0 or 3)
  G4int GetHexRingParameters(G4int i, Layout hexType, G4int& ringSize,
			     G4int& filledSlots, G4int& outerSlots);

  // Fill input buffer (slots) with coordinates of all towers in ring
  G4int GetRingPositions(G4int iRing, Layout hexType,
			 std::vector<G4ThreeVector>& slots);

  // Compute coorindates of hexagon in specified location around ring
  void PositionInHexCenter(G4int islot, G4int iRing, G4ThreeVector& pos);
  void PositionInHexCorner(G4int islot, G4int iRing, G4ThreeVector& pos);
 
  // Compute coordinates of a hexagon identified by a (rho,phiC) "corner" and
  // an offset along a different (phiS) direction
  void HexagonCoordinates(G4double phiCorner, G4double ringRadius,
			  G4double phiSide, G4double distAlong,
			  G4ThreeVector& pos);

public:
  // Must be public for use as argument to std::sort
  static G4bool CompareVectorRadii(const G4ThreeVector& a,
				   const G4ThreeVector& b);

private:
  G4int verboseLevel;
  
  Layout layoutCode;
  G4int numberOfTowers;
  G4int towerSides;
  G4double towerRadius;
  G4double towerHeight;
  G4double towerSpace;
  G4ThreeVector center;

  G4bool valid;				// Flag if position data is valid
  G4RotationMatrix rotation;
  std::vector<G4ThreeVector> position;
  std::vector<G4Transform3D> transform;

  const CDMSTowerConstruction* usingTower;	// Cache to reduce computations

  CDMSTowerPatternMessenger *messenger;
};

#endif	/* CDMSTowerPattern_hh */
