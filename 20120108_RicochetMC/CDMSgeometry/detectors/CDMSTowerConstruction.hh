#ifndef CDMSTowerConstruction_hh
#define CDMSTowerConstruction_hh 1
// $Id: CDMSTowerConstruction.hh,v 1.12 2011/07/22 21:09:17 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSTowerConstruction.hh                             //
//                                                                    //
//  Description: Construction of ZIP Tower with readout strings       //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        3 November 2010                                      //
//                                                                    //
//  20101130  M. Kelsey -- Override verbosity to pass through value.  //
//		Remove UseZipConstruction().  Add reporting.          //
//  20101203  M. Kelsey -- Make number of ZIPs per tower free param.  //
//  20101208  M. Kelsey -- Get housing parameters directly from ZIP.  //
//  20101210  M. Kelsey -- Follow base class change to dimensions     //
//  20101211  M. Kelsey -- Add full radius without clearance space    //
//  20110204  M. Kelsey -- Add accessor so clients can get ZIP data.  //
//  20110326  M. Kelsey -- Drop accessor; ZIP accessed from Manager.  //
//  20110722  M. Kelsey -- Remove materials table; use singleon in .cc//
//////////////////////////////////////////////////////////////////////// 

#include "CDMSg4base/CDMSVDetectorGeometry.hh"
#include "globals.hh"

class G4LogicalVolume;
class CDMSZipConstruction;
class CDMSTowerMessenger;


class CDMSTowerConstruction : public CDMSVDetectorGeometry {
public:
  CDMSTowerConstruction();
  virtual ~CDMSTowerConstruction();

  virtual G4double GetRadius() const { return towerRadius; }
  virtual G4double GetLength() const { return towerHeight; }

  virtual G4LogicalVolume* BuildGeometry();	// Construct for positioning
  virtual void FillExtraParameters();	// Copy values from ZipConstruction
  virtual void PrintParameters(std::ostream& os) const;

  void SetVerboseLevel(G4int verbose=0);

  void SetNTowerSides(G4int value)     { NTowerSides = value; }
  void SetNZipsPerTower(G4int value)   { NZipsPerTower = value; }
  void SetLidClearance(G4double value) { LidClearance = value; }
  void SetExtraSpace(G4double value)   { ExtraSpace = value; }
  void SetStripThick(G4double value)   { StripThick = value; }
  void SetStripWidth(G4double value)   { StripWidth = value; }

  G4int    GetNTowerSides() const      { return NTowerSides; }
  G4int    GetNZipsPerTower() const    { return NZipsPerTower; }
  G4double GetLidClearance() const     { return LidClearance; }
  G4double GetExtraSpace() const       { return ExtraSpace; }
  G4double GetStripThick() const       { return StripThick; }
  G4double GetStripWidth() const       { return StripWidth; }
  G4double GetStripRadius() const      { return stripRadius; }
  G4double GetZipDeltaHeight() const   { return zipDeltaHeight; }
  G4double GetHousingRadius() const    { return housingRinner; }
  G4double GetHousingThickness() const { return housingThickness; }
  G4double GetHousingHeight() const    { return housingHeight; }
  G4double GetTowerRadius() const      { return towerRadius; }
  G4double GetTowerHeight() const      { return towerHeight; }
  G4double GetTowerClearanceR() const  { return towerClearanceR; }

private:
  virtual G4LogicalVolume* BuildHousing();
  virtual G4LogicalVolume* BuildSideStrip(G4int sideIndex);
  virtual void AttachSideStrip(G4int sideIndex, G4LogicalVolume* zipTower);
  virtual G4double GetSideStripLength(G4int sideIndex);

  enum TowerLidSide { upperLid, lowerLid };
  virtual G4LogicalVolume* BuildLid(TowerLidSide side);

  CDMSZipConstruction* zipBuilder;

  G4int NTowerSides;
  G4int NZipsPerTower;

  G4double LidClearance;
  G4double ExtraSpace;

  G4double StripThick;
  G4double StripWidth;

  G4double zipDeltaHeight;
  G4double housingRinner;
  G4double housingRouter;
  G4double housingThickness;
  G4double housingHeight;
  G4double towerRadius;
  G4double towerClearanceR;
  G4double towerHeight;
  G4double stripRadius;

  CDMSTowerMessenger* messenger;
};

#endif	/* CDMSTowerConstruction_hh */
