#ifndef CDMSZipConstruction_hh
#define CDMSZipConstruction_hh 1
// $Id: CDMSZipConstruction.hh,v 1.11 2011/07/22 21:09:17 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSZipConstruction.hh                               //
//                                                                    //
//  Description: Construction of single ZIP with enclosure            //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        3 November 2010                                      //
//                                                                    //
//  20101203  M. Kelsey -- Add flag to suppress housing construction  //
//  20101208  M. Kelsey -- Move housing parameters to data members    //
//  20101209  M. Kelsey -- Add support for sensitive detector         //
//  20101210  M. Kelsey -- Follow base class change to dimensions     //
//  20110204  M. Kelsey -- Add accessors to get material properties   //
//  20110502  M. Kelsey -- Mother volume with housing must be solid   //
//  20110722  M. Kelsey -- Remove materials table; use singleon in .cc//
//////////////////////////////////////////////////////////////////////// 

#include "CDMSg4base/CDMSVDetectorGeometry.hh"
#include "globals.hh"

class CDMSZipMessenger;
class G4LogicalVolume;
class G4Material;
class G4VSolid;
class G4VSensitiveDetector;


class CDMSZipConstruction : public CDMSVDetectorGeometry {
public:
  CDMSZipConstruction();
  virtual ~CDMSZipConstruction();

  virtual G4double GetRadius() const;	  // Radius of maximum extent
  virtual G4double GetLength() const;	  // Z-length of maximum extent

  virtual G4LogicalVolume* BuildGeometry();	// Construct for positioning

  virtual void FillExtraParameters();		// Store housing parameters

  virtual void PrintParameters(std::ostream& os) const;

  void BuildZipWithHousing(G4bool qhousing=true) {
    makeHousingForZip = qhousing;
  }

  void SetZipMaterial(const G4String& value) { ZipMaterial = value; }
  void SetHousingMaterial(const G4String& value) { HousingMaterial = value; }
  void SetDetectorName(const G4String& value) { DetectorName = value; }

  void SetZipRad(G4double value)           { ZipRad = value; }
  void SetZipThick(G4double value)         { ZipThick = value; }
  void SetZipAxis1Len(G4double value)      { ZipAxis1Len = value; }
  void SetZipAxis2Len(G4double value)      { ZipAxis2Len = value; }
  void SetZipClearanceR(G4double value)    { ZipClearanceR = value; }
  void SetZipClearanceZ(G4double value)    { ZipClearanceZ = value; }
  void SetHousingSides(G4int value)        { HousingSides = value; }
  void SetHousingThickness(G4double value) { HousingThickness = value; }

  G4bool GetMakeHousingForZip() const { return makeHousingForZip; }

  const G4String& GetZipMaterial() const     { return ZipMaterial; }
  const G4String& GetHousingMaterial() const { return HousingMaterial; }
  const G4String& GetDetectorName() const    { return DetectorName; }

  G4double GetZipRad() const           { return ZipRad; }
  G4double GetZipThick() const         { return ZipThick; }
  G4double GetZipAxis1Len() const      { return ZipAxis1Len; }
  G4double GetZipAxis2Len() const      { return ZipAxis2Len; }
  G4double GetZipClearanceR() const    { return ZipClearanceR; }
  G4double GetZipClearanceZ() const    { return ZipClearanceZ; }
  G4int    GetHousingSides() const     { return HousingSides; }
  G4double GetHousingThickness() const { return HousingThickness; }
  G4double GetHousingHeight() const    { return housingHeight; }
  G4double GetHousingRinner() const    { return housingRinner; }
  G4double GetHousingRouter() const    { return housingRouter; }

  G4Material* GetZipG4Material() const;
  G4Material* GetHousingG4Material() const;

private:
  G4LogicalVolume* BuildZip();
  G4LogicalVolume* BuildHousing();
  G4VSolid* BuildHousingShape(G4bool solid=false);	// True fills interior
  G4VSensitiveDetector* BuildSensitiveDetector();

  G4bool makeHousingForZip;	// Flag (= false) to build "naked ZIP"

  G4String ZipMaterial;		// Should be "G4_Si" or "G4_Ge"
  G4String HousingMaterial;	// Should be "G4_Cu"
  G4String DetectorName;	// Name string for sensitive detector

  G4double ZipRad;		// Basic cylindrical crystal geometry
  G4double ZipThick;
  G4double ZipAxis1Len;		// Rectangular "box" defining crystal flats
  G4double ZipAxis2Len;
  G4double ZipClearanceR;	// Space between crystal and housing
  G4double ZipClearanceZ;

  G4int    HousingSides;	// Polygonal housing (set by enclosing Tower)
  G4double HousingThickness;
  G4double housingHeight;
  G4double housingRinner;
  G4double housingRouter;

  CDMSZipMessenger* messenger;
};

#endif	/* CDMSZipConstruction_hh */
