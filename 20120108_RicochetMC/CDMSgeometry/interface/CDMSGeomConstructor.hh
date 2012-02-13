#ifndef CDMSGeomConstructor_hh
#define CDMSGeomConstructor_hh 1
// $Id: CDMSGeomConstructor.hh,v 1.10 2011/07/22 21:09:17 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//  File:        CDMSGeomConstructor.hh                               //
//  Description: Wrapper class to build any user-defined geometry     //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        26 October 2010                                      //
//                                                                    //
//  20101208  M. Kelsey -- Drop partial cleanup actions (not safe).   //
//  20110105  M. Kelsey -- Drop "shield" -- move to MultiTower det.   //
//  20110520  M. Kelsey -- Move SetVerboseLevel to .cc.               //
//  20110722  M. Kelsey -- Remove materials table; use singleon in .cc//
//////////////////////////////////////////////////////////////////////// 

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include <vector>

class CDMSGeomMessenger;
class CDMSVDetectorGeometry;
class CDMSVLabConstruction;
class G4LogicalVolume;
class G4VPhysicalVolume;


class CDMSGeomConstructor : public G4VUserDetectorConstruction {
public:
  CDMSGeomConstructor();
  virtual ~CDMSGeomConstructor();

  virtual G4VPhysicalVolume* Construct();
  virtual void UpdateGeometry();

  virtual void SetVerboseLevel(G4int verbose=0);
  virtual G4int GetVerboseLevel() const { return verboseLevel; }

  // Specify components of CDMS experiment
  virtual void SetLab(CDMSVLabConstruction* theLab);
  virtual void SetDetector(CDMSVDetectorGeometry* theDet);

  // Users may have multiple sources simultaneously
  virtual void AddSource(CDMSVDetectorGeometry* aSource);

protected:
  virtual void ConstructWorld();
  virtual void ConstructLab();
  virtual void ConstructDetector();
  virtual void ConstructSources();
  virtual void ConstructSource(size_t iSrc);

  virtual void ClearGeometry();		// Discard any existing geometry

  virtual G4double GetWorldRadius() const;
  virtual G4double GetWorldLength() const;

private:
  G4int verboseLevel;
  CDMSGeomMessenger* messenger;

  CDMSVLabConstruction* lab;
  CDMSVDetectorGeometry* detector;
  std::vector<CDMSVDetectorGeometry*> sources;

  G4VPhysicalVolume* theWorld;		// Local cache for use in construction
  G4LogicalVolume* theOutside;		// Not owned by this class!
  G4ThreeVector detectorOffset;

  static const G4ThreeVector theOrigin;	// (0,0,0) for convenience
};

#endif	/* CDMSGeomConstructor_hh */
