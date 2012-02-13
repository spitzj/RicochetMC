#ifndef RMCGeomConstructor_hh
#define RMCGeomConstructor_hh 1
// $Id: RMCGeomConstructor.hh,v 1.10 2011/07/22 21:09:17 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//  File:        RMCGeomConstructor.hh                                //
//  Description: Wrapper class to build any user-defined geometry     //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               (adapted from Michael Kelsey (SLAC))                 //
//  Date:        8 January 2012                                       //
//////////////////////////////////////////////////////////////////////// 

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include <vector>

class RMCGeomMessenger;
class RMCVDetectorGeometry;
class RMCVLabConstruction;
class G4LogicalVolume;
class G4VPhysicalVolume;


class RMCGeomConstructor : public G4VUserDetectorConstruction {
public:
  RMCGeomConstructor();
  virtual ~RMCGeomConstructor();

  virtual G4VPhysicalVolume* Construct();
  virtual void UpdateGeometry();

  virtual void SetVerboseLevel(G4int verbose=0);
  virtual G4int GetVerboseLevel() const { return verboseLevel; }

  // Specify components of RMC experiment
  virtual void SetLab(RMCVLabConstruction* theLab);
  virtual void SetDetector(RMCVDetectorGeometry* theDet);

  // Users may have multiple sources simultaneously
  virtual void AddSource(RMCVDetectorGeometry* aSource);

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
  RMCGeomMessenger* messenger;

  RMCVLabConstruction* lab;
  RMCVDetectorGeometry* detector;
  std::vector<RMCVDetectorGeometry*> sources;

  G4VPhysicalVolume* theWorld;		// Local cache for use in construction
  G4LogicalVolume* theOutside;		// Not owned by this class!
  G4ThreeVector detectorOffset;

  static const G4ThreeVector theOrigin;	// (0,0,0) for convenience
};

#endif	/* RMCGeomConstructor_hh */
