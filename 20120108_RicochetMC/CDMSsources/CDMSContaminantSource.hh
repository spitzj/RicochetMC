////////////////////////////////////////////////////////////////////////
// $Id: CDMSContaminantSource.hh,v 1.4 2011/07/08 23:26:34 kelsey Exp $
//  File:        CDMSContaminantSource.hh                             //
//  Description: Base class (fully functional) for sources which are  //
//		 contaminants on/in detector components.              //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        5 July 2011                                          //
//                                                                    //
//  20110707  Drop list of contaminants; use MultiGenerator instead.  //
//////////////////////////////////////////////////////////////////////// 

#ifndef CDMSContaminantSource_hh
#define CDMSContaminantSource_hh 1

#include "CDMSg4base/CDMSVSourceConstruction.hh"
#include "CDMSsources/CDMSMultiGenerator.hh" 
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include <vector>

class CDMSContaminantMessenger;
class G4Event;
class G4LogicalVolume;
class G4RadioactiveDecay;
class G4VPhysicalVolume;


class CDMSContaminantSource : public CDMSVSourceConstruction {
public:
  CDMSContaminantSource(const G4String& name, const G4String& volume="");
  virtual ~CDMSContaminantSource();

  // Create event
  virtual void GeneratePrimaries(G4Event* evt);
  virtual void GeneratePrimary(G4Event* evt);

  // Associated source with detector volume (not done until first event)
  virtual void SetVolumeName(const G4String& val) { volumeName = val; }
  virtual const G4String& GetVolumeName() const { return volumeName; }

  // To generate multiple particles per event
  virtual void SetNumberOfParticles(G4int val=1) { numberOfParticles = val; }
  virtual G4int GetNumberOfParticles() const { return numberOfParticles; }

protected:
  void FindVolumes();			// Fill lists of named volumes 
  void FindLogicalVolumes();
  void FindPhysicalVolumes();

  G4VPhysicalVolume* ChooseVolume() const;	// Select volume from list

  // Determine transformation from local to global coordinates
  G4Transform3D GetWorldFrame(G4VPhysicalVolume* vol) const;

  // Get position and direction for particle, in absolute coordinates
  void ChoosePointAndDir(G4VPhysicalVolume* detector,
			 G4ThreeVector& location,
			 G4ThreeVector& direction) const;

  // Acquire pointer to decay physics process in order to set collimation
  G4RadioactiveDecay* findRadioactiveDecayProcess() const;

private:
  G4String volumeName;			// Name of detector volume for source
  CDMSMultiGenerator sources;		// Collection of sources to use
  G4int numberOfParticles;		// For multiple backgrounds per event
  G4int halfAngle;			// For use with pseudo-collimation

  std::vector<G4LogicalVolume*> lVolumes;
  std::vector<G4VPhysicalVolume*> pVolumes;

  CDMSContaminantMessenger* messenger;
};

#endif	/* CDMSContaminantSource_hh */
