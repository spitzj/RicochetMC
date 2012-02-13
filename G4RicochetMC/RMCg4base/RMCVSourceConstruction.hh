#ifndef RMCVSourceConstruction_hh
#define RMCVSourceConstruction_hh 1
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCVSourceConstruction.hh                            //
//  Description: base class for geometrically based RMC sources       //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Adapted from: Dennis Wright (SLAC)                   //
//  Date:        2 September 2010                                     //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "RMCg4base/RMCVDetectorGeometry.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"

class G4Event;
class G4ParticleGun;
class G4LogicalVolume;
class G4ParticleDefinition;


class RMCVSourceConstruction
  : public G4VUserPrimaryGeneratorAction, public RMCVDetectorGeometry {
public:
  RMCVSourceConstruction(const G4String& sourceName = "aSource");
  virtual ~RMCVSourceConstruction();

  // For GeneratorAction
  virtual void GeneratePrimaries(G4Event*) = 0;
  virtual void SetParticleType(G4ParticleDefinition*);

  virtual const G4ThreeVector& GetDirection() const { return direction; }
  virtual void SetDirection(const G4ThreeVector& val) { direction = val; }

  // Default behaviour for non-geometric sources (e.g. pure particleGun)
  // Structural sources must override these with proper implementations
  virtual G4double GetRadius() const { return 0.; }
  virtual G4double GetLength() const { return 0.; }
  virtual G4LogicalVolume* BuildGeometry() { return 0; }

protected:
  G4ThreeVector direction;
  G4ParticleGun* particleGun;
};

#endif	/* RMCVSourceConstruction_hh */
