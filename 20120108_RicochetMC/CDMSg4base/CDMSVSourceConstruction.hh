#ifndef CDMSVSourceConstruction_hh
#define CDMSVSourceConstruction_hh 1
// $Id: CDMSVSourceConstruction.hh,v 1.4 2011/06/29 22:24:54 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSVSourceConstruction.hh                           //
//  Description: base class for geometrically based CDMS sources      //
//               (radiogenic or cosmogenic)                           //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        2 September 2010                                     //
//                                                                    //
//  20101026  M. Kelsey -- Move to CDMSsources, add inheritance from  //
//		CDMSVDetectorGeometry.                                //
//  20101129  M. Kelsey -- Implement ctor and dtor in .cc file        //
//  20101210  M. Kelsey -- Follow geometry class change to dimensions //
//  20110426  M. Kelsey -- Add interface to set beam type             //
//  20110629  M. Kelsey -- Fix constness of ctor arg                  //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSg4base/CDMSVDetectorGeometry.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"

class G4Event;
class G4ParticleGun;
class G4LogicalVolume;
class G4ParticleDefinition;


class CDMSVSourceConstruction
  : public G4VUserPrimaryGeneratorAction, public CDMSVDetectorGeometry {
public:
  CDMSVSourceConstruction(const G4String& sourceName = "aSource");
  virtual ~CDMSVSourceConstruction();

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

#endif	/* CDMSVSourceConstruction_hh */
