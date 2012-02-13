#ifndef MITReactor_hh
#define MITReactor_hh 1

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        MITReactor.hh                                        //
//  Description: MIT reactor environment                              //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//  Date:        28 September 2011                                    //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSg4base/CDMSVLabConstruction.hh"
#include "G4ThreeVector.hh"
#include "CDMSgeometry/labs/MITReactorMessenger.hh"

class G4LogicalVolume;
class MITReactorMessenger;

class MITReactor: public CDMSVLabConstruction {
public:
  MITReactor();
  virtual ~MITReactor();

  virtual G4double GetRadius() const {return domeRadius;}
  virtual G4double GetLength() const {return domeHeight;}
  virtual G4LogicalVolume* BuildGeometry();

  virtual G4double GetOverBurden() const {return overBurden;}
  virtual G4double GetCavernLength() const {return cavernLength;}
  virtual G4double GetCavernWidth() const {return vaultRadius;}
  virtual G4double GetCavernHeight() const {return wallHeight + vaultRadius;}

  inline void SetOverBurden(G4double oB) {overBurden = oB;}
  
private:
  G4double overBurden;
  G4double minWallThickness;
  G4double cavernLength;
  G4double vaultRadius;
  G4double wallHeight;

  G4double domeRadius;
  G4double domeHeight;
  G4double domeThickness;

  G4double coreRadius;     // of nuclear reactor
  G4double coreHeight;     // of nuclear reactor

  G4double roomXDim;       // x-dimension of room
  G4double roomYDim;       // y-dimension of room
  G4double roomZDim;       // z-dimension of room
  G4double roomWallThickness; // thickness of room walls

  MITReactorMessenger* theMessenger;
};

#endif
