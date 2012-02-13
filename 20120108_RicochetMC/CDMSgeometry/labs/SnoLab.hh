#ifndef SnoLab_hh
#define SnoLab_hh 1
// $Id: SnoLab.hh,v 1.3 2011/07/22 21:09:18 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        SnoLab.hh                                            //
//  Description: SnoLab environment (cavern within norite rock)       //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        29 April 2011                                        //
//                                                                    //
//  20110722  M. Kelsey -- Remove materials table; use singleon in .cc//
//////////////////////////////////////////////////////////////////////// 

#include "CDMSg4base/CDMSVLabConstruction.hh"
#include "G4ThreeVector.hh"

class G4LogicalVolume;
class SnoLabMessenger;


class SnoLab: public CDMSVLabConstruction {
public:
  SnoLab();
  virtual ~SnoLab();

  virtual G4double GetRadius() const {return radius;}
  virtual G4double GetLength() const {return height;}
  virtual G4LogicalVolume* BuildGeometry();

  virtual G4double GetOverBurdenDensity() const {return ovbDens;}
  virtual G4double GetOverBurden() const {return overBurden;}
  virtual G4double GetCavernLength() const {return cavernLength;}
  virtual G4double GetCavernWidth() const {return vaultRadius;}
  virtual G4double GetCavernHeight() const {return wallHeight + vaultRadius;}

  inline void SetOverBurden(G4double oB) {overBurden = oB;}
  inline void SetMinWallThickness(G4double minWT) {minWallThickness = minWT;}
  inline void SetCavernLength(G4double cL) {cavernLength = cL;}
  
private:
  G4double overBurden;
  G4double minWallThickness;
  G4double cavernLength;
  G4double vaultRadius;
  G4double wallHeight;

  G4double radius;   // of rock cylinder
  G4double height;   // of rock cylinder
  G4double ovbDens;  // density of rock

  SnoLabMessenger* theMessenger;
};

#endif
