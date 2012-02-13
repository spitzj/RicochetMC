#ifndef NoLab_hh
#define NoLab_hh 1
// $Id: NoLab.hh,v 1.8 2011/07/22 21:09:18 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        NoLab.hh                                             //
//  Description: default lab environment (simple vacuum box)          //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        19 May 2010                                          //
//                                                                    //
//  20101026  M. Kelsey -- Drop Messenger and sub-detector args; now  //
//            handed by CDMSGeomConstructor. Change base class        //
//  20101210  M. Kelsey -- Follow base class change to dimensions     //
//  20110722  M. Kelsey -- Remove materials table; use singleon in .cc//
//////////////////////////////////////////////////////////////////////// 

#include "CDMSg4base/CDMSVLabConstruction.hh"
#include "G4ThreeVector.hh"

class G4LogicalVolume;


class NoLab: public CDMSVLabConstruction {
public:
  NoLab();
  virtual ~NoLab();

  virtual G4double GetRadius() const { return radius; }
  virtual G4double GetLength() const { return height; }
  virtual G4LogicalVolume* BuildGeometry();

  virtual G4double GetOverBurdenDensity() const {return ovbDens;}
  virtual G4double GetCavernLength() const {return 1*m;}
  virtual G4double GetCavernWidth() const {return 1*m;}
  virtual G4double GetCavernHeight() const {return 1*m;}

private:
  G4double radius;
  G4double height;
  G4double ovbDens;
};

#endif	/* NoLab_hh */
