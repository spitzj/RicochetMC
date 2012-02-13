#ifndef NoLab_hh
#define NoLab_hh 1
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        NoLab.hh                                             //
//  Description: default lab environment (simple vacuum box)          //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Dennis Wright (SLAC)                                 //
//  Date:        14 January 2012                                      //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "RMCg4base/RMCVLabConstruction.hh"
#include "G4ThreeVector.hh"

class G4LogicalVolume;

class NoLab: public RMCVLabConstruction {
public:
    NoLab();
    virtual ~NoLab();
    
    virtual G4LogicalVolume* BuildGeometry();
    
    virtual G4double GetRadius() const {return radius;}	  // Radius of maximum extent
    virtual G4double GetLength() const {return height;}	  // Z-length of maximum extent
    
private:
    G4double radius;
    G4double height;
};

#endif	/* NoLab_hh */
