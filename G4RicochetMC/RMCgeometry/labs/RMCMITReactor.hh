#ifndef RMCMITReactor_hh
#define RMCMITReactor_hh 1

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCMITReactor.hh                                     //
//  Description: MIT reactor environment                              //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//  Date:        28 September 2011                                    //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "RMCg4base/RMCVLabConstruction.hh"
#include "G4ThreeVector.hh"
#include "RMCgeometry/labs/RMCMITReactorMessenger.hh"

class G4LogicalVolume;
class RMCMITReactorMessenger;

class RMCMITReactor: public RMCVLabConstruction {
public:
    RMCMITReactor();
    virtual ~RMCMITReactor();
    
    virtual G4LogicalVolume* BuildGeometry();
    virtual G4double GetRadius() const {return domeRadius;};	  // Radius of maximum extent
    virtual G4double GetLength() const {return domeHeight;};	  // Z-length of maximum extent
    
private:
    G4double domeRadius;
    G4double domeHeight;
    G4double domeThickness;
    
    G4double coreRadius;     // of nuclear reactor
    G4double coreHeight;     // of nuclear reactor
    
    
    G4double roomXDim;       // x-dimension of room
    G4double roomYDim;       // y-dimension of room
    G4double roomZDim;       // z-dimension of room
    G4double roomWallThickness; // thickness of room walls
    
    RMCMITReactorMessenger* theMessenger;
};

#endif
