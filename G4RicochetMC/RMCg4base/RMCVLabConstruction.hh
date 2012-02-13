#ifndef RMCVLabConstruction_hh
#define RMCVLabConstruction_hh 1
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCVLabConstruction.hh                               //
//  Description: specialized base class for lab geometries            //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Adapted from: Dennis Wright (SLAC)                   //
//  Date:        13 January 2012                                      //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "RMCg4base/RMCVDetectorGeometry.hh"
#include "globals.hh"

class G4LogicalVolume;

class RMCVLabConstruction : public RMCVDetectorGeometry {
public:
    RMCVLabConstruction(const G4String& name="RMCVLab", const G4String& hall="")
    : RMCVDetectorGeometry(name), hallName(hall.empty()?name:hall) {;}
    
    virtual ~RMCVLabConstruction() {;}
    
    // Return space (e.g., interior vault) into which detectors are placed
    virtual void GetLaboratoryHall() const;
    virtual G4double GetRadius() const { return 0.; }	  // Radius of maximum extent
    virtual G4double GetLength() const { return 0.; }	  // Z-length of maximum extent
    
protected:
    G4String hallName;	// Volume name for interior vault space, if any
};

#endif
