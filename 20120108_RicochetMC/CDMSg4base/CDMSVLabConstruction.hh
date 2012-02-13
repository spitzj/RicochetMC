#ifndef CDMSVLabConstruction_hh
#define CDMSVLabConstruction_hh 1
// $Id: CDMSVLabConstruction.hh,v 1.2 2011/06/30 23:49:10 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSVLabConstruction.hh                              //
//  Description: specialized base class for lab geometries            //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        1 June 2011                                          //
//                                                                    //
//  20110630  M. Kelsey -- Add name string argument for base class    //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSg4base/CDMSVDetectorGeometry.hh"
#include "G4ThreeVector.hh"


class CDMSVLabConstruction : public CDMSVDetectorGeometry {
  public:
    CDMSVLabConstruction(const G4String& name="CDMSVLab")
      : CDMSVDetectorGeometry(name) {;}

    virtual ~CDMSVLabConstruction() {;}

    virtual G4double GetOverBurdenDensity() const {return 0.;}
    virtual G4double GetOverBurden() const {return 0.;}
    virtual G4double GetCavernLength() const {return 0.;} 
    virtual G4double GetCavernWidth() const {return 0.;} 
    virtual G4double GetCavernHeight() const {return 0.;} 
};

#endif
