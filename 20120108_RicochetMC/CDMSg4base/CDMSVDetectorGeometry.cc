// $Id: CDMSVDetectorGeometry.cc,v 1.2 2011/04/21 22:00:56 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSVDetectorGeometry.hh                             //
//                                                                    //
//  Description: virtual base class for CDMS detector geometries      //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        19 May 2010                                          //
//                                                                    //
//  20101130  M. Kelsey -- Implement reporting function.              //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSg4base/CDMSVDetectorGeometry.hh"
#include <iostream>


G4LogicalVolume* CDMSVDetectorGeometry::BuildGeometryCopy(G4int copy) {
  if (verboseLevel>1)
    G4cout << "CDMSVDetectorGeometry building copy " << copy << G4endl;

  copyNumber = copy;
  return BuildGeometry();
}


void CDMSVDetectorGeometry::PrintParameters(std::ostream& os) const {
  os << " CDMSVDetectorGeometry\n copyNumber " << copyNumber << std::endl;
}
