////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCVDetectorGeometry.cc                              //
//                                                                    //
//  Description: virtual base class for detector geometries           //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//  Date:        8 January 2012                                       //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

//#include "CDMSg4base/RMCVDetectorGeometry.hh"
#include "RMCVDetectorGeometry.hh"
#include <iostream>


G4LogicalVolume* RMCVDetectorGeometry::BuildGeometryCopy(G4int copy) {
  if (verboseLevel>1)
    G4cout << "RMCVDetectorGeometry building copy " << copy << G4endl;

  copyNumber = copy;
  return BuildGeometry();
}


void RMCVDetectorGeometry::PrintParameters(std::ostream& os) const {
  os << " RMCVDetectorGeometry\n copyNumber " << copyNumber << std::endl;
}
