////////////////////////////////////////////////////////////////////////
//  File:        RMCVLabConstruction.cc                               //
//  Description: specialized base class for lab geometries            //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Adapted from: Michael Kelsey (SLAC)                  //
//  Date:        13 January 2012                                      //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "RMCg4base/RMCVLabConstruction.hh"
#include "G4LogicalVolume.hh"
#include <vector>


// Retrieve currently built lab volume from GEANT4 Store
// This is a dummy for now...

//G4LogicalVolume* RMCVLabConstruction::GetLaboratoryHall() const {
void RMCVLabConstruction::GetLaboratoryHall() const {
  if (verboseLevel) 
    G4cout << "RMCVLabConstruction::GetLaboratoryHall " << hallName << G4endl;

  /*static RMCGeometryTools* gtool = RMCGeometryTools::GetInstance();

  // There should be exactly one lab built -- will exit on error
  vlv namedLabs = gtool->GetLogicalVolumes(hallName, 1);

  return namedLabs[0];*/
}
