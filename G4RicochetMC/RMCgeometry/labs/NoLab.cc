////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        NoLab.cc                                             //
//  Description: default lab environment (simple vacuum box)          //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Dennis Wright (SLAC)                                 //
//  Date:        14 January 2012                                      //
//                                                                    //
////////////////////////////////////////////////////////////////////////
 
#include "RMCgeometry/labs/NoLab.hh"
#include "RMCg4base/RMCMaterialTable.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"


NoLab::NoLab()
 : RMCVLabConstruction("NoLab"), radius(2.5*m), height(5.*m)
{ }

NoLab::~NoLab() {}


G4LogicalVolume* NoLab::BuildGeometry() {
  G4Tubs* labBox = new G4Tubs(GetName(), 0., radius, height/2., 0., 360*deg);
  G4Material* air = RMCMaterialTable::GetMaterial("G4_AIR");
  G4LogicalVolume* logicalLab = new G4LogicalVolume(labBox, air, GetName());
  logicalLab->SetVisAttributes(G4VisAttributes::Invisible);

  return logicalLab;
}
