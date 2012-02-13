// $Id: NoLab.cc,v 1.10 2011/07/22 21:09:17 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        NoLab.cc                                             //
//  Description: default lab environment (simple vacuum box)          //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        19 May 2010                                          //
//                                                                    //
//  20101026  M. Kelsey -- Drop Messenger and sub-detector args; now  //
//            handed by CDMSGeomConstructor.                          //
//  20101210  M. Kelsey -- Follow base class change to dimensions     //
//  20110106  M. Kelsey -- Enlarge structure to support shielding     //
//  20110630  M. Kelsey -- Set and use name string                    //
//  20110722  M. Kelsey -- Discard matTable, using singleton instead. //
////////////////////////////////////////////////////////////////////////
 
#include "CDMSgeometry/labs/NoLab.hh"
#include "CDMSg4base/CDMSMaterialTable.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"


NoLab::NoLab()
 : CDMSVLabConstruction("NoLab"), radius(2.5*m), height(5.*m), ovbDens(1*g/cm3)
{}

NoLab::~NoLab() {}


G4LogicalVolume* NoLab::BuildGeometry() {
  G4Tubs* labBox = new G4Tubs(GetName(), 0., radius, height/2., 0., 360*deg);
  G4Material* air = CDMSMaterialTable::GetMaterial("G4_AIR");
  G4LogicalVolume* logicalLab = new G4LogicalVolume(labBox, air, GetName());
  ovbDens = air->GetDensity();
  logicalLab->SetVisAttributes(G4VisAttributes::Invisible);

  return logicalLab;
}
