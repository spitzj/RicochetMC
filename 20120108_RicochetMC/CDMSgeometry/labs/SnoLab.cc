// $Id: SnoLab.cc,v 1.4 2011/07/22 21:09:18 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        SnoLab.cc                                            //
//  Description: SnoLab environment (cavern within norite rock)       //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        29 April 2011                                        //
//                                                                    //
//  20110630  M. Kelsey -- Set and use name string                    //
//  20110722  M. Kelsey -- Discard matTable, using singleton instead. //
////////////////////////////////////////////////////////////////////////
 
#include "CDMSgeometry/labs/SnoLab.hh"
#include "CDMSgeometry/labs/SnoLabMessenger.hh"
#include "CDMSg4base/CDMSMaterialTable.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4UnionSolid.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"


SnoLab::SnoLab()
 : CDMSVLabConstruction("SnoLab"), overBurden(5.*m), minWallThickness(10.*m),
   cavernLength(10.*m), vaultRadius(3.94*m), wallHeight(3.13*m),
   ovbDens(1*g/cm3), theMessenger(new SnoLabMessenger(this) )
{
  // Set these here for the world volume, then re-build if params change 
  radius = sqrt(vaultRadius*vaultRadius + cavernLength*cavernLength/4.)
         + minWallThickness;
  height = wallHeight + vaultRadius + 2.*overBurden;
}

SnoLab::~SnoLab() {}


G4LogicalVolume* SnoLab::BuildGeometry()
{
  radius = sqrt(vaultRadius*vaultRadius + cavernLength*cavernLength/4.)
         + minWallThickness;
  height = wallHeight + vaultRadius + 2.*overBurden;

  // Rock surrounding the cavern
  G4Tubs* labCyl = new G4Tubs(GetName(), 0., radius, height/2., 0., 360*deg);
  G4Material* norite = CDMSMaterialTable::GetMaterial("Norite");
  ovbDens = norite->GetDensity();
  G4LogicalVolume* logicalLab = new G4LogicalVolume(labCyl, norite, GetName());
  G4VisAttributes* rockAtt = new G4VisAttributes(G4Color(1.0,0.0,0.0));
  rockAtt->SetForceWireframe(true);
  logicalLab->SetVisAttributes(rockAtt);
  //  logicalLab->SetVisAttributes(G4VisAttributes::Invisible);

  // Build the cavern
  G4double eps = 1*mm;
  G4Tubs* cavernCyl = new G4Tubs("CavCyl", 0., vaultRadius, cavernLength/2. - eps,
                                 -95*deg, 190*deg);
  G4Box* cavernBox = new G4Box("CavBox", cavernLength/2., vaultRadius, wallHeight/2.);
  G4RotationMatrix* rot = new G4RotationMatrix;
  rot->rotateY(90*deg);
  G4ThreeVector trans(0, 0, wallHeight/2.);
  G4UnionSolid* cavern = new G4UnionSolid("Cavern", cavernBox, cavernCyl, rot, trans);
  G4Material* air = CDMSMaterialTable::GetMaterial("G4_AIR");
  G4LogicalVolume* logicalCavern = new G4LogicalVolume(cavern, air, "SnoLabCavern");
  G4VisAttributes* cavernAtt = new G4VisAttributes(G4Color(0.0,1.0,1.0));
  cavernAtt->SetForceSolid(true);
  logicalCavern->SetVisAttributes(cavernAtt);

  // Insert cavern in rock
  G4double zPos = -height/2. + overBurden + wallHeight/2.;
  new G4PVPlacement(0, G4ThreeVector(0., 0., zPos), logicalCavern,
                      "SnoLabCavern", logicalLab, false, 0);

  // Set position of rock cylinder so that the origin coincides with the center
  // of the box part of the cavern
  // G4double fudge = 1*m; // fix this when detector support stuff is ready
  position.setX(0.0);
  position.setY(0.0);
  position.setZ(0.0);

  return logicalLab;
}
