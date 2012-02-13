// $Id: SnoLab.cc,v 1.4 2011/07/22 21:09:18 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCMITReactor.cc                                     //
//  Description: MIT reactor environment                              //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//  Date:        28 Adam Anderson 2011                                //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "RMCgeometry/labs/RMCMITReactor.hh"
#include "RMCgeometry/labs/RMCMITReactorMessenger.hh"
#include "RMCg4base/RMCMaterialTable.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"


RMCMITReactor::RMCMITReactor()
: RMCVLabConstruction("RMCMITReactor"),
domeRadius(10.*m), domeHeight(10.*m), domeThickness(1.0*cm),
coreRadius(2.*m), coreHeight(5.*m),
roomXDim(2.8*m), roomYDim(2.8*m), roomZDim(2.8*m), roomWallThickness(121.92*cm),
theMessenger(new RMCMITReactorMessenger(this) )
{}

RMCMITReactor::~RMCMITReactor() {}


G4LogicalVolume* RMCMITReactor::BuildGeometry()
{
    G4cout << "RMCMITReactor::BuildGeometry()" << G4endl;
    
    // Get the necessary materials
    G4Material* iron = RMCMaterialTable::GetMaterial("G4_Fe");
    G4Material* air = RMCMaterialTable::GetMaterial("G4_AIR");
    G4Material* concrete = RMCMaterialTable::GetMaterial("G4_CONCRETE");
    G4Material* HDconcrete = new G4Material("HDconcrete", 4.0*g/cm3, 10);
    HDconcrete->AddMaterial(concrete, 1.0);
    
    
    // Build the reactor core
    G4Tubs* reactorCoreCyl = new G4Tubs("reactorCoreCyl", 0., coreRadius, coreHeight/2., 0.*deg, 360.*deg);
    G4LogicalVolume* reactorCoreLog = new G4LogicalVolume(reactorCoreCyl, concrete, "reactorCoreLog");
    G4VisAttributes* reactorCoreAttributes = new G4VisAttributes(G4Color(0.0,1.0,1.0));
    reactorCoreAttributes->SetForceSolid(true);
    reactorCoreLog->SetVisAttributes(reactorCoreAttributes);
    
    
    // Build the fission converter room
    G4Box* fissionConverterWall = new G4Box("fissionConverterWall", (roomXDim/2.0)+roomWallThickness,
                                            (roomYDim/2.0)+roomWallThickness, (roomZDim/2.0)+roomWallThickness);
    G4Box* fissionConverterRoom = new G4Box("fissionConverterRoomInt", (roomXDim+roomWallThickness)/2., roomYDim/2.,
                                            (roomXDim+roomWallThickness)/2.);
    G4SubtractionSolid* fissionConverterRoomSubtr = new G4SubtractionSolid("fissionConverterRoom",
                                                                           fissionConverterWall, fissionConverterRoom,
                                                                           G4Transform3D(G4RotationMatrix(),
                                                                                         G4ThreeVector(-1.*roomWallThickness/2., 0.,
                                                                                                       -1.*roomWallThickness/2.)));
    G4LogicalVolume* fissionConverterRoomLog = new G4LogicalVolume(fissionConverterRoomSubtr, concrete,
                                                                   "fissionConverterRoomLog");
    G4VisAttributes* fissionConverterRoomAttributes = new G4VisAttributes(G4Color(0.0, 0.0, 1.0));
    fissionConverterRoomAttributes->SetForceSolid(true);
    fissionConverterRoomLog->SetVisAttributes(fissionConverterRoomAttributes);
    
    
    // Build the reactor dome
    G4Tubs* reactorDomeCyl = new G4Tubs("reactorDomeCyl", domeRadius, domeRadius+domeThickness,
                                        domeHeight/2., 0.*deg, 360.*deg);
    G4LogicalVolume* reactorDomeLog = new G4LogicalVolume(reactorDomeCyl, iron, "reactorDomeLog");
    G4VisAttributes* domeAttributes = new G4VisAttributes(G4Color(0.0,1.0,1.0));
    domeAttributes->SetForceSolid(true);
    reactorDomeLog->SetVisAttributes(domeAttributes);
    
    
    // Build dome ceiling
    G4Tubs* reactorDomeCeiling = new G4Tubs("reactorDomeCeiling", 0., domeRadius+domeThickness,
                                            domeThickness/2., 0.*deg, 360.*deg);
    G4LogicalVolume* reactorDomeCeilingLog = new G4LogicalVolume(reactorDomeCeiling, iron, "RMCMITReactorDomeCeiling");
    reactorDomeCeilingLog->SetVisAttributes(domeAttributes);
    
    
    // Build the sky volume.  
    G4double skyRadius = sqrt(pow(domeRadius,2) + pow(domeHeight,2));
    G4Sphere* skyHemisphere = new G4Sphere("SkyHemisphere", 0.0, skyRadius, 0.0, 360.0*deg, 0.0, 90.0*deg);
    G4LogicalVolume* skyHemisphereLog = new G4LogicalVolume(skyHemisphere, air, GetName());
    G4VisAttributes* skyAttributes = new G4VisAttributes(G4Color(1., 0., 0.));
    skyHemisphereLog->SetVisAttributes(G4VisAttributes::Invisible);
    
    
    // Insert dome and ceiling in the sky volume
    new G4PVPlacement(0, G4ThreeVector(0., 0., domeHeight/2.), reactorDomeLog,
                      "reactorDomePlace", skyHemisphereLog, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0., 0., domeHeight+(domeThickness/2.)), reactorDomeCeilingLog,
                      "reactorDomeCeilingPlace", skyHemisphereLog, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0., 0., coreHeight/2.), reactorCoreLog,
                      "rectorCorePlace", skyHemisphereLog, false, 0);
    new G4PVPlacement(0, G4ThreeVector(coreRadius+roomWallThickness+(roomXDim/2.0), 0., (roomZDim/2.0)+roomWallThickness),
                      fissionConverterRoomLog, "fissionConverterRoomPlace", skyHemisphereLog, false, 0);
    
    // Set position of rock cylinder so that the origin coincides with the center
    // of the box part of the cavern
    // G4double fudge = 1*m; // fix this when detector support stuff is ready
    position.setX(0.0);
    position.setY(0.0);
    position.setZ(0.0);
    
    return skyHemisphereLog;
}
