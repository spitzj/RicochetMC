////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCZipMessenger.cc                                  //
//                                                                    //
//  Description: Messenger class to allow setting RMC single-ZIP     //
//               geometry parameters                                  //
//                                                                    //
//  Author:      Michael Kesey (SLAC)                                 //
//  Date:        24 November 2010                                     //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "RMCgeometry/interface/RMCZipMessenger.hh"
#include "RMCgeometry/detectors/RMCZipConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"


// Constructor and destructor

RMCZipMessenger::RMCZipMessenger(RMCZipConstruction* builder)
  : theBuilder(builder) {
  cmdDir = new G4UIdirectory("/RMC/Zip/");
  cmdDir->SetGuidance("Configuration of RMC individual mZIP or iZIP");

  MakeHousingCmd = new G4UIcmdWithABool("/RMC/Zip/MakeHousing", this);
  MakeHousingCmd->SetGuidance("Build housing around ZIP");
  MakeHousingCmd->SetParameterName("makeHousingForZip",true,true);
  MakeHousingCmd->SetDefaultValue(theBuilder->GetMakeHousingForZip());

  ZipMaterialCmd = new G4UIcmdWithAString("/RMC/Zip/ZipMaterial", this);
  ZipMaterialCmd->SetGuidance("ZIP crystal material (usually Ge or Si)");
  ZipMaterialCmd->SetParameterName("ZipMaterial",true,true);
  ZipMaterialCmd->SetDefaultValue(theBuilder->GetZipMaterial());

  HousingMaterialCmd = new G4UIcmdWithAString("/RMC/Zip/HousingMaterial", this);
  HousingMaterialCmd->SetGuidance("ZIP housing material");
  HousingMaterialCmd->SetParameterName("HousingMaterial",true,true);
  HousingMaterialCmd->SetDefaultValue(theBuilder->GetHousingMaterial());

  ZipRadCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Zip/Radius", this);
  ZipRadCmd->SetGuidance("Radius of ZIP crystal");
  ZipRadCmd->SetParameterName("Radius",true,true);
  ZipRadCmd->SetUnitCategory("Length");
  ZipRadCmd->SetDefaultValue(theBuilder->GetZipRad());

  ZipThickCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Zip/Thickness", this);
  ZipThickCmd->SetGuidance("Thickness of ZIP crystal");
  ZipThickCmd->SetParameterName("Thickness",true,true);
  ZipThickCmd->SetUnitCategory("Length");
  ZipThickCmd->SetDefaultValue(theBuilder->GetZipThick());

  ZipAxis1LenCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Zip/Axis1Length", this);
  ZipAxis1LenCmd->SetGuidance("X-axis distance between crystal flats");
  ZipAxis1LenCmd->SetParameterName("Axis1Length",true,true);
  ZipAxis1LenCmd->SetUnitCategory("Length");
  ZipAxis1LenCmd->SetDefaultValue(theBuilder->GetZipAxis1Len());

  ZipAxis2LenCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Zip/Axis2Length", this);
  ZipAxis2LenCmd->SetGuidance("Y-axis distance between crystal flats");
  ZipAxis2LenCmd->SetParameterName("Axis2Length",true,true);
  ZipAxis2LenCmd->SetUnitCategory("Length");
  ZipAxis2LenCmd->SetDefaultValue(theBuilder->GetZipAxis2Len());

  HousingSidesCmd = new G4UIcmdWithAnInteger("/RMC/Zip/HousingSides", this);
  HousingSidesCmd->SetGuidance("Number of sides in ZIP housing");
  HousingSidesCmd->SetParameterName("HousingSides",true,true);
  HousingSidesCmd->SetDefaultValue(theBuilder->GetHousingSides());

  HousingThicknessCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Zip/HousingThickness", this);
  HousingThicknessCmd->SetGuidance("Thickness of housing");
  HousingThicknessCmd->SetParameterName("HousingThickness",true,true);
  HousingThicknessCmd->SetUnitCategory("Length");
  HousingThicknessCmd->SetDefaultValue(theBuilder->GetHousingThickness());

  ZipClearanceRCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Zip/ClearanceR", this);
  ZipClearanceRCmd->SetGuidance("Radial gap between ZIP and housing");
  ZipClearanceRCmd->SetParameterName("ClearanceR",true,true);
  ZipClearanceRCmd->SetUnitCategory("Length");
  ZipClearanceRCmd->SetDefaultValue(theBuilder->GetZipClearanceR());

  ZipClearanceZCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Zip/ClearanceZ", this);
  ZipClearanceZCmd->SetGuidance("Vertical gap between ZIP and housing");
  ZipClearanceZCmd->SetParameterName("ClearanceZ",true,true);
  ZipClearanceZCmd->SetUnitCategory("Length");
  ZipClearanceZCmd->SetDefaultValue(theBuilder->GetZipClearanceZ());

  ZipPositionCmd = new G4UIcmdWith3VectorAndUnit("/RMC/Zip/Position",this);
  ZipPositionCmd->SetGuidance("Set ZIP position in mother volume");
  ZipPositionCmd->SetParameterName("X", "Y", "Z", false);
  ZipPositionCmd->SetUnitCategory("Length");  
  ZipPositionCmd->SetDefaultValue(theBuilder->GetPosition());
}

RMCZipMessenger::~RMCZipMessenger() {
  delete MakeHousingCmd;
  delete ZipMaterialCmd;
  delete HousingMaterialCmd;
  delete ZipRadCmd;
  delete ZipThickCmd;
  delete ZipAxis1LenCmd;
  delete ZipAxis2LenCmd;
  delete HousingSidesCmd;
  delete HousingThicknessCmd;
  delete ZipClearanceRCmd;
  delete ZipClearanceZCmd;
  delete ZipPositionCmd;
  delete cmdDir;
}


// Propagate user input to geometry model

void RMCZipMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  if (cmd == ZipMaterialCmd)     theBuilder->SetZipMaterial(value);
  if (cmd == HousingMaterialCmd) theBuilder->SetHousingMaterial(value);

  if (cmd == MakeHousingCmd)
    theBuilder->BuildZipWithHousing(MakeHousingCmd->GetNewBoolValue(value));

  if (cmd == ZipRadCmd) 
    theBuilder->SetZipRad(ZipRadCmd->GetNewDoubleValue(value));
  if (cmd == ZipThickCmd) 
    theBuilder->SetZipThick(ZipThickCmd->GetNewDoubleValue(value));
  if (cmd == ZipAxis1LenCmd) 
    theBuilder->SetZipAxis1Len(ZipAxis1LenCmd->GetNewDoubleValue(value));
  if (cmd == ZipAxis2LenCmd) 
    theBuilder->SetZipAxis2Len(ZipAxis2LenCmd->GetNewDoubleValue(value));
  if (cmd == HousingSidesCmd) 
    theBuilder->SetHousingSides(HousingSidesCmd->GetNewIntValue(value));
  if (cmd == HousingThicknessCmd) 
    theBuilder->SetHousingThickness(HousingThicknessCmd->GetNewDoubleValue(value));
  if (cmd == ZipClearanceRCmd) 
    theBuilder->SetZipClearanceR(ZipClearanceRCmd->GetNewDoubleValue(value));
  if (cmd == ZipClearanceZCmd) 
    theBuilder->SetZipClearanceZ(ZipClearanceZCmd->GetNewDoubleValue(value));

  if (cmd == ZipPositionCmd)
    theBuilder->SetPosition(ZipPositionCmd->GetNew3VectorValue(value));
}

