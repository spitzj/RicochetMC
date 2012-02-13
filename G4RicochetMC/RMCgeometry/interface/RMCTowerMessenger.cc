////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCTowerMessenger.cc                                //
//                                                                    //
//  Description: Messenger class to allow setting RMC single-tower   //
//               geometry parameters                                  //
//                                                                    //
//  Author:      Michael Kesey (SLAC)                                 //
//  Date:        24 November 2010                                     //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "RMCgeometry/interface/RMCTowerMessenger.hh"
#include "RMCgeometry/detectors/RMCTowerConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"


// Constructor and destructor

RMCTowerMessenger::RMCTowerMessenger(RMCTowerConstruction* builder)
  : theBuilder(builder) {
  cmdDir = new G4UIdirectory("/RMC/Tower/");
  cmdDir->SetGuidance("Configuration of RMC mZIP or iZIP tower");

  NTowerSidesCmd = new G4UIcmdWithAnInteger("/RMC/Tower/Sides", this);
  NTowerSidesCmd->SetGuidance("Number of sides for tower housing");
  NTowerSidesCmd->SetParameterName("Sides",true,true);
  NTowerSidesCmd->SetDefaultValue(theBuilder->GetNTowerSides());

  NTowerZipsCmd = new G4UIcmdWithAnInteger("/RMC/Tower/Zips", this);
  NTowerZipsCmd->SetGuidance("Number of ZIPs to include in tower");
  NTowerZipsCmd->SetParameterName("Zips",true,true);
  NTowerZipsCmd->SetDefaultValue(theBuilder->GetNZipsPerTower());

  LidClearanceCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Tower/LidClearance", this);
  LidClearanceCmd->SetGuidance("Gap between end ZIPs and lids");
  LidClearanceCmd->SetParameterName("LidClearance",true,true);
  LidClearanceCmd->SetUnitCategory("Length");
  LidClearanceCmd->SetDefaultValue(theBuilder->GetLidClearance());

  ExtraSpaceCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Tower/ExtraSpace", this);
  ExtraSpaceCmd->SetGuidance("Minimum gaps between towers");
  ExtraSpaceCmd->SetParameterName("ExtraSpace",true,true);
  ExtraSpaceCmd->SetUnitCategory("Length");
  ExtraSpaceCmd->SetDefaultValue(theBuilder->GetExtraSpace());

  StripThickCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Tower/StripThick", this);
  StripThickCmd->SetGuidance("Radial thickness of coaxial readout cable");
  StripThickCmd->SetParameterName("StripThick",true,true);
  StripThickCmd->SetUnitCategory("Length");
  StripThickCmd->SetDefaultValue(theBuilder->GetStripThick());

  StripWidthCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Tower/StripWidth", this);
  StripWidthCmd->SetGuidance("Width of coaxial readout cable");
  StripWidthCmd->SetParameterName("StripWidth",true,true);
  StripWidthCmd->SetUnitCategory("Length");
  StripWidthCmd->SetDefaultValue(theBuilder->GetStripWidth());
}

RMCTowerMessenger::~RMCTowerMessenger() {
  delete NTowerSidesCmd;
  delete NTowerZipsCmd;
  delete LidClearanceCmd;
  delete ExtraSpaceCmd;
  delete StripThickCmd;
  delete StripWidthCmd;
  delete cmdDir;
}


// Propagate user input to geometry model

void RMCTowerMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  if (cmd == NTowerSidesCmd)
    theBuilder->SetNTowerSides(NTowerSidesCmd->GetNewIntValue(value));
  if (cmd == NTowerZipsCmd)
    theBuilder->SetNZipsPerTower(NTowerZipsCmd->GetNewIntValue(value));
  if (cmd == LidClearanceCmd)
    theBuilder->SetLidClearance(LidClearanceCmd->GetNewDoubleValue(value));
  if (cmd == ExtraSpaceCmd)
    theBuilder->SetExtraSpace(ExtraSpaceCmd->GetNewDoubleValue(value));
  if (cmd == StripThickCmd)
    theBuilder->SetStripThick(StripThickCmd->GetNewDoubleValue(value));
  if (cmd == StripWidthCmd)
    theBuilder->SetStripWidth(StripWidthCmd->GetNewDoubleValue(value));
}

