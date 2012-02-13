////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCTowerSupportMessenger.cc                         //
//                                                                    //
//  Description: Messenger class to allow setting RMC tower support  //
//               geometry parameters                                  //
//                                                                    //
//  Author:      Michael Kesey (SLAC)                                 //
//  Date:        8 December 2010                                      //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "RMCgeometry/interface/RMCTowerSupportMessenger.hh"
#include "RMCgeometry/detectors/RMCTowerSupportConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"


// Constructor and destructor

RMCTowerSupportMessenger::RMCTowerSupportMessenger(RMCTowerSupportConstruction* builder)
  : theBuilder(builder) {
  cmdDir = new G4UIdirectory("/RMC/Support/");
  cmdDir->SetGuidance("Define dimensions and parameters for RMC tower support");

  NStagesCmd = new G4UIcmdWithAnInteger("/RMC/Support/Stages", this);
  NStagesCmd->SetGuidance("Number of nested cryostats");
  NStagesCmd->SetParameterName("Stages",true,true);
  NStagesCmd->SetRange("Stages>0");
  NStagesCmd->SetDefaultValue(theBuilder->GetNStages());

  StageHeightCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Support/StageHeight", this);
  StageHeightCmd->SetGuidance("Vertical separation between cryostats");
  StageHeightCmd->SetParameterName("StageHeight",true,true);
  StageHeightCmd->SetUnitCategory("Length");
  StageHeightCmd->SetDefaultValue(theBuilder->GetStageHeight());

  StageWallThickCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Support/StageWallThick", this);
  StageWallThickCmd->SetGuidance("Thickness of tower support structure");
  StageWallThickCmd->SetParameterName("StageWallThick",true,true);
  StageWallThickCmd->SetUnitCategory("Length");
  StageWallThickCmd->SetDefaultValue(theBuilder->GetStageWallThick());

  StageGapCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Support/StageGap", this);
  StageGapCmd->SetGuidance("Vertical spacing between tower support stages");
  StageGapCmd->SetParameterName("StageGap",true,true);
  StageGapCmd->SetUnitCategory("Length");
  StageGapCmd->SetDefaultValue(theBuilder->GetStageGap());

  VesselMountStageCmd = new G4UIcmdWithAnInteger("/RMC/Support/VesselMountStage", this);
  VesselMountStageCmd->SetGuidance("Index of stage mounted to innermost cryostat");
  VesselMountStageCmd->SetParameterName("VesselMountStage",true,true);
  VesselMountStageCmd->SetRange("VesselMountStage>=0");
  VesselMountStageCmd->SetDefaultValue(theBuilder->GetVesselMountStage());

  
  SpoolLengthCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Support/SpoolLength", this);
  SpoolLengthCmd->SetGuidance("Length of tower support column (spool)");
  SpoolLengthCmd->SetParameterName("SpoolLength",true,true);
  SpoolLengthCmd->SetUnitCategory("Length");
  SpoolLengthCmd->SetDefaultValue(theBuilder->GetSpoolLength());
 
  SpoolIRCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Support/SpoolIR", this);
  SpoolIRCmd->SetGuidance("Inner radius of tower support");
  SpoolIRCmd->SetParameterName("SpoolIR",true,true);
  SpoolIRCmd->SetUnitCategory("Length");
  SpoolIRCmd->SetDefaultValue(theBuilder->GetSpoolIR());

  SpoolThicknessCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Support/SpoolThickness", this);
  SpoolThicknessCmd->SetGuidance("Radial thickness of tower support");
  SpoolThicknessCmd->SetParameterName("SpoolThickness",true,true);
  SpoolThicknessCmd->SetUnitCategory("Length");
  SpoolThicknessCmd->SetDefaultValue(theBuilder->GetSpoolThickness());

  CTubeIRCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Support/CTubeIR", this);
  CTubeIRCmd->SetGuidance("Inner radius of carbon-fiber support");
  CTubeIRCmd->SetParameterName("CTubeIR",true,true);
  CTubeIRCmd->SetUnitCategory("Length");
  CTubeIRCmd->SetDefaultValue(theBuilder->GetCTubeIR());

  CTubeThickCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Support/CTubeThick", this);
  CTubeThickCmd->SetGuidance("Radial thickness of carbon-fiber support");
  CTubeThickCmd->SetParameterName("CTubeThick",true,true);
  CTubeThickCmd->SetUnitCategory("Length");
  CTubeThickCmd->SetDefaultValue(theBuilder->GetCTubeThick());


  SquidLenCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Support/SquidLen", this);
  SquidLenCmd->SetGuidance("Vertical height of SQUID readout card");
  SquidLenCmd->SetParameterName("SquidLen",true,true);
  SquidLenCmd->SetUnitCategory("Length");
  SquidLenCmd->SetDefaultValue(theBuilder->GetSquidLength());

  SquidThickCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Support/SquidThick", this);
  SquidThickCmd->SetGuidance("Thickness of SQUID readout card");
  SquidThickCmd->SetParameterName("SquidThick",true,true);
  SquidThickCmd->SetUnitCategory("Length");
  SquidThickCmd->SetDefaultValue(theBuilder->GetSquidThick());

  SquidStageCmd = new G4UIcmdWithAnInteger("/RMC/Support/SquidStage", this);
  SquidStageCmd->SetGuidance("Index of stage gap where SQUID card is mounted");
  SquidStageCmd->SetGuidance("Negative index counts down from outside");
  SquidStageCmd->SetParameterName("SquidStage",true,true);
  SquidStageCmd->SetDefaultValue(theBuilder->GetSquidStage());

  FETLenCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Support/FETLen", this);
  FETLenCmd->SetGuidance("Vertical height of FET readout card");
  FETLenCmd->SetGuidance("NOTE: width and thickness are same as SQUID card");
  FETLenCmd->SetParameterName("FETLen",true,true);
  FETLenCmd->SetUnitCategory("Length");
  FETLenCmd->SetDefaultValue(theBuilder->GetFETLength());

  FETStageCmd = new G4UIcmdWithAnInteger("/RMC/Support/FETStage", this);
  FETStageCmd->SetGuidance("Index of stage gap where FET card is mounted");
  FETStageCmd->SetGuidance("Negative index counts down from outside");
  FETStageCmd->SetParameterName("FETStage",true,true);
  FETStageCmd->SetDefaultValue(theBuilder->GetFETStage());
}

RMCTowerSupportMessenger::~RMCTowerSupportMessenger() {
  delete NStagesCmd;
  delete StageHeightCmd;
  delete StageWallThickCmd;
  delete StageGapCmd;
  delete VesselMountStageCmd;
  delete SpoolLengthCmd; 
  delete SpoolIRCmd;
  delete SpoolThicknessCmd;
  delete CTubeIRCmd;
  delete CTubeThickCmd;
  delete SquidLenCmd;
  delete SquidThickCmd;
  delete SquidWidthCmd;
  delete SquidStageCmd;
  delete FETLenCmd;
  delete FETStageCmd;
  delete cmdDir;
}


// Propagate user input to geometry model

void RMCTowerSupportMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  if (cmd == NStagesCmd)
    theBuilder->SetNStages(NStagesCmd->GetNewIntValue(value));
  if (cmd == StageHeightCmd)
    theBuilder->SetStageHeight(StageHeightCmd->GetNewDoubleValue(value));
  if (cmd == StageWallThickCmd)
    theBuilder->SetStageWallThick(StageWallThickCmd->GetNewDoubleValue(value));
  if (cmd == StageGapCmd)
    theBuilder->SetStageGap(StageGapCmd->GetNewDoubleValue(value));
  if (cmd == VesselMountStageCmd)
    theBuilder->SetVesselMountStage(VesselMountStageCmd->GetNewIntValue(value));
  if (cmd == SpoolLengthCmd)
    theBuilder->SetSpoolLength(SpoolLengthCmd->GetNewDoubleValue(value)); 
  if (cmd == SpoolIRCmd)
    theBuilder->SetSpoolIR(SpoolIRCmd->GetNewDoubleValue(value));
  if (cmd == SpoolThicknessCmd)
    theBuilder->SetSpoolThickness(SpoolThicknessCmd->GetNewDoubleValue(value));
  if (cmd == CTubeIRCmd)
    theBuilder->SetCTubeIR(CTubeIRCmd->GetNewDoubleValue(value));
  if (cmd == CTubeThickCmd)
    theBuilder->SetCTubeThick(CTubeThickCmd->GetNewDoubleValue(value));
  if (cmd == SquidLenCmd)
    theBuilder->SetSquidLength(SquidLenCmd->GetNewDoubleValue(value));
  if (cmd == SquidThickCmd)
    theBuilder->SetSquidThick(SquidThickCmd->GetNewDoubleValue(value));
  if (cmd == SquidWidthCmd)
    theBuilder->SetSquidWidth(SquidWidthCmd->GetNewDoubleValue(value));
  if (cmd == SquidStageCmd)
    theBuilder->SetSquidStage(SquidStageCmd->GetNewIntValue(value));
  if (cmd == FETLenCmd)
    theBuilder->SetFETLength(FETLenCmd->GetNewDoubleValue(value));
  if (cmd == FETStageCmd)
    theBuilder->SetFETStage(FETStageCmd->GetNewIntValue(value));
}

