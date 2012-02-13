////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCVesselMessenger.cc                                //
//                                                                    //
//  Description: Messenger class to allow setting RMC cryostat        //
//               geometry parameters                                  //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Adapted from: Michael Kesey (SLAC)                   //
//  Date:        14 January 2012                                      //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "RMCgeometry/interface/RMCVesselMessenger.hh"
#include "RMCgeometry/detectors/RMCVesselConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"


// Constructor and destructor

RMCVesselMessenger::RMCVesselMessenger(RMCVesselConstruction* builder)
  : theBuilder(builder) {
  cmdDir = new G4UIdirectory("/RMC/Vessel/");
  cmdDir->SetGuidance("Define dimensions and parameters for RMC cryostat");

  NStagesCmd = new G4UIcmdWithAnInteger("/RMC/Vessel/Stages", this);
  NStagesCmd->SetGuidance("Number of nested cryostats");
  NStagesCmd->SetParameterName("Stages",true,true);
  NStagesCmd->SetDefaultValue(theBuilder->GetNStages());

  NVacuumCmd = new G4UIcmdWithAnInteger("/RMC/Vessel/VacuumStages", this);
  NVacuumCmd->SetGuidance("Number of (outer) vacuum-pumped cryostats");
  NVacuumCmd->SetParameterName("VacuumStages",true,true);
  NVacuumCmd->SetDefaultValue(theBuilder->GetNVacuumStages());

  Vessel0RadCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Vessel/InnerRadius", this);
  Vessel0RadCmd->SetGuidance("Inner radius of innermost cryostat cylinder");
  Vessel0RadCmd->SetParameterName("InnerRadius",true,true);
  Vessel0RadCmd->SetUnitCategory("Length");
  Vessel0RadCmd->SetDefaultValue(theBuilder->GetVessel0Rad());

  VesselDeltaRadCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Vessel/DeltaRadius", this);
  VesselDeltaRadCmd->SetGuidance("Step in radius between adjacent cryostats");
  VesselDeltaRadCmd->SetParameterName("DeltaRadius",true,true);
  VesselDeltaRadCmd->SetUnitCategory("Length");
  VesselDeltaRadCmd->SetDefaultValue(theBuilder->GetVesselDeltaRad());

  VesselDeltaHeightCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Vessel/DeltaHeight", this);
  VesselDeltaHeightCmd->SetGuidance("Step in height between adjacent cryostats");
  VesselDeltaHeightCmd->SetParameterName("DeltaHeight",true,true);
  VesselDeltaHeightCmd->SetUnitCategory("Length");
  VesselDeltaHeightCmd->SetDefaultValue(theBuilder->GetVesselDeltaHeight());

  VesselExtraHeightCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Vessel/ExtraHeight", this);
  VesselExtraHeightCmd->SetGuidance("Additional height step between outer cryostats");
  VesselExtraHeightCmd->SetParameterName("ExtraHeight",true,true);
  VesselExtraHeightCmd->SetUnitCategory("Length");
  VesselExtraHeightCmd->SetDefaultValue(theBuilder->GetVesselExtraHeight());

  VesselGapCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Vessel/Gap", this);
  VesselGapCmd->SetGuidance("");
  VesselGapCmd->SetParameterName("Gap",true,true);
  VesselGapCmd->SetUnitCategory("Length");
  VesselGapCmd->SetDefaultValue(theBuilder->GetVesselGap());

  StemRadCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Vessel/StemRadius", this);
  StemRadCmd->SetGuidance("Cryogen pipe radius (innermost)");
  StemRadCmd->SetParameterName("StemRadius",true,true);
  StemRadCmd->SetUnitCategory("Length");
  StemRadCmd->SetDefaultValue(theBuilder->GetStemInnerRad());

  VacRadCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Vessel/VacRadius", this);
  VacRadCmd->SetGuidance("Vacuum pipe radius (innermost)");
  VacRadCmd->SetParameterName("VacRadius",true,true);
  VacRadCmd->SetUnitCategory("Length");
  VacRadCmd->SetDefaultValue(theBuilder->GetVacInnerRad());

  PipeThickCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Vessel/PipeThickness", this);
  PipeThickCmd->SetGuidance("Thickness of all pipe walls");
  PipeThickCmd->SetParameterName("PipeThick",true,true);
  PipeThickCmd->SetUnitCategory("Length");
  PipeThickCmd->SetDefaultValue(theBuilder->GetPipeThick());

  PipeBaseLenCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Vessel/PipeBaseLength", this);
  PipeBaseLenCmd->SetGuidance("Length of pipe stems outside cryostat");
  PipeBaseLenCmd->SetParameterName("PipeOutsideLength",true,true);
  PipeBaseLenCmd->SetUnitCategory("Length");
  PipeBaseLenCmd->SetDefaultValue(theBuilder->GetPipeOutsideLen());
}

RMCVesselMessenger::~RMCVesselMessenger() {
  delete NStagesCmd;
  delete NVacuumCmd;
  delete Vessel0RadCmd;
  delete VesselDeltaRadCmd;
  delete VesselDeltaHeightCmd;
  delete VesselExtraHeightCmd;
  delete VesselGapCmd;
  delete StemRadCmd;
  delete VacRadCmd;
  delete PipeThickCmd;
  delete PipeBaseLenCmd;
  delete cmdDir;
}


// Propagate user input to geometry model

void RMCVesselMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  if (cmd == NStagesCmd)
    theBuilder->SetNStages(NStagesCmd->GetNewIntValue(value));
  if (cmd == NVacuumCmd)
    theBuilder->SetNVacuumStages(NVacuumCmd->GetNewIntValue(value));

  if (cmd == Vessel0RadCmd)
    theBuilder->SetVessel0Rad(Vessel0RadCmd->GetNewDoubleValue(value));
  if (cmd == VesselDeltaRadCmd)
    theBuilder->SetVesselDeltaRad(VesselDeltaRadCmd->GetNewDoubleValue(value));
  if (cmd == VesselDeltaHeightCmd)
    theBuilder->SetVesselDeltaHeight(VesselDeltaHeightCmd->GetNewDoubleValue(value));
  if (cmd == VesselExtraHeightCmd)
    theBuilder->SetVesselExtraHeight(VesselExtraHeightCmd->GetNewDoubleValue(value));
  if (cmd == VesselGapCmd)
    theBuilder->SetVesselGap(VesselGapCmd->GetNewDoubleValue(value));

  if (cmd == StemRadCmd)
    theBuilder->SetStemInnerRad(StemRadCmd->GetNewDoubleValue(value));
  if (cmd == VacRadCmd)
    theBuilder->SetVacInnerRad(VacRadCmd->GetNewDoubleValue(value));
  if (cmd == PipeThickCmd)
    theBuilder->SetPipeThick(PipeThickCmd->GetNewDoubleValue(value));
  if (cmd == PipeBaseLenCmd)
    theBuilder->SetPipeOutsideLen(PipeBaseLenCmd->GetNewDoubleValue(value));
}

