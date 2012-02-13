// $Id: CDMSVesselMessenger.cc,v 1.8 2011/01/05 23:26:07 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSVesselMessenger.cc                               //
//                                                                    //
//  Description: Messenger class to allow setting CDMS cryostat       //
//               geometry parameters                                  //
//                                                                    //
//  Author:      Michael Kesey (SLAC)                                 //
//  Date:        24 November 2010                                     //
//                                                                    //
//  20101208  M. Kelsey -- Move tower support structure to new class  //
//  20101215  M. Kelsey -- Add number of vacuum vessels               //
//  20101222  M. Kelsey -- Add units category to all parameters       //
//  20101223  M. Kelsey -- Use new vector-input command (CDMS only)   //
////////////////////////////////////////////////////////////////////////

#include "CDMSgeometry/interface/CDMSVesselMessenger.hh"
#include "CDMSgeometry/detectors/CDMSVesselConstruction.hh"
#include "CDMSg4base/CDMS_UIcmdDoublesListAndUnit.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"


// Constructor and destructor

CDMSVesselMessenger::CDMSVesselMessenger(CDMSVesselConstruction* builder)
  : theBuilder(builder) {
  cmdDir = new G4UIdirectory("/CDMS/Vessel/");
  cmdDir->SetGuidance("Define dimensions and parameters for CDMS cryostat");

  NStagesCmd = new G4UIcmdWithAnInteger("/CDMS/Vessel/Stages", this);
  NStagesCmd->SetGuidance("Number of nested cryostats");
  NStagesCmd->SetParameterName("Stages",true,true);
  NStagesCmd->SetDefaultValue(theBuilder->GetNStages());

  NVacuumCmd = new G4UIcmdWithAnInteger("/CDMS/Vessel/VacuumStages", this);
  NVacuumCmd->SetGuidance("Number of (outer) vacuum-pumped cryostats");
  NVacuumCmd->SetParameterName("VacuumStages",true,true);
  NVacuumCmd->SetDefaultValue(theBuilder->GetNVacuumStages());

  VesselThickCmd = new CDMS_UIcmdDoublesListAndUnit("/CDMS/Vessel/Thickness", this);
  VesselThickCmd->SetGuidance("Radial thickness of cryostat cylinders");
  VesselThickCmd->SetParameterName("Thickness",true,true);
  VesselThickCmd->SetUnitCategory("Length");
  VesselThickCmd->SetDefaultValue(theBuilder->GetVesselThick());

  Vessel0RadCmd = new G4UIcmdWithADoubleAndUnit("/CDMS/Vessel/InnerRadius", this);
  Vessel0RadCmd->SetGuidance("Inner radius of innermost cryostat cylinder");
  Vessel0RadCmd->SetParameterName("InnerRadius",true,true);
  Vessel0RadCmd->SetUnitCategory("Length");
  Vessel0RadCmd->SetDefaultValue(theBuilder->GetVessel0Rad());

  VesselDeltaRadCmd = new G4UIcmdWithADoubleAndUnit("/CDMS/Vessel/DeltaRadius", this);
  VesselDeltaRadCmd->SetGuidance("Step in radius between adjacent cryostats");
  VesselDeltaRadCmd->SetParameterName("DeltaRadius",true,true);
  VesselDeltaRadCmd->SetUnitCategory("Length");
  VesselDeltaRadCmd->SetDefaultValue(theBuilder->GetVesselDeltaRad());

  VesselDeltaHeightCmd = new G4UIcmdWithADoubleAndUnit("/CDMS/Vessel/DeltaHeight", this);
  VesselDeltaHeightCmd->SetGuidance("Step in height between adjacent cryostats");
  VesselDeltaHeightCmd->SetParameterName("DeltaHeight",true,true);
  VesselDeltaHeightCmd->SetUnitCategory("Length");
  VesselDeltaHeightCmd->SetDefaultValue(theBuilder->GetVesselDeltaHeight());

  VesselExtraHeightCmd = new G4UIcmdWithADoubleAndUnit("/CDMS/Vessel/ExtraHeight", this);
  VesselExtraHeightCmd->SetGuidance("Additional height step between outer cryostats");
  VesselExtraHeightCmd->SetParameterName("ExtraHeight",true,true);
  VesselExtraHeightCmd->SetUnitCategory("Length");
  VesselExtraHeightCmd->SetDefaultValue(theBuilder->GetVesselExtraHeight());

  VesselGapCmd = new G4UIcmdWithADoubleAndUnit("/CDMS/Vessel/Gap", this);
  VesselGapCmd->SetGuidance("");
  VesselGapCmd->SetParameterName("Gap",true,true);
  VesselGapCmd->SetUnitCategory("Length");
  VesselGapCmd->SetDefaultValue(theBuilder->GetVesselGap());


  LidThicknessCmd = new CDMS_UIcmdDoublesListAndUnit("/CDMS/Vessel/LidThickness", this);
  LidThicknessCmd->SetGuidance("Thickness of cryostat lids");
  LidThicknessCmd->SetParameterName("LidThickness",true,true);
  LidThicknessCmd->SetUnitCategory("Length");
  LidThicknessCmd->SetDefaultValue(theBuilder->GetLidThickness());

  StemRadCmd = new G4UIcmdWithADoubleAndUnit("/CDMS/Vessel/StemRadius", this);
  StemRadCmd->SetGuidance("Cryogen pipe radius (innermost)");
  StemRadCmd->SetParameterName("StemRadius",true,true);
  StemRadCmd->SetUnitCategory("Length");
  StemRadCmd->SetDefaultValue(theBuilder->GetStemInnerRad());

  VacRadCmd = new G4UIcmdWithADoubleAndUnit("/CDMS/Vessel/VacRadius", this);
  VacRadCmd->SetGuidance("Vacuum pipe radius (innermost)");
  VacRadCmd->SetParameterName("VacRadius",true,true);
  VacRadCmd->SetUnitCategory("Length");
  VacRadCmd->SetDefaultValue(theBuilder->GetVacInnerRad());

  PipeThickCmd = new G4UIcmdWithADoubleAndUnit("/CDMS/Vessel/PipeThickness", this);
  PipeThickCmd->SetGuidance("Thickness of all pipe walls");
  PipeThickCmd->SetParameterName("PipeThick",true,true);
  PipeThickCmd->SetUnitCategory("Length");
  PipeThickCmd->SetDefaultValue(theBuilder->GetPipeThick());

  PipeBaseLenCmd = new G4UIcmdWithADoubleAndUnit("/CDMS/Vessel/PipeBaseLength", this);
  PipeBaseLenCmd->SetGuidance("Length of pipe stems outside cryostat");
  PipeBaseLenCmd->SetParameterName("PipeOutsideLength",true,true);
  PipeBaseLenCmd->SetUnitCategory("Length");
  PipeBaseLenCmd->SetDefaultValue(theBuilder->GetPipeOutsideLen());
}

CDMSVesselMessenger::~CDMSVesselMessenger() {
  delete NStagesCmd;
  delete NVacuumCmd;
  delete VesselThickCmd;
  delete Vessel0RadCmd;
  delete VesselDeltaRadCmd;
  delete VesselDeltaHeightCmd;
  delete VesselExtraHeightCmd;
  delete VesselGapCmd;
  delete LidThicknessCmd;
  delete StemRadCmd;
  delete VacRadCmd;
  delete PipeThickCmd;
  delete PipeBaseLenCmd;
  delete cmdDir;
}


// Propagate user input to geometry model

void CDMSVesselMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  if (cmd == NStagesCmd)
    theBuilder->SetNStages(NStagesCmd->GetNewIntValue(value));
  if (cmd == NVacuumCmd)
    theBuilder->SetNVacuumStages(NVacuumCmd->GetNewIntValue(value));

  if (cmd == VesselThickCmd)
    theBuilder->SetVesselThick(VesselThickCmd->GetNewListValue(value));

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

  if (cmd == LidThicknessCmd)
    theBuilder->SetLidThickness(LidThicknessCmd->GetNewListValue(value));

  if (cmd == StemRadCmd)
    theBuilder->SetStemInnerRad(StemRadCmd->GetNewDoubleValue(value));
  if (cmd == VacRadCmd)
    theBuilder->SetVacInnerRad(VacRadCmd->GetNewDoubleValue(value));
  if (cmd == PipeThickCmd)
    theBuilder->SetPipeThick(PipeThickCmd->GetNewDoubleValue(value));
  if (cmd == PipeBaseLenCmd)
    theBuilder->SetPipeOutsideLen(PipeBaseLenCmd->GetNewDoubleValue(value));
}

