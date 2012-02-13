//
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        SnoLabMessenger.cc                                   //
//                                                                    //
//  Description: Messenger class to allow setting SnoLab geometry     //
//               parameters                                           //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        3 May 2011                                           //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "SnoLabMessenger.hh"

#include "SnoLab.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"


SnoLabMessenger::SnoLabMessenger(SnoLab* theSetup)
 :labSetup(theSetup)
{ 
  labDir = new G4UIdirectory("/CDMS/SnoLabGeom/");
  labDir->SetGuidance("UI commands to configure SnoLab geometry");

  OverBurdenCmd = new G4UIcmdWithADoubleAndUnit("/CDMS/SnoLabGeom/setOverBurden", this);
  OverBurdenCmd->SetGuidance("Set thickness of rock above cavern");
  OverBurdenCmd->SetParameterName("Size", false);
  OverBurdenCmd->SetRange("Size>=0.");
  OverBurdenCmd->SetUnitCategory("Length");
  OverBurdenCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  MinWallThicknessCmd = new G4UIcmdWithADoubleAndUnit("/CDMS/SnoLabGeom/setWallThickness", this);
  MinWallThicknessCmd->SetGuidance("Set minimum thickness of rock wall around cavern");
  MinWallThicknessCmd->SetParameterName("Size", false);
  MinWallThicknessCmd->SetRange("Size>=0.");
  MinWallThicknessCmd->SetUnitCategory("Length");
  MinWallThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  CavernLengthCmd = new G4UIcmdWithADoubleAndUnit("/CDMS/SnoLabGeom/setCavernLength", this);
  CavernLengthCmd->SetGuidance("Set length of cavern");
  CavernLengthCmd->SetParameterName("Size", false);
  CavernLengthCmd->SetRange("Size>=0.");
  CavernLengthCmd->SetUnitCategory("Length");
  CavernLengthCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}


SnoLabMessenger::~SnoLabMessenger()
{
  delete OverBurdenCmd;
  delete MinWallThicknessCmd;
  delete CavernLengthCmd;
  delete labDir;
}


void SnoLabMessenger::SetNewValue(G4UIcommand* command, G4String value)
{ 
  if (command == OverBurdenCmd) labSetup->SetOverBurden(OverBurdenCmd->GetNewDoubleValue(value));

  if (command == MinWallThicknessCmd) labSetup->SetMinWallThickness(MinWallThicknessCmd->GetNewDoubleValue(value));

  if (command == CavernLengthCmd) labSetup->SetCavernLength(CavernLengthCmd->GetNewDoubleValue(value));
}
