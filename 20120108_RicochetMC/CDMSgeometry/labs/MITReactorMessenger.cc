//
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        MITReactorMessenger.cc                               //
//                                                                    //
//  Description: Messenger class to allow setting MITReactor geometry //
//               parameters                                           //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//  Date:        29 September 2011                                    //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "MITReactorMessenger.hh"

#include "MITReactor.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"


MITReactorMessenger::MITReactorMessenger(MITReactor* theSetup)
 :labSetup(theSetup)
{ 
  labDir = new G4UIdirectory("/CDMS/MITReactorGeom/");
  labDir->SetGuidance("UI commands to configure MITReactor geometry");

  OverBurdenCmd = new G4UIcmdWithADoubleAndUnit("/CDMS/MITReactorGeom/setOverBurden", this);
  OverBurdenCmd->SetGuidance("Set thickness of rock above cavern");
  OverBurdenCmd->SetParameterName("Size", false);
  OverBurdenCmd->SetRange("Size>=0.");
  OverBurdenCmd->SetUnitCategory("Length");
  OverBurdenCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}


MITReactorMessenger::~MITReactorMessenger()
{
  delete OverBurdenCmd;
  delete labDir;
}


void MITReactorMessenger::SetNewValue(G4UIcommand* command, G4String value)
{ 
  if (command == OverBurdenCmd) labSetup->SetOverBurden(OverBurdenCmd->GetNewDoubleValue(value));
}
