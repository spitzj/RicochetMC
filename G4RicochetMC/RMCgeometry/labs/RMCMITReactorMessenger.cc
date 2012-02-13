////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCMITReactorMessenger.cc                               //
//                                                                    //
//  Description: Messenger class to allow setting MITReactor geometry //
//               parameters                                           //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//  Date:        29 September 2011                                    //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "RMCMITReactorMessenger.hh"
#include "RMCMITReactor.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"


RMCMITReactorMessenger::RMCMITReactorMessenger(RMCMITReactor* theSetup)
 :labSetup(theSetup)
{ 

}


RMCMITReactorMessenger::~RMCMITReactorMessenger()
{

}


void RMCMITReactorMessenger::SetNewValue(G4UIcommand* command, G4String)
{}