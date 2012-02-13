////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        NoLabMessenger.cc                                    //
//                                                                    //
//  Description: Messenger class to allow setting NoLab geometry      //
//               parameters                                           //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Dennis Wright (SLAC)                                 //
//  Date:        14 January 2012                                      //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "NoLabMessenger.hh"
#include "NoLab.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"


NoLabMessenger::NoLabMessenger(NoLab* theSetup)
 :labSetup(theSetup)
{ 
  labDir = new G4UIdirectory("/RMC/NoLabGeom/");
  labDir->SetGuidance("UI commands to configure NoLab geometry");
}


NoLabMessenger::~NoLabMessenger()
{
  delete labDir;
}


void NoLabMessenger::SetNewValue(G4UIcommand* command, G4String)
{}
