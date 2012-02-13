//
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        NoLabMessenger.cc                                    //
//                                                                    //
//  Description: Messenger class to allow setting NoLab geometry      //
//               parameters                                           //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        23 June 2010                                         //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "NoLabMessenger.hh"

#include "NoLab.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"


NoLabMessenger::NoLabMessenger(NoLab* theSetup)
 :labSetup(theSetup)
{ 
  labDir = new G4UIdirectory("/CDMS/NoLabGeom/");
  labDir->SetGuidance("UI commands to configure NoLab geometry");
}


NoLabMessenger::~NoLabMessenger()
{
  delete labDir;
}


void NoLabMessenger::SetNewValue(G4UIcommand* command, G4String)
{}
