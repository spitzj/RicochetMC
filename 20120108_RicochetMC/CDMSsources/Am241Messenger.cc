// $Id: Am241Messenger.cc,v 1.4 2010/12/14 07:41:00 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        Am241Messenger.cc                                    //
//                                                                    //
//  Description: messenger class for Am241 source                     //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        2 September 2010                                     //
//                                                                    //
//  20101210  M. Kelsey -- Follow interface changes to action class.  //
////////////////////////////////////////////////////////////////////////

#include "CDMSsources/Am241Messenger.hh"

#include "CDMSsources/Am241Source.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"


Am241Messenger::Am241Messenger(Am241Source* theSource)
 :source(theSource)
{ 
  Am241Dir = new G4UIdirectory("/CDMS/Am241/");
  Am241Dir->SetGuidance("UI commands to configure Am241 source geometry");
  
  SrcRadCmd = new G4UIcmdWithADoubleAndUnit("/CDMS/Am241/Radius", this);
  SrcRadCmd->SetGuidance("Set active radius of source");
  SrcRadCmd->SetParameterName("Size",false);
  SrcRadCmd->SetRange("Size>=0.");
  SrcRadCmd->SetUnitCategory("Length");  
  SrcRadCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
  SrcLenCmd = new G4UIcmdWithADoubleAndUnit("/CDMS/Am241/Length",this);
  SrcLenCmd->SetGuidance("Set length of source");
  SrcLenCmd->SetParameterName("Size",false);
  SrcLenCmd->SetRange("Size>=0.");
  SrcLenCmd->SetUnitCategory("Length");
  SrcLenCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SrcPosCmd = new G4UIcmdWith3VectorAndUnit("/CDMS/Am241/Position",this);
  SrcPosCmd->SetGuidance("Set source position");
  SrcPosCmd->SetParameterName("X", "Y", "Z", false);
  SrcPosCmd->SetUnitCategory("Length");  
  SrcPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SrcDirCmd = new G4UIcmdWith3Vector("/CDMS/Am241/Direction",this);
  SrcDirCmd->SetGuidance("Set source direction");
  SrcDirCmd->SetParameterName("VX", "VY", "VZ", false);
  SrcDirCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}


Am241Messenger::~Am241Messenger()
{
  delete SrcRadCmd;
  delete SrcLenCmd;
  delete SrcPosCmd;
  delete SrcDirCmd;
  delete Am241Dir;
}


void Am241Messenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if (command == SrcRadCmd) source->SetRadius(SrcRadCmd->GetNewDoubleValue(newValue) );
  if (command == SrcLenCmd) source->SetLength(SrcLenCmd->GetNewDoubleValue(newValue) );

  if (command == SrcPosCmd) source->SetPosition(SrcPosCmd->GetNew3VectorValue(newValue) );
  if (command == SrcDirCmd) source->SetDirection(SrcDirCmd->GetNew3VectorValue(newValue) );
}
