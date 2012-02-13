// $Id: CDMSMiniGeomMessenger.cc,v 1.3 2010/10/26 00:03:43 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSMiniGeomMessenger.cc                             //
//                                                                    //
//  Description: Messenger class to allow setting CDMS mini           //
//               geometry parameters                                  //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        23 June 2010                                         //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "CDMSgeometry/detectors/CDMSMiniGeomMessenger.hh"

#include "CDMSgeometry/detectors/CDMSMiniDetConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"


CDMSMiniGeomMessenger::
CDMSMiniGeomMessenger(CDMSMiniDetConstruction* theSetup)
 :miniSetup(theSetup)
{ 
  miniGeomDir = new G4UIdirectory("/CDMS/miniGeom/");
  miniGeomDir->SetGuidance("UI commands to configure CDMS mini geometry");
  
  ZipRadCmd = new G4UIcmdWithADoubleAndUnit("/CDMS/miniGeom/setZipRad", this);
  ZipRadCmd->SetGuidance("Set radius of zip crystal");
  ZipRadCmd->SetParameterName("Size",false);
  ZipRadCmd->SetRange("Size>=0.");
  ZipRadCmd->SetUnitCategory("Length");  
  ZipRadCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
  ZipLenCmd = new G4UIcmdWithADoubleAndUnit("/CDMS/miniGeom/setZipLen",this);
  ZipLenCmd->SetGuidance("Set length of zip crystal");
  ZipLenCmd->SetParameterName("Size",false);
  ZipLenCmd->SetRange("Size>=0.");
  ZipLenCmd->SetUnitCategory("Length");
  ZipLenCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  DetBoxShimCmd = new G4UIcmdWithADoubleAndUnit("/CDMS/miniGeom/setDetBoxShimThick",this);
  DetBoxShimCmd->SetGuidance("Set detector box shim thickness");
  DetBoxShimCmd->SetParameterName("Size",false);
  DetBoxShimCmd->SetRange("Size>=0.");
  DetBoxShimCmd->SetUnitCategory("Length");  
  DetBoxShimCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
}


CDMSMiniGeomMessenger::~CDMSMiniGeomMessenger()
{
  delete ZipRadCmd;
  delete ZipLenCmd;
  delete DetBoxShimCmd;
  delete miniGeomDir;
}


void CDMSMiniGeomMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if (command == ZipRadCmd) miniSetup->SetZipRad(ZipRadCmd->GetNewDoubleValue(newValue) );

  if (command == ZipLenCmd) miniSetup->SetZipLen(ZipLenCmd->GetNewDoubleValue(newValue) );

  if (command == DetBoxShimCmd) miniSetup->SetDetBoxShimThick(DetBoxShimCmd->GetNewDoubleValue(newValue) );
}
