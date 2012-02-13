// $Id: CDMSMultiTowerMessenger.cc,v 1.3 2011/05/25 05:10:27 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//  File:        CDMSMultiTowerMessenger.cc                           //
//                                                                    //
//  Description: Messenger class to construct full 100kg detector     //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        5 January 2011                                       //
//                                                                    //
//  20110524  M. Kelsey -- Use new CreateCommand interface.           //
////////////////////////////////////////////////////////////////////////

#include "CDMSgeometry/interface/CDMSMultiTowerMessenger.hh"
#include "CDMSgeometry/detectors/CDMSMultiTowerConstruction.hh"
#include "G4UIcmdWithABool.hh"


// NOTE:  Putting this command at the top level, parallel to "/CDMS/Detector"

CDMSMultiTowerMessenger::
CDMSMultiTowerMessenger(CDMSMultiTowerConstruction* builder)
  : CDMSMessengerBase("/CDMS/","UI commands to configure CDMS geometry"),
    theBuilder(builder) {
  UseShieldCmd = CreateCommand<G4UIcmdWithABool>("UseShield", 
	"Build veto shield around cryostat");
  UseShieldCmd->SetParameterName("useShield", true, false);
  UseShieldCmd->SetDefaultValue(true);
    }

CDMSMultiTowerMessenger::~CDMSMultiTowerMessenger() {
  delete UseShieldCmd;
}

void CDMSMultiTowerMessenger::ActionAfterSetVerbose() {
  if (verboseLevel)
    G4cout << " CDMSMultiTowerMessenger::ActionAfterSetVerbose" << G4endl;

  theBuilder->SetVerboseLevel(verboseLevel);
}


void CDMSMultiTowerMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  CDMSMessengerBase::SetNewValue(cmd,value);	// Check base class first!

  if (cmd == UseShieldCmd)
    theBuilder->UseShield(UseShieldCmd->GetNewBoolValue(value));
}



