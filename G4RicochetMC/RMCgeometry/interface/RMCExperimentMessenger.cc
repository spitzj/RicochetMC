////////////////////////////////////////////////////////////////////////
//  File:        RMCExperimentMessenger.cc                            //
//                                                                    //
//  Description: Messenger class to construct detector with modular   //
//               components.                                          //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Adapted from: Michael Kelsey (SLAC)                  //
//  Date:        12 February 2012                                     //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "RMCgeometry/interface/RMCExperimentMessenger.hh"
#include "RMCgeometry/detectors/RMCExperimentConstruction.hh"
#include "G4UIcmdWithABool.hh"

RMCExperimentMessenger::RMCExperimentMessenger(RMCExperimentConstruction* builder)
: RMCMessengerBase("/RMC/","UI commands to configure RMC geometry"),
theBuilder(builder) {
    UseShieldCmd = CreateCommand<G4UIcmdWithABool>("UseShield", 
                                                   "Build veto shield around cryostat");
    UseShieldCmd->SetParameterName("UseShield", true, false);
    UseShieldCmd->SetDefaultValue(true);
    
    UseCryostatCmd = CreateCommand<G4UIcmdWithABool>("UseCryostat", 
                                                   "Build cryostat");   //
    UseCryostatCmd->SetParameterName("UseCryostat", true, false);
    UseCryostatCmd->SetDefaultValue(true);
    
    UseTowerCmd = CreateCommand<G4UIcmdWithABool>("UseTower", 
                                                     "Build detector tower");
    UseTowerCmd->SetParameterName("UseTower", true, false);
    UseTowerCmd->SetDefaultValue(true);
}

RMCExperimentMessenger::~RMCExperimentMessenger() {
    delete UseShieldCmd;
    delete UseCryostatCmd;
    delete UseTowerCmd;
}

void RMCExperimentMessenger::ActionAfterSetVerbose() {
    if (verboseLevel)
        G4cout << " RMCExperimentMessenger::ActionAfterSetVerbose" << G4endl;
    
    theBuilder->SetVerboseLevel(verboseLevel);
}


void RMCExperimentMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
    RMCMessengerBase::SetNewValue(cmd,value);	// Check base class first!
    
    if (cmd == UseShieldCmd)
        theBuilder->UseShield(UseShieldCmd->GetNewBoolValue(value));
    if (cmd == UseCryostatCmd)
        theBuilder->UseCryostat(UseCryostatCmd->GetNewBoolValue(value));
    if (cmd == UseTowerCmd)
        theBuilder->UseTower(UseTowerCmd->GetNewBoolValue(value));
}



