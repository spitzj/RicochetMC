#ifndef RMCGeomMessenger_hh
#define RMCGeomMessenger_hh 1
////////////////////////////////////////////////////////////////////////
//  File:        RMCExperimentMessenger.hh                            //
//                                                                    //
//  Description: Messenger class to construct detector with modular   //
//               components.                                          //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Adapted from: Michael Kelsey (SLAC)                  //
//  Date:        13 January 2012                                       //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "RMCg4base/RMCMessengerBase.hh"
#include "globals.hh"

class G4UIcmdWithABool;
class RMCExperimentConstruction;


class RMCExperimentMessenger : public RMCMessengerBase {
public:
    RMCExperimentMessenger(RMCExperimentConstruction* builder);
    virtual ~RMCExperimentMessenger();
    
    void SetNewValue(G4UIcommand* cmd, G4String value);
    
    void ActionAfterSetVerbose();		// Used by base class for pass-through
    
private:
    RMCExperimentConstruction* theBuilder;
    G4UIcmdWithABool* UseShieldCmd;
    G4UIcmdWithABool* UseCryostatCmd;
    G4UIcmdWithABool* UseTowerCmd;
};

#endif	/* RMCExperimentMessenger_hh */
