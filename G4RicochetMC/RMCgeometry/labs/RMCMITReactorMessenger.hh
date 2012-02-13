////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCMITReactorMessenger.hh                               //
//                                                                    //
//  Description: Messenger class to allow setting MITReactor geometry //
//               parameters                                           //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//  Date:        28 September 2011                                    //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#ifndef RMCMITReactorMessenger_h
#define RMCMITReactorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class RMCMITReactor;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;


class RMCMITReactorMessenger: public G4UImessenger
{
public:
    RMCMITReactorMessenger(RMCMITReactor*);
    ~RMCMITReactorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
private:
    RMCMITReactor* labSetup;
    G4UIdirectory* labDir;
    G4UIcmdWithADoubleAndUnit* OverBurdenCmd;
};

#endif

