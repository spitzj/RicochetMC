//
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        MITReactorMessenger.hh                               //
//                                                                    //
//  Description: Messenger class to allow setting MITReactor geometry //
//               parameters                                           //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//  Date:        28 September 2011                                    //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#ifndef MITReactorMessenger_h
#define MITReactorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class MITReactor;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;


class MITReactorMessenger: public G4UImessenger
{
  public:
    MITReactorMessenger(MITReactor*);
   ~MITReactorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    MITReactor* labSetup;
    G4UIdirectory* labDir;
    G4UIcmdWithADoubleAndUnit* OverBurdenCmd;
};

#endif

