//
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        SnoLabMessenger.hh                                   //
//                                                                    //
//  Description: Messenger class to allow setting SnoLab geometry     //
//               parameters                                           //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        3 May 2011                                           //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#ifndef SnoLabMessenger_h
#define SnoLabMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class SnoLab;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;


class SnoLabMessenger: public G4UImessenger
{
  public:
    SnoLabMessenger(SnoLab*);
   ~SnoLabMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    SnoLab* labSetup;
    G4UIdirectory* labDir;
    G4UIcmdWithADoubleAndUnit* OverBurdenCmd;
    G4UIcmdWithADoubleAndUnit* MinWallThicknessCmd;
    G4UIcmdWithADoubleAndUnit* CavernLengthCmd;
};

#endif

