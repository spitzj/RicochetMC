#ifndef RMCTowerMessenger_h
#define RMCTowerMessenger_h 1
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCTowerMessenger.hh                                //
//                                                                    //
//  Description: Messenger class to allow setting RMC single tower   //
//               geometry parameters                                  //
//                                                                    //
//  Author:      Michael Kesey (SLAC)                                 //
//  Date:        24 November 2010                                     //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "G4UImessenger.hh"
#include "globals.hh"

class RMCTowerConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;


class RMCTowerMessenger: public G4UImessenger {
public:
  RMCTowerMessenger(RMCTowerConstruction* builder);
  virtual ~RMCTowerMessenger();
  
  void SetNewValue(G4UIcommand* cmd, G4String value);

private:
  RMCTowerConstruction* theBuilder;
  G4UIdirectory* cmdDir;
  G4UIcmdWithAnInteger*      NTowerSidesCmd;
  G4UIcmdWithAnInteger*      NTowerZipsCmd;
  G4UIcmdWithADoubleAndUnit* LidClearanceCmd;
  G4UIcmdWithADoubleAndUnit* ExtraSpaceCmd;
  G4UIcmdWithADoubleAndUnit* StripThickCmd;
  G4UIcmdWithADoubleAndUnit* StripWidthCmd;
};

#endif	/* RMCTowerMessenger_hh */
