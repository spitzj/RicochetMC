#ifndef RMCShieldMessenger_h
#define RMCShieldMessenger_h 1
////////////////////////////////////////////////////////////////////////
//  File:        RMCShieldMessenger.hh                               //
//                                                                    //
//  Description: Messenger class to allow setting of RMC shielding   //
//               geometry parameters                                  //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        22 December 2010                                     //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "RMCg4base/RMCMessengerBase.hh"

class RMCShieldConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;


class RMCShieldMessenger: public RMCMessengerBase {
public:
  RMCShieldMessenger(RMCShieldConstruction* builder);
  virtual ~RMCShieldMessenger();
  
  void SetNewValue(G4UIcommand* cmd, G4String value);

private:
  RMCShieldConstruction* theBuilder;

  G4UIcmdWithADoubleAndUnit* ScintThickCmd;
  G4UIcmdWithAnInteger*      TopPanelsXCmd;
  G4UIcmdWithAnInteger*      TopPanelsYCmd;
  G4UIcmdWithADoubleAndUnit* OverhangCmd;
  G4UIcmdWithADoubleAndUnit* OverlapXCmd;
  G4UIcmdWithADoubleAndUnit* OverlapYCmd;
  G4UIcmdWithADoubleAndUnit* TopVetoSHCmd;
  G4UIcmdWithADoubleAndUnit* BottomVetoSHCmd;
  G4UIcmdWithAnInteger*      SidePanelsCmd;
  G4UIcmdWithADoubleAndUnit* SideClearanceCmd;
  G4UIcmdWithADoubleAndUnit* SideOverlapCmd;
  G4UIcmdWithADoubleAndUnit* SideCornerCmd;
};

#endif
