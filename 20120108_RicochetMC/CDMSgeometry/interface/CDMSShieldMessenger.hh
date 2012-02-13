#ifndef CDMSShieldMessenger_h
#define CDMSShieldMessenger_h 1
// $Id: CDMSShieldMessenger.hh,v 1.2 2011/01/13 17:51:49 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//  File:        CDMSShieldMessenger.hh                               //
//                                                                    //
//  Description: Messenger class to allow setting of CDMS shielding   //
//               geometry parameters                                  //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        22 December 2010                                     //
//                                                                    //
//  20110112  M. Kelsey -- Use new ShieldLayer command interface.     //
////////////////////////////////////////////////////////////////////////

#include "CDMSg4base/CDMSMessengerBase.hh"

class CDMSShieldConstruction;
class CDMS_UIcmdShieldLayer;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;


class CDMSShieldMessenger: public CDMSMessengerBase {
public:
  CDMSShieldMessenger(CDMSShieldConstruction* builder);
  virtual ~CDMSShieldMessenger();
  
  void SetNewValue(G4UIcommand* cmd, G4String value);

private:
  CDMSShieldConstruction* theBuilder;

  CDMS_UIcmdShieldLayer*     MuMetalCmd;
  CDMS_UIcmdShieldLayer*     InnerPolyCmd;
  CDMS_UIcmdShieldLayer*     InnerLeadCmd;
  CDMS_UIcmdShieldLayer*     OuterPolyCmd;
  CDMS_UIcmdShieldLayer*     OuterLeadCmd;
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
