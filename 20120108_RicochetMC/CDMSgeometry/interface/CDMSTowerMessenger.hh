#ifndef CDMSTowerMessenger_h
#define CDMSTowerMessenger_h 1
// $Id: CDMSTowerMessenger.hh,v 1.2 2010/12/07 19:23:49 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSTowerMessenger.hh                                //
//                                                                    //
//  Description: Messenger class to allow setting CDMS single tower   //
//               geometry parameters                                  //
//                                                                    //
//  Author:      Michael Kesey (SLAC)                                 //
//  Date:        24 November 2010                                     //
//                                                                    //
//  20101203  M. Kelsey -- Make number of ZIPs per tower free param.  //
////////////////////////////////////////////////////////////////////////

#include "G4UImessenger.hh"
#include "globals.hh"

class CDMSTowerConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;


class CDMSTowerMessenger: public G4UImessenger {
public:
  CDMSTowerMessenger(CDMSTowerConstruction* builder);
  virtual ~CDMSTowerMessenger();
  
  void SetNewValue(G4UIcommand* cmd, G4String value);

private:
  CDMSTowerConstruction* theBuilder;
  G4UIdirectory* cmdDir;
  G4UIcmdWithAnInteger*      NTowerSidesCmd;
  G4UIcmdWithAnInteger*      NTowerZipsCmd;
  G4UIcmdWithADoubleAndUnit* LidClearanceCmd;
  G4UIcmdWithADoubleAndUnit* ExtraSpaceCmd;
  G4UIcmdWithADoubleAndUnit* StripThickCmd;
  G4UIcmdWithADoubleAndUnit* StripWidthCmd;
};

#endif	/* CDMSTowerMessenger_hh */
