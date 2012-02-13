// $Id: CDMSRunActionMessenger.hh,v 1.5 2011/07/21 21:19:00 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
// File:        CDMSRunActionMessenger.hh                             //
// Description: run action messenger for CDMS mini detector           //
//                                                                    //
// Author:      adapted from existing CDMS mini code by D.H. Wright   //
//              (SLAC)                                                //
// Date:        25 June 2010                                          //
//                                                                    //
// 20101026  M. Kelsey -- Add verbosity                               //
// 20110630  M. Kelsey -- Drop unused commands                        //
// 20110706  M. Kelsey -- Add command to select event filter.         //
// 20110721  M. Kelsey -- Add command to initialize run number.       //
////////////////////////////////////////////////////////////////////////

#ifndef CDMSRunActionMessenger_h
#define CDMSRunActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class CDMSEventSelector;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class CDMSRunAction;


class CDMSRunActionMessenger: public G4UImessenger {
public:
  CDMSRunActionMessenger(CDMSRunAction* action);
  ~CDMSRunActionMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);

protected:
  CDMSEventSelector* MakeEventSelector(const G4String& name) const;

private:
  CDMSRunAction* runAction;
  G4UIdirectory* CDMSRunDir;
  
  G4UIcmdWithAString* setRunFileName;
  G4UIcmdWithAnInteger* RunNumberCmd;
  G4UIcmdWithABool* OutputDataToFileCmd;
  G4UIcmdWithABool* OutputTreesCmd;
  G4UIcmdWithABool* setAutoSeedCmd;
  G4UIcmdWithAString* EventFilterCmd;
  G4UIcmdWithAnInteger* setVerboseCmd;
};

#endif	/* CDMSRunActionMessenger_h */
