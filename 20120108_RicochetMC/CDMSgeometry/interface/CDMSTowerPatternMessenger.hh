#ifndef CDMSTowerPatternMessenger_h
#define CDMSTowerPatternMessenger_h 1
// $Id: CDMSTowerPatternMessenger.hh,v 1.3 2011/02/05 05:54:26 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//  File:        CDMSTowerPatternMessenger.hh                         //
//                                                                    //
//  Description: Messenger class to allow setting CDMS multi-tower    //
//               layout parameters                                    //
//                                                                    //
//  Author:      Michael Kesey (SLAC)                                 //
//  Date:        24 November 2010                                     //
//                                                                    //
//  20110121  M. Kelsey -- Use string input for layouts               //
//  20110204  M. Kelsey -- Add new command for Optimize(mass)         //
////////////////////////////////////////////////////////////////////////

#include "G4UImessenger.hh"
#include "globals.hh"

class CDMSTowerPattern;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;


class CDMSTowerPatternMessenger: public G4UImessenger {
public:
  CDMSTowerPatternMessenger(CDMSTowerPattern* pattern);
  virtual ~CDMSTowerPatternMessenger();
  
  void SetNewValue(G4UIcommand* cmd, G4String value);

private:
  CDMSTowerPattern* towerPattern;
  G4UIdirectory* cmdDir;
  G4UIcmdWithAString*   LayoutCodeCmd;
  G4UIcmdWithAnInteger* NumberOfTowersCmd;
  G4UIcmdWithADoubleAndUnit* OptimizeCmd;
};

#endif	/* CDMSTowerPatternMessenger_hh */
