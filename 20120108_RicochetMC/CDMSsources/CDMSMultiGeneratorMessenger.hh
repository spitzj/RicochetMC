#ifndef CDMSMultiGeneratorMessenger_hh
#define CDMSMultiGeneratorMessenger_hh 1
// $Id: CDMSMultiGeneratorMessenger.hh,v 1.1 2011/07/19 18:44:28 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//  File:        CDMSMultiGeneratorMessenger.hh                       //
//  Description: messenger class for multiple-particle generator      //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        6 July 2011                                          //
//                                                                    //
//  20110707  M. Kelsey -- Extract source additions from GammaSphere  //
////////////////////////////////////////////////////////////////////////

#include "CDMSg4base/CDMSMessengerBase.hh"
#include "globals.hh"

class CDMSMultiGenerator;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

class CDMSMultiGeneratorMessenger: public CDMSMessengerBase {
public:
  CDMSMultiGeneratorMessenger(CDMSMultiGenerator* theGenerator);
  ~CDMSMultiGeneratorMessenger();
    
  void SetNewValue(G4UIcommand* cmd, G4String value);

  void ActionAfterSetVerbose();		// To pass verbosity through

private:
  void AddSpectrum(const G4String& value);
  void AddIsotope(const G4String& value);

  CDMSMultiGenerator* generator;
  
  G4UIcmdWithAString* AddSpectrumCmd;
  G4UIcmdWithAString* AddIsotopeCmd;
  G4UIcmdWithoutParameter* PrintCmd;
};

#endif	/* CDMSMultiGeneratorMessenger_hh */
