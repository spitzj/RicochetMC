#ifndef CDMSGammaSphereMessenger_hh
#define CDMSGammaSphereMessenger_hh 1
// $Id: CDMSGammaSphereMessenger.hh,v 1.3 2011/07/07 05:04:47 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//  File:        CDMSGammaSphereMessenger.hh                          //
//  Description: messenger class for spherical gamma source           //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        24 May 2011                                          //
//                                                                    //
//  20110706  M. Kelsey -- Drop "cmdWithList" to allow optional dbl.  //
//  20110707  M. Kelsey -- Move source addition to MultiGenMessenger  //
////////////////////////////////////////////////////////////////////////

#include "CDMSg4base/CDMSMessengerBase.hh"
#include "globals.hh"

class CDMSGammaSphere;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;


class CDMSGammaSphereMessenger: public CDMSMessengerBase {
public:
  CDMSGammaSphereMessenger(CDMSGammaSphere* theSource);
  ~CDMSGammaSphereMessenger();
    
  void SetNewValue(G4UIcommand* cmd, G4String value);

  void ActionAfterSetVerbose();		// To pass verbosity through

private:
  CDMSGammaSphere* source;
  
  G4UIcmdWithADoubleAndUnit* RadiusCmd;
  G4UIcmdWith3VectorAndUnit* CenterCmd;
  G4UIcmdWithABool* InwardCmd;
  G4UIcmdWithADoubleAndUnit* HalfAngleCmd;
  G4UIcmdWithAnInteger* NParticlesCmd;
};

#endif	/* CDMSGammaSphereMessenger_hh */
