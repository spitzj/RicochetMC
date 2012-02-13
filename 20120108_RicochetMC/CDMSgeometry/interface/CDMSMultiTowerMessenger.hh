#ifndef CDMSGeomMessenger_hh
#define CDMSGeomMessenger_hh 1
// $Id: CDMSMultiTowerMessenger.hh,v 1.1 2011/01/05 06:56:45 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//  File:        CDMSMultiTowerMessenger.hh                           //
//                                                                    //
//  Description: Messenger class to construct full 100kg detector     //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        5 January 2011                                       //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "CDMSg4base/CDMSMessengerBase.hh"
#include "globals.hh"

class G4UIcmdWithABool;
class CDMSMultiTowerConstruction;


class CDMSMultiTowerMessenger : public CDMSMessengerBase {
public:
  CDMSMultiTowerMessenger(CDMSMultiTowerConstruction* builder);
  virtual ~CDMSMultiTowerMessenger();
  
  void SetNewValue(G4UIcommand* cmd, G4String value);
  
  void ActionAfterSetVerbose();		// Used by base class for pass-through
  
private:
  CDMSMultiTowerConstruction* theBuilder;
  G4UIcmdWithABool* UseShieldCmd;
};

#endif	/* CDMSMultiTowerMessenger_hh */
