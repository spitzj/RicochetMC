#ifndef CDMSGeomMessenger_hh
#define CDMSGeomMessenger_hh 1
// $Id: CDMSGeomMessenger.hh,v 1.6 2011/04/27 17:04:11 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//  File:        CDMSGeomMessenger.hh                                 //
//                                                                    //
//  Description: Messenger class to construct complete CDMS models    //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        26 October 2010                                      //
//                                                                    //
//  20101203  M. Kelsey -- Discard "icebox" as separate activity.     //
//  20101226  M. Kelsey -- Use new CDMSMessengerBase                  //
//  20110105  M. Kelsey -- Drop "shield" -- move to MultiTower det.   //
//  20110426  M. Kelsey -- Make GeometryManager a data member.        //
////////////////////////////////////////////////////////////////////////

#include "CDMSg4base/CDMSMessengerBase.hh"
#include "globals.hh"

class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class CDMSGeomConstructor;
class CDMSGeometryManager;

class CDMSGeomMessenger: public CDMSMessengerBase {
public:
  CDMSGeomMessenger(CDMSGeomConstructor* theSetup);
  virtual ~CDMSGeomMessenger();
  
  void SetNewValue(G4UIcommand* cmd, G4String value);

  void ActionAfterSetVerbose();		// Used by base class for pass-through

protected:
  void assignLab(const G4String& name);
  void assignDetector(const G4String& name);
  void assignSource(const G4String& name);

private:
  CDMSGeomConstructor* setup;
  CDMSGeometryManager* geomManager;

  G4UIcmdWithoutParameter* UpdateGeomCmd;
  G4UIcmdWithAString* SelectLabCmd;
  G4UIcmdWithAString* SelectDetCmd;
  G4UIcmdWithAString* SelectSrcCmd;
};

#endif	/* CDMSGeomMessenger_hh */
