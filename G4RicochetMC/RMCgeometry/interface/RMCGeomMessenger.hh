#ifndef RMCGeomMessenger_hh
#define RMCGeomMessenger_hh 1
////////////////////////////////////////////////////////////////////////
//  File:        RMCGeomMessenger.hh                                  //
//                                                                    //
//  Description: Messenger class to construct complete RMC models     //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Adapted from: Michael Kelsey (SLAC)                  //
//  Date:        14 January 2012                                      //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "RMCg4base/RMCMessengerBase.hh"
#include "globals.hh"

class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class RMCGeomConstructor;
class RMCGeometryManager;

class RMCGeomMessenger: public RMCMessengerBase {
public:
  RMCGeomMessenger(RMCGeomConstructor* theSetup);
  virtual ~RMCGeomMessenger();
  
  void SetNewValue(G4UIcommand* cmd, G4String value);

  void ActionAfterSetVerbose();		// Used by base class for pass-through

protected:
  void assignLab(const G4String& name);
  void assignDetector(const G4String& name);
  void assignSource(const G4String& name);

private:
  RMCGeomConstructor* setup;
  RMCGeometryManager* geomManager;

  G4UIcmdWithoutParameter* UpdateGeomCmd;
  G4UIcmdWithAString* SelectLabCmd;
  G4UIcmdWithAString* SelectDetCmd;
  G4UIcmdWithAString* SelectSrcCmd;
};

#endif	/* RMCGeomMessenger_hh */
