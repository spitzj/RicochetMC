#ifndef RMCMessengerBase_hh
#define RMCMessengerBase_hh 1
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCMessengerBase.hh                                  //     
//  Description: Base class for all RMC GEANT4 Messengers. Provides   //
//		 common functionality and some diagnostic utilities.          //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Adapted from: Michael Kelsey (SLAC)                  //
//  Date:        14 January 2012                                      //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "G4UImessenger.hh"
#include "globals.hh"
#include <iosfwd>

class G4UIdirectory;
class G4UIcommandTree;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;


class RMCMessengerBase : public G4UImessenger {
public:
  RMCMessengerBase(const char* path, const char* desc);
  virtual ~RMCMessengerBase();

  // Interface command needed by G4UImanager -- subclasses should call back!
  virtual void SetNewValue(G4UIcommand* command, G4String newValue);

  // Optional subclass action to pass verbosity to other classes
  virtual void ActionAfterSetVerbose() {}	// Subclasses may implement

protected:
  // Create G4UIcommand (arbitrary subclass) within current command path
  template <class T>
  T* CreateCommand(const G4String& commandName, const G4String& description);

  G4int verboseLevel;		// Subclasses may access this directly

  // Set verbosity, and call optional user pass-along action
  void SetVerboseLevel(G4int verbose=0) {
    verboseLevel = verbose;
    ActionAfterSetVerbose();
  }

  G4int GetVerboseLevel() const { return verboseLevel; }

private:
  // Print names and current values of all parametric commands
  void PrintCurrentValues(std::ostream& os, G4bool printAllTrees=true) const;
  void PrintCurrentValues(std::ostream& os, G4UIcommandTree* tree,
			  G4bool printAllTrees=true) const;
  void PrintCurrentValue(std::ostream& os, const G4UIcommand* cmd) const;

  // Configuration functions used by constructor
  void CreateDirectory(const char* path, const char* desc);
  void CreateCommands(const char* path);

private:
  G4bool localCmdDir;		// Flag if directory was created or found
  G4UIdirectory* cmdDir;
  G4UIcommandTree* cmdTree;

  G4UIcmdWithAnInteger* verboseCmd;
  G4UIcmdWithAString* reportCmd;
};

#include "RMCg4base/RMCMessengerBase.icc"

#endif /* RMCMessengerBase_hh */
