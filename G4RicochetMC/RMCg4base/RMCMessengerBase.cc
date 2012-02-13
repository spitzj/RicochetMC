////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCMessengerBase.cc                                  //     
//  Description: Base class for all RMC GEANT4 Messengers. Provides   //
//		 common functionality and some diagnostic utilities.          //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Adapted from: Michael Kelsey (SLAC)                  //
//  Date:        14 January 2012                                      //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "RMCg4base/RMCMessengerBase.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcommand.hh"
#include "G4UIcommandTree.hh"
#include "G4UIdirectory.hh"
#include "G4UImanager.hh"


// Constructor looks to see if directory exists before creating it

RMCMessengerBase::RMCMessengerBase(const char* path, const char* desc) :
  verboseLevel(0), localCmdDir(false), cmdDir(0), cmdTree(0),
  verboseCmd(0), reportCmd(0) {
  CreateDirectory(path, desc);
  CreateCommands(path);
}

// Destructor only deletes directory if it was created here

RMCMessengerBase::~RMCMessengerBase() {
  delete verboseCmd;
  delete reportCmd;
  if (localCmdDir) delete cmdDir;
}


// Create new command directory, or reference existing version

void RMCMessengerBase::CreateDirectory(const char* path, const char* desc) {
  G4UImanager* UIman = G4UImanager::GetUIpointer();
  if (!UIman) return;

  // See if input path has already been registered
  G4UIcommand* foundPath = UIman->GetTree()->FindPath(path);
  if (foundPath) cmdDir = dynamic_cast<G4UIdirectory*>(foundPath);

  if (!cmdDir) {		// Create local deletable directory
    localCmdDir = true;
    cmdDir = new G4UIdirectory(path);
    cmdDir->SetGuidance(desc);
  }

  cmdTree = UIman->GetTree()->FindCommandTree(path);	// For printing later
}


// Define local commands, available to all subclasses

void RMCMessengerBase::CreateCommands(const char* path) {
  if (!cmdDir) return;		// No valid directory to contain commands

  G4String basePath(path);

  verboseCmd =
    CreateCommand<G4UIcmdWithAnInteger>("verbose",
					"Report information during processing");
  verboseCmd->SetParameterName("verbose",true,true);
  verboseCmd->SetGuidance(" 0 : Silent (default)");
  verboseCmd->SetGuidance(" 1 : Display main actions");
  verboseCmd->SetGuidance(" 2 : Display detailed activities");
  verboseCmd->SetDefaultValue(0);

  reportCmd =
    CreateCommand<G4UIcmdWithAString>("printValues",
				      "Dump commands with current values");
  reportCmd->SetParameterName("recursive",true,true);
  reportCmd->SetGuidance("all : Include nested command trees (optional)");
  reportCmd->SetCandidates("all self");
  reportCmd->SetDefaultValue("self");
}


// Process base-class commands (subclasses should call back here!)

void RMCMessengerBase::SetNewValue(G4UIcommand* command, G4String val) {
  if (command == verboseCmd) {
    verboseLevel = StoI(val);
    ActionAfterSetVerbose();		// For subclass propagation
  }

  if (command == reportCmd) {
    // FIXME:  Direct indexing gives "ambiguous overload" error w/GCC 4.0.1
    G4bool allTrees = ('a'==val.c_str()[0] || 'A'==val.c_str()[0]);
    PrintCurrentValues(G4cout, allTrees);
  }
}

// Print names and current values of all parametric commands

void RMCMessengerBase::PrintCurrentValues(std::ostream& os,
					   G4bool printAllTrees) const {
  PrintCurrentValues(os, cmdTree, printAllTrees);
}


// DANGER:  G4UIcommandTree indexes its lists from 1, not from 0!
void RMCMessengerBase::PrintCurrentValues(std::ostream& os,
					   G4UIcommandTree* tree,
					   G4bool printAllTrees) const {
  // Report commands in current tree first
  G4int Ncommands = tree->GetCommandEntry();

  // DANGER:  G4UIcommandTree indexes its lists from 1, not from 0!
  for (G4int i=0; i<Ncommands; i++) {
    G4UIcommand* cmdI = tree->GetCommand(i+1);
    if (!cmdI || cmdI->GetParameterEntries() == 0) continue;

    // NOTE:  Cannot use misnamed "GetCurrentValue()" to print current value
    PrintCurrentValue(os, cmdI);
  }

  // If user requested recursion, do the same for nested trees
  if (printAllTrees) {
    G4int Ntrees = tree->GetTreeEntry();

    // DANGER:  G4UIcommandTree indexes its lists from 1, not from 0!
    for (G4int i=0; i<Ntrees; i++) {
      PrintCurrentValues(os, tree->GetTree(i+1), printAllTrees);
    }
  }
}

// NOTE:  Cannot use "GetCurrentValue()" to print most recently entered value
void RMCMessengerBase::PrintCurrentValue(std::ostream& os,
					  const G4UIcommand* cmd) const {
  G4int Nparam = cmd->GetParameterEntries();
  if (Nparam < 1) return;

  // DANGER:  Cannot actually get "current value" without "useCurrentAsDefault"
  os << cmd->GetCommandPath();
  for (G4int i=0; i<Nparam; i++)
    os << " " << cmd->GetParameter(i)->GetDefaultValue();
    
  os << G4endl;
}
