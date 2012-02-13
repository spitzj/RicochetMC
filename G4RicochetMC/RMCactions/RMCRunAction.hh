////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCRunAction.hh                                      //
//  Description: run action class                                     //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               (adapted from Dennis Wright (SLAC))                  //
//  Date:        8 January 2012                                       //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#ifndef RMCRunAction_h
#define RMCRunAction_h 1

//#include "RMCactions/CDMSSimRoot.hh"
//#include "RMCactions/CDMSSimData.hh"
#include "G4UserRunAction.hh"
#include "globals.hh"
 
class G4Run;
class G4Event;
//class RMCEventSelector;
//class RMCRunActionMessenger;


class RMCRunAction : public G4UserRunAction {
public:
  RMCRunAction();
  virtual ~RMCRunAction();
  
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);
 
  //void transferEvent(const G4Event*); 

  void SetVerboseLevel(const G4int val);

  //void SetAutoSeed(const G4bool val) { AutoSeed = val; }
  //void SetOutputDataToFile(const G4bool val) { OutputDataToFile = val; }
  //void SetOutputTrees(const G4bool val) { OutputTrees = val; }
  //void SetDataFileNamePrefix(const G4String& val) { filePrefix = val; }

  // NOTE:  Transfers ownership here, so Messenger doesn't have to save
  //void SetEventSelector(RMCEventSelector* val=0);

protected:
  G4int GetLaboratoryCode() const;	// Determine lab code from geometry

private:
  G4int zipCollID;
  G4int vetoCollID;
  G4String filePrefix;
  //CDMSSimRoot rootOut;
  //CDMSSimData dataOut;

  G4int verboseLevel;		// Will be passed to other entities
  G4bool AutoSeed;
  G4bool OutputDataToFile;
  G4bool OutputTrees;

  //RMCEventSelector* EventSelector;

  //RMCRunActionMessenger* runActionMessenger;
};

#endif	/* RMCRunAction_h */

