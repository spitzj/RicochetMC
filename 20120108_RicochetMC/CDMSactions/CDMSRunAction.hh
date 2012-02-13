// $Id: CDMSRunAction.hh,v 1.9 2011/07/06 21:24:46 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSRunAction.hh                                     //
//  Description: run action class for CDMSMini                        //
//                                                                    //
//  Author:      Dennis Wright                                        //
//  Date:        13 May 2010                                          //
//                                                                    //
//  20110628  M. Kelsey -- Add hook for user-defined event rejection. //
//  20110630  M. Kelsey -- Move SetVerboseLevel implementation to .cc //
//		Add lab-code query function for file creation.        //
//  20110706  M. Kelsey -- Move SetEventSelector to .cc to do delete  //
////////////////////////////////////////////////////////////////////////

#ifndef CDMSRunAction_h
#define CDMSRunAction_h 1

#include "CDMSactions/CDMSSimRoot.hh"
#include "CDMSactions/CDMSSimData.hh"
#include "G4UserRunAction.hh"
#include "globals.hh"
 
class G4Run;
class G4Event;
class CDMSEventSelector;
class CDMSRunActionMessenger;


class CDMSRunAction : public G4UserRunAction {
public:
  CDMSRunAction();
  virtual ~CDMSRunAction();
  
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);
 
  void transferEvent(const G4Event*); 

  void SetVerboseLevel(const G4int val);

  void SetAutoSeed(const G4bool val) { AutoSeed = val; }
  void SetOutputDataToFile(const G4bool val) { OutputDataToFile = val; }
  void SetOutputTrees(const G4bool val) { OutputTrees = val; }
  void SetDataFileNamePrefix(const G4String& val) { filePrefix = val; }

  // NOTE:  Transfers ownership here, so Messenger doesn't have to save
  void SetEventSelector(CDMSEventSelector* val=0);

protected:
  G4int GetLaboratoryCode() const;	// Determine lab code from geometry

private:
  G4int zipCollID;
  G4int vetoCollID;
  G4String filePrefix;
  CDMSSimRoot rootOut;
  CDMSSimData dataOut;

  G4int verboseLevel;		// Will be passed to other entities
  G4bool AutoSeed;
  G4bool OutputDataToFile;
  G4bool OutputTrees;

  CDMSEventSelector* EventSelector;

  CDMSRunActionMessenger* runActionMessenger;
};

#endif	/* CDMSRunAction_h */

