// $Id: CDMSRunActionMessenger.cc,v 1.6 2011/07/21 21:19:00 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
// File:        CDMSRunActionMessenger.hh                             //
// Description: run action for CDMS mini detector                     //
//                                                                    //
// Author:      adapted from existing CDMS mini code by D.H. Wright   //
//              (SLAC)                                                //
// Date:        25 June 2010                                          //
//                                                                    //
// 20110630  M. Kelsey -- Drop unused commands                        //
// 20110706  M. Kelsey -- Add command to select event filter.         //
// 20110721  M. Kelsey -- Add command to initialize run number.       //
////////////////////////////////////////////////////////////////////////

#include "CDMSactions/CDMSRunActionMessenger.hh"
#include "CDMSactions/CDMSRunAction.hh"
#include "CDMSactions/CDMSEventSelector.hh"
#include "CDMSactions/CDMS_SelectMISS.hh"

#include "G4RunManager.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"


CDMSRunActionMessenger::CDMSRunActionMessenger(CDMSRunAction* action)
 : runAction(action) { 
  // Create the run directory
  CDMSRunDir = new G4UIdirectory("/CDMS/");
  CDMSRunDir->SetGuidance("CDMS run action controls.");

  // run directory already exists
  // Set autoSeed command
  setAutoSeedCmd = new G4UIcmdWithABool("/run/autoSeed",this);
  setAutoSeedCmd->SetGuidance("Switch on/off time-based random seeds");
  setAutoSeedCmd->SetGuidance(" true: run seeds determined by system time");
  setAutoSeedCmd->SetGuidance("false: use command 'random/resetEngineFrom'");
  setAutoSeedCmd->SetGuidance("Default = false");
  setAutoSeedCmd->SetParameterName("autoSeed", false);
  setAutoSeedCmd->AvailableForStates(G4State_Idle);

  RunNumberCmd = new G4UIcmdWithAnInteger("/run/setRunID",this);
  RunNumberCmd->SetGuidance("Initialize run ID counter before beamOn");
  RunNumberCmd->SetParameterName("runID", false);
  RunNumberCmd->AvailableForStates(G4State_Idle);


  // Set OutputDataToFile
  OutputDataToFileCmd = new G4UIcmdWithABool("/CDMS/writeData",this);
  OutputDataToFileCmd->SetGuidance("Output results to ASCII file.");
  OutputDataToFileCmd->SetParameterName("dataout",true,true); 
  OutputDataToFileCmd->SetDefaultValue(true);
  OutputDataToFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  // Set OutputDataToFile
  OutputTreesCmd = new G4UIcmdWithABool("/CDMS/writeTrees",this);
  OutputTreesCmd->SetGuidance("Output results to root trees");
  OutputTreesCmd->SetParameterName("treesout",true,true);
  OutputTreesCmd->SetDefaultValue(true);
  OutputTreesCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  // Set file name
  setRunFileName = new G4UIcmdWithAString("/CDMS/writeFilePrefix",this);
  setRunFileName->SetGuidance("Set the prefix-name of the output files.");
  setRunFileName->SetParameterName("prefix",true);
  setRunFileName->SetDefaultValue("G4MC");
  setRunFileName->AvailableForStates(G4State_PreInit,G4State_Idle);

  // Assign event filter to restrict output
  EventFilterCmd = new G4UIcmdWithAString("/CDMS/eventFilter",this);
  EventFilterCmd->SetGuidance("Use specified selector to filter event output.");
  EventFilterCmd->SetGuidance("MISS :  Events with multiple gamma scatters");
  EventFilterCmd->SetParameterName("name",false);
  EventFilterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  EventFilterCmd->SetCandidates("MISS miss");
}


CDMSRunActionMessenger::~CDMSRunActionMessenger() {
  delete setAutoSeedCmd;
  delete RunNumberCmd;
  delete OutputDataToFileCmd;
  delete OutputTreesCmd;
  delete setRunFileName;
  delete EventFilterCmd;
}


void 
CDMSRunActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue) { 
  if (command == OutputDataToFileCmd)
    runAction->SetOutputDataToFile(OutputDataToFileCmd->GetNewBoolValue(newValue) );

  if (command == OutputTreesCmd)
    runAction->SetOutputTrees(OutputTreesCmd->GetNewBoolValue(newValue) );

  if (command == setAutoSeedCmd)
    runAction->SetAutoSeed(setAutoSeedCmd->GetNewBoolValue(newValue) );

  if (command == setRunFileName) runAction->SetDataFileNamePrefix(newValue);

  if (command == EventFilterCmd)
    runAction->SetEventSelector(MakeEventSelector(newValue));

  if (command == RunNumberCmd)
    G4RunManager::GetRunManager()->SetRunIDCounter(RunNumberCmd->GetNewIntValue(newValue));
}


// Create event-filter corresponding to name string

CDMSEventSelector* 
CDMSRunActionMessenger::MakeEventSelector(const G4String& name) const {
  CDMSEventSelector* filter = 0;

  // NOTE:  G4String::compareTo() like strcmp(): -1, 0, +1 lexical order

  //***********************************************************************
  // KM 07/20/11: This const_cast is not optimal, but the "compareTo" function
  // is not const as it should be.  This should be fixed in the future.
  //***********************************************************************
  if (0 == const_cast<G4String&>(name).compareTo("MISS", G4String::ignoreCase))
    filter = new CDMS_SelectMISS();

  return filter;
}
