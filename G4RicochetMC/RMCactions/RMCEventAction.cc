////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCEventAction.cc                                   //
//  Description: user event action class for CDMSMini                 //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        14 May 2010                                          //
//                                                                    //
////////////////////////////////////////////////////////////////////////

//#include "RMCactions/RMCEventAction.hh"
//#include "RMCactions/RMCRunAction.hh"
#include "RMCEventAction.hh"
#include "RMCRunAction.hh"
#include "G4Event.hh"


RMCEventAction::
RMCEventAction(RMCRunAction* runAct)
 :runAction(runAct)
{}

RMCEventAction::~RMCEventAction()
{}


void RMCEventAction::BeginOfEventAction(const G4Event* /*evt*/)
{}

 
void RMCEventAction::EndOfEventAction(const G4Event* evt)
{
  //runAction->transferEvent(evt);
}  
