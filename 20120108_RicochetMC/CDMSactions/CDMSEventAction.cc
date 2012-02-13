// $Id: CDMSEventAction.cc,v 1.2 2010/08/02 04:56:59 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSEventAction.cc                                   //
//  Description: user event action class for CDMSMini                 //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        14 May 2010                                          //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "CDMSactions/CDMSEventAction.hh"
#include "CDMSactions/CDMSRunAction.hh"
#include "G4Event.hh"


CDMSEventAction::
CDMSEventAction(CDMSRunAction* runAct)
 :runAction(runAct)
{}

CDMSEventAction::~CDMSEventAction()
{}


void CDMSEventAction::BeginOfEventAction(const G4Event* /*evt*/)
{}

 
void CDMSEventAction::EndOfEventAction(const G4Event* evt)
{
  runAction->transferEvent(evt);
}  
